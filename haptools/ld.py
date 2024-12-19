from __future__ import annotations
import logging
from pathlib import Path
from dataclasses import dataclass, field

import numpy as np
import numpy.typing as npt

from . import data
from .logging import getLogger
from .data import Haplotype as HaplotypeBase


@dataclass
class Haplotype(HaplotypeBase):
    """
    A haplotype with sufficient fields for the ld command

    Properties and functions are shared with the base Haplotype object,
    :py:attr:`~.HaplotypeBase`
    """

    ld: float
    _extras: tuple = field(
        repr=False,
        init=False,
        default=(data.Extra("ld", ".3f", "Linkage-disequilibrium"),),
    )


def pearson_corr_ld(arrA: npt.NDArray, arrB: npt.NDArray) -> float | npt.NDArray:
    """
    Compute the Pearson correlation coefficient (LD) between two vectors (1D arrays)

    If 2D array(s) are given instead, the rows will be treated as samples and the columns
    as variants

    Parameters
    ----------
    arrA: npt.NDArray
        The first 1D (or 2D) numpy array
    arrB: npt.NDArray
        The second 1D (or 2D) numpy array

    Returns
    -------
    The signed LD between the genotypes in arrA and the genotypes in arrB

    If either arrA or arrB is 2D, then the return value will be an array instead. If
    both are 2D, the shape of the array will correspond with the second dimensions of
    arrA and arrB, respectively. Otherwise, it will be 1D.
    """
    if arrA.ndim == 1:
        dim = 1
    elif arrA.ndim == 2:
        dim = arrA.shape[1]
    else:
        raise ValueError("We can only compute LD for 2D genotype arrays")
    ld_mat = np.corrcoef(arrA, arrB, rowvar=False)[:dim:, dim:]
    if arrA.ndim == 1:
        ld_mat = ld_mat[0, :]
        if arrB.ndim == 1:
            ld_mat = ld_mat[0]
    elif arrB.ndim == 1:
        ld_mat = ld_mat[:, 0]
    return ld_mat


def calc_ld(
    target: str,
    genotypes: Path,
    haplotypes: Path,
    region: str = None,
    samples: set[str] = None,
    ids: tuple[str] = None,
    chunk_size: int = None,
    discard_missing: bool = False,
    from_gts: bool = False,
    output: Path = Path("/dev/stdout"),
    log: logging.Logger = None,
):
    """
    Creates a VCF composed of haplotypes

    Parameters
    ----------
    target : str
        The ID of the haplotype or variant with which we will calculate LD
    genotypes : Path
        The path to the genotypes
    haplotypes : Path
        The path to the haplotypes in a .hap file
    region : str, optional
        See documentation for :py:meth:`~.data.Genotypes.read`
        and :py:meth:`~.data.Haplotypes.read`
    samples : set[str], optional
        See documentation for :py:meth:`~.data.Genotypes.read`
    ids: set[str], optional
        A subset of haplotype IDs to obtain from the .hap file. All others
        are ignored.

        Alternatively, if the --from-gts switch is specified, this will be interpreted
        as a subset of variant IDs to obtain from the genotypes file.

        Defaults to loading all haplotypes or variants if not specified
    chunk_size: int, optional
        The max number of variants to fetch from the PGEN file at any given time

        If this value is provided, variants from the PGEN file will be loaded in
        chunks so as to use less memory. This argument is ignored if the genotypes are
        not in PGEN format.
    discard_missing : bool, optional
        Discard any samples that are missing any of the required genotypes

        The default is simply to complain about it
    output : Path, optional
        The location to which to write output
    log : Logger, optional
        A logging module to which to write messages about progress and any errors
    """
    if log is None:
        log = getLogger(name="ld", level="ERROR")

    # convert IDs to set but save the tuple
    ids_tup, ids = ids, (set(ids) if ids is not None else None)

    log.info("Loading haplotypes")
    hp = data.Haplotypes(haplotypes, log=log)
    haplotype_ids = None
    if not from_gts:
        haplotype_ids = ids
        if haplotype_ids is not None:
            haplotype_ids.add(target)
    hp.read(region=region, haplotypes=haplotype_ids)

    # remove all repeats from the haplotypes object since we don't yet support them
    for repeat_id in hp.type_ids["R"]:
        del hp.data[repeat_id]
    num_repeats = len(hp.type_ids["R"])
    if num_repeats:
        log.info(f"Ignoring {num_repeats} repeats in .hap file")
        hp.type_ids["R"] = []

    if from_gts:
        variants = None
        if target in hp.data and ids:
            variants = ids.copy()
            log.info("Extracting variants from haplotypes")
            variants.update(var.id for var in hp.data[target].variants)
    else:
        log.info("Extracting variants from haplotypes")
        variants = {var.id for h in hp.type_ids["H"] for var in hp.data[h].variants}

    # check to see whether the target was a haplotype
    try:
        target = hp.data.pop(target)
    except:
        # the target is a variant, instead
        pass
    else:
        log.info(f"Identified target '{target}' as a haplotype")
        hp.index(force=True)
        if len(hp.data) == 0 and not from_gts:
            log.error(
                "There must be at least one more haplotype in the .hap file "
                "than the TARGET haplotype specified."
            )

    # check that all of the haplotypes were loaded successfully and warn otherwise
    if ids is not None and not from_gts and len(ids) > len(hp.data):
        diff = list(ids.difference(hp.data.keys()))
        first_few = 5 if len(diff) > 5 else len(diff)
        log.warning(
            f"{len(diff)} haplotypes could not be found in the .hap file. Check "
            "that the IDs in your .hap file correspond with those you provided. "
            f"Here are the first few missing haplotypes: {diff[:first_few]}"
        )

    # if the target was not a haplotype, then it must be a variant
    # so we must load it from the genotype file
    if not isinstance(target, data.Haplotype):
        variants.add(target)

    if genotypes.suffix == ".pgen":
        log.info("Loading genotypes from PGEN file")
        gt = data.GenotypesPLINK(fname=genotypes, log=log, chunk_size=chunk_size)
    else:
        log.info("Loading genotypes from VCF/BCF file")
        gt = data.GenotypesVCF(fname=genotypes, log=log)
    # gt._prephased = True
    gt.read(region=region, samples=samples, variants=variants)
    gt.check_missing(discard_also=discard_missing)
    gt.check_biallelic()
    gt.check_phase()

    # check that all of the variants were loaded successfully and warn otherwise
    if variants and len(variants) < len(gt.variants):
        # check to see whether the target got loaded if it wasn't a haplotype
        if not (isinstance(target, data.Haplotype) or (target in gt.variants["id"])):
            raise ValueError(
                "Could not find the provided target ID among either the haplotypes "
                "in the .hap file or the variants in the genotype file. Check that "
                f"'{target}' appears in either the .hap file or the genotype file."
            )
        log.info(f"Identified target '{target}' as a variant")
        # report the missing variants
        diff = list(variants.difference(gt.variants["id"]))
        first_few = 5 if len(diff) > 5 else len(diff)
        log.warning(
            f"{len(diff)} variants could not be found in the genotypes file. Check "
            "that the IDs in your .hap file correspond with those in the genotypes "
            f"file. Here are the first few missing variants: {diff[:first_few]}"
        )

    if not from_gts:
        log.info("Transforming genotypes via haplotypes")
        hp_gt = data.GenotypesVCF(fname=None, log=log)
        hp.transform(gt, hp_gt)

    log.info("Obtaining target genotypes")
    if isinstance(target, data.Haplotype):
        target_gts = target.transform(gt).sum(axis=1)
        if from_gts and ids is not None:
            gt.subset(variants=ids_tup, inplace=True)
    else:
        target_gts = gt.subset(variants=(target,)).data[:, 0, :2].sum(axis=1)

    if from_gts:
        log.info("Computing LD between genotypes and the target")
        with data.Data.hook_compressed(output, mode="w") as ld_file:
            log.info("Outputting .ld file with LD values")
            ld_file.write("CHR\tBP\tSNP\tR\n")
            for idx, variant in enumerate(gt.variants[["chrom", "pos", "id"]]):
                var_chr, var_bp, var_snp = variant
                variant_gts = gt.data[:, idx, :2].sum(axis=1)
                variant_ld = pearson_corr_ld(target_gts, variant_gts)
                ld_file.write(f"{var_chr}\t{var_bp}\t{var_snp}\t{variant_ld:.3f}\n")
    else:
        log.info("Computing LD between haplotypes and the target")
        # construct a new Haplotypes object that also stores the LD values
        hp_out = data.Haplotypes(fname=output, haplotype=Haplotype, log=log)
        hp_out.data = {}
        for hap_id in hp.type_ids["H"]:
            # break the BaseHaplotype instance up into its properties
            hapd = hp.data[hap_id].__dict__
            hapd_variants = hapd.pop("variants")
            # obtain genotypes for the current haplotype
            hap_gts = hp_gt.subset(variants=(hap_id,)).data[:, 0, :2].sum(axis=1)
            # compute the LD between the current haplotype and the target
            hapd["ld"] = pearson_corr_ld(target_gts, hap_gts)
            # create the Haplotype instance and add the variants in
            hp_out.data[hap_id] = Haplotype(**hapd)
            hp_out.data[hap_id].variants = hapd_variants
        log.info("Outputting .hap file with LD values")
        hp_out.write()
