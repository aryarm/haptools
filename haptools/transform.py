from __future__ import annotations
import logging
from pathlib import Path
from collections import namedtuple
from dataclasses import dataclass, field

import numpy as np
from cyvcf2 import VCF
import numpy.typing as npt
from pysam import VariantFile

from . import data
from .logging import getLogger


@dataclass
class HaplotypeAncestry(data.Haplotype):
    """
    A haplotype with an ancestry field for the transform subcommand

    Properties and functions are shared with the base "Haplotype" object
    """

    ancestry: str
    _extras: tuple = field(
        repr=False,
        init=False,
        default=(data.Extra("ancestry", "s", "Local ancestry"),),
    )

    def transform(self, genotypes: data.GenotypesVCF) -> npt.NDArray:
        """
        Transform a genotypes matrix via the current haplotype and its ancestral
        population

        See documentation for :py:meth:`~.Haplotype.transform` for more details
        """
        var_IDs = self.varIDs
        gts = genotypes.subset(variants=var_IDs)
        # check: were any of the variants absent from the genotypes?
        if len(gts.variants) < len(var_IDs):
            missing_IDs = set(var_IDs) - set(gts.variants["id"])
            raise ValueError(
                f"Variants {missing_IDs} are present in haplotype '{self.id}' but "
                "absent in the provided genotypes"
            )
        # create a np array denoting the alleles that we want
        # note: the excessive use of square-brackets gives us shape (1, p, 1)
        # where p denotes the number of alleles in this haplotype
        # That shape is broadcastable with gts.data which has shape (n, p, 2)
        allele_arr = np.array(
            [
                [
                    [int(var.allele != gts.variants[i]["alleles"][0])]
                    for i, var in enumerate(self.variants)
                ]
            ]
        )
        # look for the presence of each allele in each chromosomal strand
        # and then just AND them together
        hap_gts = np.all(allele_arr == gts.data, axis=1)
        # first, obtain the encoding of this haplotype's ancestry within the genotype
        # matrix. This will be an integer like 0, 1, 2, or 3 (or -1 if not found)
        ancestry_label = gts.ancestry_labels.get(self.ancestry, -1)
        # look for the presence of the desired ancestry in each chromosomal strand
        # and then just AND across all of the variants in the haplotype
        ancestry_arr = np.all(gts.ancestry == ancestry_label, axis=1)
        return np.logical_and(hap_gts, ancestry_arr)


class HaplotypesAncestry(data.Haplotypes):
    """
    A set of haplotypes with an ancestry field for the transform subcommand

    Properties and functions are shared with the base "Haplotypes" object
    """

    def __init__(
        self,
        fname: Path | str,
        haplotype: type[HaplotypeAncestry] = HaplotypeAncestry,
        variant: type[data.Variant] = data.Variant,
        log: logging.Logger = None,
    ):
        """
        Contrasting with the base Haplotypes class: this class uses HaplotypeAncestry
        as its default Haplotype class
        """
        super().__init__(fname, haplotype=haplotype, variant=variant, log=log)

    def transform(
        self,
        gts: data.GenotypesAncestry,
        hap_gts: data.GenotypesVCF = None,
    ) -> data.GenotypesVCF:
        self.index()
        haps = [self.data[hap] for hap in self.type_ids["H"]]
        # Initialize GenotypesVCF return value
        if hap_gts is None:
            hap_gts = data.GenotypesVCF(fname=None, log=self.log)
        hap_gts.samples = gts.samples
        hap_gts.variants = np.array(
            [(hap.id, hap.chrom, hap.start, ("A", "T")) for hap in haps],
            dtype=hap_gts.variants.dtype,
        )
        # build a fast data structure for querying the alleles in each haplotype:
        # a dict mapping (variant ID, allele) -> a unique index
        alleles = {}
        # and a list of arrays containing the indices of each hap's alleles
        idxs = [None] * len(haps)
        # and lastly, a list of ancestral population labels for each hap
        ancestries = np.empty(len(haps), dtype=np.uint8)
        count = 0
        for i, hap in enumerate(haps):
            ancestries[i] = gts.ancestry_labels.get(hap.ancestry, -1)
            idxs[i] = np.empty(len(hap.variants), dtype=np.uintc)
            for j, variant in enumerate(hap.variants):
                key = (variant.id, variant.allele)
                if key not in alleles:
                    alleles[key] = count
                    count += 1
                idxs[i][j] = alleles[key]
        self.log.debug(f"Copying genotypes for {len(alleles)} distinct alleles")
        gts = gts.subset(variants=tuple(k[0] for k in alleles))
        self.log.debug(f"Creating array denoting alt allele status")
        # initialize a np array denoting the allele integer in each haplotype
        # with shape (1, gts.data.shape[1], 1) for broadcasting later
        allele_arr = np.array(
            [
                int(allele != gts.variants[i]["alleles"][0])
                for i, (vID, allele) in enumerate(alleles)
            ],
            dtype=gts.data.dtype,
        )[np.newaxis, :, np.newaxis]
        # finally, obtain and merge the haplotype genotypes
        self.log.info(f"Transforming genotypes for {len(haps)} haplotypes")
        equality_arr = np.equal(allele_arr, gts.data)
        self.log.debug(
            f"Allocating array with dtype {gts.data.dtype} and size "
            f"{(len(gts.samples), len(haps), 2)}"
        )
        hap_gts.data = np.empty((gts.data.shape[0], len(haps), 2), dtype=np.bool_)
        self.log.debug("Computing haplotype genotypes. This may take a while")
        for i in range(len(haps)):
            hap_gts.data[:, i] = np.logical_and(
                np.all(gts.ancestry[:, idxs[i]] == ancestries[i], axis=1),
                np.all(equality_arr[:, idxs[i]], axis=1),
            )
        return hap_gts


class GenotypesAncestry(data.GenotypesVCF):
    """
    Extends the GenotypesVCF class for ancestry data

    The ancestry information is stored within the FORMAT field of the VCF

    Attributes
    ----------
    data : np.array
        See documentation for :py:attr:`~.Genotypes.data`
    fname : Path | str
        See documentation for :py:attr:`~.Genotypes.fname`
    samples : tuple[str]
        See documentation for :py:attr:`~.Genotypes.samples`
    variants : np.array
        See documentation for :py:attr:`~.GenotypesVCF.variants`
    valid_labels: np.array
        Reference VCF sample and respective variant grabbed for
        each sample.
    ancestry : np.array
        The ancestral population of each allele in each sample of
        :py:attr:`~.GenotypesAncestry.data`
    log: logging.Logger
        See documentation for :py:attr:`~.Genotypes.log`
    """

    def __init__(self, fname: Path | str, log: logging.Logger = None):
        super().__init__(fname, log)
        self.ancestry = None
        self.valid_labels = None
        # goes from population code to encoding number
        self.ancestry_labels = {}
        # goes from encoding number to population code
        self.popnum_ancestry = {}

    def _iterate(self, vcf: VCF, region: str = None, variants: set[str] = None):
        """
        See documentation for :py:meth:`~.Genotypes._iterate`
        """
        self.log.info(f"Loading genotypes from {len(self.samples)} samples")
        Record = namedtuple("Record", "data ancestry variants")
        num_seen = 0
        pop_count = 0
        # iterate over each line in the VCF
        # note, this can take a lot of time if there are many samples
        for variant in vcf(region):
            if variants is not None and variant.ID not in variants:
                if num_seen >= len(variants):
                    # exit early if we've already found all the variants
                    break
                continue
            # save meta information about each variant
            variant_arr = self._variant_arr(variant)
            # extract the genotypes to a matrix of size n x 3
            # the last dimension has three items:
            # 1) presence of REF in strand one
            # 2) presence of REF in strand two
            # 3) whether the genotype is phased (if self._prephased is False)
            data = np.array(variant.genotypes, dtype=np.uint8)
            data = data[:, : (2 + (not self._prephased))]
            # also extract the ancestral population of each variant in each individual
            ancestry = np.empty((data.shape[0], 2), dtype=np.uint8)
            for i, sample in enumerate(variant.format("POP")):
                pops = sample.split(",")
                for pop in pops:
                    if pop not in self.ancestry_labels:
                        self.ancestry_labels[pop] = pop_count
                        self.popnum_ancestry[pop_count] = pop
                        pop_count += 1
                ancestry[i] = tuple(map(self.ancestry_labels.get, pops))
            # finally, output everything
            yield Record(data, ancestry, variant_arr)
            num_seen += 1
        vcf.close()

    def read(
        self,
        region: str = None,
        samples: set[str] = None,
        variants: set[str] = None,
        max_variants: int = None,
    ):
        """
        See documentation for :py:meth:`~.Genotypes.read`
        """
        super(data.Genotypes, self).read()
        records = self.__iter__(region=region, samples=samples, variants=variants)
        if variants is not None:
            max_variants = len(variants)
        # check whether we can preallocate memory instead of making copies
        if max_variants is None:
            self.log.warning(
                "The max_variants parameter was not specified. We have no choice but to"
                " append to an ever-growing array, which can lead to memory overuse!"
            )
            variants_arr = []
            data_arr = []
            ancestry_arr = []
            for rec in records:
                variants_arr.append(rec.variants)
                data_arr.append(rec.data)
                ancestry_arr.append(rec.ancestry)
            self.log.info(f"Copying {len(variants_arr)} variants into np arrays.")
            # convert to np array for speedy operations later on
            self.variants = np.array(variants_arr, dtype=self.variants.dtype)
            self.data = np.array(data_arr, dtype=np.uint8)
            self.ancestry = np.array(ancestry_arr, dtype=np.uint8)
        else:
            # preallocate arrays! this will save us lots of memory and speed b/c
            # appends can sometimes make copies
            self.variants = np.empty((max_variants,), dtype=self.variants.dtype)
            # in order to check_phase() later, we must store the phase info, as well
            self.data = np.empty(
                (max_variants, len(self.samples), (2 + (not self._prephased))),
                dtype=np.uint8,
            )
            self.ancestry = np.empty(
                (max_variants, len(self.samples), 2),
                dtype=np.uint8,
            )
            self.valid_labels = np.empty(
                (max_variants, len(self.samples), 2),
                dtype=object,
            )
            num_seen = 0
            for rec in records:
                if num_seen >= max_variants:
                    break
                self.variants[num_seen] = rec.variants
                self.data[num_seen] = rec.data
                self.ancestry[num_seen] = rec.ancestry
                num_seen += 1
            if max_variants > num_seen:
                self.log.info(
                    f"Removing {max_variants-num_seen} unneeded variant records that "
                    "were preallocated b/c max_variants was specified."
                )
                self.variants = self.variants[:num_seen]
                self.data = self.data[:num_seen]
                self.ancestry = self.ancestry[:num_seen]
        if 0 in self.data.shape:
            self.log.warning(
                "Failed to load genotypes. If you specified a region, check that the"
                " contig name matches! For example, double-check the 'chr' prefix."
            )
        # transpose the GT matrix so that samples are rows and variants are columns
        self.log.info(f"Transposing genotype matrix of size {self.data.shape}.")
        self.data = self.data.transpose((1, 0, 2))
        self.ancestry = self.ancestry.transpose((1, 0, 2))

    def subset(
        self,
        samples: tuple[str] = None,
        variants: tuple[str] = None,
        inplace: bool = False,
    ):
        """
        See documentation for :py:meth:`~.Genotypes.subset`
        """
        # First, initialize variables
        gts = self
        if not inplace:
            gts = self.__class__(self.fname, self.log)
        gts.samples = self.samples
        gts.variants = self.variants
        gts.data = self.data
        gts.ancestry = self.ancestry
        gts.ancestry_labels = self.ancestry_labels
        # Index the current set of samples and variants so we can have fast look-up
        self.index(samples=(samples is not None), variants=(variants is not None))
        # Subset the samples
        if samples is not None:
            gts.samples = tuple(samp for samp in samples if samp in self._samp_idx)
            if len(gts.samples) < len(samples):
                diff = len(samples) - len(gts.samples)
                self.log.warning(
                    f"Saw {diff} fewer samples than requested. Proceeding with "
                    f"{len(gts.samples)} samples."
                )
            samp_idx = tuple(self._samp_idx[samp] for samp in gts.samples)
            if inplace:
                self._samp_idx = None
            gts.data = gts.data[samp_idx, :]
            gts.ancestry = gts.ancestry[samp_idx, :]
        # Subset the variants
        if variants is not None:
            var_idx = [self._var_idx[var] for var in variants if var in self._var_idx]
            if len(var_idx) < len(variants):
                diff = len(variants) - len(var_idx)
                self.log.warning(
                    f"Saw {diff} fewer variants than requested. Proceeding with "
                    f"{len(var_idx)} variants."
                )
            gts.variants = self.variants[var_idx]
            if inplace:
                self._var_idx = None
            gts.data = gts.data[:, var_idx]
            gts.ancestry = gts.ancestry[:, var_idx]
        if not inplace:
            return gts

    def check_missing(self, discard_also=False):
        """
        See documentation for :py:meth:`~.Genotypes.check_missing`
        """
        # check: are there any samples that have genotype values that are empty?
        # A genotype value equal to the max for uint8 indicates the value was missing
        missing = np.any(self.data[:, :, :2] == np.iinfo(np.uint8).max, axis=2)
        if np.any(missing):
            samp_idx, variant_idx = np.nonzero(missing)
            if discard_also:
                original_num_samples = len(self.samples)
                self.data = np.delete(self.data, samp_idx, axis=0)
                self.ancestry = np.delete(self.ancestry, samp_idx, axis=0)
                self.samples = tuple(np.delete(self.samples, samp_idx))
                self.log.info(
                    "Ignoring missing genotypes from "
                    f"{original_num_samples - len(self.samples)} samples"
                )
                self._samp_idx = None
            else:
                raise ValueError(
                    "Genotype with ID {} at POS {}:{} is missing for sample {}".format(
                        *tuple(self.variants[variant_idx[0]])[:3],
                        self.samples[samp_idx[0]],
                    )
                )
        if discard_also and not self.data.shape[0]:
            self.log.warning(
                "All samples were discarded! Check that that none of your variants are"
                " missing genotypes (GT: '.|.')."
            )

    def check_biallelic(self, discard_also=False):
        """
        See documentation for :py:meth:`~.Genotypes.check_biallelic`
        """
        if self.data.dtype == np.bool_:
            self.log.warning("All genotypes are already biallelic")
            return
        # check: are there any variants that have genotype values above 1?
        # A genotype value above 1 would imply the variant has more than one ALT allele
        multiallelic = np.any(self.data[:, :, :2] > 1, axis=2)
        if np.any(multiallelic):
            samp_idx, variant_idx = np.nonzero(multiallelic)
            if discard_also:
                self.log.info(f"Ignoring {len(variant_idx)} multiallelic variants")
                self.data = np.delete(self.data, variant_idx, axis=1)
                self.ancestry = np.delete(self.ancestry, variant_idx, axis=1)
                self.variants = np.delete(self.variants, variant_idx)
                self._var_idx = None
            else:
                raise ValueError(
                    "Variant with ID {} at POS {}:{} is multiallelic for sample {}"
                    .format(
                        *tuple(self.variants[variant_idx[0]])[:3],
                        self.samples[samp_idx[0]],
                    )
                )
        if discard_also and not self.data.shape[1]:
            self.log.warning(
                "All variants were discarded! Check that there are biallelic variants "
                "in your dataset."
            )
        self.data = self.data.astype(np.bool_)

    def write(self, chroms=None):
        # Assumption is the data must be phased
        """
        Write the variants in this class to a VCF at :py:attr:`~.GenotypesAncestry.fname`
        """
        vcf = VariantFile(str(self.fname), mode="w")

        # make sure the header is properly structured
        if not chroms:
            for contig in set(self.variants["chrom"]):
                vcf.header.contigs.add(contig)
        else:
            # make sure the header is properly structured with contig names from ref VCF
            for contig in set(self.variants["chrom"]):
                # remove chr in front of seqname if present and compare
                if contig.startswith("chr"):
                    if contig[3:] in chroms:
                        vcf.header.contigs.add(contig)
                if contig in chroms:
                    vcf.header.contigs.add(contig)

        vcf.header.add_meta(
            "FORMAT",
            items=[
                ("ID", "GT"),
                ("Number", 1),
                ("Type", "String"),
                ("Description", "Genotype"),
            ],
        )
        if not self.ancestry is None:
            vcf.header.add_meta(
                "FORMAT",
                items=[
                    ("ID", "POP"),
                    ("Number", 2),
                    ("Type", "String"),
                    (
                        "Description",
                        "Origin Population of each respective allele in GT",
                    ),
                ],
            )
        if not self.valid_labels is None:
            vcf.header.add_meta(
                "FORMAT",
                items=[
                    ("ID", "SAMPLE"),
                    ("Number", 2),
                    ("Type", "String"),
                    (
                        "Description",
                        "Origin sample and haplotype of each respective allele in GT",
                    ),
                ],
            )
        try:
            vcf.header.add_samples(self.samples)
        except AttributeError:
            self.log.warning(
                "Upgrade to pysam >=0.19.1 to reduce the time required to create "
                "VCFs. See https://github.com/pysam-developers/pysam/issues/1104"
            )
            for sample in self.samples:
                vcf.header.add_sample(sample)
        self.log.info("Writing VCF records")
        phased = self._prephased or (self.data.shape[2] < 3)
        for var_idx, var in enumerate(self.variants):
            rec = {
                "contig": var["chrom"],
                "start": var["pos"],
                "stop": var["pos"] + len(var["alleles"][0]) - 1,
                "qual": None,
                "alleles": var["alleles"],
                "id": var["id"],
                "filter": None,
            }
            # handle pysam increasing the start site by 1
            rec["start"] -= 1
            # parse the record into a pysam.VariantRecord
            record = vcf.new_record(**rec)
            for samp_idx, sample in enumerate(self.samples):
                # TODO: make this work when there are missing values
                record.samples[sample]["GT"] = tuple(self.data[samp_idx, var_idx, :2])
                if not self.ancestry is None:
                    record.samples[sample]["POP"] = tuple(
                        map(
                            self.popnum_ancestry.get,
                            self.ancestry[samp_idx, var_idx, :],
                        )
                    )
                if not self.valid_labels is None:
                    record.samples[sample]["SAMPLE"] = tuple(
                        self.valid_labels[samp_idx, var_idx, :]
                    )
                # add proper phasing info
                if phased:
                    record.samples[sample].phased = True
                else:
                    record.samples[sample].phased = self.data[samp_idx, var_idx, 2]
            # write the record to a file
            vcf.write(record)
        vcf.close()

    def merge_variants(
        cls, objs: tuple[data.Genotypes], check_samples: bool = True, **kwargs
    ) -> data.Genotypes:
        """
        See documentation for :py:meth:`~.data.Genotypes.merge_variants`
        """
        raise NotImplementedError


def transform_haps(
    genotypes: Path,
    haplotypes: Path,
    region: str = None,
    samples: set[str] = None,
    haplotype_ids: set[str] = None,
    chunk_size: int = None,
    discard_missing: bool = False,
    ancestry: bool = False,
    output: Path = Path("-"),
    log: logging.Logger = None,
):
    """
    Creates a VCF composed of haplotypes

    Parameters
    ----------
    genotypes : Path
        The path to the genotypes in VCF or PGEN format
    haplotypes : Path
        The path to the haplotypes in a .hap file
    region : str, optional
        See documentation for :py:meth:`~.data.Genotypes.read`
        and :py:meth:`~.data.Haplotypes.read`
    samples : set[str], optional
        See documentation for :py:meth:`~.data.Genotypes.read`
    haplotype_ids: set[str], optional
        A set of haplotype IDs to obtain from the .hap file. All others are ignored.

        If not provided, all haplotypes will be used.
    chunk_size: int, optional
        The max number of variants to fetch from the PGEN file at any given time

        If this value is provided, variants from the PGEN file will be loaded in
        chunks so as to use less memory. This argument is ignored if the genotypes are
        not in PGEN format.
    discard_missing : bool, optional
        Discard any samples that are missing any of the required genotypes

        The default is simply to complain about it
    ancestry : bool, optional
        Whether to also match ancestral population labels from the VCF against those in
        the .hap file
    output : Path, optional
        The location to which to write output
    log : Logger, optional
        A logging module to which to write messages about progress and any errors
    """
    if log is None:
        log = getLogger(name="transform", level="ERROR")

    haps_class = HaplotypesAncestry if ancestry else data.Haplotypes
    log.info("Loading haplotypes")
    hp = haps_class(haplotypes, log=log)
    hp.read(region=region, haplotypes=haplotype_ids)

    # check that all of the haplotypes were loaded successfully and warn otherwise
    if haplotype_ids is not None and len(haplotype_ids) > len(hp.data):
        diff = list(haplotype_ids.difference(hp.data.keys()))
        first_few = 5 if len(diff) > 5 else len(diff)
        log.warning(
            f"{len(diff)} haplotypes could not be found in the .hap file. Check "
            "that the IDs in your .hap file correspond with those you provided. "
            f"Here are the first few missing haplotypes: {diff[:first_few]}"
        )
    if len(hp.data) == 0:
        raise ValueError("Didn't load any haplotypes from the .hap file")

    log.info("Extracting variants from haplotypes")
    variants = {vr.id for id in hp.type_ids["H"] for vr in hp.data[id].variants}

    # load the genotypes, but first get the path to the breakpoints file
    if genotypes.suffix == ".gz":
        bps_file = genotypes.with_suffix("").with_suffix(".bp")
    else:
        bps_file = genotypes.with_suffix(".bp")
    # now, get the genotypes
    if genotypes.suffix == ".pgen":
        log.info("Loading genotypes from PGEN file")
        gt = data.GenotypesPLINK(fname=genotypes, log=log, chunk_size=chunk_size)
    else:
        log.info("Loading genotypes from VCF/BCF file")
        if ancestry and not bps_file.exists():
            gt = GenotypesAncestry(fname=genotypes, log=log)
        else:
            gt = data.GenotypesVCF(fname=genotypes, log=log)
    # gt._prephased = True
    gt.read(region=region, samples=samples, variants=variants)
    gt.check_missing(discard_also=discard_missing)
    gt.check_phase()

    # check that all of the variants were loaded successfully and warn otherwise
    if len(variants) > len(gt.variants):
        diff = list(variants.difference(gt.variants["id"]))
        first_few = 5 if len(diff) > 5 else len(diff)
        log.warning(
            f"{len(diff)} variant(s) could not be found in the genotypes file. Check "
            "that the IDs in your .hap file correspond with those in the genotypes "
            f"file. Here are the first few missing variants: {diff[:first_few]}"
        )
        # subset the set of haplotypes so that we keep only those that we can transform
        gt_variants = set(gt.variants["id"])
        original_num_haps = len(hp.data)
        haplotype_ids = tuple(
            hap_id
            for hap_id, hap in hp.data.items()
            if gt_variants.issuperset(hap.varIDs)
        )
        hp.subset(haplotypes=haplotype_ids, inplace=True)
        log.info(f"Proceeding with {len(hp.data)} of {original_num_haps} haplotypes")

    if ancestry and not isinstance(gt, GenotypesAncestry):
        log.info("Loading ancestry info from .bp file")
        if not bps_file.exists():
            raise ValueError("A .bp file is needed when using --ancestry")
        bps = data.Breakpoints(fname=bps_file, log=log)
        bps.read(samples=set(gt.samples))
        bps.encode()
        # convert the GenotypesVCF object to a GenotypesAncestry object
        # TODO: figure out a better solution for this
        # this is just a temp hack to get output from simgenotype to load a bit faster
        gta = GenotypesAncestry(fname=None, log=log)
        gta.data = gt.data
        gta.samples = gt.samples
        gta.variants = gt.variants
        gta.ancestry_labels = bps.labels
        gta.ancestry = bps.population_array(gt.variants[["chrom", "pos"]])
        gt = gta

    if output.suffix == ".pgen":
        out_file_type = "PGEN"
        hp_gt = data.GenotypesPLINK(fname=output, log=log, chunk_size=chunk_size)
    else:
        out_file_type = "VCF/BCF"
        hp_gt = data.GenotypesVCF(fname=output, log=log)
    log.info("Transforming genotypes via haplotypes")
    hp.transform(gt, hp_gt)

    log.info(f"Writing haplotypes to {out_file_type} file")
    hp_gt.write()

    log.debug("Done!")
    return hp_gt
