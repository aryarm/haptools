from __future__ import annotations
import logging
from pathlib import Path
from dataclasses import dataclass, field

import numpy as np
import numpy.typing as npt

from .logging import getLogger
from .data import (
    Extra,
    Genotypes,
    Phenotypes,
    Haplotypes,
    GenotypesTR,
    GenotypesPLINK,
    GenotypesPLINKTR,
    Repeat as RepeatBase,
    Haplotype as HaplotypeBase,
)


@dataclass
class Effect:
    """
    A variable in the simphenotype linear model

    Attributes
    ----------
    id : str
        The ID of the variable; corresponds to a variant in a Genotypes object
    beta : float
        The effect size of the variable
    """

    id: str
    beta: float

    @classmethod
    def from_hap_spec(cls: Effect, line: str) -> Effect:
        """
        Convert a .snplist line into an Effect object

        Parameters
        ----------
        line: str
            A line from a .snplist file

        Returns
        -------
        Effect
            The converted Effect instance
        """
        ID, beta = line.split("\t")[:2]
        return cls(id=ID, beta=float(beta))


@dataclass
class Haplotype(HaplotypeBase):
    """
    A haplotype with sufficient fields for simphenotype

    Properties and functions are shared with the base Haplotype object, "HaplotypeBase"
    """

    beta: float
    _extras: tuple = field(
        repr=False,
        init=False,
        default=(Extra("beta", ".2f", "Effect size in linear model"),),
    )


@dataclass
class Repeat(RepeatBase):
    """
    A repeat with sufficient fields for simphenotype

    Properties and functions are shared with the base Repeat object, "RepeatBeta"
    """

    beta: float
    _extras: tuple = field(
        repr=False,
        init=False,
        default=(Extra("beta", ".2f", "Effect size in linear model"),),
    )


class PhenoSimulator:
    """
    Simulate phenotypes from genotypes

    Attributes
    ----------
    gens: dict[Genotypes]
        Genotypes to simulate
    phens: Phenotypes
        Simulated phenotypes; filled by :py:meth:`~.PhenoSimular.run`
    rng: np.random.Generator, optional
        A numpy random number generator
    log: logging.Logger
        A logging instance for recording debug statements

    Examples
    --------
    >>> gens = Genotypes.load("tests/data/example.vcf.gz")
    >>> tr_gens = GenotypesTR.load("tests/data/simple_tr.vcf")
    >>> haps = Haplotypes.load("tests/data/basic.hap")
    >>> haps_gts = haps.transform(gens)
    >>> phenosim = PhenoSimulator(haps_gts)
    >>> phenosim.run(haps.data.values())
    >>> phenotypes = phenosim.phens
    """

    def __init__(
        self,
        genotypes: Genotypes,
        output: Path = Path("/dev/stdout"),
        seed: int = None,
        log: logging.Logger = None,
    ):
        """
        Initialize a PhenoSimulator object

        Parameters
        ----------
        genotypes: Genotypes
            Genotypes for each haplotype
        output: Path, optional
            Path to a '.pheno' file to which the generated phenotypes could be written

            Defaults to stdout if not provided
        seed: int, optional
            A seed to initialize the random number generator

            This is useful if you want the generated phenotypes to be the same across
            multiple PhenoSimulator instances. If not provided, it will be random.
        log: logging.Logger, optional
            A logging instance for recording debug statements
        """
        self.gens = genotypes
        self.phens = Phenotypes(fname=output)
        self.phens.data = None
        self.phens.samples = self.gens.samples
        self.rng = np.random.default_rng(seed)
        self.log = log or logging.getLogger(self.__class__.__name__)

    def run(
        self,
        effects: list[Effect | Haplotype | Repeat],
        heritability: float = None,
        prevalence: float = None,
        normalize: bool = True,
        environment: float = None,
    ) -> npt.NDArray:
        """
        Simulate phenotypes for an entry in the Genotypes object

        The generated phenotypes will also be added to
        :py:attr:`~.PhenoSimulator.output`

        Parameters
        ----------
        effects: list[Effect|Haplotype|Repeat]
            A list of Haplotypes to use in an additive fashion within the simulations
        heritability: float, optional
            The simulated heritability of the trait

            If not provided, this will default to the sum of the squared effect sizes
        prevalence: float, optional
            How common should the disease be within the population?

            If this value is specified, case/control phenotypes will be generated
            instead of quantitative traits.
        normalize: bool, optional
            If True, normalize the genotypes before using them to simulate the
            phenotypes. Otherwise, use the raw values.
        environment: float, optional
            The variance (aka strength) of the environmental contribution to the trait.
            This is inferred from the betas if it isn't specified.

        Returns
        -------
        npt.NDArray
            The simulated phenotypes, as a np array of shape num_samples x 1
        """
        # extract the relevant haplotype info from the Haplotype objects
        ids = [effect.id for effect in effects]
        betas = np.array([effect.beta for effect in effects])
        self.log.debug(f"Beta values are {betas}")
        self.log.debug(f"Extracting haplotype genotypes for haps: {ids}")
        gts = self.gens.subset(variants=ids).data[:, :, :2].sum(axis=2)

        if normalize:
            gts = self.normalize_gts(gts, ids)

        self.log.info(f"Computing genetic component w/ {gts.shape[1]} causal effects")

        # generate the genetic component
        pt = (betas * gts).sum(axis=1)
        # compute the heritability
        if heritability is None and environment is None:
            self.log.debug("Computing heritability as the sum of the squared betas")
            heritability = np.power(betas, 2).sum()
            if heritability > 1:
                self.log.warning(
                    "Variance of error term exceeds 1. Check your betas! Capping at 1 "
                    "for now."
                )
                heritability = 1
            # compute the environmental effect
            noise = 1 - heritability
        else:
            noise = environment
            if environment is None:
                # compute the environmental effect
                noise = np.var(pt)
                if noise == 0:
                    self.log.warning(
                        "Your genotypes have 0 variance. Creating artificial noise..."
                    )
                    noise = 1
            elif heritability is None:
                heritability = 0.5
            # TODO: handle a heritability of 0 somehow?
            noise *= np.reciprocal(heritability) - 1
        self.log.info(f"Adding environmental component {noise} for h^2 {heritability}")
        # finally, add everything together to get the simulated phenotypes
        pt_noise = self.rng.normal(0, np.sqrt(noise), size=pt.shape)
        if self.log.getEffectiveLevel() == logging.DEBUG:
            # if we're in debug mode, compute the pearson correlation and report it
            # but don't do this otherwise to keep things fast
            corr = np.corrcoef(pt, pt + pt_noise)[1, 0]
            self.log.debug(f"Estimated heritability is {corr}")
        pt += pt_noise
        # now, handle case/control
        name_suffix = ""
        if prevalence is not None:
            self.log.info(f"Converting to case/control with prevalence {prevalence}")
            # first, find the number of desired positives
            k = int(prevalence * len(pt))
            # choose the top k values and label them positive
            bool_pt = np.zeros(len(pt), dtype=np.bool_)
            if k == len(pt):
                bool_pt = ~bool_pt
            else:
                max_indices = np.argpartition(-pt, k)[:k]
                bool_pt[max_indices] = True
            pt = bool_pt
            name_suffix = "-cc"
        # now, save the archived phenotypes for later
        self.phens.append(name="-".join(ids) + name_suffix, data=pt.astype(np.float64))
        return pt

    def normalize_gts(self, gts, ids):
        """
        Normalize variant or repeats genotypes

        Parameters
        ----------
        gts: np.array
            Genotypes variant array stored in data.
        ids: list[str]
            IDs for variants

        Returns
        -------
        normalized_gts: np.array
            Normalized Genotypes variant array.
        """
        std = gts.std(axis=0)
        gts = (gts - gts.mean(axis=0)) / std
        # when the stdev is 0, just set all values to 0 instead of nan
        zero_elements = std == 0
        num_zero_elements = np.sum(zero_elements)
        if num_zero_elements:
            # get the first five causal variables with variances == 0
            zero_elements_ids = np.array(ids)[zero_elements]
            if len(zero_elements_ids) > 5:
                zero_elements_ids = zero_elements_ids[:5]
            self.log.warning(
                "Some of your causal variables have genotypes with variance 0. "
                f"Here are the first few: {zero_elements_ids}"
            )
        gts[:, zero_elements] = np.zeros((gts.shape[0], num_zero_elements))
        return gts

    def write(self):
        """
        Write the generated phenotypes to the file specified in
        :py:meth:`~.PhenoSimular.__init__`
        """
        self.phens.write()


def simulate_pt(
    genotypes: Path,
    haplotypes: Path,
    num_replications: int = 1,
    environment: float = None,
    heritability: float = None,
    prevalence: float = None,
    normalize: bool = True,
    region: str = None,
    samples: set[str] = None,
    haplotype_ids: set[str] = None,
    chunk_size: int = None,
    repeats: Path = None,
    seed: int = None,
    output: Path = Path("-"),
    log: logging.Logger = None,
):
    """
    Haplotype-aware phenotype simulation. Create a set of simulated phenotypes from a
    set of haplotypes.

    GENOTYPES must be formatted as a VCF or PGEN file and HAPLOTYPES must be formatted
    according to the .hap format spec

    Note: GENOTYPES must be the output from the the transform subcommand.

    \f
    Examples
    --------
    >>> haptools simphenotype tests/data/example.vcf.gz tests/data/example.hap.gz > simu_phens.tsv

    Parameters
    ----------
    genotypes : Path
        The path to the transformed genotypes in VCF or PGEN format
    haplotypes : Path
        The path to the haplotypes in a .hap file
    replications : int, optional
        The number of rounds of simulation to perform
    heritability : int, optional
        The heritability of the simulated trait; must be a float between 0 and 1

        If not provided, it will be computed from the sum of the squared effect sizes
    prevalence : int, optional
        The prevalence of the disease if the trait should be simulated as case/control;
        must be a float between 0 and 1

        If not provided, a quantitative trait will be simulated, instead
    normalize: bool, optional
        If True, normalize the genotypes before using them to simulate the phenotypes.
        Otherwise, use the raw values.
    region : str, optional
        The region from which to extract haplotypes; ex: 'chr1:1234-34566' or 'chr7'

        For this to work, the VCF and .hap file must be indexed and the seqname must
        match!

        Defaults to loading all haplotypes
    samples : set[str], optional
        A subset of the samples from which to extract genotypes

        Defaults to loading genotypes from all samples
    haplotype_ids: set[str], optional
        A list of haplotype IDs to obtain from the .hap file. All others are ignored.

        If not provided, all haplotypes will be used.
    chunk_size: int, optional
        The max number of variants to fetch from the PGEN file at any given time

        If this value is provided, variants from the PGEN file will be loaded in
        chunks so as to use less memory. This argument is ignored if the genotypes are
        not in PGEN format.
    repeats: Path, optional
        The path to a genotypes file containing tandem repeats. This is only necessary
        when simulating both haplotypes *and* repeats as causal effects
    environment: float, optional
        The variance (aka strength) of the environmental term. This will be inferred if
        it isn't specified.
    seed: int, optional
        Seed for random processes
    output : Path, optional
        The location to which to write the simulated phenotypes
    log : logging.Logger, optional
        The logging module for this task
    """
    if log is None:
        log = getLogger(name="simphenotype", level="ERROR")

    load_as_haps = True
    # either load SNPs from the snplist file or load haps/repeats from the hap file
    if haplotypes.suffix == ".snplist":
        log.info("Loading from .snplist")
        with open(haplotypes) as snplist_file:
            effects = map(Effect.from_hap_spec, snplist_file.readlines())
        if haplotype_ids is None:
            effects = list(effects)
            haplotype_ids = set(effect.id for effect in effects)
        else:
            effects = list(filter(lambda e: e.id in haplotype_ids, effects))
    else:
        log.info("Loading from .hap")
        hp = Haplotypes(haplotypes, haplotype=Haplotype, repeat=Repeat, log=log)
        hp.read(region=region, haplotypes=haplotype_ids)
        effects = hp.data.values()

        if haplotype_ids is None:
            haplotype_ids = set(hp.data.keys())

        # check if these are all repeat IDs, haplotype IDs, or a mix of them
        if len(hp.type_ids["R"]) >= len(haplotype_ids) and repeats is None:
            # if they're all repeat IDs or --repeats was specified
            if genotypes.suffix == ".pgen":
                log.info("Loading TR genotypes from PGEN file")
                gt = GenotypesPLINKTR(fname=genotypes, log=log, chunk_size=chunk_size)
            else:
                log.info("Loading TR genotypes from VCF/BCF file")
                gt = GenotypesTR(fname=genotypes, log=log)
            load_as_haps = False
        else:
            # the genotypes variable must contain haplotype genotypes
            # but first, check if they're a mix but --repeats wasn't specified
            if len(hp.type_ids["H"]) < len(haplotype_ids) and repeats is None:
                raise ValueError(
                    "The --repeats option must be specified when simulating a mix of"
                    " both haplotypes and repeats as causal effects."
                )

    if load_as_haps:
        # load these as haplotype pseudo-genotypes
        if genotypes.suffix == ".pgen":
            log.info("Loading haplotype genotypes from PGEN file")
            gt = GenotypesPLINK(fname=genotypes, log=log, chunk_size=chunk_size)
        else:
            log.info("Loading haplotype genotypes from VCF/BCF file")
            gt = Genotypes(fname=genotypes, log=log)

    # gt._prephased = True
    gt.read(region=region, samples=samples, variants=haplotype_ids)
    log.info("QC-ing genotypes")
    gt.check_missing()

    if repeats:
        log.info("Merging with TR genotypes")
        if repeats.suffix == ".pgen":
            log.info("Loading repeat genotypes from PGEN file")
            tr_gt = GenotypesPLINKTR(fname=repeats, log=log, chunk_size=chunk_size)
        else:
            log.info("Loading repeat genotypes from VCF/BCF file")
            tr_gt = GenotypesTR(fname=repeats, log=log)
        tr_gt.read(region=region, samples=samples, variants=haplotype_ids)
        tr_gt.check_missing()
        gt = Genotypes.merge_variants((gt, tr_gt), fname=None)

    # check that all of the genotypes were loaded successfully and warn otherwise
    if len(haplotype_ids) > len(gt.variants):
        diff = list(haplotype_ids.difference(gt.variants["id"]))
        first_few = 5 if len(diff) > 5 else len(diff)
        log.error(
            f"{len(diff)} effects could not be found in the genotypes file. Check "
            "that the IDs in your .snplist or .hap file correspond with those in the "
            "genotypes file. Here are the first few missing variants: "
            f"{diff[:first_few]}"
        )

    # Initialize phenotype simulator (haptools simphenotype)
    log.info("Simulating phenotypes")
    pt_sim = PhenoSimulator(gt, output=output, seed=seed, log=log)
    for i in range(num_replications):
        pt_sim.run(effects, heritability, prevalence, normalize, environment)
    log.info("Writing phenotypes")
    pt_sim.write()
