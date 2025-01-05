#!/usr/bin/env python

from __future__ import annotations
import sys
from pathlib import Path

import click

# AVOID IMPORTING ANYTHING ABOVE
# any imports we put here will make it slower to use the command line client
# a basic "haptools --help" should be quick and require very few imports, for example


################### Haptools ##################
@click.group()
@click.version_option(message="%(version)s")
def main():
    """
    haptools: A toolkit for simulating and analyzing genotypes and
    phenotypes while taking into account haplotype information
    """
    pass


@main.command()
@click.option(
    "--bp",
    type=str,
    required=True,
    help="Path to .bp file with breakpoints",
)
@click.option(
    "--sample",
    type=str,
    required=True,
    help="Sample ID to plot",
)
@click.option(
    "--out",
    type=str,
    required=True,
    help="Name of output file",
)
@click.option(
    "--title",
    type=str,
    required=False,
    help="Optional plot title",
)
@click.option(
    "--centromeres",
    type=str,
    required=False,
    help="Optional file with telomere/centromere cM positions",
)
@click.option(
    "--colors",
    type=str,
    required=False,
    help=(
        "Optional color dictionary. Input can be from the matplotlib list of colors "
        "or in hexcode. Format is e.g. 'YRI:blue,CEU:green'"
    ),
)
@click.option(
    "-v",
    "--verbosity",
    type=click.Choice(["CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG", "NOTSET"]),
    default="INFO",
    show_default=True,
    help="The level of verbosity desired",
)
def karyogram(bp, sample, out, title, centromeres, colors, verbosity):
    """
    Visualize a karyogram of local ancestry tracks
    """
    from .karyogram import PlotKaryogram
    from .logging import getLogger

    log = getLogger(name="karyogram", level=verbosity)

    if colors is not None:
        colors = dict([item.split(":") for item in colors.split(",")])
    log.info("Generating Karyogram...")
    PlotKaryogram(
        bp, sample, out, log, centromeres_file=centromeres, title=title, colors=colors
    )


@main.command()
@click.option(
    "--model",
    type=str,
    required=True,
    help=(
        "Admixture model in .dat format. See File Formats under simgenotype in the "
        "docs for complete info."
    ),
)
@click.option(
    "--mapdir",
    required=True,
    type=click.Path(exists=True, file_okay=False, dir_okay=True, readable=True),
    help=(
        "Directory containing files with chr{1-22,X} and ending in .map in the file "
        "name with genetic map coords."
    ),
)
@click.option(
    "--out",
    type=str,
    required=True,
    help=(
        "Path to desired output file. E.g. /path/to/output.vcf.gz "
        "Possible outputs are vcf|bcf|vcf.gz|pgen and there will be an "
        "additional breakpoints output with extension bp e.g. /path/to/output.bp."
    ),
)
@click.option(
    "--chroms",
    type=str,
    required=False,
    default="1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X",
    help="Sorted and comma delimited list of chromosomes to simulate",
)
@click.option(
    "--seed",
    type=int,
    required=False,
    default=None,
    help="Random seed. Set to make simulations reproducible",
)
@click.option(
    "--popsize",
    type=int,
    hidden=True,
    default=10000,
    required=False,
    help="Number of samples to simulate each generation",
)
@click.option(
    "--ref_vcf",
    required=True,
    help=(
        "VCF or PGEN file used as reference for creation of simulated samples"
        " respective genotypes."
    ),
)
@click.option(
    "--sample_info",
    required=True,
    help=(
        "File that maps samples from the reference VCF (--invcf) to population"
        " codes describing the populations in the header of the model file."
    ),
)
@click.option(
    "--region",
    required=False,
    default=None,
    help=(
        "Subset the simulation to a specific region in a chromosome using the"
        " form chrom:start-end. Example 2:1000-2000"
    ),
)
@click.option(
    "--pop_field",
    required=False,
    is_flag=True,
    default=False,
    help=(
        "Flag for outputting the population field in your VCF output. NOTE this"
        " flag does not work when your output file is in PGEN format."
    ),
)
@click.option(
    "--sample_field",
    required=False,
    is_flag=True,
    default=False,
    help=(
        "Flag for outputting the sample field in your VCF output. NOTE this"
        " flag does not work when your output file is in PGEN format."
    ),
)
@click.option(
    "--no_replacement",
    is_flag=True,
    required=False,
    default=False,
    help=(
        "Flag used to determine whether to sample reference haplotypes"
        " with or without replacement. (Default = Replacement)"
    ),
)
@click.option(
    "--only_breakpoint",
    is_flag=True,
    required=False,
    help=(
        "Flag used to determine whether to only output breakpoints or"
        " continue to simulate a vcf file."
    ),
)
@click.option(
    "-c",
    "--chunk-size",
    type=int,
    default=None,
    show_default="all variants",
    help=(
        "If requesting a PGEN output file, write genotypes in chunks of X variants; "
        "reduces memory"
    ),
)
@click.option(
    "-v",
    "--verbosity",
    type=click.Choice(["CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG", "NOTSET"]),
    default="INFO",
    show_default=True,
    help="The level of verbosity desired",
)
def simgenotype(
    ref_vcf,
    sample_info,
    model,
    mapdir,
    out,
    popsize,
    seed,
    chroms,
    region,
    pop_field,
    sample_field,
    no_replacement,
    only_breakpoint,
    chunk_size,
    verbosity,
):
    """
    Simulate admixed genomes under a pre-defined model.
    """
    import re
    import time
    from .sim_genotype import (
        output_vcf,
        simulate_gt,
        validate_params,
        write_breakpoints,
    )
    from .logging import getLogger

    log = getLogger(name="simgenotype", level=verbosity)
    start = time.time()

    # immediately set pop_filed and sample_field flags to false if pgen file
    if out.endswith(".pgen"):
        pop_field = False
        sample_field = False

    # parse region and chroms parameters
    if not (chroms or region):
        raise Exception("Either chroms or region must be specified.")
    if region:
        region_info = re.split(":|-", region)
        try:
            region = {
                "chr": region_info[0],
                "start": int(region_info[1]),
                "end": int(region_info[2]),
            }
            chroms = [region["chr"]]
        except:
            raise Exception(
                "Unable to parse region. Please ensure it has the correct format "
                "<chr>:<start>-<end> eg. 1:1-2000"
            )
    else:
        chroms = chroms.split(",")

    # Handle if mapdir has a '/' at the end
    if mapdir[-1] == "/":
        mapdir = mapdir[:-1]

    # grab prefix from --out for outputting breakpoint
    out_prefix = re.split(r"(\.vcf|\.bcf|\.vcf\.gz|\.pgen)$", out)[0]

    # simulate breakpoints
    popsize = validate_params(
        model,
        mapdir,
        chroms,
        popsize,
        ref_vcf,
        sample_info,
        no_replacement,
        region,
        only_breakpoint,
    )
    samples, pop_dict, breakpoints = simulate_gt(
        model, mapdir, chroms, region, popsize, log, seed
    )
    breakpoints = write_breakpoints(samples, pop_dict, breakpoints, out_prefix, log)
    bp_end = time.time()

    # simulate vcfs
    vcf_start = time.time()
    if not only_breakpoint:
        output_vcf(
            breakpoints,
            chroms,
            model,
            ref_vcf,
            sample_info,
            region,
            pop_field,
            sample_field,
            no_replacement,
            out,
            log,
            chunk_size,
        )
    end = time.time()

    log.debug(f"Time elapsed for breakpoint simulation: {bp_end - start}")
    log.debug(f"Time elapsed for creating vcf: {end - vcf_start}")
    log.debug(f"Total time elapsed for simgenotype execution: {end - start}")


@main.command()
@click.argument("genotypes", type=click.Path(exists=True, path_type=Path))
@click.argument("haplotypes", type=click.Path(exists=True, path_type=Path))
@click.option(
    "-r",
    "--replications",
    type=click.IntRange(min=1),
    default=1,
    show_default=True,
    help="Number of rounds of simulation to perform",
)
@click.option(
    "--environment",
    type=click.FloatRange(min=0),
    default=None,
    show_default=True,
    help="Variance of environmental term; inferred if not specified",
)
@click.option(
    "-h",
    "--heritability",
    type=click.FloatRange(min=0, max=1),
    default=None,
    show_default="0.5",
    help="Trait heritability",
)
@click.option(
    "-p",
    "--prevalence",
    type=click.FloatRange(min=0, max=1, min_open=False, max_open=True),
    show_default="quantitative trait",
    help="Disease prevalence if simulating a case-control trait",
)
@click.option(
    "--normalize/--no-normalize",
    is_flag=True,
    default=True,
    show_default=True,
    help="Whether to normalize the genotypes before using them for simulation",
)
@click.option(
    "--region",
    type=str,
    default=None,
    show_default="all haplotypes",
    help=(
        "The region from which to extract haplotypes; ex: 'chr1:1234-34566' or 'chr7'."
        "\nFor this to work, the VCF and .hap file must be indexed and the seqname "
        "provided must correspond with one in the files"
    ),
)
@click.option(
    "-s",
    "--sample",
    "samples",
    type=str,
    multiple=True,
    show_default="all samples",
    help=(
        "A list of the samples to subset from the genotypes file (ex: '-s sample1 -s"
        " sample2')"
    ),
)
@click.option(
    "-S",
    "--samples-file",
    type=click.File("r"),
    show_default="all samples",
    help=(
        "A single column txt file containing a list of the samples (one per line) to"
        " subset from the genotypes file"
    ),
)
@click.option(
    "-i",
    "--id",
    "ids",
    type=str,
    multiple=True,
    show_default="all haplotypes",
    help=(
        "A list of the haplotype IDs from the .hap file to use as causal variables "
        "(ex: '-i H1 -i H2')."
    ),
)
@click.option(
    "-I",
    "--ids-file",
    type=str,
    multiple=True,
    show_default="all haplotypes",
    help=(
        "A single column txt file containing a list of the haplotype IDs "
        "(one per line) to subset from the .hap file"
    ),
)
@click.option(
    "-c",
    "--chunk-size",
    type=int,
    default=None,
    show_default="all variants",
    help="If using a PGEN file, read genotypes in chunks of X variants; reduces memory",
)
@click.option(
    "--repeats",
    type=click.Path(exists=True, path_type=Path),
    default=None,
    help=(
        "Path to a genotypes file containing tandem repeats. This is only necessary "
        "when simulating both haplotypes *and* repeats as causal effects"
    ),
)
@click.option(
    "--seed",
    type=int,
    default=None,
    show_default="chosen randomly",
    help="Use this option across executions to make the output reproducible",
)
@click.option(
    "-o",
    "--output",
    type=click.Path(path_type=Path),
    default=Path("/dev/stdout"),
    show_default="stdout",
    help="A TSV file containing simulated phenotypes",
)
@click.option(
    "-v",
    "--verbosity",
    type=click.Choice(["CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG", "NOTSET"]),
    default="INFO",
    show_default=True,
    help="The level of verbosity desired",
)
def simphenotype(
    genotypes: Path,
    haplotypes: Path,
    replications: int = 1,
    environment: float = None,
    heritability: float = None,
    prevalence: float = None,
    normalize: bool = True,
    region: str = None,
    samples: tuple[str] = tuple(),
    samples_file: Path = None,
    ids: tuple[str] = tuple(),
    ids_file: Path = None,
    chunk_size: int = None,
    repeats: Path = None,
    seed: int = None,
    output: Path = Path("-"),
    verbosity: str = "INFO",
):
    """
    Haplotype-aware phenotype simulation. Create a set of simulated phenotypes from a
    set of haplotypes.

    GENOTYPES must be formatted as a VCF or PGEN file and HAPLOTYPES must be formatted
    according to the .hap format spec

    Note: GENOTYPES must be the output from the transform subcommand.
    """
    from .logging import getLogger
    from .sim_phenotype import simulate_pt

    log = getLogger(name="simphenotype", level=verbosity)

    # handle samples
    if samples and samples_file:
        raise click.UsageError(
            "You may only use one of --sample or --samples-file but not both."
        )
    if samples_file:
        with samples_file as samps_file:
            samples = set(samps_file.read().splitlines())
    elif samples:
        # needs to be converted from tuple to set
        samples = set(samples)
    else:
        samples = None

    if ids_file:
        with ids_file as id_file:
            ids = set(id_file.read().splitlines())
    elif ids:
        ids = set(ids)
    else:
        ids = None

    if heritability is None and environment is None and not normalize:
        # in this case, the assumptions of the model break
        # The variances of both sides of the equation will no longer properly sum to 1
        log.error("A --heritability value should be specified with --no-normalize")

    # Run simulation
    simulate_pt(
        genotypes,
        haplotypes,
        replications,
        environment,
        heritability,
        prevalence,
        normalize,
        region,
        samples,
        ids,
        chunk_size,
        repeats,
        seed,
        output,
        log,
    )


@main.command(short_help="Transform a genotypes matrix via a set of haplotypes")
@click.argument("genotypes", type=click.Path(exists=True, path_type=Path))
@click.argument("haplotypes", type=click.Path(exists=True, path_type=Path))
@click.option(
    "--region",
    type=str,
    default=None,
    show_default="all haplotypes",
    help=(
        "The region from which to extract haplotypes; ex: 'chr1:1234-34566' or 'chr7'."
        "\nFor this to work, the VCF and .hap file must be indexed and the seqname "
        "provided must correspond with one in the files"
    ),
)
@click.option(
    "-s",
    "--sample",
    "samples",
    type=str,
    multiple=True,
    show_default="all samples",
    help=(
        "A list of the samples to subset from the genotypes file (ex: '-s sample1 -s"
        " sample2')"
    ),
)
@click.option(
    "-S",
    "--samples-file",
    type=click.File("r"),
    show_default="all samples",
    help=(
        "A single column txt file containing a list of the samples (one per line) to"
        " subset from the genotypes file"
    ),
)
@click.option(
    "-i",
    "--id",
    "ids",
    type=str,
    multiple=True,
    show_default="all haplotypes",
    help="A list of the haplotype IDs to use from the .hap file (ex: '-i H1 -i H2').",
)
@click.option(
    "-I",
    "--ids-file",
    type=str,
    multiple=True,
    show_default="all haplotypes",
    help=(
        "A single column txt file containing a list of the haplotype IDs "
        "(one per line) to subset from the .hap file"
    ),
)
@click.option(
    "-c",
    "--chunk-size",
    type=int,
    default=None,
    show_default="all variants",
    help="If using a PGEN file, read genotypes in chunks of X variants; reduces memory",
)
@click.option(
    "--discard-missing",
    is_flag=True,
    show_default=True,
    default=False,
    help="Ignore any samples that are missing genotypes for the required variants",
)
@click.option(
    "--ancestry",
    is_flag=True,
    show_default=True,
    default=False,
    help="Also transform using VCF 'POP' FORMAT field and 'ancestry' .hap extra field",
)
@click.option(
    "-o",
    "--output",
    type=click.Path(path_type=Path),
    default=Path("-"),
    show_default="stdout",
    help="A VCF file containing haplotype 'genotypes'",
)
@click.option(
    "-v",
    "--verbosity",
    type=click.Choice(["CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG", "NOTSET"]),
    default="INFO",
    show_default=True,
    help="The level of verbosity desired",
)
def transform(
    genotypes: Path,
    haplotypes: Path,
    region: str = None,
    samples: tuple[str] = tuple(),
    samples_file: Path = None,
    ids: tuple[str] = tuple(),
    ids_file: Path = None,
    chunk_size: int = None,
    discard_missing: bool = False,
    ancestry: bool = False,
    output: Path = Path("-"),
    verbosity: str = "INFO",
):
    """
    Creates a VCF composed of haplotypes

    GENOTYPES must be formatted as a VCF or PGEN and HAPLOTYPES must be formatted
    according to the .hap format spec
    """
    from .logging import getLogger
    from .transform import transform_haps

    log = getLogger(name="transform", level=verbosity)

    # handle samples
    if samples and samples_file:
        raise click.UsageError(
            "You may only use one of --sample or --samples-file but not both."
        )
    if samples_file:
        with samples_file as samps_file:
            samples = set(samps_file.read().splitlines())
    elif samples:
        # needs to be converted from tuple to set
        samples = set(samples)
    else:
        samples = None

    if ids_file:
        with ids_file as id_file:
            ids = set(id_file.read().splitlines())
    elif ids:
        ids = set(ids)
    else:
        ids = None

    transform_haps(
        genotypes,
        haplotypes,
        region,
        samples,
        ids,
        chunk_size,
        discard_missing,
        ancestry,
        output,
        log,
    )


@main.command(short_help="Compute pair-wise LD")
@click.argument("target", type=str)
@click.argument("genotypes", type=click.Path(exists=True, path_type=Path))
@click.argument("haplotypes", type=click.Path(exists=True, path_type=Path))
@click.option(
    "--region",
    type=str,
    default=None,
    show_default="all haplotypes",
    help=(
        "The region from which to extract haplotypes; ex: 'chr1:1234-34566' or 'chr7'."
        "\nFor this to work, the VCF and .hap file must be indexed and the seqname "
        "provided must correspond with one in the files"
    ),
)
@click.option(
    "-s",
    "--sample",
    "samples",
    type=str,
    multiple=True,
    show_default="all samples",
    help=(
        "A list of the samples to subset from the genotypes file (ex: '-s sample1 -s"
        " sample2')"
    ),
)
@click.option(
    "-S",
    "--samples-file",
    type=click.File("r"),
    show_default="all samples",
    help=(
        "A single column txt file containing a list of the samples (one per line) to"
        " subset from the genotypes file"
    ),
)
@click.option(
    "-i",
    "--id",
    "ids",
    type=str,
    multiple=True,
    show_default="all haplotypes",
    help=(
        "A list of the haplotype IDs to use from the .hap file (ex: '-i H1 -i H2'). "
        "Or, if --from-gts, a list of the variant IDs to use from the genotypes file."
        "\nFor this to work, the .hap file must be indexed"
    ),
)
@click.option(
    "-I",
    "--ids-file",
    type=str,
    multiple=True,
    show_default="all haplotypes",
    help=(
        "A single column txt file containing a list of the haplotype (or variant) IDs "
        "(one per line) to subset from the .hap (or genotype) file"
    ),
)
@click.option(
    "-c",
    "--chunk-size",
    type=int,
    default=None,
    show_default="all variants",
    help="If using a PGEN file, read genotypes in chunks of X variants; reduces memory",
)
@click.option(
    "--discard-missing",
    is_flag=True,
    show_default=True,
    default=False,
    help="Ignore any samples that are missing genotypes for the required variants",
)
@click.option(
    "--from-gts",
    is_flag=True,
    show_default=True,
    default=False,
    help=(
        "By default, LD is computed with the haplotypes in the .hap file. Use this "
        "switch to compute LD with the genotypes in the genotypes file, instead."
    ),
)
@click.option(
    "-o",
    "--output",
    type=click.Path(path_type=Path),
    default=Path("/dev/stdout"),
    show_default="stdout",
    help="A .hap file containing haplotypes and their LD with TARGET",
)
@click.option(
    "-v",
    "--verbosity",
    type=click.Choice(["CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG", "NOTSET"]),
    default="INFO",
    show_default=True,
    help="The level of verbosity desired",
)
def ld(
    target: str,
    genotypes: Path,
    haplotypes: Path,
    region: str = None,
    samples: tuple[str] = tuple(),
    samples_file: Path = None,
    ids: tuple[str] = tuple(),
    ids_file: Path = None,
    chunk_size: int = None,
    discard_missing: bool = False,
    from_gts: bool = False,
    output: Path = Path("/dev/stdout"),
    verbosity: str = "INFO",
):
    """
    Compute the pair-wise LD (Pearson's correlation) between haplotypes (or variants)
    and a single TARGET haplotype (or variant)

    GENOTYPES must be formatted as a VCF or PGEN and HAPLOTYPES must be formatted
    according to the .hap format spec

    TARGET refers to the ID of a variant or haplotype. LD is computed pair-wise between
    TARGET and all of the other haplotypes in the .hap (or genotype) file

    If TARGET is a variant ID, the ID must appear in GENOTYPES. Otherwise, it must
    be present in the .hap file
    """
    from .ld import calc_ld
    from .logging import getLogger

    log = getLogger(name="ld", level=verbosity)

    # handle samples
    if samples and samples_file:
        raise click.UsageError(
            "You may only use one of --sample or --samples-file but not both."
        )
    if samples_file:
        with samples_file as samps_file:
            samples = set(samps_file.read().splitlines())
    elif samples:
        # needs to be converted from tuple to set
        samples = set(samples)
    else:
        samples = None

    if ids_file:
        with ids_file as id_file:
            ids = tuple(id_file.read().splitlines())
    elif ids:
        ids = tuple(ids)
    else:
        ids = None

    calc_ld(
        target,
        genotypes,
        haplotypes,
        region,
        samples,
        ids,
        chunk_size,
        discard_missing,
        from_gts,
        output,
        log,
    )


@main.command(short_help="Sort and index .hap files")
@click.argument("haplotypes", type=click.Path(exists=True, path_type=Path))
@click.option(
    "--sort/--no-sort",
    is_flag=True,
    default=True,
    show_default=True,
    help="Sorting of the file will not be performed",
)
@click.option(
    "-o",
    "--output",
    type=click.Path(path_type=Path),
    default=None,
    show_default="input file",
    help="A .hap file containing sorted and indexed haplotypes and variants",
)
@click.option(
    "-v",
    "--verbosity",
    type=click.Choice(["CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG", "NOTSET"]),
    default="INFO",
    show_default=True,
    help="The level of verbosity desired",
)
def index(
    haplotypes: Path,
    sort: bool = True,
    output: Path = None,
    verbosity: str = "INFO",
):
    """
    Takes in an unsorted .hap file and outputs it as a .gz and a .tbi file
    """

    from .index import index_haps
    from .logging import getLogger

    log = getLogger(name="index", level=verbosity)

    index_haps(haplotypes, sort, output, log)


@main.command(short_help="Clump summary stat files.")
@click.option(
    "--summstats-snps",
    type=click.Path(path_type=Path),
    help="File to load snps summary statistics",
)
@click.option(
    "--summstats-strs",
    type=click.Path(path_type=Path),
    help="File to load strs summary statistics",
)
@click.option(
    "--gts-snps", type=click.Path(path_type=Path), help="SNP genotypes (VCF or PGEN)"
)
@click.option("--gts-strs", type=click.Path(path_type=Path), help="STR genotypes (VCF)")
@click.option(
    "--clump-p1", type=float, default=0.0001, help="Max pval to start a new clump"
)
@click.option(
    "--clump-p2", type=float, default=0.01, help="Filter for pvalue less than"
)
@click.option(
    "--clump-id-field", type=str, default="SNP", help="Column header of the variant ID"
)
@click.option(
    "--clump-field", type=str, default="P", help="Column header of the p-values"
)
@click.option(
    "--clump-chrom-field",
    type=str,
    default="CHR",
    help="Column header of the chromosome",
)
@click.option(
    "--clump-pos-field", type=str, default="POS", help="Column header of the position"
)
@click.option(
    "--clump-kb",
    type=float,
    default=250,
    help="clump kb radius",
)
@click.option(
    "--clump-r2",
    type=float,
    default=0.5,
    help="r^2 threshold",
)
@click.option(
    "--ld",
    type=click.Choice(["Exact", "Pearson"]),
    default="Pearson",
    show_default=True,
    help=(
        "Calculation type to infer LD, Exact Solution or "
        "Pearson R. (Exact|Pearson). Note the Exact Solution "
        "works best when all three genotypes are present (0,1,2) in "
        "the variants being compared."
    ),
)
@click.option(
    "--out",
    type=click.Path(path_type=Path),
    required=True,
    help="Output filename",
)
@click.option(
    "-v",
    "--verbosity",
    type=click.Choice(["CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG", "NOTSET"]),
    default="INFO",
    show_default=True,
    help="The level of verbosity desired",
)
def clump(
    summstats_snps: Path,
    summstats_strs: Path,
    gts_snps: Path,
    gts_strs: Path,
    clump_p1: float,
    clump_p2: float,
    clump_id_field: str,
    clump_field: str,
    clump_chrom_field: str,
    clump_pos_field: str,
    clump_kb: float,
    clump_r2: float,
    ld: str,
    out: Path,
    verbosity: str = "INFO",
):
    """
    Performs clumping on datasets with SNPs, SNPs and STRs, and STRs.
    Clumping is the process of identifying SNPs or STRs that are highly
    correlated with one another and concatenating them all together into
    a single "clump" in order to not repeat the same effect size due to
    LD.
    """
    from .logging import getLogger
    from .clump import clumpstr

    log = getLogger(name="clump", level=verbosity)
    log.debug(f"Loading SNPs from {summstats_snps} {gts_snps}")
    log.debug(f"Loading STRs from {summstats_strs} {gts_strs}")

    clumpstr(
        summstats_snps,
        summstats_strs,
        gts_snps,
        gts_strs,
        clump_p1,
        clump_p2,
        clump_id_field,
        clump_field,
        clump_chrom_field,
        clump_pos_field,
        clump_kb,
        clump_r2,
        ld,
        out,
        log,
    )


@main.command(short_help="Validate the structure of a .hap file")
@click.argument("filename", type=click.Path(exists=True, path_type=Path))
@click.option(
    "--sorted/--not-sorted",
    is_flag=True,
    default=False,
    show_default=True,
    help="Has the file been sorted already?",
)
@click.option(
    "--genotypes",
    type=click.Path(path_type=Path),
    default=None,
    show_default="optional .pvar file to compare against",
    help=(
        "A .pvar file containing variant IDs in order to compare them to the .hap file"
    ),
)
@click.option(
    "-v",
    "--verbosity",
    type=click.Choice(["CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG", "NOTSET"]),
    default="INFO",
    show_default=True,
    help="The level of verbosity desired",
)
def validate(
    filename: Path,
    sorted: bool = False,
    genotypes: Path | None = None,
    verbosity: str = "INFO",
):
    """
    Validate the formatting of a .hap file

    Output warnings/errors explaining how the formatting of your .hap file may
    be improved.
    """
    from .logging import getLogger
    from .validate import is_hapfile_valid

    log = getLogger(name="validate", level=verbosity)

    # if the hap file is compressed and a .tbi index exists for it, assume it is sorted
    if filename.suffix == ".gz" and filename.with_suffix(".gz.tbi").exists():
        sorted = True

    is_valid = is_hapfile_valid(filename, sorted=sorted, log=log, pvar=genotypes)

    if not is_valid:
        raise click.ClickException(
            "Found several warnings and / or errors in the .hap file"
        )


if __name__ == "__main__":
    # run the CLI if someone tries 'python -m haptools' on the command line
    main(prog_name="haptools")
