#!/usr/bin/env python

import sys
import pickle
import shutil
import subprocess
from pathlib import Path
from time import process_time

import click
import matplotlib
import numpy as np
from pysam import VariantFile
import matplotlib.pyplot as plt

# Force matplotlib to not use any Xwindows backend.
matplotlib.use("Agg")

from haptools import logging
from haptools.data import (
    GenotypesTR,
    GenotypesVCF,
    GenotypesPLINK,
    GenotypesPLINKTR,
)


# COMMAND FOR GENERATING UKB PLOT:
# tests/bench_genotypes.py \
# --default-variants 18472 --default-samples 487409 --intervals-variants 1 80 4 \
# --intervals-samples 1 80 4 -o plot.pdf -a results.pickle
# (allocate 30 mins, 4 CPUs, and 4GB of RAM)

DATADIR = Path(__file__).parent.joinpath("data")


class GenotypesVCFTR(GenotypesVCF):
    def write(self):
        """
        Write the variants in this class to a VCF at :py:attr:`~.GenotypesTR.fname`

        Output the VCF in HipSTR format
        """
        vcf = VariantFile(str(self.fname), mode="w")
        # make sure the header is properly structured
        for contig in set(self.variants["chrom"]):
            vcf.header.contigs.add(contig)
        d = "Inclusive {} coodinate for the repetitive portion of the reference allele"
        vcf.header.add_meta(
            "INFO",
            items=[
                ("ID", "START"),
                ("Number", 1),
                ("Type", "Integer"),
                ("Description", d.format("start")),
            ],
        )
        vcf.header.add_meta(
            "INFO",
            items=[
                ("ID", "END"),
                ("Number", 1),
                ("Type", "Integer"),
                ("Description", d.format("end")),
            ],
        )
        vcf.header.add_meta(
            "INFO",
            items=[
                ("ID", "PERIOD"),
                ("Number", 1),
                ("Type", "Integer"),
                ("Description", "Length of STR motif"),
            ],
        )
        vcf.header.add_meta(
            "FORMAT",
            items=[
                ("ID", "GT"),
                ("Number", 1),
                ("Type", "String"),
                ("Description", "Genotype"),
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
        missing_val = np.iinfo(np.uint8).max
        for var_idx, var in enumerate(self.variants):
            rec = {
                "contig": var["chrom"],
                "start": var["pos"],
                "stop": var["pos"] + len(var["alleles"][0]),
                "qual": None,
                "alleles": var["alleles"],
                "id": var["id"],
                "filter": None,
            }
            # handle pysam increasing the start site by 1
            rec["start"] -= 1
            # parse the record into a pysam.VariantRecord
            record = vcf.new_record(**rec)
            # add INFO flags expected of HipSTR
            # Note: this is only possible because we guarantee that the REF allele
            # is a single copy of a dinucleotide repeat
            record.info["START"] = int(rec["start"])
            record.stop = int(rec["stop"])
            record.info["PERIOD"] = len(var["alleles"][0])
            for samp_idx, sample in enumerate(self.samples):
                record.samples[sample]["GT"] = tuple(
                    None if val == missing_val else val
                    for val in self.data[samp_idx, var_idx, :2]
                )
                # add proper phasing info
                if phased:
                    record.samples[sample].phased = True
                else:
                    record.samples[sample].phased = self.data[samp_idx, var_idx, 2]
            # write the record to a file
            vcf.write(record)
        vcf.close()


def write_tr_files(plink2: Path, vcf: Path, pgen: Path):
    # first, add "HipSTR" to the header so that we can trick TRTools
    subprocess.call(["sed", "-i", "2i##command=HipSTR-v0.7 --test/", vcf])
    subprocess.call(
        [
            plink2,
            "--vcf-half-call",
            "m",
            "--make-pgen",
            "pvar-cols=vcfheader,qual,filter,info",
            "--vcf",
            vcf,
            "--out",
            pgen.with_suffix(""),
        ],
        stderr=subprocess.DEVNULL,
        stdout=subprocess.DEVNULL,
    )


def create_variant_files(gts, intervals, num_samps, plink2: Path = None):
    samples = gts.samples[:num_samps]
    variant_dir = gts.fname / "variant"
    variant_dir.mkdir(parents=True, exist_ok=True)
    for val in intervals:
        variants = tuple(gts.variants["id"][:val])
        sub = gts.subset(samples=samples, variants=variants)
        sub.fname = variant_dir / f"{val}.pgen"
        # write VCF files
        if plink2:
            vcf = GenotypesVCFTR(sub.fname.with_suffix(".vcf"))
        else:
            vcf = GenotypesVCF(sub.fname.with_suffix(".vcf"))
        vcf.data = sub.data
        vcf.samples = sub.samples
        vcf.variants = sub.variants
        vcf.write()
        # write PLINK2 files
        if plink2:
            write_tr_files(plink2, vcf.fname, sub.fname)
        else:
            sub.write()
    return variant_dir


def create_sample_files(gts, intervals, num_vars, plink2: Path = None):
    variants = tuple(gts.variants["id"][:num_vars])
    sample_dir = gts.fname / "sample"
    sample_dir.mkdir(parents=True, exist_ok=True)
    for val in intervals:
        samples = gts.samples[:val]
        sub = gts.subset(samples=samples, variants=variants)
        sub.fname = sample_dir / f"{val}.pgen"
        # write VCF files
        if plink2:
            vcf = GenotypesVCFTR(sub.fname.with_suffix(".vcf"))
        else:
            vcf = GenotypesVCF(sub.fname.with_suffix(".vcf"))
        vcf.data = sub.data
        vcf.samples = sub.samples
        vcf.variants = sub.variants
        vcf.write()
        # write PLINK2 files
        if plink2:
            write_tr_files(plink2, vcf.fname, sub.fname)
        else:
            sub.write()
    return sample_dir


def time_vcf(vcf, max_variants, chunk_size=None):
    GenotypesVCF(vcf).read(max_variants=max_variants)


def time_vcf_tr(vcf, max_variants, chunk_size=None):
    GenotypesTR(vcf, vcftype="hipstr").read(max_variants=max_variants)


def time_plink(pgen, max_variants, chunk_size=None):
    GenotypesPLINK(pgen).read(max_variants=max_variants)


def time_plink_tr(pgen, max_variants, chunk_size=None):
    GenotypesPLINKTR(pgen, vcftype="hipstr").read(max_variants=max_variants)


def time_plink_chunk(pgen, max_variants, chunk_size=None):
    GenotypesPLINK(pgen, chunk_size=chunk_size).read(max_variants=max_variants)


def time_plink_chunk_tr(pgen, max_variants, chunk_size=None):
    GenotypesPLINKTR(pgen, chunk_size=chunk_size, vcftype="hipstr").read(
        max_variants=max_variants
    )


def progressbar(it, prefix="", size=60, out=sys.stdout):  # Python3.6+
    count = len(it)

    def show(j):
        x = int(size * j / count)
        print(
            f"{prefix}[{u'█'*x}{('.'*(size-x))}] {j}/{count}",
            end="\r",
            file=out,
            flush=True,
        )

    show(0)
    for i, item in enumerate(it):
        yield item
        show(i + 1)
    print("\n", flush=True, file=out)


def create_tr_alleles():
    return tuple("AT" * j for j in range(1, np.random.randint(2, 11) + 1))


@click.command()
@click.option(
    "-t",
    "--temp",
    type=click.Path(path_type=Path, file_okay=False, dir_okay=True),
    default=Path("bench_temp_dir"),
    show_default=True,
    help="A temp directory into which to place generated temporary files",
)
@click.option(
    "--default-variants",
    type=int,
    default=4,
    show_default=True,
    help="The number of variants to use when we vary the number of samples",
)
@click.option(
    "--default-samples",
    type=int,
    default=5,
    show_default=True,
    help="The number of samples to use when we vary the number of variants",
)
@click.option(
    "--intervals-variants",
    type=int,
    nargs=3,
    default=(1, 4, 1),
    show_default=True,
    help="The start, end, and step values for the x-axis of the variants plot",
)
@click.option(
    "--intervals-samples",
    type=int,
    nargs=3,
    default=(1, 5, 1),
    show_default=True,
    help="The start, end, and step values for the x-axis of the samples plot",
)
@click.option(
    "--reps",
    type=int,
    default=3,
    show_default=True,
    help="For each benchmark value, we take the mean of --reps X replicates",
)
@click.option(
    "-o",
    "--output",
    type=click.Path(path_type=Path),
    default=Path("bench-plink.png"),
    show_default=True,
    help="A PNG file containing the results of the benchmark",
)
@click.option(
    "-p",
    "--progress",
    is_flag=True,
    default=False,
    show_default=True,
    help="Whether to display a progress bar. Useful for large benchmarks",
)
@click.option(
    "-s",
    "--silent",
    is_flag=True,
    default=False,
    show_default=True,
    help="Whether to be entirely silent",
)
@click.option(
    "-a",
    "--archive",
    type=click.Path(path_type=Path),
    default=None,
    show_default="do not save generated results",
    help="A python pickle file into which to store results",
)
@click.option(
    "--plink2",
    type=click.Path(path_type=Path, exists=True),
    default=None,
    show_default="do not benchmark repeats",
    help="Benchmark repeats in addition to SNPs. Create PGENs with this plink2 binary",
)
def main(
    temp,
    default_variants,
    default_samples,
    intervals_variants,
    intervals_samples,
    reps,
    output,
    progress=False,
    silent=False,
    archive=None,
    plink2=None,
):
    """
    Benchmarks classes in the data.genotypes module
    """
    DEFAULT_VARIANTS, DEFAULT_SAMPLES, INTERVALS_VARIANTS, INTERVALS_SAMPLES, REPS = (
        default_variants,
        default_samples,
        range(*intervals_variants),
        range(*intervals_samples),
        reps,
    )
    LOG = logging.getLogger("run", 0 if silent else "DEBUG")
    gts = GenotypesPLINK(None, log=LOG)
    num_variants = max(DEFAULT_VARIANTS, INTERVALS_VARIANTS.stop)
    sample_size = max(DEFAULT_SAMPLES, INTERVALS_SAMPLES.stop)
    gts.samples = tuple(f"sample{i}" for i in range(sample_size))
    if plink2:
        gts.variants = np.array(
            [
                (f"id{i}", "chr0", i, create_tr_alleles())
                for i in range(1, num_variants + 1)
            ],
            dtype=gts.variants.dtype,
        )
    else:
        gts.variants = np.array(
            [(f"id{i}", "chr0", i, ("A", "T")) for i in range(1, num_variants + 1)],
            dtype=gts.variants.dtype,
        )
    np.random.seed(12345)
    if not temp.exists():
        LOG.info("Generating fake genotypes")
        gts.data = np.empty(
            (sample_size, num_variants, 2), dtype=(np.uint8 if plink2 else np.bool_)
        )
        for i in range(num_variants):
            gts.data[:, i] = np.random.choice(
                range(len(gts.variants[i]["alleles"])), size=sample_size * 2
            ).reshape((sample_size, 2))
    gts.fname = temp
    # set initial variables
    SAMPLES = gts.samples
    VARIANTS = gts.variants["id"]
    DEFAULT_SAMPLES = min(DEFAULT_SAMPLES, len(SAMPLES))
    DEFAULT_VARIANTS = min(DEFAULT_VARIANTS, len(VARIANTS))
    if INTERVALS_VARIANTS.stop > len(VARIANTS):
        INTERVALS_VARIANTS = range(
            INTERVALS_VARIANTS.start, len(VARIANTS), INTERVALS_VARIANTS.step
        )
    INTERVALS_VARIANTS = list(INTERVALS_VARIANTS)
    if INTERVALS_SAMPLES.stop > len(SAMPLES):
        INTERVALS_SAMPLES = range(
            INTERVALS_SAMPLES.start, len(SAMPLES), INTERVALS_SAMPLES.step
        )
    INTERVALS_SAMPLES = list(INTERVALS_SAMPLES)
    FILE_TYPES = {
        "vcf": "VCF",
        "pgen": "PLINK2",
        "chunked3000": "PLINK2 chunk_size: 3000",
        "chunked5000": "PLINK2 chunk_size: 5000",
        "chunked7000": "PLINK2 chunk_size: 7000",
        "chunked9000": "PLINK2 chunk_size: 9000",
    }

    # create the files we will try to load if they haven't been created already
    if not gts.fname.exists():
        LOG.info("Creating VCF and PGEN files that we can load")
        variant_dir = create_variant_files(
            gts, INTERVALS_VARIANTS, DEFAULT_SAMPLES, plink2
        )
        sample_dir = create_sample_files(
            gts, INTERVALS_SAMPLES, DEFAULT_VARIANTS, plink2
        )
    else:
        variant_dir = gts.fname / "variant"
        sample_dir = gts.fname / "sample"
        LOG.info(
            "Temp directory already exists. Assuming files have already been created",
        )

    LOG.setLevel(level="ERROR")

    # run each test
    LOG.info("Benchmarking the loading of each file")
    results = {}
    len_variants = len(VARIANTS)
    for arg in ("samples", "variants"):
        genotype_dir = variant_dir
        if arg == "samples":
            genotype_dir = sample_dir
        results[arg] = {}
        intervals = INTERVALS_SAMPLES if arg == "samples" else INTERVALS_VARIANTS
        for file_type in FILE_TYPES.keys():
            item_iter = (
                progressbar(intervals, prefix=f"{arg}, {file_type}: ", out=sys.stderr)
                if progress
                else intervals
            )
            results[arg][file_type] = []
            for val in item_iter:
                chunk_size = None
                if file_type == "vcf":
                    func = time_vcf_tr if plink2 else time_vcf
                    file = genotype_dir / f"{val}.vcf"
                elif file_type == "pgen" or file_type.startswith("chunked"):
                    if file_type == "pgen":
                        func = time_plink_tr if plink2 else time_plink
                    else:
                        func = time_plink_chunk_tr if plink2 else time_plink_chunk
                    file = genotype_dir / f"{val}.pgen"
                    if file_type.startswith("chunked"):
                        chunk_size = int(file_type[len("chunked") :])
                        if arg == "variants" and chunk_size > val:
                            continue
                else:
                    continue
                times = np.empty(REPS, dtype=np.float64)
                for rep in range(REPS):
                    start = process_time()
                    func(file, max_variants=len_variants, chunk_size=chunk_size)
                    end = process_time()
                    times[rep] = end - start
                results[arg][file_type].append(times.mean())

    # plot the results
    LOG.info("Generating plot of results")
    fig, (ax_samples, ax_variants) = plt.subplots(1, 2, figsize=(10, 5))
    lab_font = {"fontsize": "xx-small"}
    for file_type in FILE_TYPES.keys():
        x_vals = INTERVALS_SAMPLES
        y_vals = results["samples"][file_type]
        # fit a line to each so that we can report the slope
        slope = np.polyfit(x_vals, y_vals, 1)[0]
        ax_samples.plot(
            x_vals,
            y_vals,
            marker="o",
            label=FILE_TYPES[file_type],
        )
        ax_samples.text(
            x_vals[-1],
            y_vals[-1] + (y_vals[-1] / 16),
            f"m = {slope:.3E}",
            fontdict=lab_font,
        )
    ax_samples.set_xlabel(f"Number of samples\nnum_variants = {DEFAULT_VARIANTS}")
    ax_samples.set_ylim(ymin=0)
    for file_type in FILE_TYPES.keys():
        x_vals = INTERVALS_VARIANTS
        y_vals = results["variants"][file_type]
        if not len(y_vals):
            continue
        # fit a line to each so that we can report the slope
        slope = np.polyfit(x_vals, y_vals, 1)[0]
        ax_variants.plot(
            x_vals,
            y_vals,
            marker="o",
        )
        ax_variants.text(
            x_vals[-1],
            y_vals[-1] + (y_vals[-1] / 16),
            f"m = {slope:.3E}",
            fontdict=lab_font,
        )
    ax_variants.set_xlabel(f"Number of variants\nnum_samples = {DEFAULT_SAMPLES}")
    ax_variants.set_ylim(ymin=0)
    fig.supylabel("CPU Time (s)")
    fig.legend(loc="lower left", fontsize="xx-small")
    fig.tight_layout()
    fig.savefig(str(output), bbox_inches="tight", dpi=400)

    # saving results
    if archive is not None:
        with open(archive, "wb") as fh:
            pickle.dump(results, fh)


def test_bench_genotypes():
    tmp_dir = Path("bench_temp_dir")
    tmp_plot = Path("bench-plink.png")

    try:
        main(["-t", tmp_dir, "-o", tmp_plot, "--reps", 1, "-s"])
    except SystemExit as err:
        # re-raise unless main() finished without an error
        if err.code:
            raise

    shutil.rmtree(tmp_dir)
    tmp_plot.unlink()


if __name__ == "__main__":
    main()
