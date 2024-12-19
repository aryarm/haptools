#!/usr/bin/env python

import sys
import pickle
from pathlib import Path
from time import process_time
from datetime import datetime

import click
import matplotlib
import numpy as np
import matplotlib.pyplot as plt

# Force matplotlib to not use any Xwindows backend.
matplotlib.use("Agg")

from haptools.logging import getLogger
from haptools.data import GenotypesVCF, Haplotypes, Haplotype, Variant


# ---------------USAGE----------------
# COMMAND FOR GENERATING THE MAIN PLOT:
# tests/bench_transform.py --name 'hap-loop' --reps 5 --archive archive.pickle \
# --default-variants 115 --default-samples 5000 --default-haplotypes 20 \
# --intervals-variants 1 80 4 --intervals-samples 1 3200 160 --intervals-haplotypes 10 91 4 -o plot.png


def create_genotypes(log, samples, variants, with_phase=False):
    gts = GenotypesVCF(fname=None, log=log)
    shape = (samples, variants, 2 + with_phase)
    # create a GT matrix with shape: samples x SNPs x (strands+phase)
    gts.data = np.random.choice([0, 1], size=np.prod(shape))
    gts.data = gts.data.reshape(shape).astype(np.uint8)
    gts.variants = np.array(
        [(f"SNP{i}", "chr1", i, 0, "T", "C") for i in range(gts.data.shape[1])],
        dtype=[
            ("id", "U50"),
            ("chrom", "U10"),
            ("pos", np.uint32),
            ("aaf", np.float64),
            ("ref", "U100"),
            ("alt", "U100"),
        ],
    )
    gts.samples = tuple(f"sample{i}" for i in range(gts.data.shape[0]))
    return gts


def create_haplotypes(gts, log, haplotypes, variants):
    haps = Haplotypes(fname=None, log=log)
    haps.data = {}
    for i in range(haplotypes):
        hap_id = f"H{i}"
        haps.data[hap_id] = Haplotype(chrom="chr1", start=0, end=0, id=f"H{i}")
        num_alleles = np.random.randint(low=1, high=(variants + 1))
        hap_variants = np.random.choice(gts.variants, size=num_alleles, replace=False)
        haps.data[hap_id].variants = tuple(
            Variant(0, 0, variant["id"], np.random.choice(["T", "C"]))
            for variant in hap_variants
        )
    return haps


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


@click.command()
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
    "--default-haplotypes",
    type=int,
    default=5,
    show_default=True,
    help="The number of haps to use when we vary the number of variants and samples",
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
    "--intervals-haplotypes",
    type=int,
    nargs=3,
    default=(1, 5, 1),
    show_default=True,
    help="The start, end, and step values for the x-axis of the haplotypes plot",
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
    "-a",
    "--archive",
    type=click.Path(path_type=Path),
    default=None,
    show_default="do not load generated results",
    help="A python pickle file into which to store results",
)
@click.option(
    "-n",
    "--name",
    type=str,
    default=None,
    show_default="do not save generated results",
    help="A unique name for the result when it gets saved in the pickle file",
)
@click.option(
    "--skip-bench",
    is_flag=True,
    default=False,
    show_default=True,
    help="Load directly from the archive instead of benchmarking",
)
def main(
    default_variants,
    default_samples,
    default_haplotypes,
    intervals_variants,
    intervals_samples,
    intervals_haplotypes,
    reps,
    output,
    archive,
    name,
    skip_bench,
):
    """
    Benchmarks transform() in the data.haplotypes module
    """
    INTERVALS_VALS = {
        "vars": range(*intervals_variants),
        "samps": range(*intervals_samples),
        "haps": range(*intervals_haplotypes),
    }
    REPS = reps
    DEFAULT_VAL = {
        "vars": default_variants,
        "samps": default_samples,
        "haps": default_haplotypes,
    }
    NAME = name
    VARIABLES = {"samples": "samps", "alleles_max": "vars", "haplotypes": "haps"}
    LOG = getLogger("run", "ERROR")
    ALG_NAME = {"hap-loop": "loop over haplotypes", "allele-loop": "loop over alleles"}

    # loading results
    results = {}
    if archive is not None:
        if archive.exists():
            with open(archive, "rb") as fh:
                for name, val in pickle.load(fh).items():
                    results[name] = val
                print(f"Loaded items from pickle: {tuple(results.keys())}")

    # run each test
    if not skip_bench:
        print("Benchmarking the loading of each file", file=sys.stderr)
        results[NAME] = {}
        for arg, short in VARIABLES.items():
            results[NAME][arg] = []
            vals = DEFAULT_VAL.copy()
            intervals = INTERVALS_VALS[short]
            for val in progressbar(
                intervals, prefix=f"{arg}, {NAME}: ", out=sys.stderr
            ):
                vals[short] = val
                gts = create_genotypes(LOG, vals["samps"], vals["vars"])
                hps = create_haplotypes(gts, LOG, vals["haps"], vals["vars"])
                times = np.empty(REPS, dtype=np.float64)
                for rep in range(REPS):
                    start = process_time()
                    hps.transform(gts)
                    end = process_time()
                    times[rep] = end - start
                results[NAME][arg].append(times.mean())

    # plot the results
    print("Generating plot of results", file=sys.stderr)
    fig, axs = plt.subplots(1, 3, figsize=(17, 5))
    lab_font = {"fontsize": "xx-small"}
    for name in results:
        with_label = True
        for i, (arg, short) in enumerate(VARIABLES.items()):
            x_vals = INTERVALS_VALS[short]
            y_vals = results[name][arg]
            # fit a line to each so that we can report the slope
            slope = np.polyfit(x_vals, y_vals, 1)[0]
            axs[i].plot(
                x_vals,
                y_vals,
                marker="o",
                label=(ALG_NAME[name] if with_label else None),
            )
            with_label = False
            axs[i].text(
                x_vals[-1],
                y_vals[-1] + (y_vals[-1] / 16),
                f"m = {slope:.3E}",
                fontdict=lab_font,
            )
    for i, (arg, short) in enumerate(VARIABLES.items()):
        label = f"Number of {arg}"
        for other, other_short in VARIABLES.items():
            if other != arg:
                label += f"\nnum_{other} = {DEFAULT_VAL[other_short]}"
        axs[i].set_xlabel(label)
        axs[i].set_ylim(ymin=0)
    fig.supylabel("CPU Time (s)")
    fig.legend(loc="lower left", fontsize="x-small")
    fig.tight_layout()
    fig.savefig(str(output), bbox_inches="tight", dpi=400)

    # saving results
    if NAME is not None and archive is not None:
        with open(archive, "wb") as fh:
            pickle.dump(results, fh)


if __name__ == "__main__":
    main()
