import os
from pathlib import Path

import pytest
import numpy as np
import numpy.lib.recfunctions as rfn
from click.testing import CliRunner

from haptools.transform import (
    HaplotypeAncestry,
    HaplotypesAncestry,
    GenotypesAncestry,
)
from haptools.__main__ import main
from .test_data import TestGenotypesVCF, TestHaplotypes
from haptools.data import Variant, GenotypesVCF, GenotypesPLINK


DATADIR = Path(__file__).parent.joinpath("data")


class TestGenotypesAncestry:
    file = DATADIR / "simple-ancestry.vcf"

    def _get_expected_ancestry(self):
        # create a matrix with shape: samples x SNPs x strands
        expected = np.zeros(40).reshape((5, 4, 2)).astype(np.uint8)
        expected[0, 0, 1] = 1
        expected[1, 1, 1] = 1
        expected[2, 2, 0] = 1
        expected[3, 3, 0] = 1
        expected[4, 1, 1] = 1
        expected[3, 1, :] = 1
        expected[3, 0, 0] = 2
        expected[0, 2, 1] = 2
        expected[-1, -1, :] = 2
        return expected

    def _get_fake_genotypes(self):
        base_gts = TestGenotypesVCF()._get_fake_genotypes_refalt(with_phase=True)
        # copy all of the fields
        gts = GenotypesAncestry(fname=None)
        gts.data = base_gts.data
        gts.samples = base_gts.samples
        gts.variants = base_gts.variants
        # add additional ancestry data
        gts.ancestry = self._get_expected_ancestry()
        gts.ancestry_labels = {"YRI": 0, "CEU": 1, "ASW": 2}
        return gts

    def test_load_genotypes_ancestry(self, caplog):
        expected = self._get_fake_genotypes()

        # can we load the data from the VCF?
        gts = GenotypesAncestry(self.file)
        gts.read()
        np.testing.assert_allclose(gts.data, expected.data)
        assert gts.ancestry_labels == expected.ancestry_labels
        np.testing.assert_allclose(gts.ancestry, expected.ancestry)
        assert gts.samples == expected.samples

    def test_load_genotypes_iterate(self, caplog):
        expected = self._get_expected_ancestry().transpose((1, 0, 2))

        # can we load the data from the VCF?
        gts = GenotypesAncestry(self.file)
        for idx, line in enumerate(gts):
            np.testing.assert_allclose(line.ancestry, expected[idx])

    def test_load_genotypes_discard_multiallelic(self):
        gts = self._get_fake_genotypes()

        # make a copy for later
        data_copy = gts.data.copy().astype(np.bool_)
        ancestry_copy = gts.ancestry.copy()
        variant_shape = list(gts.variants.shape)
        variant_shape[0] -= 1

        # force one of the SNPs to have more than one allele and check that it gets dicarded
        gts.data[1, 1, 1] = 2
        gts.check_biallelic(discard_also=True)

        data_copy_without_biallelic = np.delete(data_copy, [1], axis=1)
        ancestry_copy_without_biallelic = np.delete(ancestry_copy, [1], axis=1)
        np.testing.assert_equal(gts.data, data_copy_without_biallelic)
        np.testing.assert_equal(gts.ancestry, ancestry_copy_without_biallelic)
        assert gts.variants.shape == tuple(variant_shape)

    def test_subset_genotypes(self):
        gts = self._get_fake_genotypes()

        # subset to just the samples we want
        expected_data = gts.data[:3]
        expected_ancestry = gts.ancestry[:3]
        expected_variants = gts.variants
        samples = ("HG00096", "HG00097", "HG00099")
        gts_sub = gts.subset(samples=samples)
        assert gts_sub.samples == samples
        np.testing.assert_allclose(gts_sub.data, expected_data)
        assert gts_sub.ancestry_labels == gts.ancestry_labels
        np.testing.assert_allclose(gts_sub.ancestry, expected_ancestry)
        assert np.array_equal(gts_sub.variants, expected_variants)

        # subset to just the variants we want
        expected_data = gts.data[:, [1, 2]]
        expected_ancestry = gts.ancestry[:, [1, 2]]
        expected_variants = gts.variants[[1, 2]]
        variants = ("1:10116:A:G", "1:10117:C:A")
        gts_sub = gts.subset(variants=variants)
        assert gts_sub.samples == gts.samples
        np.testing.assert_allclose(gts_sub.data, expected_data)
        assert gts_sub.ancestry_labels == gts.ancestry_labels
        np.testing.assert_allclose(gts_sub.ancestry, expected_ancestry)
        assert np.array_equal(gts_sub.variants, expected_variants)

        # subset both: samples and variants
        expected_data = gts.data[[3, 4], :2]
        expected_ancestry = gts.ancestry[[3, 4], :2]
        expected_variants = gts.variants[:2]
        samples = ("HG00100", "HG00101")
        variants = ("1:10114:T:C", "1:10116:A:G")
        gts_sub = gts.subset(samples=samples, variants=variants)
        assert gts_sub.samples == samples
        np.testing.assert_allclose(gts_sub.data, expected_data)
        assert gts_sub.ancestry_labels == gts.ancestry_labels
        np.testing.assert_allclose(gts_sub.ancestry, expected_ancestry)
        assert np.array_equal(gts_sub.variants, expected_variants)

    @pytest.mark.xfail(reason="not implemented yet")
    def test_write_genotypes(self):
        assert False

    @pytest.mark.xfail(reason="not implemented yet")
    def test_write_genotypes_phase(self, prephased=True):
        assert False


class TestHaplotypesAncestry:
    def _get_dummy_haps(self):
        # create three haplotypes
        haplotypes = {
            "H1": HaplotypeAncestry(
                chrom="1", start=10114, end=8, id="H1", ancestry="YRI"
            ),
            "H2": HaplotypeAncestry(
                chrom="1", start=10114, end=10119, id="H2", ancestry="YRI"
            ),
            "H3": HaplotypeAncestry(
                chrom="1", start=10116, end=10119, id="H3", ancestry="ASW"
            ),
        }
        haplotypes["H1"].variants = (
            Variant(start=10114, end=10115, id="1:10114:T:C", allele="T"),
            Variant(start=10116, end=10117, id="1:10116:A:G", allele="G"),
        )
        haplotypes["H2"].variants = (
            Variant(start=10114, end=10115, id="1:10114:T:C", allele="C"),
            Variant(start=10117, end=10118, id="1:10117:C:A", allele="C"),
        )
        haplotypes["H3"].variants = (
            Variant(start=10116, end=10117, id="1:10116:A:G", allele="A"),
            Variant(start=10117, end=10118, id="1:10117:C:A", allele="A"),
        )
        haps = HaplotypesAncestry(fname=None)
        haps.data = haplotypes
        return haps

    def test_hap_transform(self):
        expected = np.array(
            [
                [0, 0],
                [0, 0],
                [1, 1],
                [0, 0],
                [0, 0],
            ],
            dtype=np.uint8,
        )

        hap = list(self._get_dummy_haps().data.values())[0]
        gens = TestGenotypesAncestry()._get_fake_genotypes()
        gens.check_phase()
        hap_gt = hap.transform(gens)
        np.testing.assert_allclose(hap_gt, expected)

    def test_haps_transform(self):
        expected = np.array(
            [
                [[0, 0], [0, 0], [0, 0]],
                [[0, 0], [0, 0], [0, 0]],
                [[1, 1], [0, 0], [0, 0]],
                [[0, 0], [0, 0], [0, 0]],
                [[0, 0], [0, 0], [0, 0]],
            ],
            dtype=np.uint8,
        )

        haps = self._get_dummy_haps()
        gens = TestGenotypesAncestry()._get_fake_genotypes()
        gens.check_phase()

        hap_gt = GenotypesVCF(fname=None)
        haps.transform(gens, hap_gt)
        np.testing.assert_allclose(hap_gt.data, expected)

        expected[2, 0, 1] = 0
        expected[[2, 4], 1, 1] = 1
        gens.data[[2, 4], 0, 1] = 1
        gens.data[[1, 4], 2, 0] = 1

        hap_gt = GenotypesVCF(fname=None)
        haps.transform(gens, hap_gt)
        np.testing.assert_allclose(hap_gt.data, expected)


def test_basic(capfd):
    expected = """##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##contig=<ID=1>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG00096\tHG00097\tHG00099\tHG00100\tHG00101
1\t10114\tH1\tA\tT\t.\t.\t.\tGT\t0|1\t0|1\t1|1\t1|1\t0|0
1\t10114\tH2\tA\tT\t.\t.\t.\tGT\t0|0\t0|0\t0|0\t0|0\t0|0
1\t10116\tH3\tA\tT\t.\t.\t.\tGT\t0|0\t0|0\t0|0\t0|0\t0|0
"""
    gt_file = DATADIR / "simple.vcf.gz"
    hp_file = DATADIR / "simple.hap"

    cmd = f"transform {gt_file} {hp_file}"
    runner = CliRunner()
    result = runner.invoke(main, cmd.split(" "), catch_exceptions=False)
    captured = capfd.readouterr()
    assert captured.out == expected
    assert result.exit_code == 0


def test_basic_multiallelic(capfd):
    expected = """##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##contig=<ID=1>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG00096\tHG00097\tHG00099\tHG00100\tHG00101
1\t10114\tH1\tA\tT\t.\t.\t.\tGT\t0|0\t0|0\t1|0\t0|0\t0|1
1\t10114\tH2\tA\tT\t.\t.\t.\tGT\t0|0\t0|0\t0|0\t0|0\t0|0
1\t10116\tH3\tA\tT\t.\t.\t.\tGT\t0|0\t0|0\t0|0\t0|0\t0|0
"""
    gt_file = DATADIR / "simple-multiallelic.vcf"
    hp_file = DATADIR / "simple-multiallelic.hap"

    # first, create the multiallelic hap file
    hp = TestHaplotypes()._get_dummy_haps_multiallelic()
    hp.fname = hp_file
    hp.write()

    cmd = f"transform {gt_file} {hp_file}"
    runner = CliRunner()
    result = runner.invoke(main, cmd.split(" "), catch_exceptions=False)
    captured = capfd.readouterr()
    assert captured.out == expected
    assert result.exit_code == 0

    hp.fname.unlink()


def test_basic_subset(capfd):
    expected = """##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##contig=<ID=1>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG00096\tHG00097\tHG00099\tHG00100\tHG00101
1\t10114\tH1\tA\tT\t.\t.\t.\tGT\t0|1\t0|1\t1|1\t1|1\t0|0
"""
    hp_file = DATADIR / "simple.hap"

    # first, remove the last two variants in the genotypes file
    gts = GenotypesVCF.load(DATADIR / "simple.vcf")
    gts.fname = Path("simple_minus_two.vcf")
    gts.subset(variants=tuple(gts.variants["id"][:-2]), inplace=True)
    gts.write()

    cmd = f"transform simple_minus_two.vcf {hp_file}"
    runner = CliRunner()
    result = runner.invoke(main, cmd.split(" "), catch_exceptions=False)
    captured = capfd.readouterr()
    assert captured.out == expected
    assert result.exit_code == 0

    gts.fname.unlink()


def test_basic_pgen_input(capfd):
    expected = """##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##contig=<ID=1>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG00096\tHG00097\tHG00099\tHG00100\tHG00101
1\t10114\tH1\tA\tT\t.\t.\t.\tGT\t0|1\t0|1\t1|1\t1|1\t0|0
1\t10114\tH2\tA\tT\t.\t.\t.\tGT\t0|0\t0|0\t0|0\t0|0\t0|0
1\t10116\tH3\tA\tT\t.\t.\t.\tGT\t0|0\t0|0\t0|0\t0|0\t0|0
"""
    gt_file = DATADIR / "simple.pgen"
    hp_file = DATADIR / "simple.hap"

    cmd = f"transform {gt_file} {hp_file}"
    runner = CliRunner()
    result = runner.invoke(main, cmd.split(" "), catch_exceptions=False)
    captured = capfd.readouterr()
    assert captured.out == expected
    assert result.exit_code == 0


def test_pgen_two_samples(capfd):
    pytest.importorskip("pgenlib")
    expected = np.array(
        [
            [[0, 1, 1]],
            [[0, 0, 1]],
        ],
        dtype=np.uint8,
    )
    gt_file = DATADIR / "apoe.vcf.gz"
    hp_file = DATADIR / "apoe4.hap"

    cmd = f"transform -o output.pgen -s HG00097 -s NA12878 {gt_file} {hp_file}"
    runner = CliRunner()
    result = runner.invoke(main, cmd.split(" "), catch_exceptions=False)
    captured = capfd.readouterr()
    assert captured.out == ""
    output = Path("output.pgen")
    gt = GenotypesPLINK(output)
    gt.read()
    np.testing.assert_allclose(gt.data, expected)
    assert tuple(gt.variants[["id", "chrom", "pos"]][0]) == ("APOe4", "19", 45411941)
    assert gt.samples == ("HG00097", "NA12878")
    output.unlink()
    output.with_suffix(".pvar").unlink()
    output.with_suffix(".psam").unlink()
    assert result.exit_code == 0


ancestry_results = """##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##contig=<ID=1>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG00096\tHG00097\tHG00099\tHG00100\tHG00101
1\t10114\tH1\tA\tT\t.\t.\t.\tGT\t0|0\t0|0\t1|1\t0|0\t0|0
1\t10114\tH2\tA\tT\t.\t.\t.\tGT\t0|0\t0|0\t0|0\t0|0\t0|0
1\t10116\tH3\tA\tT\t.\t.\t.\tGT\t0|0\t0|0\t0|0\t0|0\t0|0
"""


def test_ancestry_from_vcf(capfd):
    gt_file = DATADIR / "simple-ancestry.vcf"
    hp_file = DATADIR / "simple.hap"

    cmd = f"transform --ancestry {gt_file} {hp_file}"
    runner = CliRunner()
    result = runner.invoke(main, cmd.split(" "), catch_exceptions=False)
    captured = capfd.readouterr()
    assert captured.out == ancestry_results
    assert result.exit_code == 0


def test_ancestry_from_bp(capfd):
    gt_file = DATADIR / "simple.vcf"
    hp_file = DATADIR / "simple.hap"

    cmd = f"transform --ancestry {gt_file} {hp_file}"
    runner = CliRunner()
    result = runner.invoke(main, cmd.split(" "), catch_exceptions=False)
    captured = capfd.readouterr()
    assert captured.out == ancestry_results
    assert result.exit_code == 0


def test_transform_empty_hap(capfd):
    gt_file = DATADIR / "simple.vcf.gz"
    hp_file = Path("empty.hap")
    hp_file_gz = Path("empty.hap.gz")
    hp_file_idx = Path("empty.hap.gz.tbi")

    # create an empty .hap file
    with open(hp_file, "w") as f:
        f.write("")

    # can we run transform with the empty hap file?
    cmd = f"transform --region 1:10116-10122 {gt_file} {hp_file}"
    runner = CliRunner()
    result = runner.invoke(main, cmd.split(" "))
    captured = capfd.readouterr()
    assert all(line for line in captured.out.split("\n") if line.startswith("#"))
    assert result.exit_code != 0

    # now, index the empty hap file and try again
    cmd = f"index {hp_file}"
    runner = CliRunner()
    result = runner.invoke(main, cmd.split(" "), catch_exceptions=False)
    captured = capfd.readouterr()
    assert result.exit_code == 0
    assert hp_file_gz.exists()
    assert hp_file_idx.exists()

    # what about now? does it still fail?
    cmd = f"transform --region 1:10116-10122 {gt_file} {hp_file_gz}"
    runner = CliRunner()
    result = runner.invoke(main, cmd.split(" "))
    captured = capfd.readouterr()
    assert result.exit_code != 0

    hp_file.unlink()
    hp_file_gz.unlink()
    hp_file_idx.unlink()
