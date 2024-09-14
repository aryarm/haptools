import os
from pathlib import Path
from dataclasses import dataclass, field

import pytest
import numpy as np
import numpy.lib.recfunctions as rfn
from haptools.sim_phenotype import Haplotype as HaptoolsHaplotype
from haptools.sim_phenotype import Repeat as HaptoolsRepeat
from haptools.data import (
    Extra,
    Repeat,
    Variant,
    HapBlock,
    Haplotype,
    Genotypes,
    Phenotypes,
    Covariates,
    Haplotypes,
    Breakpoints,
    GenotypesTR,
    GenotypesVCF,
    GenotypesPLINK,
    GenotypesPLINKTR,
)


DATADIR = Path(__file__).parent.joinpath("data")


class TestGenotypes:
    def _get_expected_genotypes(self):
        # create a GT matrix with shape: samples x SNPs x (strands+phase)
        expected = np.zeros(60).reshape((5, 4, 3)).astype(np.uint8)
        expected[:4, 1, 1] = 1
        expected[2:4, 1, 0] = 1
        expected[:, :, 2] = 1
        return expected

    def _get_fake_genotypes(self):
        gts = Genotypes(fname=None)
        gts.data = self._get_expected_genotypes()
        gts.variants = np.array(
            [
                ("1:10114:T:C", "1", 10114),
                ("1:10116:A:G", "1", 10116),
                ("1:10117:C:A", "1", 10117),
                ("1:10122:A:G", "1", 10122),
            ],
            dtype=gts.variants.dtype,
        )
        gts.samples = ("HG00096", "HG00097", "HG00099", "HG00100", "HG00101")
        gts.check_phase()
        return gts

    def test_load_genotypes(self, caplog):
        expected = self._get_expected_genotypes()

        # can we load the data from the VCF?
        gts = Genotypes(DATADIR / "simple.vcf")
        gts.read()
        np.testing.assert_allclose(gts.data, expected)
        assert gts.samples == ("HG00096", "HG00097", "HG00099", "HG00100", "HG00101")

        # try loading the data again - it should warn b/c we've already done it
        caplog.clear()
        gts.read()
        assert len(caplog.records) > 0 and caplog.records[0].levelname == "WARNING"

        # force one of the samples to have a missing GT and check that we get an error
        # we convert to an int8 and then convert back
        gts.data = gts.data.astype(np.int8)
        gts.data[1, 1, 1] = -1
        gts.data = gts.data.astype(np.uint8)
        with pytest.raises(ValueError) as info:
            gts.check_missing()
        assert (
            str(info.value)
            == "Genotype with ID 1:10116:A:G at POS 1:10116 is missing for sample"
            " HG00097"
        )
        gts.data[1, 1, 1] = 1

        # force one of the SNPs to have more than one allele and check that we get an error
        gts.data[1, 1, 1] = 2
        with pytest.raises(ValueError) as info:
            gts.check_biallelic()
        assert (
            str(info.value)
            == "Variant with ID 1:10116:A:G at POS 1:10116 is multiallelic for sample"
            " HG00097"
        )
        gts.data[1, 1, 1] = 1

        # check biallelic-ness and convert to bool_
        gts.check_biallelic()
        expected = expected.astype(np.bool_)
        np.testing.assert_allclose(gts.data, expected)

        # force one of the het SNPs to be unphased and check that we get an error message
        gts.data[1, 1, 2] = 0
        with pytest.raises(ValueError) as info:
            gts.check_phase()
        assert (
            str(info.value)
            == "Variant with ID 1:10116:A:G at POS 1:10116 is unphased for sample"
            " HG00097"
        )
        gts.data[1, 1, 2] = 1

        # check phase and remove the phase axis
        gts.check_phase()
        expected = expected[:, :, :2]
        np.testing.assert_allclose(gts.data, expected)

        # try to check phase again - it should warn b/c we've already done it before
        caplog.clear()
        gts.check_phase()
        assert len(caplog.records) > 0 and caplog.records[0].levelname == "WARNING"

    def test_load_genotypes_example(self):
        samples = (
            "HG00096",
            "HG00097",
            "HG00099",
            "HG00100",
            "HG00101",
            "HG00102",
            "HG00103",
            "HG00105",
            "HG00106",
            "HG00107",
            "HG00108",
            "HG00109",
            "HG00110",
            "HG00111",
            "HG00112",
            "HG00113",
            "HG00114",
            "HG00115",
            "HG00116",
            "HG00117",
            "HG00118",
            "HG00119",
            "HG00120",
            "HG00121",
            "HG00122",
        )
        gts = Genotypes(DATADIR / "example.vcf.gz")
        gts.read()
        assert gts.samples[:25] == samples

    def test_load_genotypes_iterate(self, caplog):
        expected = self._get_fake_genotypes()

        # can we load the data from the VCF?
        gts = Genotypes(DATADIR / "simple.vcf")
        for idx, line in enumerate(gts):
            np.testing.assert_allclose(line.data[:, :2], expected.data[:, idx])
            for col in ("chrom", "pos", "id"):
                assert line.variants[col] == expected.variants[col][idx]
        assert gts.samples == expected.samples

    def test_load_genotypes_discard_multiallelic(self):
        gts = self._get_fake_genotypes()

        # make a copy for later
        data_copy = gts.data.copy().astype(np.bool_)
        variant_shape = list(gts.variants.shape)
        variant_shape[0] -= 1

        # force one of the SNPs to have more than one allele and check that it gets dicarded
        gts.data[1, 1, 1] = 2
        gts.check_biallelic(discard_also=True)

        data_copy_without_biallelic = np.delete(data_copy, [1], axis=1)
        np.testing.assert_equal(gts.data, data_copy_without_biallelic)
        assert gts.variants.shape == tuple(variant_shape)

    def test_load_genotypes_subset(self):
        expected = self._get_expected_genotypes()

        # subset for the region we want
        expected = expected[:, 1:3]

        # can we load the data from the VCF?
        gts = Genotypes(DATADIR / "simple.vcf.gz")
        gts.read(region="1:10115-10117")
        np.testing.assert_allclose(gts.data, expected)
        assert gts.samples == ("HG00096", "HG00097", "HG00099", "HG00100", "HG00101")

        # subset for just the samples we want
        expected = expected[[1, 3]]

        gts = Genotypes(DATADIR / "simple.vcf.gz")
        samples = ["HG00097", "HG00100"]
        samples_set = set(samples)
        gts.read(region="1:10115-10117", samples=samples_set)
        np.testing.assert_allclose(gts.data, expected)
        assert gts.samples == tuple(samples)

        # subset to just one of the variants
        expected = expected[:, [1]]

        gts = Genotypes(DATADIR / "simple.vcf.gz")
        variants = {"1:10117:C:A"}
        gts.read(region="1:10115-10117", samples=samples_set, variants=variants)
        np.testing.assert_allclose(gts.data, expected)
        assert gts.samples == tuple(samples)

    def test_subset_genotypes(self):
        gts = self._get_fake_genotypes()

        # subset to just the samples we want
        expected_data = gts.data[:3]
        expected_variants = gts.variants
        samples = ("HG00096", "HG00097", "HG00099")
        gts_sub = gts.subset(samples=samples)
        assert gts_sub.samples == samples
        np.testing.assert_allclose(gts_sub.data, expected_data)
        assert np.array_equal(gts_sub.variants, expected_variants)

        # subset to just the variants we want
        expected_data = gts.data[:, [1, 2]]
        expected_variants = gts.variants[[1, 2]]
        variants = ("1:10116:A:G", "1:10117:C:A")
        gts_sub = gts.subset(variants=variants)
        assert gts_sub.samples == gts.samples
        np.testing.assert_allclose(gts_sub.data, expected_data)
        assert np.array_equal(gts_sub.variants, expected_variants)

        # subset both: samples and variants
        expected_data = gts.data[[3, 4], :2]
        expected_variants = gts.variants[:2]
        samples = ("HG00100", "HG00101")
        variants = ("1:10114:T:C", "1:10116:A:G")
        gts_sub = gts.subset(samples=samples, variants=variants)
        assert gts_sub.samples == samples
        np.testing.assert_allclose(gts_sub.data, expected_data)
        assert np.array_equal(gts_sub.variants, expected_variants)

    def test_check_maf(self, caplog):
        gts = self._get_fake_genotypes()
        expected_maf = np.array([0, 0.4, 0, 0])

        maf = gts.check_maf()
        np.testing.assert_allclose(maf, expected_maf)

        msg = "Variant with ID 1:10114:T:C at POS 1:10114 has MAF 0.0 < 0.01"
        with pytest.raises(ValueError) as info:
            gts.check_maf(threshold=0.01)
        assert str(info.value) == msg

        # test just the warning system
        caplog.clear()
        maf = gts.check_maf(threshold=0.01, warn_only=True)
        assert len(caplog.records) > 0 and caplog.records[0].levelname == "WARNING"

        maf = gts.check_maf(threshold=0, discard_also=True)
        np.testing.assert_allclose(maf, expected_maf)
        assert len(gts.variants) == 4
        assert gts.data.shape[1] == 4

        maf = gts.check_maf(threshold=0.01, discard_also=True)
        np.testing.assert_allclose(maf, expected_maf[1])
        assert len(gts.variants) == 1
        assert gts.data.shape[1] == 1

        maf = gts.check_maf(threshold=0.5, discard_also=True)
        np.testing.assert_allclose(maf, expected_maf[:0])
        assert len(gts.variants) == 0
        assert gts.data.shape[1] == 0

    def test_check_sorted(self, caplog):
        gts = self._get_fake_genotypes()
        gts.check_sorted()

        # swap two of the rows
        gts.variants[[0, 2]] = gts.variants[[2, 0]]

        with pytest.raises(ValueError) as info:
            gts.check_sorted()

    def test_merge_variants(self):
        gts1 = self._get_fake_genotypes()
        gts2 = Genotypes(DATADIR / "example.vcf.gz")
        gts2.read()
        gts2.check_phase()
        gts2.subset(samples=gts1.samples, inplace=True)

        gts = Genotypes.merge_variants((gts1, gts2), fname=None)

        assert gts.samples == gts1.samples
        assert len(gts.variants) == len(gts1.variants) + len(gts2.variants)
        assert gts.data.shape[0] == gts1.data.shape[0]
        assert gts.data.shape[2] == gts1.data.shape[2]
        assert gts.data.shape[1] == (gts1.data.shape[1] + gts2.data.shape[1])

        # also try it when the phases are different
        gts2 = Genotypes(DATADIR / "example.vcf.gz")
        gts2.read()
        gts2.subset(samples=gts1.samples, inplace=True)

        gts = Genotypes.merge_variants((gts1, gts2), fname=None)

        assert gts.samples == gts1.samples
        assert len(gts.variants) == len(gts1.variants) + len(gts2.variants)
        assert gts.data.shape[0] == gts1.data.shape[0]
        assert gts.data.shape[2] == gts2.data.shape[2]
        assert gts.data.shape[1] == (gts1.data.shape[1] + gts2.data.shape[1])


class TestGenotypesPLINK:
    def _get_fake_genotypes_plink(self):
        gts_ref_alt = TestGenotypesVCF()._get_fake_genotypes_refalt()
        gts = GenotypesPLINK(gts_ref_alt.fname)
        gts.data = gts_ref_alt.data
        gts.samples = gts_ref_alt.samples
        gts.variants = gts_ref_alt.variants
        return gts

    def _get_fake_genotypes_multiallelic(self):
        gts_ref_alt = TestGenotypesVCF()._get_fake_genotypes_multiallelic()
        gts = GenotypesPLINK(gts_ref_alt.fname)
        gts.data = gts_ref_alt.data
        gts.samples = gts_ref_alt.samples
        gts.variants = gts_ref_alt.variants
        return gts

    def _get_fake_genotypes_multiallelic_tr(self):
        gts_tr = GenotypesTR(DATADIR / "simple-tr.vcf")
        gts_tr.read()

        gts = GenotypesPLINK(fname="")
        gts.samples = ("HG00096", "HG00097", "HG00099", "HG00100", "HG00101")

        na = np.iinfo(np.uint8).max
        data = np.array(
            [
                [[0, 1], [2, 3], [4, 5], [6, 7], [8, 8]],
                [[0, 1], [0, 1], [1, 1], [1, 1], [0, 0]],
                [[0, 0], [0, 0], [0, 0], [0, 0], [0, 0]],
                [[4, 4], [na, na], [3, 0], [1, 2], [4, na]],
                [[1, na], [na, na], [0, na], [na, na], [na, na]],
            ],
            dtype=np.uint8,
        )
        # append the phasing info, transpose, and clean up
        phasing = np.ones(data.shape[:2] + (1,)).astype(np.uint8)
        phasing[4, [1, 3, 4]] = 0
        phasing[3, [1, 4]] = 0
        phasing[:2, 1:3] = 0
        phasing[1, 2] = 1  # since homozygous unphased variants are marked as phased
        gts.data = np.concatenate((data, phasing), axis=-1).transpose((1, 0, 2))
        # handle half-calls and chrX according to the flag:
        # --vcf-half-call m
        gts.data[4, 3, 0] = na
        gts.data[:, 4, 1] = gts.data[:, 4, 0]
        # TODO: figure out whether this is the correct way to handle half-calls and try
        # to see if we can change chrX representations to properly simulate from it

        # add the appropriate alleles
        ref_alt = [
            tuple("GTT" * i for i in range(1, 10)),
            ("ACACAC", "AC"),
            ("AAA", "AAAA"),
            (
                "GTGT",
                "GTGTGT",
                "GTGTGTGT",
                "GTGTGTGTGTGTGT",
                "GTGTGTGTGTGTGTGTGTGTGT",
                "GT",
            ),
            ("GGG", "GGGGG"),
        ]
        gts.variants = np.array(
            [tuple(rec) + (ref_alt[idx],) for idx, rec in enumerate(gts_tr.variants)],
            dtype=gts.variants.dtype,
        )

        return gts

    def test_load_genotypes(self):
        expected = self._get_fake_genotypes_plink()

        gts = GenotypesPLINK(DATADIR / "simple.pgen")
        gts.read()
        gts.check_phase()

        # check that everything matches what we expected
        np.testing.assert_allclose(gts.data, expected.data)
        assert gts.samples == expected.samples
        for i, x in enumerate(expected.variants):
            for col in ("chrom", "pos", "id", "alleles"):
                assert gts.variants[col][i] == expected.variants[col][i]

    def test_load_genotypes_multiallelic(self):
        expected = self._get_fake_genotypes_multiallelic()

        gts = GenotypesPLINK(DATADIR / "simple-multiallelic.pgen")
        gts.read()
        gts.check_phase()

        # check that everything matches what we expected
        np.testing.assert_allclose(gts.data, expected.data)
        assert gts.samples == expected.samples
        for i, x in enumerate(expected.variants):
            for col in ("chrom", "pos", "id", "alleles"):
                assert gts.variants[col][i] == expected.variants[col][i]

    def test_load_genotypes_multiallelic_tr(self):
        expected = self._get_fake_genotypes_multiallelic_tr()

        gts = GenotypesPLINK(DATADIR / "simple-tr.pgen")
        gts.read()

        # check that everything matches what we expected
        np.testing.assert_array_equal(gts.data, expected.data)
        assert gts.samples == expected.samples
        for i, x in enumerate(expected.variants):
            for col in ("chrom", "pos", "id", "alleles"):
                assert gts.variants[col][i] == expected.variants[col][i]

    def test_load_genotypes_chunked(self):
        expected = self._get_fake_genotypes_plink()

        gts = GenotypesPLINK(DATADIR / "simple.pgen", chunk_size=1)
        gts.read()
        gts.check_phase()

        # check that everything matches what we expected
        np.testing.assert_allclose(gts.data, expected.data)
        assert gts.samples == expected.samples
        for i, x in enumerate(expected.variants):
            for col in ("chrom", "pos", "id", "alleles"):
                assert gts.variants[col][i] == expected.variants[col][i]

    def test_load_genotypes_chunked_cpus(self):
        gts = GenotypesPLINK(DATADIR / "example.pgen", num_cpus=1)
        gts.read()
        # TODO: test with missing values
        # and with different combinations of num_cpus and chunk_size
        breakpoint()

        # check that everything matches what we expected
        np.testing.assert_allclose(gts.data, expected.data)
        assert gts.samples == expected.samples
        for i, x in enumerate(expected.variants):
            for col in ("chrom", "pos", "id", "alleles"):
                assert gts.variants[col][i] == expected.variants[col][i]

    def test_load_genotypes_prephased(self):
        expected = self._get_fake_genotypes_plink()

        gts = GenotypesPLINK(DATADIR / "simple.pgen")
        gts._prephased = True
        gts.read()

        # check that everything matches what we expected
        np.testing.assert_allclose(gts.data, expected.data)
        assert gts.samples == expected.samples
        for i, x in enumerate(expected.variants):
            for col in ("chrom", "pos", "id", "alleles"):
                assert gts.variants[col][i] == expected.variants[col][i]

    def test_load_genotypes_iterate(self):
        expected = self._get_fake_genotypes_plink()

        gts = GenotypesPLINK(DATADIR / "simple.pgen")

        # check that everything matches what we expected
        for idx, line in enumerate(gts):
            np.testing.assert_allclose(line.data[:, :2], expected.data[:, idx])
            for col in ("chrom", "pos", "id"):
                assert line.variants[col] == expected.variants[col][idx]
            assert (
                line.variants["alleles"].tolist() == expected.variants["alleles"][idx]
            )
        assert gts.samples == expected.samples

    def test_iter_multiallelic_tr(self):
        expected = self._get_fake_genotypes_multiallelic_tr()

        gts = GenotypesPLINK(DATADIR / "simple-tr.pgen")

        # Check that everything matches what we expected
        for idx, line in enumerate(gts):
            np.testing.assert_allclose(line.data[:, :3], expected.data[:, idx])
            for col in ("chrom", "pos", "id"):
                assert line.variants[col] == expected.variants[col][idx]
            assert (
                line.variants["alleles"].tolist() == expected.variants["alleles"][idx]
            )

        # Check samples
        assert gts.samples == expected.samples

    def test_load_genotypes_subset(self):
        expected = self._get_fake_genotypes_plink()

        # subset for the region we want
        expected_data = expected.data[:, 1:3]

        # can we load the data from the VCF?
        gts = GenotypesPLINK(DATADIR / "simple.pgen")
        gts.read(region="1:10115-10117")
        gts.check_phase()
        np.testing.assert_allclose(gts.data, expected_data)
        assert gts.samples == expected.samples

        # subset for just the samples we want
        expected_data = expected_data[[1, 3]]

        gts = GenotypesPLINK(DATADIR / "simple.pgen")
        samples = [expected.samples[1], expected.samples[3]]
        samples_set = set(samples)
        gts.read(region="1:10115-10117", samples=samples_set)
        gts.check_phase()
        np.testing.assert_allclose(gts.data, expected_data)
        assert gts.samples == tuple(samples)

        # subset to just one of the variants
        expected_data = expected_data[:, [1]]

        gts = GenotypesPLINK(DATADIR / "simple.pgen")
        variants = {"1:10117:C:A"}
        gts.read(region="1:10115-10117", samples=samples_set, variants=variants)
        gts.check_phase()
        np.testing.assert_allclose(gts.data, expected_data)
        assert gts.samples == tuple(samples)

    def test_write_genotypes_chunked(self):
        gts = self._get_fake_genotypes_plink()

        fname = DATADIR / "test_write_chunked.pgen"
        gts.fname = fname
        gts.write()

        new_gts = GenotypesPLINK(fname, chunk_size=1)
        new_gts.read()
        new_gts.check_phase()

        # check that everything matches what we expected
        np.testing.assert_allclose(gts.data, new_gts.data)
        assert gts.samples == new_gts.samples
        for i in range(len(new_gts.variants)):
            for col in ("chrom", "pos", "id", "alleles"):
                assert gts.variants[col][i] == new_gts.variants[col][i]

        # clean up afterwards: delete the files we created
        fname.with_suffix(".psam").unlink()
        fname.with_suffix(".pvar").unlink()
        fname.unlink()

    def test_write_genotypes(self):
        gts = self._get_fake_genotypes_plink()

        fname = DATADIR / "test_write.pgen"
        gts.fname = fname
        gts.write()

        new_gts = GenotypesPLINK(fname)
        new_gts.read()
        new_gts.check_phase()

        # check that everything matches what we expected
        np.testing.assert_allclose(gts.data, new_gts.data)
        assert gts.samples == new_gts.samples
        for i in range(len(new_gts.variants)):
            for col in ("chrom", "pos", "id", "alleles"):
                assert gts.variants[col][i] == new_gts.variants[col][i]

        # clean up afterwards: delete the files we created
        fname.with_suffix(".psam").unlink()
        fname.with_suffix(".pvar").unlink()
        fname.unlink()

    def test_write_genotypes_prephased(self):
        gts = self._get_fake_genotypes_plink()

        fname = DATADIR / "test_write.pgen"
        gts.fname = fname
        gts._prephased = True
        gts.write()

        new_gts = GenotypesPLINK(fname)
        new_gts._prephased = True
        new_gts.read()

        # check that everything matches what we expected
        np.testing.assert_allclose(gts.data, new_gts.data)
        assert gts.samples == new_gts.samples
        for i in range(len(new_gts.variants)):
            for col in ("chrom", "pos", "id", "alleles"):
                assert gts.variants[col][i] == new_gts.variants[col][i]

        # clean up afterwards: delete the files we created
        fname.with_suffix(".psam").unlink()
        fname.with_suffix(".pvar").unlink()
        fname.unlink()

    def test_write_genotypes_unphased(self):
        gts = self._get_fake_genotypes_plink()
        # add phasing information back
        gts.data = np.dstack((gts.data, np.ones(gts.data.shape[:2], dtype=np.uint8)))
        gts.data[:2, 1, 2] = 0

        fname = DATADIR / "test_unphased.pgen"
        gts.fname = fname
        gts.write()

        new_gts = GenotypesPLINK(fname)
        new_gts.read()

        # check that everything matches what we expected
        np.testing.assert_allclose(gts.data, new_gts.data)
        assert gts.samples == new_gts.samples
        for i in range(len(new_gts.variants)):
            for col in ("chrom", "pos", "id", "alleles"):
                assert gts.variants[col][i] == new_gts.variants[col][i]

        # clean up afterwards: delete the files we created
        fname.with_suffix(".psam").unlink()
        fname.with_suffix(".pvar").unlink()
        fname.unlink()

    def test_write_genotypes_biallelic(self):
        gts = self._get_fake_genotypes_plink()

        # convert to bool - we should be able to handle both
        gts.data = gts.data.astype(np.bool_)

        fname = DATADIR / "test_write.pgen"
        gts.fname = fname
        gts.write()

        new_gts = GenotypesPLINK(fname)
        new_gts.read()
        new_gts.check_phase()

        # check that everything matches what we expected
        np.testing.assert_allclose(gts.data, new_gts.data)
        assert gts.samples == new_gts.samples
        for i in range(len(new_gts.variants)):
            for col in ("chrom", "pos", "id", "alleles"):
                assert gts.variants[col][i] == new_gts.variants[col][i]

        # clean up afterwards: delete the files we created
        fname.with_suffix(".psam").unlink()
        fname.with_suffix(".pvar").unlink()
        fname.unlink()

    def test_write_multiallelic(self):
        # Create fake multi-allelic genotype data to write to the PLINK file
        gts_multiallelic = self._get_fake_genotypes_multiallelic()
        #  Name of the PLINK file where the data will be written
        fname = DATADIR / "test_write_multiallelic.pgen"

        # Write multiallelic genotype data to the PLINK file
        gts_multiallelic.fname = fname
        gts_multiallelic.write()

        # Read the data from the newly created PLINK file
        new_gts = GenotypesPLINK(fname)
        new_gts.read()
        new_gts.check_phase()

        # Verify that the uploaded data matches the original data
        np.testing.assert_allclose(gts_multiallelic.data, new_gts.data)
        assert gts_multiallelic.samples == new_gts.samples
        for i in range(len(new_gts.variants)):
            for col in ("chrom", "pos", "id", "alleles"):
                assert gts_multiallelic.variants[col][i] == new_gts.variants[col][i]

        # clean the files created after the test
        fname.with_suffix(".psam").unlink()
        fname.with_suffix(".pvar").unlink()
        fname.unlink()

    def test_write_multiallelic_tr(self):
        # Create fake multi-allelic genotype data to write to the PLINK file
        gts_multiallelic = self._get_fake_genotypes_multiallelic_tr()
        #  Name of the PLINK file where the data will be written
        fname = DATADIR / "test_write_multiallelic_tr.pgen"

        # Write multiallelic genotype data to the PLINK file
        gts_multiallelic.fname = fname
        gts_multiallelic.write()

        # Read the data from the newly created PLINK file
        new_gts = GenotypesPLINK(fname)
        new_gts.read()

        # Verify that the uploaded data matches the original data
        np.testing.assert_allclose(gts_multiallelic.data, new_gts.data)
        assert gts_multiallelic.samples == new_gts.samples
        for i in range(len(new_gts.variants)):
            for col in ("chrom", "pos", "id", "alleles"):
                assert gts_multiallelic.variants[col][i] == new_gts.variants[col][i]

        # clean the files created after the test
        fname.with_suffix(".psam").unlink()
        fname.with_suffix(".pvar").unlink()
        fname.unlink()


class TestGenotypesPLINKTR:
    def _get_fake_genotypes_multiallelic(self):
        # define the code for an empty allele
        na = np.iinfo(np.uint8).max

        # fill out the class with the same variants and samples as GenotypesPLINK
        gts_plink = TestGenotypesPLINK()._get_fake_genotypes_multiallelic_tr()
        gts = GenotypesPLINKTR(gts_plink.fname)
        gts.variants = gts_plink.variants
        gts.samples = gts_plink.samples

        # fill out the data array with the repeat unit counts
        data = np.array(
            [
                [[1, 2], [3, 4], [5, 6], [7, 8], [9, 9]],
                [[3, 1], [3, 1], [1, 1], [1, 1], [3, 3]],
                [[3, 3], [3, 3], [3, 3], [3, 3], [3, 3]],
                [[11, 11], [na, na], [7, 2], [3, 4], [11, na]],
                [[5, na], [na, na], [3, na], [na, na], [na, na]],
            ],
            dtype=np.uint8,
        ).transpose((1, 0, 2))

        # append the phasing info, transpose, and clean up
        phasing = gts_plink.data[:, :, 2][:, :, np.newaxis]
        gts.data = np.concatenate((data, phasing), axis=-1)

        # handle half-calls and chrX according to the flag:
        # --vcf-half-call m
        gts.data[4, 3, 0] = na
        gts.data[:, 4, 1] = gts.data[:, 4, 0]
        # TODO: figure out whether this is the correct way to handle half-calls and try
        # to see if we can change chrX representations to properly simulate from it

        return gts

    def test_iter(self):
        # Get the expected data
        expected = self._get_fake_genotypes_multiallelic()
        gts = GenotypesPLINKTR(DATADIR / "simple-tr.pgen")
        # Check that everything matches what we expected
        for idx, line in enumerate(gts):
            np.testing.assert_allclose(line.data[:, :3], expected.data[:, idx])
            for col in ("chrom", "pos", "id"):
                assert line.variants[col] == expected.variants[col][idx]
            assert (
                line.variants["alleles"].tolist() == expected.variants["alleles"][idx]
            )

        # Check samples
        assert gts.samples == expected.samples

    def test_read(self):
        expected_alleles = self._get_fake_genotypes_multiallelic().data
        gts = GenotypesPLINKTR(DATADIR / "simple-tr.pgen")
        gts.read()
        # check genotypes
        np.testing.assert_allclose(expected_alleles, gts.data)


class TestPhenotypes:
    def _get_expected_phenotypes(self):
        # create a phen matrix with shape: samples x phenotypes
        expected = np.array(
            [
                [1, 1, 2, 2, 0],
                [3, 6, 8, 1, 4],
            ],
            dtype="float64",
        ).T
        return expected

    def _get_fake_phenotypes(self):
        gts = Phenotypes(fname=None)
        gts.data = self._get_expected_phenotypes()
        gts.samples = ("HG00096", "HG00097", "HG00099", "HG00100", "HG00101")
        gts.names = ("height", "bmi")
        return gts

    def test_load_phenotypes(self, caplog):
        expected_phen = self._get_fake_phenotypes()
        expected = expected_phen.data

        # can we load the data from the phenotype file?
        phens = Phenotypes(DATADIR / "simple.pheno")
        phens.read()
        np.testing.assert_allclose(phens.data, expected)
        assert phens.samples == expected_phen.samples

        # try loading the data again - it should warn b/c we've already done it
        caplog.clear()
        phens.read()
        assert len(caplog.records) > 0 and caplog.records[0].levelname == "WARNING"

    def test_check_missing(self, caplog):
        phens = Phenotypes(DATADIR / "simple.na.pheno")
        caplog.clear()
        phens.read()

        # we should have logged two ERRORS respectively for the 'NA' and 'na' values
        assert len(caplog.records) > 0
        assert sum(msg.levelname == "ERROR" for msg in caplog.records) == 2
        # any samples with 'NA' or 'na' should have been dropped
        assert len(phens.data) == 3
        assert len(phens.samples) == 3

        with pytest.raises(ValueError) as info:
            phens.check_missing()
        assert "HG00097" in str(info.value) and "bmi" in str(info.value)

        # now: check that it works when we discard the sample
        phens.check_missing(discard_also=True)
        assert len(phens.data) == 2
        assert len(phens.samples) == 2

    def test_load_phenotypes_iterate(self):
        expected_phen = self._get_fake_phenotypes()
        expected = expected_phen.data
        samples = expected_phen.samples

        # can we load the data from the phenotype file?
        phens = Phenotypes(DATADIR / "simple.pheno")
        for idx, line in enumerate(phens):
            np.testing.assert_allclose(line.data, expected[idx])
            assert line.samples == samples[idx]

    def test_load_phenotypes_subset(self):
        expected_phen = self._get_fake_phenotypes()
        expected = expected_phen.data
        samples = ["HG00097", "HG00100"]

        # subset for just the samples we want
        expected = expected[[1, 3]]

        # can we load the data from the phenotype file?
        phens = Phenotypes(DATADIR / "simple.pheno")
        phens.read(samples=set(samples))
        np.testing.assert_allclose(phens.data, expected)
        assert phens.samples == tuple(samples)

    def test_standardize(self):
        expected_phen = self._get_fake_phenotypes()
        exp_data = expected_phen.data

        exp_data = (exp_data - np.mean(exp_data, axis=0)) / np.std(exp_data, axis=0)
        expected_phen.standardize()
        np.testing.assert_allclose(expected_phen.data, exp_data)

        # also test case where the stdev is 0
        zero_phen = np.zeros((exp_data.shape[0], 1), dtype=exp_data.dtype)
        exp_data = np.concatenate((exp_data, zero_phen), axis=1)
        expected_phen.data = np.concatenate(
            (expected_phen.data, zero_phen + 5),
            dtype=exp_data.dtype,
            axis=1,
        )

        expected_phen.standardize()
        np.testing.assert_allclose(expected_phen.data, exp_data)

    def test_write(self):
        exp_phen = self._get_fake_phenotypes()

        # first, we write the data
        exp_phen.fname = DATADIR / "test.pheno"
        exp_phen.write()

        # now, let's load the data and check that it's what we wrote
        result = Phenotypes(exp_phen.fname)
        result.read()
        np.testing.assert_allclose(exp_phen.data, result.data, rtol=0, atol=0)
        assert exp_phen.names == result.names
        assert exp_phen.samples == result.samples

        # try to standardize the data and see if it's still close
        # to validate that our code can handle phenotypes with arbitrary precision
        exp_phen.standardize()
        exp_phen.write()

        # now, let's load the data and check that it's what we wrote
        result = Phenotypes(exp_phen.fname)
        result.read()
        np.testing.assert_allclose(exp_phen.data, result.data, rtol=0, atol=0)
        assert exp_phen.names == result.names
        assert exp_phen.samples == result.samples

        # try to make some of the data negative/positive and check that, as well
        exp_phen.data[[1, 4], [0, 1]] = -439.58
        exp_phen.data[[2, 3], [0, 1]] = 439.58
        exp_phen.write()

        # now, let's load the data and check that it's what we wrote
        result = Phenotypes(exp_phen.fname)
        result.read()
        np.testing.assert_allclose(exp_phen.data, result.data, rtol=0, atol=0)
        assert exp_phen.names == result.names
        assert exp_phen.samples == result.samples

        # let's just try to create random, fake phenotypes
        shape = exp_phen.data.shape[0] * exp_phen.data.shape[1]
        exp_phen.data = np.random.normal(size=shape).reshape(exp_phen.data.shape)
        exp_phen.write()

        # now, let's load the data and check that it's what we wrote
        result = Phenotypes(exp_phen.fname)
        result.read()
        np.testing.assert_allclose(exp_phen.data, result.data, rtol=0, atol=0)
        assert exp_phen.names == result.names
        assert exp_phen.samples == result.samples

        # let's clean up after ourselves and delete the file
        exp_phen.fname.unlink()

    def test_write_casecontrol(self):
        exp_phen = self._get_fake_phenotypes()

        # first, we write the data
        exp_phen.fname = DATADIR / "test_cc.pheno"

        # also check that plain integers like 0 and 1 get written in a way that PLINK2
        # will recognize them as case/control
        shape = exp_phen.data.shape[0] * exp_phen.data.shape[1]
        exp_phen.data = (
            np.random.choice([True, False], size=shape)
            .reshape(exp_phen.data.shape)
            .astype(exp_phen.data.dtype)
        )
        exp_phen.write()

        # now, let's load the data and check that it's what we wrote
        result = Phenotypes(exp_phen.fname)
        result.read()
        np.testing.assert_allclose(exp_phen.data, result.data, rtol=0, atol=0)
        assert exp_phen.names == result.names
        assert exp_phen.samples == result.samples

        # let's also load it using a simple python reader
        with open(exp_phen.fname, "r") as result_file:
            for res in result_file.read().splitlines()[1:]:
                height, bmi = res.split("\t")[1:]
                assert len(height) <= 2
                assert len(bmi) <= 2
                assert int(float(height)) <= 1
                assert int(float(bmi)) <= 1

        # let's clean up after ourselves and delete the file
        exp_phen.fname.unlink()

    def test_append_phenotype(self):
        expected1 = self._get_fake_phenotypes()
        expected2 = self._get_fake_phenotypes()

        # add a phenotype called "test" to the set of phenotypes
        new_phen = np.array([5, 2, 8, 0, 3], dtype=expected2.data.dtype)
        expected2.append("test", new_phen)

        # did it get added correctly?
        assert expected1.data.shape[1] == expected2.data.shape[1] - 1
        assert len(expected1.names) == len(expected2.names) - 1
        assert expected2.names[2] == "test"
        assert len(expected1.samples) == len(expected2.samples)
        np.testing.assert_allclose(expected1.data, expected2.data[:, :2])
        np.testing.assert_allclose(expected2.data[:, 2], new_phen)

    def test_subset_phenotypes(self):
        pts = self._get_fake_phenotypes()

        # subset to just the samples we want
        expected_data = pts.data[:3]
        expected_names = pts.names
        samples = ("HG00096", "HG00097", "HG00099")
        pts_sub = pts.subset(samples=samples)
        assert pts_sub.samples == samples
        np.testing.assert_allclose(pts_sub.data, expected_data)
        assert np.array_equal(pts_sub.names, expected_names)

        # subset to just the names we want
        expected_data = pts.data[:, [1]]
        assert len(expected_data.shape) == 2
        expected_names = (pts.names[1],)
        names = ("bmi",)
        pts_sub = pts.subset(names=names)
        assert pts_sub.samples == pts.samples
        np.testing.assert_allclose(pts_sub.data, expected_data)
        assert np.array_equal(pts_sub.names, expected_names)

        # subset both: samples and names
        expected_data = pts.data[[3, 4], [1]]
        expected_data = expected_data[:, np.newaxis]
        assert len(expected_data.shape) == 2
        expected_names = (pts.names[1],)
        samples = ("HG00100", "HG00101")
        names = ("bmi",)
        pts_sub = pts.subset(samples=samples, names=names)
        assert pts_sub.samples == samples
        np.testing.assert_allclose(pts_sub.data, expected_data)
        assert np.array_equal(pts_sub.names, expected_names)

        # also try in-place
        expected_data = pts.data[[3, 4], [1]]
        expected_data = expected_data[:, np.newaxis]
        assert len(expected_data.shape) == 2
        expected_names = (pts.names[1],)
        samples = ("HG00100", "HG00101")
        names = ("bmi",)
        pts.subset(samples=samples, names=names, inplace=True)
        assert pts_sub.samples == samples
        np.testing.assert_allclose(pts_sub.data, expected_data)
        assert np.array_equal(pts_sub.names, expected_names)


class TestCovariates:
    def _get_expected_covariates(self):
        # create a covar matrix with shape: samples x covariates
        expected = np.array([(0, 4), (1, 20), (1, 33), (0, 15), (0, 78)])
        return expected

    def _get_fake_covariates(self):
        gts = Covariates(fname=None)
        gts.data = self._get_expected_covariates()
        gts.samples = ("HG00096", "HG00097", "HG00099", "HG00100", "HG00101")
        gts.names = ("sex", "age")
        return gts

    def test_load_covariates(self, caplog):
        # create a covariate vector with shape: num_samples x num_covars
        expected_covar = self._get_fake_covariates()
        expected = expected_covar.data
        samples = expected_covar.samples

        # can we load the data from the covariates file?
        covars = Covariates(DATADIR / "simple.covar")
        covars.read()
        np.testing.assert_allclose(covars.data, expected)
        assert covars.samples == samples
        assert covars.names == ("sex", "age")

        # try loading the data again - it should warn b/c we've already done it
        caplog.clear()
        covars.read()
        assert len(caplog.records) > 0 and caplog.records[0].levelname == "WARNING"

    def test_load_covariates_iterate(self):
        # create a covariate vector with shape: num_samples x num_covars
        expected_covar = self._get_fake_covariates()
        expected = expected_covar.data
        samples = expected_covar.samples

        # can we load the data from the covariates file?
        covars = Covariates(DATADIR / "simple.covar")
        for idx, line in enumerate(covars):
            np.testing.assert_allclose(line.data, expected[idx])
            assert line.samples == samples[idx]
        assert covars.names == ("sex", "age")

    def test_load_covariates_subset(self):
        # create a covariate vector with shape: num_samples x num_covars
        expected_covar = self._get_fake_covariates()
        expected = expected_covar.data
        samples = ["HG00097", "HG00100"]

        # subset for just the samples we want
        expected = expected[[1, 3]]

        # can we load the data from the covariate file?
        covars = Covariates(DATADIR / "simple.covar")
        covars.read(samples=set(samples))
        np.testing.assert_allclose(covars.data, expected)
        assert covars.samples == tuple(samples)

    def test_subset_covariates(self):
        cvs = self._get_fake_covariates()

        # subset to just the samples we want
        expected_data = cvs.data[:3]
        expected_names = cvs.names
        samples = ("HG00096", "HG00097", "HG00099")
        cvs_sub = cvs.subset(samples=samples)
        assert cvs_sub.samples == samples
        np.testing.assert_allclose(cvs_sub.data, expected_data)
        assert np.array_equal(cvs_sub.names, expected_names)

        # subset to just the names we want
        expected_data = cvs.data[:, [1]]
        assert len(expected_data.shape) == 2
        expected_names = (cvs.names[1],)
        names = ("age",)
        cvs_sub = cvs.subset(names=names)
        assert cvs_sub.samples == cvs.samples
        np.testing.assert_allclose(cvs_sub.data, expected_data)
        assert np.array_equal(cvs_sub.names, expected_names)

        # subset both: samples and names
        expected_data = cvs.data[[3, 4], [1]]
        expected_data = expected_data[:, np.newaxis]
        assert len(expected_data.shape) == 2
        expected_names = (cvs.names[1],)
        samples = ("HG00100", "HG00101")
        names = ("age",)
        cvs_sub = cvs.subset(samples=samples, names=names)
        assert cvs_sub.samples == samples
        np.testing.assert_allclose(cvs_sub.data, expected_data)
        assert np.array_equal(cvs_sub.names, expected_names)


class TestHaplotypes:
    def _basic_haps(self):
        # what do we expect to see from the basic.hap file?
        expected = {
            "chr21.q.3365*1": Haplotype("21", 26928472, 26941960, "chr21.q.3365*1"),
            "chr21.q.3365*10": Haplotype("21", 26938989, 26941960, "chr21.q.3365*10"),
            "chr21.q.3365*11": Haplotype("21", 26938353, 26938989, "chr21.q.3365*11"),
            "21_26941880_STR": Repeat("21", 26941880, 26941900, "21_26941880_STR"),
            "21_26938989_STR": Repeat("21", 26939000, 26939010, "21_26938989_STR"),
            "21_26938353_STR": Repeat("21", 26938353, 26938400, "21_26938353_STR"),
        }
        expected["chr21.q.3365*1"].variants = (
            Variant(26928472, 26928472, "21_26928472_C_A", "C"),
            Variant(26938353, 26938353, "21_26938353_T_C", "T"),
            Variant(26940815, 26940815, "21_26940815_T_C", "C"),
            Variant(26941960, 26941960, "21_26941960_A_G", "G"),
        )
        expected["chr21.q.3365*10"].variants = (
            Variant(26938989, 26938989, "21_26938989_G_A", "A"),
            Variant(26940815, 26940815, "21_26940815_T_C", "T"),
            Variant(26941960, 26941960, "21_26941960_A_G", "A"),
        )
        expected["chr21.q.3365*11"].variants = (
            Variant(26938353, 26938353, "21_26938353_T_C", "T"),
            Variant(26938989, 26938989, "21_26938989_G_A", "A"),
        )
        return expected

    def _basic_unordered_first_field_haps(self):
        # what do we expect to see from the basic.hap file?
        expected = {
            "21.q.3365*1": Haplotype("21", 26928472, 26941960, "21.q.3365*1"),
            "22.q.3365*10": Haplotype("22", 26938989, 26941960, "22.q.3365*10"),
            "22.q.3365*11": Haplotype("22", 26938353, 26938989, "22.q.3365*11"),
            "21_26938353_STR": Repeat("21", 26938353, 26938400, "21_26938353_STR"),
        }
        expected["21.q.3365*1"].variants = (
            Variant(26928472, 26928472, "21_26928472_C_A", "C"),
            Variant(26938353, 26938353, "21_26938353_T_C", "T"),
            Variant(26940815, 26940815, "21_26940815_T_C", "C"),
            Variant(26941960, 26941960, "21_26941960_A_G", "G"),
        )
        expected["22.q.3365*10"].variants = (
            Variant(26938989, 26938989, "21_26938989_G_A", "A"),
            Variant(26940815, 26940815, "21_26940815_T_C", "T"),
            Variant(26941960, 26941960, "21_26941960_A_G", "A"),
        )
        expected["22.q.3365*11"].variants = (
            Variant(26938353, 26938353, "21_26938353_T_C", "T"),
            Variant(26938989, 26938989, "21_26938989_G_A", "A"),
        )
        return expected

    def _get_dummy_haps(self):
        # create three haplotypes
        haplotypes = {
            "H1": Haplotype(chrom="1", start=10114, end=8, id="H1"),
            "H2": Haplotype(chrom="1", start=10114, end=10119, id="H2"),
            "H3": Haplotype(chrom="1", start=10116, end=10119, id="H3"),
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
        haps = Haplotypes(fname=None)
        haps.data = haplotypes
        return haps

    def _get_dummy_haps_multiallelic(self):
        haps = self._get_dummy_haps()
        h1_vars = list(haps.data["H1"].variants)
        h1_vars[1] = Variant(start=10116, end=10117, id="1:10116:A:G", allele="T")
        haps.data["H1"].variants = tuple(h1_vars)
        h3_vars = list(haps.data["H3"].variants)
        h3_vars[1] = Variant(start=10122, end=10123, id="1:10122:A:G", allele="C")
        haps.data["H3"].variants = tuple(h3_vars)
        return haps

    def test_load(self):
        expected = self._basic_haps()

        # can we load this data from the hap file?
        haps = Haplotypes.load(DATADIR / "basic.hap")
        assert expected == haps.data

        # also check the indexed file
        # it should be the same
        haps = Haplotypes.load(DATADIR / "basic.hap.gz")
        assert expected == haps.data

    def test_load_no_header(self):
        expected = self._basic_haps()

        # what if we remove the header line, can we still load it?
        # let's try to make a version without the header
        no_header_file = DATADIR / "basic_no_header.hap"
        with open(DATADIR / "basic.hap", "r") as fr, open(no_header_file, "w") as fw:
            fw.writelines([ln for ln in fr.readlines() if not ln.startswith("#\t")])

        haps = Haplotypes.load(no_header_file)
        assert expected == haps.data

        no_header_file.unlink()

    def test_iterate(self):
        exp_full = self._basic_haps()

        exp_single_hap = [exp_full["chr21.q.3365*1"]]
        exp_single_hap += exp_single_hap[0].variants
        exp_single_hap2 = [exp_full["chr21.q.3365*11"]]
        exp_single_hap2 += exp_single_hap2[0].variants
        exp_single_hap2 += [exp_full["21_26938353_STR"]]
        exp_single_hap3 = [exp_full["21_26938353_STR"]]

        expected = [hp for hp in exp_full.values() if isinstance(hp, Haplotype)]
        all_haps = [hp for hp in exp_full.values() if isinstance(hp, Haplotype)]
        repeats = [hap for hap in exp_full.values() if isinstance(hap, Repeat)]
        for i, hap in enumerate(tuple(all_haps)):
            if isinstance(hap, Haplotype):
                expected += hap.variants
                expected += [repeats[i]]
                hap.variants = ()

        # can we load this data from the hap file?
        haps = Haplotypes(DATADIR / "basic.hap")

        # currently read like repeats after haplotypes which isnt
        for exp_hap, line in zip(expected, haps):
            assert exp_hap == line

        # also check whether it works when we pass function params
        haps = Haplotypes(DATADIR / "basic.hap.gz")
        haps_iter = list(haps.__iter__(region="21:26928472-26941960"))
        assert len(haps_iter) == len(expected)
        assert all(line in expected for line in haps_iter)

        haps_iter = list(haps.__iter__(region="21"))
        assert len(haps_iter) == len(expected)
        assert all(line in expected for line in haps_iter)

        haps_iter = list(haps.__iter__(region="21:"))
        assert len(haps_iter) == len(expected)
        assert all(line in expected for line in haps_iter)

        haps_iter = list(haps.__iter__(region="21:26928472-"))
        assert len(haps_iter) == len(expected)
        assert all(line in expected for line in haps_iter)

        haps_iter = list(haps.__iter__(region="21:26928472-26938989"))
        assert len(haps_iter) == len(exp_single_hap2)
        assert all(line in exp_single_hap2 for line in haps_iter)

        # also, try adding the hap ID
        i = haps.__iter__(region="21:26928472-26941960", haplotypes={"chr21.q.3365*1"})
        for exp_hap, line in zip(exp_single_hap, i):
            assert exp_hap == line

        # also, try adding the hap ID
        haps_iter = haps.__iter__(haplotypes={"chr21.q.3365*1"})
        for exp_hap, line in zip(exp_single_hap, haps_iter):
            assert exp_hap == line

        # try adding repeat ID
        i = haps.__iter__(region="21:26928472-26938989", haplotypes={"21_26938353_STR"})
        for exp_hap, line in zip(exp_single_hap3, i):
            assert exp_hap == line

    def test_iterate_unordered_first_field(self):
        exp_full = self._basic_unordered_first_field_haps()

        exp_single_hap = [exp_full["21.q.3365*1"]]
        exp_single_hap += exp_single_hap[0].variants
        exp_single_hap += [exp_full["21_26938353_STR"]]
        exp_single_hap2 = [exp_full["22.q.3365*11"]]
        exp_single_hap2 += exp_single_hap2[0].variants

        expected = [hap for hap in exp_full.values()]
        for hap in tuple(expected):
            if isinstance(hap, Haplotype):
                expected += hap.variants
                hap.variants = ()

        # also check whether it works when we pass function params
        haps = Haplotypes(DATADIR / "unordered_first_field.hap.gz")
        haps_iter = list(haps.__iter__(region="21:26928472-26941960"))
        assert len(haps_iter) == len(exp_single_hap)
        assert all(line in exp_single_hap for line in haps_iter)

        haps_iter = list(haps.__iter__(region="21"))
        assert len(haps_iter) == len(exp_single_hap)
        assert all(line in exp_single_hap for line in haps_iter)

        haps_iter = list(haps.__iter__(region="21:"))
        assert len(haps_iter) == len(exp_single_hap)
        assert all(line in exp_single_hap for line in haps_iter)

        haps_iter = list(haps.__iter__(region="21:26928472-"))
        assert len(haps_iter) == len(exp_single_hap)
        assert all(line in exp_single_hap for line in haps_iter)

        haps_iter = list(haps.__iter__(region="22:26928472-26938989"))
        assert len(haps_iter) == len(exp_single_hap2)
        assert all(line in exp_single_hap2 for line in haps_iter)

        # also, try adding the hap ID
        i = haps.__iter__(region="21:26928472-26941960", haplotypes={"21.q.3365*1"})
        for exp_hap, line in zip(exp_single_hap, i):
            assert exp_hap == line

        # also, try adding the hap ID
        haps_iter = haps.__iter__(haplotypes={"22.q.3365*11"})
        for exp_hap, line in zip(exp_single_hap2, haps_iter):
            assert exp_hap == line

    def test_read_subset(self):
        expected = {}
        expected["chr21.q.3365*1"] = self._basic_haps()["chr21.q.3365*1"]
        tr_expected = {}
        tr_expected["21_26938989_STR"] = self._basic_haps()["21_26938989_STR"]

        haps = Haplotypes(DATADIR / "basic.hap")
        # this shouldn't fail anymore as of version 0.1.0
        haps.read(haplotypes={"chr21.q.3365*1"})
        assert expected == haps.data

        haps = Haplotypes(DATADIR / "basic.hap")
        haps.read(haplotypes={"21_26938989_STR"})
        assert tr_expected == haps.data

        haps = Haplotypes(DATADIR / "basic.hap.gz")
        haps.read(haplotypes={"chr21.q.3365*1"})
        assert expected == haps.data

        haps = Haplotypes(DATADIR / "basic.hap.gz")
        haps.read(haplotypes={"21_26938989_STR"})
        assert tr_expected == haps.data

        haps = Haplotypes(DATADIR / "basic.hap.gz")
        haps.read(region="21:26928472-26941960", haplotypes={"chr21.q.3365*1"})
        assert expected == haps.data

        haps = Haplotypes(DATADIR / "basic.hap.gz")
        haps.read(region="21:26928472-26941960", haplotypes={"21_26938989_STR"})
        assert tr_expected == haps.data

        # check that haplotypes that overlap but don't fit perfectly are excluded!
        haps = Haplotypes(DATADIR / "basic.hap.gz")
        haps.read(region="21:26928473-26941960", haplotypes={"chr21.q.3365*1"})
        assert {} == haps.data
        haps = Haplotypes(DATADIR / "basic.hap.gz")
        haps.read(region="21:26928472-26941959", haplotypes={"chr21.q.3365*1"})
        assert {} == haps.data

        # test repeats overlap but don't fit perfectly
        haps = Haplotypes(DATADIR / "basic.hap.gz")
        haps.read(region="21:26938273-26938360", haplotypes={"21_26938353_STR"})
        assert {} == haps.data
        haps = Haplotypes(DATADIR / "basic.hap.gz")
        haps.read(region="21:26941890-26941959", haplotypes={"21_26941880_STR"})
        assert {} == haps.data

        expected = self._basic_haps()

        haps = Haplotypes(DATADIR / "basic.hap.gz")
        haps.read(region="21:26928472-26941960")
        assert len(expected) == len(haps.data)
        assert expected == haps.data

    def test_subset(self):
        expected = Haplotypes(DATADIR / "basic.hap")
        expected.read(haplotypes={"chr21.q.3365*1"})

        haps = Haplotypes(DATADIR / "basic.hap")
        haps.read()
        haps = haps.subset(haplotypes=("chr21.q.3365*1",))

        assert len(expected.data) == len(haps.data)
        assert expected.data == haps.data

    def test_read_extras(self):
        # what do we expect to see from the simphenotype.hap file?
        expected = {
            "chr21.q.3365*1": HaptoolsHaplotype(
                "21", 26928472, 26941960, "chr21.q.3365*1", 0.73
            ),
            "chr21.q.3365*10": HaptoolsHaplotype(
                "21", 26938989, 26941960, "chr21.q.3365*10", 0.30
            ),
            "chr21.q.3365*11": HaptoolsHaplotype(
                "21", 26938353, 26938989, "chr21.q.3365*11", 0.49
            ),
            "21_26938353_STR": HaptoolsRepeat(
                "21", 26938353, 26938400, "21_26938353_STR", 0.45
            ),
        }
        for hap_id, hap in self._basic_haps().items():
            if isinstance(hap, Haplotype):
                expected[hap_id].variants = hap.variants

        # can we load this data from the hap file?
        haps = Haplotypes(
            DATADIR / "simphenotype.hap",
            haplotype=HaptoolsHaplotype,
            repeat=HaptoolsRepeat,
        )
        haps.read()
        assert expected == haps.data

    def test_read_extras_large(self):
        """
        try reading a large-ish file
        """
        haps = Haplotypes(DATADIR / "example.hap.gz", HaptoolsHaplotype)
        haps.read()
        assert len(self._basic_haps().keys() & haps.data.keys())

    def test_write(self):
        haps = Haplotypes(DATADIR / "test.hap")
        haps.data = self._basic_haps()
        haps.write()

        haps.data = None
        haps.read()
        assert self._basic_haps() == haps.data

        # remove the file
        haps.fname.unlink()

    def _get_writable_haplotypes(self):
        expected = {
            "chr21.q.3365*1": HaptoolsHaplotype(
                "21", 26928472, 26941960, "chr21.q.3365*1", 0.73
            ),
            "chr21.q.3365*10": HaptoolsHaplotype(
                "21", 26938989, 26941960, "chr21.q.3365*10", 0.30
            ),
            "chr21.q.3365*11": HaptoolsHaplotype(
                "21", 26938353, 26938989, "chr21.q.3365*11", 0.49
            ),
            "21_26938353_STR": HaptoolsRepeat(
                "21", 26938353, 26938400, "21_26938353_STR", 0.45
            ),
        }
        for hap_id, hap in self._basic_haps().items():
            if isinstance(hap, Haplotype):
                expected[hap_id].variants = hap.variants
        return expected

    def test_write_extras(self):
        # what do we expect to see from the test.hap file?
        expected = self._get_writable_haplotypes()

        haps = Haplotypes(
            DATADIR / "test.hap",
            haplotype=HaptoolsHaplotype,
            repeat=HaptoolsRepeat,
        )
        haps.data = expected
        haps.write()

        haps.data = None
        haps.read()
        assert expected == haps.data

        # remove the file
        haps.fname.unlink()

    def test_write_plus_extra(self):
        @dataclass
        class HaplotypePlusExtra(HaptoolsHaplotype):
            """
            A haplotype with an additional, unnecessary extra field

            Properties and functions are shared with the HaptoolsHaplotype object
            """

            score: float
            _extras: tuple = field(
                repr=False,
                init=False,
                default=(
                    Extra("score", ".2f", "Score for a thing"),
                    Extra("beta", ".2f", "Effect size in linear model"),
                ),
            )

        @dataclass
        class RepeatPlusExtra(HaptoolsRepeat):
            """
            A repeat with an additional, unnecessary extra field

            Properties and functions are shared with the HaptoolsRepeat object
            """

            score: float
            _extras: tuple = field(
                repr=False,
                init=False,
                default=(
                    Extra("score", ".2f", "Score for a thing"),
                    Extra("beta", ".2f", "Effect size in linear model"),
                ),
            )

        # what do we want to write to the test.hap file?
        expected = {
            "chr21.q.3365*1": HaplotypePlusExtra(
                "21", 26928472, 26941960, "chr21.q.3365*1", 0.73, 0.40
            ),
            "chr21.q.3365*10": HaplotypePlusExtra(
                "21", 26938989, 26941960, "chr21.q.3365*10", 0.30, 0.28
            ),
            "chr21.q.3365*11": HaplotypePlusExtra(
                "21", 26938353, 26938989, "chr21.q.3365*11", 0.49, 0.84
            ),
            "21_26938353_STR": RepeatPlusExtra(
                "21", 26938353, 26938400, "21_26938353_STR", 0.45, 0.0
            ),
        }
        for hap_id, hap in self._basic_haps().items():
            if isinstance(hap, Haplotype):
                expected[hap_id].variants = hap.variants

        haps = Haplotypes(
            DATADIR / "test.hap",
            haplotype=HaplotypePlusExtra,
            repeat=RepeatPlusExtra,
        )
        haps.data = expected
        haps.write()

        haps = Haplotypes(
            DATADIR / "test.hap",
            haplotype=HaptoolsHaplotype,
            repeat=HaptoolsRepeat,
        )
        haps.read()
        assert haps.data == self._get_writable_haplotypes()

        # remove the file
        haps.fname.unlink()

    def test_hap_transform(self):
        expected = np.array(
            [
                [0, 1],
                [0, 1],
                [1, 1],
                [1, 1],
                [0, 0],
            ],
            dtype=np.uint8,
        )

        hap = list(self._get_dummy_haps().data.values())[0]
        gens = TestGenotypesVCF()._get_fake_genotypes_refalt()
        hap_gt = hap.transform(gens)
        np.testing.assert_allclose(hap_gt, expected)

    def test_hap_transform_multiallelic(self):
        expected = np.array(
            [
                [0, 0],
                [0, 0],
                [1, 0],
                [0, 0],
                [0, 1],
            ],
            dtype=np.uint8,
        )

        hap = list(self._get_dummy_haps_multiallelic().data.values())[0]
        gens = TestGenotypesVCF()._get_fake_genotypes_multiallelic()
        hap_gt = hap.transform(gens)
        np.testing.assert_allclose(hap_gt, expected)

    def test_haps_transform(self, return_also=False):
        expected = np.array(
            [
                [[0, 1], [0, 0], [0, 0]],
                [[0, 1], [0, 0], [1, 0]],
                [[1, 0], [0, 1], [0, 0]],
                [[1, 1], [0, 0], [0, 0]],
                [[0, 0], [0, 1], [1, 0]],
            ],
            dtype=np.uint8,
        )

        haps = self._get_dummy_haps()
        gens = TestGenotypesVCF()._get_fake_genotypes_refalt()
        gens.data[[2, 4], 0, 1] = 1
        gens.data[[1, 4], 2, 0] = 1
        hap_gt = GenotypesVCF(fname=None)
        haps.transform(gens, hap_gt)
        np.testing.assert_allclose(hap_gt.data, expected)

        if return_also:
            return hap_gt

    def test_haps_transform_multiallelic(self, return_also=False):
        expected = np.array(
            [
                [[0, 0], [0, 0], [0, 0]],
                [[0, 0], [0, 0], [1, 0]],
                [[1, 0], [0, 1], [0, 0]],
                [[0, 0], [0, 0], [0, 0]],
                [[0, 0], [0, 1], [1, 0]],
            ],
            dtype=np.uint8,
        )

        haps = self._get_dummy_haps_multiallelic()
        gens = TestGenotypesVCF()._get_fake_genotypes_multiallelic()
        gens.data[[2, 4], 0, 1] = 1
        gens.data[[1, 4], 2, 0] = 1
        gens.data[[1, 4], 3, 0] = 2
        hap_gt = GenotypesVCF(fname=None)
        haps.transform(gens, hap_gt)
        np.testing.assert_allclose(hap_gt.data, expected)

        if return_also:
            return hap_gt

    def test_hap_gt_write(self):
        fname = DATADIR / "simple_haps.vcf"

        hap_gt = self.test_haps_transform(return_also=True)
        hap_gt.fname = fname
        expected_data = hap_gt.data
        expected_samples = hap_gt.samples
        hap_gt.write()

        hap_gt.data = None
        hap_gt.read()
        hap_gt.check_phase()
        np.testing.assert_allclose(hap_gt.data, expected_data)
        assert hap_gt.samples == expected_samples
        assert len(hap_gt.variants) == hap_gt.data.shape[1]

        # remove the file
        os.remove(str(fname))

    def test_lt_repeat_hap(self):
        # test comparison of haplotypes and repeats
        haps = self._get_dummy_haps()
        less_repeat = Repeat(chrom="1", start=10113, end=10118, id="U2")
        equal_repeat = Repeat(chrom="1", start=10114, end=10119, id="H2")
        greater_repeat = Repeat(chrom="2", start=10114, end=10119, id="H2")
        assert less_repeat < haps.data["H2"]
        assert not (equal_repeat < haps.data["H2"])
        assert not (greater_repeat < haps.data["H2"])

    def test_lt_haps(self):
        hap1 = Haplotype(chrom="A", start=3, end=1000, id="test1")
        hap2 = Haplotype(chrom="B", start=2, end=1000, id="test2")
        assert hap1 < hap2

    def test_gt_haps(self):
        hap1 = Haplotype(chrom="A", start=3, end=1000, id="test1")
        hap2 = Haplotype(chrom="A", start=2, end=1000, id="test2")
        assert hap1 > hap2

    def test_lt_var(self):
        var1 = Variant(start=1, end=1000, id="test1", allele="test1")
        var2 = Variant(start=1, end=1001, id="test2", allele="test2")
        assert var1 < var2

    def test_gt_var(self):
        var1 = Variant(start=7, end=1000, id="test1", allele="test1")
        var2 = Variant(start=1, end=1001, id="test2", allele="test2")
        assert var1 > var2

    def test_gt_equal_var(self):
        var1 = Variant(start=1, end=1000, id="test", allele="test")
        var2 = Variant(start=1, end=1000, id="test", allele="test")
        assert var1 >= var2

    def test_gt_equal_haps(self):
        hap1 = Haplotype(chrom="A", start=2, end=1000, id="test")
        hap2 = Haplotype(chrom="A", start=2, end=1000, id="test")
        assert hap1 >= hap2

    def test_lt_equal_var(self):
        var1 = Variant(start=1, end=1000, id="test1", allele="test1")
        var2 = Variant(start=1, end=1001, id="test2", allele="test2")
        assert var1 <= var2

    def test_lt_equal_haps(self):
        hap1 = Haplotype(chrom="A", start=2, end=1000, id="test")
        hap2 = Haplotype(chrom="A", start=2, end=1000, id="test")
        assert hap1 <= hap2

    def test_sort(self):
        test_hap1 = Haplotypes(DATADIR / "test_sort_unordered.hap")
        test_hap2 = Haplotypes(DATADIR / "test_sort_ordered.hap")
        test_hap1.read()
        test_hap1.sort()
        test_hap1.fname = Path("test_temp_sort.hap")
        test_hap1.write()

        with open(test_hap2.fname) as f1, open(test_hap1.fname) as f2:
            for line1, line2 in zip(f1, f2):
                assert line1 == line2

        test_hap1.fname.unlink()


class TestGenotypesVCF:
    def _get_fake_genotypes_refalt(self, with_phase=False):
        base_gts = TestGenotypes()._get_fake_genotypes()
        # copy all of the fields
        gts = GenotypesVCF(fname=None)
        gts.data = base_gts.data
        if with_phase:
            data_shape = (gts.data.shape[0], gts.data.shape[1], 1)
            # add phase info back
            gts.data = np.concatenate(
                (gts.data, np.ones(data_shape, dtype=gts.data.dtype)), axis=2
            )
        gts.samples = base_gts.samples
        base_dtype = {k: v[0] for k, v in base_gts.variants.dtype.fields.items()}
        ref_alt = [
            ("T", "C"),
            ("A", "G"),
            ("C", "A"),
            ("A", "G"),
        ]
        gts.variants = np.array(
            [tuple(rec) + (ref_alt[idx],) for idx, rec in enumerate(base_gts.variants)],
            dtype=(list(base_dtype.items()) + [("alleles", object)]),
        )
        return gts

    def _get_fake_genotypes_multiallelic(self, with_phase=False):
        gts = self._get_fake_genotypes_refalt(with_phase=with_phase)
        # replace the necessary properties
        gts.variants["alleles"] = np.array(
            [
                ("T", "C"),
                ("A", "G", "T"),
                ("C", "A"),
                ("A", "G", "C"),
            ],
            dtype=gts.variants["alleles"].dtype,
        )
        gts.data[[2, 4], 1, [0, 1]] = 2
        return gts

    def test_read_ref_alt(self):
        # simple.vcf
        expected = self._get_fake_genotypes_refalt()
        gts = GenotypesVCF(DATADIR / "simple.vcf")
        gts.read()
        for i, x in enumerate(expected.variants["alleles"]):
            assert gts.variants["alleles"][i] == x

        # example.vcf.gz
        gts = GenotypesVCF(DATADIR / "example.vcf.gz")
        gts.read()
        expected = np.array(
            [
                ("C", "A"),
                ("T", "C"),
                ("G", "A"),
                ("T", "C"),
                ("A", "G"),
                ("A", "G"),
                ("T", "G"),
                ("T", "A"),
                ("G", "A"),
                ("T", "C"),
                ("C", "G"),
                ("A", "G"),
                ("T", "C"),
                ("T", "C"),
                ("C", "T"),
            ],
            dtype=gts.variants["alleles"].dtype,
        )
        for i, x in enumerate(expected):
            assert gts.variants["alleles"][i] == tuple(x.tolist())

    def test_read_multiallelic(self):
        # simple-multiallelic.vcf
        expected = self._get_fake_genotypes_multiallelic(with_phase=True)

        gts = GenotypesVCF(DATADIR / "simple-multiallelic.vcf")
        gts.read()
        for i, x in enumerate(expected.variants["alleles"]):
            assert gts.variants["alleles"][i] == x
        np.testing.assert_allclose(gts.data, expected.data)

    def test_write_ref_alt(self, multiallelic=False):
        # strategy is to read in the file, write it, and then read again
        if multiallelic:
            expected = self._get_fake_genotypes_multiallelic()
            # read genotypes
            gts = GenotypesVCF(DATADIR / "simple-multiallelic.vcf")
        else:
            expected = self._get_fake_genotypes_refalt()
            # read genotypes
            gts = GenotypesVCF(DATADIR / "simple.vcf")
        gts.read()
        gts.check_phase()
        # write file to new file
        fname = DATADIR / "test.vcf"
        gts.fname = fname
        gts.write()
        # read again
        gts.read()
        gts.check_phase()
        # compare samples, data, variants, and ref/alt for equality
        assert gts.samples == expected.samples
        np.testing.assert_allclose(gts.data, expected.data)

        for i, x in enumerate(expected.variants):
            # dtype.names gives us names of columns in the array
            for col in expected.variants.dtype.names:
                # each row is a variant
                # index into col, i gets specific variant, x iterates thru
                assert gts.variants[col][i] == x[col]

        fname.unlink()

    def test_write_multiallelic(self):
        self.test_write_ref_alt(multiallelic=True)

    def test_write_phase(self, prephased=True):
        gts = self._get_fake_genotypes_refalt()

        fname = DATADIR / "test_write_phase.vcf"
        gts.fname = fname
        if prephased:
            gts._prephased = True
        else:
            # add phasing information back
            gts.data = np.dstack(
                (gts.data, np.ones(gts.data.shape[:2], dtype=np.uint8))
            )
            gts.data[:2, 1, 2] = 0
        gts.write()

        new_gts = GenotypesVCF(fname)
        if prephased:
            new_gts._prephased = True
        new_gts.read()

        # check that everything matches what we expected
        np.testing.assert_allclose(gts.data, new_gts.data)
        assert gts.samples == new_gts.samples
        for i in range(len(new_gts.variants)):
            for col in ("chrom", "pos", "id", "alleles"):
                assert gts.variants[col][i] == new_gts.variants[col][i]

        fname.unlink()

    def test_write_unphased(self):
        self.test_write_phase(prephased=False)

    def test_write_missing(self):
        gts = self._get_fake_genotypes_refalt()
        gts.fname = DATADIR / "test_write_missing.vcf"
        vals = len(np.unique(gts.data))

        gts.write()

        new_gts = GenotypesVCF(gts.fname)
        new_gts.read()

        new_gts.check_missing()

        # force two of the samples to have a missing GT
        gts.data = gts.data.astype(np.int8)
        gts.data[1, 1, 1] = -1
        gts.data = gts.data.astype(np.uint8)

        gts.write()

        new_gts = GenotypesVCF(gts.fname)
        new_gts.read()

        assert len(np.unique(gts.data)) == (vals + 1)

        with pytest.raises(ValueError) as info:
            gts.check_missing()
        assert (
            str(info.value)
            == "Genotype with ID 1:10116:A:G at POS 1:10116 is missing for sample"
            " HG00097"
        )

        gts.fname.unlink()

    def test_merge_variants_vcf(self):
        gts1 = Genotypes(DATADIR / "example.vcf.gz")
        gts2 = self._get_fake_genotypes_refalt()
        gts1.read()
        gts1.check_phase()
        gts1.subset(samples=gts2.samples, inplace=True)

        gts = Genotypes.merge_variants((gts1, gts2), fname=None)

        assert isinstance(gts, Genotypes)
        assert gts.samples == gts1.samples
        assert len(gts.variants) == len(gts1.variants) + len(gts2.variants)
        assert gts.data.shape[0] == gts1.data.shape[0]
        assert gts.data.shape[2] == gts1.data.shape[2]
        assert gts.data.shape[1] == (gts1.data.shape[1] + gts2.data.shape[1])


class TestGenotypesTR:
    def _get_fake_tr_alleles(self):
        # max uint8 codes for empty ("em") genotype
        em = np.iinfo(np.uint8).max
        return np.array(
            [
                [[1, 2, 1], [3, 4, 1], [5, 6, 1], [7, 8, 1], [9, 0, 1]],
                [[3, 1, 1], [3, 1, 1], [1, 1, 1], [1, 1, 1], [3, 3, 1]],
                [[3, 3, 1], [3, 3, 1], [3, 3, 1], [3, 3, 1], [3, 3, 1]],
                [[0, 11, 1], [em, em - 1, 1], [7, 2, 1], [3, 4, 1], [0, em, 1]],
                [[5, em, 0], [em, em, 0], [3, em, 0], [em, em, 0], [em, em, 0]],
            ],
            dtype=np.uint8,
        ).transpose((1, 0, 2))

    def test_read(self):
        expected_alleles = self._get_fake_tr_alleles()

        gts = GenotypesTR(DATADIR / "simple_tr.vcf")
        gts.read()

        # check genotypes
        np.testing.assert_allclose(expected_alleles, gts.data)

    def test_iter(self):
        # Get the expected data
        expected = self._get_fake_tr_alleles()

        gts = GenotypesTR(DATADIR / "simple_tr.vcf")
        # Check that everything matches what we expected
        for idx, line in enumerate(gts):
            np.testing.assert_allclose(line.data[:, :3], expected[:, idx])


class TestBreakpoints:
    def _get_expected_breakpoints(self):
        bps = Breakpoints(fname=None)
        create_arr = lambda *arr_list: np.array(list(arr_list), dtype=HapBlock)
        bps.data = {
            "Sample_1": [
                create_arr(
                    ("YRI", "1", 59423086, 85.107755),
                    ("CEU", "1", 239403765, 266.495714),
                    ("YRI", "2", 229668157, 244.341689),
                ),
                create_arr(
                    ("YRI", "1", 59423086, 85.107755),
                    ("YRI", "1", 239403765, 266.495714),
                    ("CEU", "2", 229668157, 244.341689),
                ),
            ],
            "Sample_2": [
                create_arr(
                    ("CEU", "1", 59423086, 85.107755),
                    ("YRI", "1", 239403765, 266.495714),
                    ("CEU", "2", 229668157, 244.341689),
                ),
                create_arr(
                    ("CEU", "1", 59423086, 85.107755),
                    ("CEU", "1", 239403765, 266.495714),
                    ("YRI", "2", 229668157, 244.341689),
                ),
            ],
        }
        return bps

    def _compare_bkpt_data(self, obs, exp):
        count = 0
        for samp, blocks in obs:
            for obs_block, exp_block in zip(blocks, exp[samp]):
                assert obs_block.tolist() == exp_block.tolist()
            count += 1
        assert count == len(exp)

    def test_load_breakpoints_iterate(self):
        expected = self._get_expected_breakpoints()

        # can we load the data from the VCF?
        bps = Breakpoints(DATADIR / "outvcf_test.bp")
        self._compare_bkpt_data(bps, expected.data)

    def test_load_breakpoints(self):
        expected = self._get_expected_breakpoints()

        # can we load the data from the VCF?
        bps = Breakpoints(DATADIR / "outvcf_test.bp")
        bps.read()

        # first, check that the samples appear in the proper order
        assert tuple(bps.data.keys()) == tuple(expected.data.keys())
        # now, check that each sample is the same
        self._compare_bkpt_data(bps.data.items(), expected.data)

    def test_write(self):
        expected = self._get_expected_breakpoints()
        expected.fname = Path("test.bp")
        expected.write()
        observed = Breakpoints(expected.fname)
        observed.read()

        # first, check that the samples appear in the proper order
        assert tuple(observed.data.keys()) == tuple(expected.data.keys())
        # now, check that each sample is the same
        self._compare_bkpt_data(observed.data.items(), expected.data)
        expected.fname.unlink()

    def test_load_underscore(self):
        """check if we can load samples with extra underscores in their IDs"""
        expected = self._get_expected_breakpoints()
        expected.fname = Path("test.bp")
        expected.data["Sam_ple_2"] = expected.data.pop("Sample_2")
        expected.write()
        observed = Breakpoints(expected.fname)
        observed.read()

        # first, check that the samples appear in the proper order
        assert tuple(observed.data.keys()) == tuple(expected.data.keys())
        # now, check that each sample is the same
        self._compare_bkpt_data(observed.data.items(), expected.data)
        expected.fname.unlink()

    def test_encode(self):
        expected = self._get_expected_breakpoints()
        expected.labels = {"YRI": 0, "CEU": 1}

        observed = self._get_expected_breakpoints()
        observed.encode()

        assert observed.labels == expected.labels
        assert len(expected.data) == len(observed.data)
        for sample in expected.data:
            for strand in range(len(expected.data[sample])):
                exp_strand = expected.data[sample][strand]
                obs_strand = observed.data[sample][strand]
                assert len(exp_strand) == len(observed.data[sample][strand])
                for obs, exp in zip(obs_strand["pop"], exp_strand["pop"]):
                    assert expected.labels[exp] == obs

    def test_recode(self):
        expected = self._get_expected_breakpoints()

        observed = self._get_expected_breakpoints()
        observed.encode()
        observed.recode()

        self._compare_bkpt_data(observed.data.items(), expected.data)

    def test_breakpoints_to_pop_array(self):
        variants = np.array(
            [("1", 59423086), ("1", 59423090), ("1", 239403770), ("2", 229668150)],
            dtype=[("chrom", "U10"), ("pos", np.uint32)],
        )
        expected_pop_arr = np.array(
            [
                [
                    [0, 0],
                    [1, 0],
                    [1, 0],
                    [0, 1],
                ],
                [
                    [1, 1],
                    [0, 1],
                    [0, 1],
                    [1, 0],
                ],
            ],
            dtype=np.uint8,
        )
        labels = {"YRI": 0, "CEU": 1}
        labels_map = np.vectorize({v: k for k, v in labels.items()}.get)

        expected = self._get_expected_breakpoints()

        pop_arr = expected.population_array(variants[[0, 1, 3]])
        exp_arr_with_labels = labels_map(expected_pop_arr[:, [0, 1, 3]])
        np.testing.assert_array_equal(pop_arr, exp_arr_with_labels)

        expected.encode()
        pop_arr = expected.population_array(variants[[0, 1, 3]])
        assert expected.labels == labels
        np.testing.assert_allclose(pop_arr, expected_pop_arr[:, [0, 1, 3]])

        with pytest.raises(ValueError) as info:
            labels, pop_arr = expected.population_array(variants)

    def test_breakpoints_to_pop_array_chrom_no_match(self):
        variants = np.array(
            [("chr1", 59423086), ("1", 59423090), ("1", 239403770), ("2", 229668150)],
            dtype=[("chrom", "U10"), ("pos", np.uint32)],
        )

        expected = self._get_expected_breakpoints()

        with pytest.raises(ValueError) as info:
            pop_arr = expected.population_array(variants[[0, 1, 3]])
        assert str(info.value).startswith("Chromosome ")


class TestDocExamples:
    def test_blocks2hap(self):
        # load the genotypes file
        # you can use either a VCF or PGEN file
        gt = GenotypesVCF.load(DATADIR / "simple.vcf.gz")
        gt = GenotypesPLINK.load(DATADIR / "simple.pgen")

        # load the haplotypes
        hp = Haplotypes("output.hap")
        hp.data = {}

        # iterate through lines of the .blocks.det file
        with open(DATADIR / "simple.blocks.det") as blocks_file:
            for idx, line in enumerate(blocks_file.read().splitlines()[1:]):
                # initialize variables and parse line from the blocks file
                hap_id = f"H{idx}"
                chrom, bp1, bp2, kb, nsnps, snps = line.strip().split()

                # create a haplotype line in the .hap file
                hp.data[hap_id] = Haplotype(
                    chrom=chrom, start=int(bp1), end=int(bp2), id=hap_id
                )

                # fetch alleles from the genotypes file
                snp_gts = gt.subset(variants=tuple(snps.split("|")))

                # create variant lines for each haplotype
                # Note that the .blocks.det file doesn't specify an allele, so
                # we simply choose the first allele (ie the REF allele) for this example
                hp.data[hap_id].variants = tuple(
                    Variant(
                        start=v["pos"],
                        end=v["pos"] + len(v["alleles"][0]),
                        id=v["id"],
                        allele=v["alleles"][0],
                    )
                    for v in snp_gts.variants
                )

        hp.write()

        # validate the output and clean up afterwards
        with open("output.hap") as hp_file:
            with open(DATADIR / "simple.blocks.det.hap") as expected:
                assert hp_file.read() == expected.read()
        hp.fname.unlink()

    def test_gts2hap(self):
        # load variants from the snplist file
        variants = {}
        with open(DATADIR / "apoe.snplist") as snplist_file:
            for line in snplist_file.readlines():
                # parse variant ID and beta from file
                ID, beta = line.split("\t")
                variants[ID] = float(beta)

        # load the genotypes file
        gt = GenotypesVCF(DATADIR / "apoe.vcf.gz")
        gt.read(variants=variants.keys())

        # initialize an empty haplotype file
        hp = Haplotypes("output.hap", haplotype=HaptoolsHaplotype)
        hp.data = {}

        for variant in gt.variants:
            ID, chrom, pos, alleles = variant[["id", "chrom", "pos", "alleles"]]
            # we arbitrarily choose to use the ALT allele but alleles[0] will give you REF
            allele = alleles[1]
            end = pos + len(allele)

            # create a haplotype line in the .hap file
            hp.data[ID] = HaptoolsHaplotype(
                chrom=chrom, start=pos, end=end, id=ID, beta=variants[ID]
            )

            # create a variant line for each haplotype
            hp.data[ID].variants = (Variant(start=pos, end=end, id=ID, allele=allele),)

        hp.write()

        # validate the output and clean up afterwards
        with open("output.hap") as hp_file:
            with open(DATADIR / "apoe.hap") as expected:
                assert hp_file.read() == expected.read()
        hp.fname.unlink()
