import math
from pathlib import Path

import numpy as np
from click.testing import CliRunner

from haptools.data import Data
from haptools.__main__ import main
from haptools.ld import pearson_corr_ld

DATADIR = Path(__file__).parent.joinpath("data")


def test_ld(seed=42):
    rng = np.random.default_rng(seed)

    # Different input shapes will give different output shapes
    # (25,) x (25,) --> float
    arrA = rng.choice((0, 1, 2), size=(25,))
    arrB = rng.choice((0, 1, 2), size=(25,))
    ld = pearson_corr_ld(arrA, arrB)
    assert isinstance(ld, float)
    assert math.isclose(ld, -0.1148198316929615)

    # (25,) x (25,1) --> (1,)
    arrB = arrB[:, np.newaxis]
    old_ld = ld
    ld = pearson_corr_ld(arrA, arrB)
    assert isinstance(ld, np.ndarray)
    assert ld.shape == (1,)
    assert old_ld == ld[0]

    # (25,1) x (25,1) --> (1,1)
    arrA = arrA[:, np.newaxis]
    ld = pearson_corr_ld(arrA, arrB)
    assert isinstance(ld, np.ndarray)
    assert ld.shape == (1, 1)
    assert old_ld == ld[0, 0]

    # (25,3) x (25,) --> (3,)
    arrA = np.hstack((np.random.choice((0, 1, 2), size=(25, 2)), arrA))
    arrB = np.squeeze(arrB)
    ld = pearson_corr_ld(arrA, arrB)
    assert isinstance(ld, np.ndarray)
    assert ld.shape == (3,)
    assert old_ld == ld[2]

    # (25,) x (25,3) --> (3,)
    arrA = arrB
    arrB = np.random.choice((0, 1, 2), size=(25, 3))
    ld = pearson_corr_ld(arrA, arrB)
    assert isinstance(ld, np.ndarray)
    assert ld.shape == (3,)

    # (25,1) x (25,3) --> (1,3)
    arrA = arrA[:, np.newaxis]
    old_ld = ld
    ld = pearson_corr_ld(arrA, arrB)
    assert isinstance(ld, np.ndarray)
    assert ld.shape == (1, 3)
    np.testing.assert_allclose(old_ld[np.newaxis, :], ld)

    # (25,2) x (25,3) --> (2,3)
    arrA = np.hstack((arrA, np.random.choice((0, 1, 2), size=(25, 1))))
    ld = pearson_corr_ld(arrA, arrB)
    assert isinstance(ld, np.ndarray)
    assert ld.shape == (2, 3)
    np.testing.assert_allclose(old_ld, ld[0])


def test_basic(capfd):
    expected = """#\torderH\tld
#\tversion\t0.2.0
#H\tld\t.3f\tLinkage-disequilibrium
H\t21\t26938353\t26938989\tchr21.q.3365*11\t0.995
H\t21\t26938989\t26941960\tchr21.q.3365*10\t-0.012
V\tchr21.q.3365*10\t26938989\t26938989\t21_26938989_G_A\tA
V\tchr21.q.3365*10\t26940815\t26940815\t21_26940815_T_C\tT
V\tchr21.q.3365*10\t26941960\t26941960\t21_26941960_A_G\tA
V\tchr21.q.3365*11\t26938353\t26938353\t21_26938353_T_C\tT
V\tchr21.q.3365*11\t26938989\t26938989\t21_26938989_G_A\tA
"""
    gt_file = DATADIR / "example.vcf.gz"
    hp_file = DATADIR / "basic.hap.gz"

    cmd = f"ld chr21.q.3365*1 {gt_file} {hp_file}"
    runner = CliRunner()
    result = runner.invoke(main, cmd.split(" "), catch_exceptions=False)
    captured = capfd.readouterr()
    assert captured.out == expected
    assert result.exit_code == 0


def test_simple_with_repeat(capfd):
    gt_file = DATADIR / "simple.vcf"
    hp_file = DATADIR / "simple.hap"

    cmd = f"ld H1 {gt_file} {hp_file}"
    runner = CliRunner()
    result = runner.invoke(main, cmd.split(" "), catch_exceptions=False)
    captured = capfd.readouterr()
    assert captured.out
    expected = captured.out
    assert result.exit_code == 0

    # now try again. The result should be the same because it ignores repeats
    hp_file = DATADIR / "simple_tr.hap"

    cmd = f"ld H1 {gt_file} {hp_file}"
    runner = CliRunner()
    result = runner.invoke(main, cmd.split(" "), catch_exceptions=False)
    captured = capfd.readouterr()
    assert expected == captured.out
    assert result.exit_code == 0


def test_basic_variant(capfd):
    expected = """#\torderH\tld
#\tversion\t0.2.0
#H\tld\t.3f\tLinkage-disequilibrium
H\t19\t45411941\t45412079\tAPOe4\t0.999
V\tAPOe4\t45411941\t45411941\trs429358\tC
V\tAPOe4\t45412079\t45412079\trs7412\tC
"""
    tmp_file = Path("apoe4_ld.hap")
    gt_file = DATADIR / "apoe.vcf.gz"
    hp_file = DATADIR / "apoe4.hap"

    cmd = f"ld -o {tmp_file} rs429358 {gt_file} {hp_file}"
    runner = CliRunner()
    result = runner.invoke(main, cmd.split(" "), catch_exceptions=False)
    captured = capfd.readouterr()
    assert captured.out == ""
    assert result.exit_code == 0

    with Data.hook_compressed(tmp_file, mode="r") as haps:
        assert haps.read() == expected

    tmp_file.unlink()


def test_from_gts(capfd):
    expected = """CHR\tBP\tSNP\tR
19\t45411941\trs429358\t0.999
19\t45411947\trs11542041\t0.027
19\t45411962\trs573658040\t-0.012
19\t45411965\trs543363163\t-0.012
19\t45412006\trs563140413\t-0.012
19\t45412007\trs531939919\t-0.012
19\t45412040\trs769455\t0.006
19\t45412079\trs7412\t-0.098
"""
    tmp_file = Path("apoe4.ld")
    gt_file = DATADIR / "apoe.vcf.gz"
    hp_file = DATADIR / "apoe4.hap"

    cmd = f"ld --from-gts -o apoe4.ld APOe4 {gt_file} {hp_file}"
    runner = CliRunner()
    result = runner.invoke(main, cmd.split(" "), catch_exceptions=False)
    captured = capfd.readouterr()
    assert captured.out == ""
    assert result.exit_code == 0

    with Data.hook_compressed(tmp_file, mode="r") as snps:
        assert snps.read() == expected

    tmp_file.unlink()


def test_from_gts_ids(capfd):
    expected = """CHR\tBP\tSNP\tR
19\t45411965\trs543363163\t-0.012
19\t45412079\trs7412\t-0.098
"""
    gt_file = DATADIR / "apoe.vcf.gz"
    hp_file = DATADIR / "apoe4.hap"

    cmd = f"ld --from-gts -i rs543363163 -i rs7412 APOe4 {gt_file} {hp_file}"
    runner = CliRunner()
    result = runner.invoke(main, cmd.split(" "), catch_exceptions=False)
    captured = capfd.readouterr()
    assert captured.out == expected
    assert result.exit_code == 0
