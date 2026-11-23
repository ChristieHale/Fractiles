"""Reference tests of the fractile program (run with pytest from scripts/).

Each test folder has a Hazard case (HAZ inputs and the HAZ outputs they give)
and a Fractiles case (a fractile run file, the HAZ run, fault and out1 files it
reads, and the output it gives).

- Every Fractiles case is run and its output compared, number by number, with
  Fractiles/Output.
- Cluster_Test (a clustered sequence with two GMMs weighted 0.5 / 0.5 and no
  other epistemic branch) and Source_Types_Test (one branch per source and one
  GMM) are also checked against the exact answers, from the HAZ out1 and out3
  files in Hazard/Output.
- If HAZ_BIN is set to the HAZ executable, the Hazard cases are rerun and their
  out1 and out3 files compared with Hazard/Output, and the Fractiles inputs are
  checked to be copies of them. The PEER case takes a few minutes.

Set FRACT_BIN to the fractile executable (default: ../build/Fractiles).
"""

import glob
import math
import os
import re
import shutil
import subprocess

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..")
TESTS = ["Reference_Test", "PEER_Verification_Test", "Cluster_Test", "Source_Types_Test"]

FRACT_BIN = os.environ.get("FRACT_BIN", os.path.join(ROOT, "build", "Fractiles"))
HAZ_BIN = os.environ.get("HAZ_BIN")

NUMBER = re.compile(r"^[-+]?(\d+\.?\d*|\.\d+)([eE][-+]?\d+)?$")

needs_fract = pytest.mark.skipif(not os.path.isfile(FRACT_BIN), reason="set FRACT_BIN to the fractile executable")


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def compare(actual_path, expected_path, rtol=1e-4, atol=1e-12):
    """Differences between two text files, number by number (empty if they agree)."""
    with open(actual_path, errors="replace") as f:
        actual = f.read().split()
    with open(expected_path, errors="replace") as f:
        expected = f.read().split()
    if len(actual) != len(expected):
        return [f"{len(actual)} tokens, expected {len(expected)}"]
    errors = []
    for i, (a, e) in enumerate(zip(actual, expected)):
        if NUMBER.match(a) and NUMBER.match(e):
            if abs(float(a) - float(e)) > atol + rtol * abs(float(e)):
                errors.append(f"token {i}: {a} != {e}")
        elif a.lower() != e.lower():
            errors.append(f"token {i}: {a!r} != {e!r}")
    return errors


def run(binary, folder, stdin):
    p = subprocess.run([os.path.abspath(binary)], cwd=folder, input=stdin.encode(),
                       capture_output=True, timeout=3600)
    out = p.stdout.decode(errors="replace") + p.stderr.decode(errors="replace")
    assert p.returncode == 0, out
    return out


def fractile_output_name(folder):
    with open(os.path.join(folder, "Run_Fractiles.txt")) as f:
        return f.read().splitlines()[4].split()[0]


def read_fractiles(path):
    """Periods (levels, {fractile: values}, mean) and UHS tables of a fractile output."""
    periods, uhs = [], []
    lines = open(path).read().splitlines()
    i = 0
    while i < len(lines):
        line = lines[i]
        if line.strip().endswith("Period (sec)"):
            per = dict(period=float(line.split()[0]), levels=np.array([float(x) for x in lines[i + 1].split()]),
                       rows={}, mean=None)
            periods.append(per)
            i += 2
            continue
        t = line.split()
        if periods and len(t) > 2 and t[0] == "mean":
            periods[-1]["mean"] = np.array([float(x) for x in t[1:]])
        elif periods and len(t) > 2 and re.match(r"^0\.\d\d$", t[0]) and periods[-1]["mean"] is None:
            periods[-1]["rows"][round(float(t[0]), 2)] = np.array([float(x) for x in t[1:]])
        elif "UHS Fractiles for hazard level" in line:
            table = dict(level=float(t[-1]), rows=[])
            uhs.append(table)
        elif uhs and len(t) == 6 and NUMBER.match(t[0]):
            uhs[-1]["rows"].append([float(x) for x in t])
        i += 1
    return periods, uhs


def read_out1(path):
    """Blocks of a HAZ out1 file: (source, width, period) -> curves (GMM, mechanism, branch, level)."""
    with open(path) as f:
        tokens = f.read().split("\n", 1)[1].split()
    blocks = {}
    i = 0
    while i < len(tokens):
        iflt, iwidth, iprob, n_atten, _, n_ftype, n_param, n_inten = (int(t) for t in tokens[i:i + 8])
        i += 8
        n = n_atten * n_ftype * n_param * n_inten
        blocks[iflt, iwidth, iprob] = np.array([float(t) for t in tokens[i:i + n]]).reshape(
            n_atten, n_ftype, n_param, n_inten)
        i += n
    return blocks


def read_out3_totals(path):
    """The Wt_Total_Events/yr line of each period block of a HAZ out3 file."""
    totals = []
    for line in open(path):
        if line.strip().startswith("Wt_Total_Events/yr"):
            totals.append(np.array([float(x) for x in line.split()[4:]]))
    return totals


def loglog_level(levels, curve, haz):
    """Ground motion where a hazard curve falls to haz (log-log interpolation)."""
    k = int(np.argmax(curve < haz))
    assert k > 0 and curve[k] > 0
    x = (math.log(haz) - math.log(curve[k - 1])) / (math.log(curve[k]) - math.log(curve[k - 1]))
    return math.exp(math.log(levels[k - 1]) + x * (math.log(levels[k]) - math.log(levels[k - 1])))


def fractiles_case(test):
    return os.path.join(ROOT, test, "Fractiles")


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def fractile_runs(tmp_path_factory):
    """Run the fractile program once on every Fractiles case."""
    if not os.path.isfile(FRACT_BIN):
        pytest.skip("set FRACT_BIN to the fractile executable")
    root = tmp_path_factory.mktemp("fractiles")
    out = {}
    for test in TESTS:
        folder = root / test
        shutil.copytree(os.path.join(fractiles_case(test), "Input"), folder)
        screen = run(FRACT_BIN, folder, "Run_Fractiles.txt\n")
        assert "Normal Termination" in screen
        out[test] = folder / fractile_output_name(folder)
    return out


@needs_fract
@pytest.mark.parametrize("test", TESTS)
def test_matches_reference_output(fractile_runs, test):
    path = fractile_runs[test]
    errors = compare(path, os.path.join(fractiles_case(test), "Output", path.name))
    assert not errors, "\n".join(errors[:20])


@needs_fract
def test_cluster_two_point_distribution(fractile_runs):
    """Two GMMs at 0.5 / 0.5 and no other branch: each sample is one of two curves."""
    periods, uhs = read_fractiles(fractile_runs["Cluster_Test"])
    out1 = read_out1(os.path.join(ROOT, "Cluster_Test", "Hazard", "Output", "TB_Cluster_Site1.out1"))
    with open(os.path.join(fractiles_case("Cluster_Test"), "Input", "Run_Fractiles.txt")) as f:
        n_sample = int(f.read().splitlines()[1].split()[0])
    assert len(periods) == 2
    assert all(len(per["rows"]) == 99 for per in periods)
    assert len(uhs) == 2 and all(len(table["rows"]) == 2 for table in uhs)
    curves = []
    for iprob, per in enumerate(periods, 1):
        # The cluster's curves are in the first member's block; the others are zero
        c = out1[1, 1, iprob][:, 0, 0, :]
        assert c.shape[0] == 2
        assert all(np.all(out1[k, 1, iprob] == 0.0) for k in (2, 3))
        lo, hi = c.min(axis=0), c.max(axis=0)
        curves.append((lo, hi))
        for p, row in per["rows"].items():
            if p <= 0.45:
                np.testing.assert_allclose(row, lo, rtol=1e-3, err_msg=f"fractile {p}")
            elif p >= 0.55:
                np.testing.assert_allclose(row, hi, rtol=1e-3, err_msg=f"fractile {p}")
        # Mean: a binomial mixture of the two curves (4 standard deviations)
        exact = 0.5 * (c[0] + c[1])
        sd = 0.5 * np.abs(c[0] - c[1]) / math.sqrt(n_sample)
        assert np.all(np.abs(per["mean"] - exact) <= 4.0 * sd + 6e-4 * exact)
    # UHS: the 5th and 95th fractiles are the two curves at each hazard level
    for table in uhs:
        for iprob, row in enumerate(table["rows"]):
            lo, hi = curves[iprob]
            levels = periods[iprob]["levels"]
            assert row[1] == pytest.approx(loglog_level(levels, lo, table["level"]), rel=1e-3)
            assert row[5] == pytest.approx(loglog_level(levels, hi, table["level"]), rel=1e-3)


@needs_fract
def test_source_types_read_like_haz(fractile_runs):
    """Source types 1, 6 and 7, one branch each and one GMM: every sample is HAZ's total."""
    periods, uhs = read_fractiles(fractile_runs["Source_Types_Test"])
    total = read_out3_totals(os.path.join(ROOT, "Source_Types_Test", "Hazard", "Output", "SourceTypes_Site1.out3"))[0]
    per = periods[0]
    assert len(per["rows"]) == 99
    for p, row in per["rows"].items():
        np.testing.assert_allclose(row, total, rtol=1e-3, err_msg=f"fractile {p}")
    np.testing.assert_allclose(per["mean"], total, rtol=1e-3)


@pytest.mark.skipif(not HAZ_BIN, reason="set HAZ_BIN to also rerun the hazard cases")
@pytest.mark.parametrize("test", TESTS)
def test_hazard_outputs(tmp_path, test):
    folder = tmp_path / "hazard"
    shutil.copytree(os.path.join(ROOT, test, "Hazard", "Input"), folder)
    run_name = os.path.basename(glob.glob(os.path.join(folder, "Run_*.txt"))[0])
    run(HAZ_BIN, folder, f"0\n0\n{run_name}\n")
    reference = os.path.join(ROOT, test, "Hazard", "Output")
    checked = 0
    for path in sorted(glob.glob(os.path.join(reference, "*.out[13]"))):
        name = os.path.basename(path)
        errors = compare(folder / name, path)
        assert not errors, name + "\n" + "\n".join(errors[:20])
        checked += 1
    assert checked >= 2
    # The Fractiles case reads copies of the hazard inputs and out1
    fractile_input = os.path.join(fractiles_case(test), "Input")
    for name in os.listdir(fractile_input):
        if name == "Run_Fractiles.txt":
            continue
        source = os.path.join(reference if name.endswith(".out1") else os.path.join(ROOT, test, "Hazard", "Input"), name)
        assert os.path.exists(source), name
        assert not compare(os.path.join(fractile_input, name), source), name
