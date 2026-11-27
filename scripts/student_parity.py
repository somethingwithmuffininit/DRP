"""Compare the original students' implementations with the MixBench reproduction preset."""
import contextlib
import csv
import io
import json
from pathlib import Path
import runpy
import subprocess
from tempfile import TemporaryDirectory

import numpy as np

from notebook_runs import BUILD, ROOT, read


def verify_students(folder, output):
    manifest = read(folder / "manifest.json")
    runs = [read(folder / f"run-{i}.json") for i in range(manifest["run_count"])]
    with contextlib.redirect_stdout(io.StringIO()):
        original = runpy.run_path(str(ROOT.parent / "Danielle/DRP_FMMC_Final.py"))
    optimized = next(r for r in runs if r["method"] == "fmmc")
    actual = np.zeros((5, 5))
    for i, j, p in optimized["transitions"]:
        actual[i, j] = p
    error = float(np.max(np.abs(actual - original["P"].value)))
    np.testing.assert_allclose(actual, original["P"].value, atol=1e-6, rtol=0)
    exe_dir = BUILD / "caroline"
    subprocess.run(["cmake", "-S", ROOT.parent / "Caroline/cycle-mixing-project", "-B", exe_dir,
                    "-DCMAKE_BUILD_TYPE=Release"], check=True, stdout=subprocess.DEVNULL)
    subprocess.run(["cmake", "--build", exe_dir, "--target", "cycle_mixing", "-j", "2"],
                   check=True, stdout=subprocess.DEVNULL)
    expected = {(r["edge"][0]+1 if r["edge"] else None, s["start"]): s["first_crossings"][0]
                for r in runs if r["method"] == "clockwise" for s in r["analysis"]["per_start"]}
    with TemporaryDirectory() as tmp:
        subprocess.run([exe_dir / "cycle_mixing", "--mode", "placement", "--n", "32",
                        "--epsilon", "0.25", "--output", tmp], check=True, capture_output=True)
        with (Path(tmp) / "starting_vertices.csv").open() as f:
            observed = {(int(r["k"]) if r["k"] else None, int(r["start_index"])): int(r["first_crossing"])
                        for r in csv.DictReader(f)}
    if observed != expected or len(expected) != 960:
        raise AssertionError("Caroline's per-start crossings disagree")
    result = dict(caroline_matched_crossings=len(expected), danielle_max_matrix_error=error)
    (output / "student-parity.json").write_text(json.dumps(result, indent=2) + "\n")
    return result
