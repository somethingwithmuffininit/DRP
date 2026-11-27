"""Compare the MixBench CLI with Dana's Theorem 1 computations at every recorded step."""
import argparse
import csv
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
TOL = 1e-12


def cli(exe, directory, method, n, diagnostic, all_starts, max_steps, eps, snapshots, jobs):
    config = dict(studies=[dict(name=method, graph=dict(family="path", sizes=[n]),
                                methods=[method], policies=[dict(policy="none")])],
                  max_steps=max_steps, diagnostic_start=diagnostic, stride=1,
                  all_starts_limit=2*n if all_starts else 0, sample_starts=1,
                  snapshots=snapshots, epsilons=eps, export_kernel=True)
    path = directory / "config.json"
    path.write_text(json.dumps(config))
    out = directory / f"{method}-{diagnostic}-{n}"
    subprocess.run([exe, "run", path, "--jobs", str(jobs), "--output", out],
                   stdout=subprocess.DEVNULL, check=True)
    return json.loads((out / "run-0.json").read_text())


def compare_curve(cpp, reference):
    np.testing.assert_allclose(cpp, reference, rtol=0, atol=TOL)
    return float(np.max(np.abs(cpp - reference)))


def curve(run, space, metric):
    rows = [c for c in run["analysis"]["curves"] if c["space"] == space]
    np.testing.assert_array_equal([c["step"] for c in rows], np.arange(len(rows)))
    return np.array([c[metric] for c in rows])


def check_kernel(run, matrix, projection):
    actual = np.zeros_like(matrix)
    for i, j, probability in run["transitions"]:
        actual[i, j] = probability
    np.testing.assert_array_equal(actual, matrix)
    np.testing.assert_array_equal(run["state_to_position"], projection)
    np.testing.assert_array_equal(run["stationary"], np.full(len(matrix), 1 / len(matrix)))


def verify(dana, n, exe, directory, jobs):
    eps = dana.TV_THRESHOLDS
    steps = 12 * n
    lifted = cli(exe, directory, "dhn", n, 0, True, steps, eps, [4*n], jobs)
    check_kernel(lifted, dana.build_matrix(n), dana.position_of_index(n))
    np.testing.assert_array_equal(lifted["analysis"]["starts"], np.arange(2*n))
    ref, block = dana.run_chain(n, steps)
    worst = curve(lifted, "full", "maximum")
    least = curve(lifted, "full", "minimum")
    zero = curve(lifted, "full", "diagnostic")
    error = max(compare_curve(worst, ref.max(axis=1)), compare_curve(least, ref.min(axis=1)),
                compare_curve(zero, ref[:, 0]))
    crossings = [c["steps"] for c in lifted["analysis"]["crossings"] if c["space"] == "full"]
    np.testing.assert_array_equal(crossings, [dana.first_time_below(ref.max(axis=1), e) for e in eps])
    # Check the crossing from each start as well as the maximum curve.
    for s in lifted["analysis"]["per_start"]:
        np.testing.assert_array_equal(s["first_crossings"][:len(eps)],
                                      [dana.first_time_below(ref[:, s["start"]], e) for e in eps])
    snap = next(s for s in lifted["analysis"]["snapshots"] if s["step"] == 4*n)
    constant = 2*n*snap["min_probability_over_starts"]
    np.testing.assert_allclose(constant, 2*n*block.min(), rtol=0, atol=TOL)
    compare_curve(np.array(snap["distribution"]), block[0])
    if zero[1:n+1].min() < dana.LOWER_BOUND or constant < dana.C_COUPLING:
        raise AssertionError("converse or minorization failed")
    for bound in (dana.C_COUPLING, constant):
        if np.any(worst > (1-bound)**(np.arange(steps+1)//(4*n)) + TOL):
            raise AssertionError("upper bound failed")

    rev, nonrev = dana.compare(n)
    position_times = []
    for method, initial, space, reference in (("dhn", 1, "position", nonrev),
                                              ("max_degree", 0, "full", rev)):
        # One physical starting position, with (+,1) represented by DHN index 1.
        run = cli(exe, directory, method, n, initial, False, len(reference)-1, eps, [], jobs)
        if method == "max_degree":
            check_kernel(run, dana.build_reversible(n), np.arange(n))
        recorded = curve(run, space, "diagnostic")
        error = max(error, compare_curve(recorded, reference))
        hits = [c["steps"] for c in run["analysis"]["crossings"] if c["space"] == space]
        np.testing.assert_array_equal(hits, [dana.first_time_below(reference, e) for e in eps])
        position_times.append(dana.first_time_below(recorded, dana.COMPARE_THRESHOLD))
    return dict(n=n, **{f"full_tv_{e}": t for e, t in zip(eps, crossings)},
                position_dhn=position_times[0], position_reversible=position_times[1],
                measured_C=constant, min_early_tv=float(zero[1:n+1].min()),
                max_start_spread=float((worst-least).max()), max_curve_error=error)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--exe", type=Path, default=ROOT/"build/mixbench")
    parser.add_argument("--sizes", type=int, nargs="+", help="default: every size in Dana/theorem1.py")
    parser.add_argument("--jobs", type=int, default=4)
    args = parser.parse_args()
    spec = importlib.util.spec_from_file_location("dana", ROOT.parent/"Dana/theorem1.py")
    dana = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(dana)
    writer = None
    for n in args.sizes or dana.N_VALUES:
        with TemporaryDirectory() as temp:
            row = verify(dana, n, args.exe.resolve(), Path(temp), args.jobs)
        if writer is None:
            writer = csv.DictWriter(sys.stdout, row.keys())
            writer.writeheader()
        writer.writerow(row)
        sys.stdout.flush()


if __name__ == "__main__":
    main()
