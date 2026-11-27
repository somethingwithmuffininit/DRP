"""Check the shared runner against archived paired single-link trials."""
import argparse
from collections import defaultdict
import csv
import json
from math import isclose
from pathlib import Path
import subprocess
from tempfile import TemporaryDirectory

ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = ROOT / "results" / "aditi"
PRESET = ROOT / "experiments" / "link-policies.json"


def read(path):
    return json.loads(path.read_text())


def check(folder, fixtures):
    expected_graphs = {(r["family"], r["trial"]): r for r in fixtures}
    with (ARCHIVE / "trials.csv").open() as handle:
        expected = {
            (r["family"], int(r["trial"]), r["policy"]): r
            for r in csv.DictReader(handle)
            if (r["family"], int(r["trial"])) in expected_graphs
        }
    manifest = read(folder / "manifest.json")
    assert manifest["status"] == "complete"
    assert manifest["run_count"] == len(expected)
    graphs = {g["id"]: g for g in manifest["graphs"]}
    assert len(graphs) == len(fixtures)
    for graph in graphs.values():
        fixture = expected_graphs[graph["generator"]["family"], graph["replicate"]]
        assert graph["vertices"] == fixture["n"]
        assert graph["edges"] == fixture["edges"]

    starts = defaultdict(set)
    for index in range(manifest["run_count"]):
        run = read(folder / f"run-{index}.json")
        graph = graphs[run["graph"]]
        family, trial = graph["generator"]["family"], graph["replicate"]
        row = expected.pop((family, trial, run["policy"]))
        edge = [int(row["u"]), int(row["v"])] if row["u"] else None
        assert run["edge"] == edge
        assert run["start_seed"] == expected_graphs[family, trial]["start_seed"]
        analysis = run["analysis"]
        assert len(analysis["starts"]) == int(row["starts"]) == 16
        starts[run["graph"]].add(tuple(analysis["starts"]))
        assert analysis["crossings"][0]["steps"] == int(row["steps_at_tv_025"])
        for key, column in [
            ("lambda2_before", "lambda2"), ("lambda3_before", "lambda3"),
            ("lambda2_after", "edited_lambda2"), ("normalized_gain", "mu_star"),
        ]:
            assert isclose(run[key], float(row[column]), rel_tol=0, abs_tol=1e-9), (family, trial, key)
    assert not expected and all(len(samples) == 1 for samples in starts.values())
    print(f"Matched {manifest['run_count']} trials across {len(graphs)} graph realizations.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    target = parser.add_mutually_exclusive_group(required=True)
    target.add_argument("--exe", type=Path)
    target.add_argument("--results", type=Path, help="verify a completed run of the full preset")
    args = parser.parse_args()
    fixtures = read(PRESET.parent / "fixtures" / "link-policy-graphs.json")
    if args.results:
        check(args.results, fixtures)
        return

    # One realization of each family checks both policy sets. The full preset
    # can be verified separately without making the regression suite lengthy.
    sample = [r for r in fixtures if r["trial"] == 0]
    with TemporaryDirectory() as temp:
        folder = Path(temp)
        (folder / "fixtures.json").write_text(json.dumps(sample))
        cfg = read(PRESET)
        for study in cfg["studies"]:
            study["graph"]["fixtures"] = "fixtures.json"
        config = folder / "preset.json"
        config.write_text(json.dumps(cfg))
        output = folder / "runs"
        subprocess.run([args.exe, "run", config, "--output", output, "--jobs", "2"], check=True)
        check(output, sample)


if __name__ == "__main__":
    main()
