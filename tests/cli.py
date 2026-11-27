"""End-to-end regression checks through the configured executable."""
import hashlib
import json
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory

def load(path):
    return json.loads(path.read_text())


def measured(value):
    if isinstance(value, dict):
        return {k: measured(v) for k, v in value.items() if not k.endswith("seconds") and k != "ess_per_second"}
    if isinstance(value, list):
        return [measured(v) for v in value]
    return value


def invoke(config, output, jobs, success=True):
    p = subprocess.run([exe, "run", config, "--jobs", str(jobs), "--output", output],
                       capture_output=True, text=True)
    assert (p.returncode == 0) == success, p.stderr
    return p


def main(temp):
    for dataset in (root / "data").glob("*.json"):
        meta = load(dataset)
        assert hashlib.sha256((dataset.parent / meta["edges"]).read_bytes()).hexdigest() == meta["sha256"]
    cfg = {"all_starts_limit": 8, "sample_starts": 4, "sdp_limit": 10, "max_steps": 200,
           "snapshots": [0, 10], "epsilons": [.25, .1], "mcmc": {"replicates": 4, "draws": 500},
           "studies": [{"name": "cycle", "graph": {"family": "cycle", "sizes": [8, 12]},
                        "methods": ["clockwise", "max_degree", "fmmc"],
                        "policies": [{"policy": "none"}, {"policy": "random", "candidates": "anchored", "samples": 3}]}]}
    config = temp / "config.json"
    config.write_text(json.dumps(cfg))
    one, four = temp / "serial", temp / "parallel"
    invoke(config, one, 1)
    invoke(config, four, 4)
    assert load(one / "manifest.json")["status"] == "complete"
    assert len(list(one.glob("run-*.json"))) == 24
    for file in one.glob("run-*.json"):
        run = load(file)
        assert measured(run) == measured(load(four / file.name)), file.name
        if run["status"] == "skipped_sdp_limit":
            assert run["vertices"] == 12 and run["method"] == "fmmc"
            continue
        expected = 4 if run["vertices"] == 12 else 8
        assert len(run["analysis"]["starts"]) == expected
        assert run["analysis"]["scope"] == ("sampled" if expected == 4 else "all")
        if run["method"] == "max_degree" and run["edge"] is None:
            assert run["analysis"]["period"] == 2
            assert all(c["status"] == "periodic" and c["steps"] is None for c in run["analysis"]["crossings"])
        assert run["analysis"]["max_mass_error"] < 1e-10
        assert run["mcmc"]["replicates"] == 4
    invoke(config, one, 1, False)
    meta = load(root / "data" / "karate.json")
    meta["edges"] = str(root / "data" / "karate.edges")
    meta["sha256"] = "incorrect"
    (temp / "bad.json").write_text(json.dumps(meta))
    cfg["studies"][0]["graph"] = {"dataset": "bad.json"}
    config.write_text(json.dumps(cfg))
    assert "checksum mismatch" in invoke(config, temp / "bad-run", 1, False).stderr


if __name__ == "__main__":
    exe, root = Path(sys.argv[1]), Path(sys.argv[2])
    with TemporaryDirectory() as folder:
        main(Path(folder))

