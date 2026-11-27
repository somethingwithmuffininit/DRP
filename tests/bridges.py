"""Check nested bridge supports, a spectral oracle, and transition probabilities."""
import json
from math import isclose
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory

def main(folder):
    root = Path(__file__).resolve().parents[1]
    cfg = json.loads((root / "experiments/bridge-count.json").read_text())
    cfg["export_kernel"] = True
    path = folder / "config.json"
    path.write_text(json.dumps(cfg))
    out = folder / "runs"
    cmd = [sys.argv[1], "run", path, "--jobs", "2", "--output"]
    subprocess.run([*cmd, out], check=True)
    manifest = json.loads((out / "manifest.json").read_text())
    previous = set()
    assert manifest["run_count"] == 12
    for g in manifest["graphs"]:
        b = g["generator"]["bridges"]
        edges = set(map(tuple, g["edges"]))
        cross = [(i,j) for i,j in edges if i < 8 <= j]
        assert len(edges) == 56+b and len(set(i for i,j in cross)) == len(set(j for i,j in cross)) == b
        assert previous <= edges
        previous = edges
    for f in out.glob("run-*.json"):
        r = json.loads(f.read_text())
        b = manifest["graphs"][int(r["graph"][1:])]["generator"]["bridges"]
        # Antisymmetric clique mode obeys x² - (m+2)x + 2b = 0, m=8.
        assert isclose(r["lambda2_after"], 5 - (25-2*b)**.5, abs_tol=1e-10)
        if r["method"] == "max_degree":
            assert all(i == j or p == 1/8 for i,j,p in r["transitions"])
        assert all(c["status"] == "resolved" for c in r["analysis"]["crossings"])
    for b in (0, 9):
        cfg["studies"][0]["graph"]["bridges"] = b
        path.write_text(json.dumps(cfg))
        err = subprocess.run([*cmd, folder / "invalid"], capture_output=True, text=True)
        assert err.returncode and "bridge count" in err.stderr

if __name__ == "__main__":
    with TemporaryDirectory() as folder:
        main(Path(folder))
