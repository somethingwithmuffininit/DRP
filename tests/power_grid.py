"""Independent SciPy reference for the full grid's two sparse transition rules."""
import argparse
from pathlib import Path
import sys
import tempfile

import numpy as np
from scipy import sparse

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from notebook_runs import ROOT, EXE, read, single_thread_env
import subprocess


def verify(folder):
    edges = np.loadtxt(ROOT/"data/power.edges", dtype=int)
    i, j = edges.T
    n = read(ROOT/"data/power.json")["vertices"]
    degrees = np.bincount(edges.ravel(), minlength=n)
    runs = [read(f) for f in folder.glob("run-*.json")]
    assert len(runs) == 3
    times = {}
    common_starts = None
    for r in runs:
        if r["method"] == "fmmc":
            assert r["status"] == "skipped_sdp_limit"
            continue
        a = r["analysis"]
        assert a["scope"] == "sampled" and len(a["starts"]) == 32
        if common_starts is not None:
            assert a["starts"] == common_starts
        common_starts = a["starts"]
        weights = (np.full(len(edges), 1/degrees.max()) if r["method"] == "max_degree"
                   else 1/np.maximum(degrees[i],degrees[j]))
        p = sparse.csr_matrix((np.tile(weights,2), (np.r_[i,j],np.r_[j,i])), shape=(n,n))
        p += sparse.diags(1-np.asarray(p.sum(axis=1)).ravel())
        x = np.zeros(n)
        x[a["diagnostic_start"]] = 1
        target = 1/n
        curves = {c["step"]: c["diagnostic"] for c in a["curves"]}
        snaps = {s["step"]: s["distribution"] for s in a["snapshots"]}
        first = None
        for t in range(max(curves)+1):
            tv = .5*np.abs(x-target).sum()
            if first is None and tv <= .25:
                first = t
            if t in curves:
                assert abs(tv-curves[t]) < 2e-11
            if t in snaps:
                np.testing.assert_allclose(x, snaps[t], rtol=0, atol=2e-12)
            x = p @ x  # symmetric P
        index = a["starts"].index(a["diagnostic_start"])
        assert first == a["per_start"][index]["first_crossings"][0]
        times[r["method"]] = a["crossings"][0]["steps"]
    assert times == {"max_degree":28515, "metropolis":7721}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path)
    parser.add_argument("--exe", type=Path, default=EXE)
    args = parser.parse_args()
    if args.results:
        verify(args.results)
    else:
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp)/"grid"
            subprocess.run([args.exe, "run", ROOT/"experiments/power-grid.json", "--output", out, "--jobs", "4"],
                           check=True, env=single_thread_env())
            verify(out)
