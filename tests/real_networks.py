"""Cross-check every configured Glasgow and karate-club case from the original edge files."""
import json
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory

import numpy as np


def verify(directory, root):
    manifest = json.loads((directory/"manifest.json").read_text())
    assert manifest["status"] == "complete" and manifest["run_count"] == 17
    graphs = {g["id"]: g for g in manifest["graphs"]}
    observed = {}
    for path in directory.glob("run-*.json"):
        r = json.loads(path.read_text())
        g = graphs[r["graph"]]
        n, method = r["vertices"], r["method"]
        a = np.zeros((n, n))
        edges = np.loadtxt(root/"data"/g["source"]["edges"], dtype=int)
        a[tuple(edges.T)] = 1
        a = np.maximum(a, a.T)
        base = np.diag(a.sum(axis=1)) - a
        if r["edge"]:
            i, j = r["edge"]
            assert a[i, j] == 0
            a[i, j] = a[j, i] = 1
        degrees = a.sum(axis=1)
        l = np.diag(degrees) - a
        before, after = (np.linalg.eigvalsh(m)[1] for m in (base, l))
        np.testing.assert_allclose([r["lambda2_before"], r["lambda2_after"], r["delta_lambda2"]],
                                   [before, after, after-before], atol=1e-11, rtol=0)
        if method == "clockwise":
            p = .5*np.eye(n) + .5*np.roll(np.eye(n), 1, axis=1)
            if r["edge"]:
                p[i, i] = p[j, j] = p[i, j] = p[j, i] = .25
        elif method == "fmmc":
            p = np.zeros_like(a)
            for i, j, value in r["transitions"]:
                p[i, j] = value
            np.testing.assert_array_equal(p, p.T)
            assert np.all(p >= 0) and (p * (1 - a - np.eye(n)) == 0).all()
            np.testing.assert_allclose(p.sum(axis=1), 1, atol=1e-12, rtol=0)
        else:
            assert method in ("max_degree", "metropolis")
            p = a / (degrees.max() if method == "max_degree" else np.maximum.outer(degrees, degrees))
            p += np.diag(1-p.sum(axis=1))
        slem = max(abs(np.linalg.eigvals(p - 1/n)))
        np.testing.assert_allclose(r["slem"], slem, atol=1e-9, rtol=0)
        np.testing.assert_allclose(r["absolute_spectral_gap"], 1 - slem, atol=1e-10)
        np.testing.assert_allclose(r["asymptotic_rate"], -np.log(slem), atol=1e-10)
        np.testing.assert_allclose(r["lambda3_after"], np.linalg.eigvalsh(l)[2], atol=1e-10)
        if r["reversible"]:
            eig = np.linalg.eigvalsh(p)
            np.testing.assert_allclose([r["transition_lambda2"], r["transition_lambda_min"], r["spectral_gap"]],
                                       [eig[-2], eig[0], 1-eig[-2]], atol=1e-10)
        if method == "fmmc":
            assert abs(slem - r["solver"]["objective"]) < 1e-5
        analysis = r["analysis"]
        assert analysis["scope"] == "all" and analysis["starts"] == list(range(n))
        terminal = analysis["crossings"][0]["checked_through"]
        mass = np.eye(n)
        history = []
        for _ in range(terminal+1):
            history.append(np.abs(mass - 1/n).sum(axis=1)*.5)
            mass = mass @ p
        history = np.array(history)
        for row in analysis["curves"]:
            values = history[row["step"]]
            np.testing.assert_allclose([row["minimum"], row["maximum"], row["diagnostic"]],
                                       [values.min(), values.max(), values[analysis["diagnostic_start"]]],
                                       atol=2e-12, rtol=0)
        for k, c in enumerate(analysis["crossings"]):
            hits = history <= c["epsilon"]
            first = np.where(hits.any(axis=0), hits.argmax(axis=0), -1)
            np.testing.assert_array_equal([s["first_crossings"][k] for s in analysis["per_start"]], first)
            common = np.flatnonzero(hits.all(axis=1))
            assert c["steps"] == (int(common[0]) if len(common) else None)
            if c["epsilon"] == .25:
                observed[g["name"], method, r["policy"]] = c["steps"]
    assert [observed["Glasgow ring", "clockwise", p] for p in ("none", "fixed")] == [43, 20]
    assert [observed["Karate club", m, "none"] for m in ("max_degree", "metropolis", "fmmc")] == [48, 35, 26]


if __name__ == "__main__":
    root = Path(__file__).resolve().parents[1]
    with TemporaryDirectory() as temp:
        out = Path(temp)/"real"
        subprocess.run([sys.argv[1], "run", root/"experiments/real-networks.json",
                        "--output", out, "--jobs", "4"], check=True)
        verify(out, root)
