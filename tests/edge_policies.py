"""Check exhaustive and heuristic single-edge edits against dense Laplacians."""
import json
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory

import numpy as np


def check_graph(graph, runs):
    n = graph["vertices"]
    adj = np.zeros((n, n))
    for i, j in graph["edges"]:
        adj[i, j] = adj[j, i] = 1
    degrees = adj.sum(axis=1)
    laplacian = np.diag(degrees) - adj
    eigenvalues, vectors = np.linalg.eigh(laplacian)
    fiedler = vectors[:, 1]
    connectivity, separation = {}, {}
    for i in range(n):
        for j in range(i+1, n):
            if adj[i, j]:
                continue
            incidence = np.zeros(n)
            incidence[i], incidence[j] = 1, -1
            edited = laplacian + np.outer(incidence, incidence)
            connectivity[i, j] = np.linalg.eigvalsh(edited)[1]
            separation[i, j] = abs(fiedler[i] - fiedler[j])

    exhaustive = [tuple(r["edge"]) for r in runs if r["policy"] == "all"]
    assert sorted(exhaustive) == sorted(connectivity)
    assert {r["policy"] for r in runs} == {"all", "none", "random", "min_degree", "fiedler", "optimal_connectivity"}
    for r in runs:
        if r["policy"] == "none":
            assert r["edge"] is None
            continue
        edge = tuple(r["edge"])
        np.testing.assert_allclose(r["lambda2_after"], connectivity[edge], rtol=0, atol=2e-12)
        if r["policy"] == "optimal_connectivity":
            assert max(connectivity.values()) - connectivity[edge] < 2e-12
        if r["policy"] == "min_degree":
            assert min(degrees[list(edge)]) == min(degrees[degrees < n-1])
        # A repeated eigenvalue permits multiple valid Fiedler vectors.
        if r["policy"] == "fiedler" and eigenvalues[2] - eigenvalues[1] > 1e-8:
            assert max(separation.values()) - separation[edge] < 1e-10


if __name__ == "__main__":
    root = Path(__file__).resolve().parents[1]
    with TemporaryDirectory() as temp:
        out = Path(temp)/"edges"
        subprocess.run([sys.argv[1], "run", root/"experiments/edge-comparison.json",
                        "--output", out, "--jobs", "4"], check=True)
        manifest = json.loads((out/"manifest.json").read_text())
        runs = [json.loads((out/f"run-{i}.json").read_text()) for i in range(manifest["run_count"])]
        for graph in manifest["graphs"]:
            check_graph(graph, [r for r in runs if r["graph"] == graph["id"]])
