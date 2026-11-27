"""Audit complete notebook outputs without rerunning the simulations."""
import argparse
from collections import Counter, defaultdict
import csv
import hashlib
from itertools import groupby, combinations
import json
from math import isclose, ulp
from pathlib import Path

from notebook_runs import ROOT, OUT, read, validate
from comparison_stats import write_csv


def audit_budget(manifest, runs):
    cfg = manifest["config"]
    settings = cfg["wall_budget"]
    assert manifest["jobs"] == 1
    group = {g["id"]: s for g, s in zip(manifest["graphs"], cfg["studies"], strict=True)}
    expected = Counter((gid, method, policy["policy"], rep) for gid, s in group.items()
                       for method in s["methods"] for policy in s["policies"]
                       for rep in range(settings["repeats"]))
    assert Counter((r["graph"], r["method"], r["policy"], r["trial"]) for r in runs) == expected
    for r in runs:
        assert r["initial_position"] == 0
        assert [c["budget_seconds"] for c in r["checkpoints"]] == settings["seconds"]
        counts = [c["draws"] for c in r["checkpoints"]]
        assert counts == sorted(counts)
        for c in r["checkpoints"]:
            if c["draws"]:
                assert c["status"] == "ok" and all(0 <= c[k] <= 1 for k in ("rare_estimate","partition_estimate"))
            else:
                assert c["partition_estimate"] is None and c["rare_estimate"] is None
    def pairing(r):
        return r["graph"], r["policy"], r["trial"]
    for _, paired in groupby(sorted(runs, key=pairing), pairing):
        assert len({(r["seed"], tuple(r["edge"] or [])) for r in paired}) == 1


def audit_design(manifest, runs):
    """Derive graph, case and exhaustive-edge counts directly from the preset."""
    cfg = manifest["config"]
    graphs = manifest["graphs"]
    expected = Counter()
    all_edges = defaultdict(set)
    for study in cfg["studies"]:
        definition = study["graph"]
        matches = [g for g in graphs if g["name"] == study["name"] and g["generator"] == definition]
        if "fixtures" in definition:
            fixtures = [r for r in read(ROOT / "experiments" / definition["fixtures"])
                        if r["family"] == definition["family"]]
            for rep, (g, fixture) in enumerate(zip(matches, fixtures, strict=True)):
                edges = sorted(sorted(edge) for edge in fixture["edges"])
                assert (g["vertices"], g["replicate"], g["edges"]) == (fixture["n"], rep, edges)
        else:
            sizes = definition.get("sizes")
            if sizes is None:
                sizes = [read(ROOT/"experiments"/definition["dataset"])["vertices"]]
            assert sorted((g["vertices"],g["replicate"]) for g in matches) == sorted(
                (n,rep) for n in sizes for rep in range(study.get("replicates",1)))
        for g in matches:
            n = g["vertices"]
            for policy in study["policies"]:
                name = policy["policy"]
                if name == "halves" and n % 2:
                    continue
                count = policy.get("samples", 1)
                if name == "all":
                    existing = set(map(tuple,g["edges"]))
                    potential = ((i,n-1) for i in range(1,n-2)) if policy.get("candidates") == "anchored" else combinations(range(n),2)
                    candidates = {edge for edge in potential if edge not in existing}
                    count = len(candidates)
                for method in study["methods"]:
                    key = (g["id"],method,name)
                    expected[key] += count
                    if name == "all":
                        all_edges[key] = candidates
    actual = Counter((r["graph"],r["method"],r["policy"]) for r in runs)
    assert actual == expected, (actual-expected, expected-actual)
    for key, candidates in all_edges.items():
        assert {tuple(r["edge"]) for r in runs if (r["graph"],r["method"],r["policy"])==key} == candidates


def audit_run(r, graph, config):
    n = graph["vertices"]
    assert n == r["vertices"] and r["lambda2_after"] >= r["lambda2_before"] - 1e-9
    if r["status"] == "skipped_sdp_limit":
        assert r["method"] == "fmmc" and n > config.get("sdp_limit", 64)
        return
    assert r["status"] == "ok"
    if r["method"] == "fmmc":
        diagnostics = r["solver"]
        assert abs(diagnostics["objective"] - r["slem"]) < 1e-5
        assert max(abs(diagnostics[k]) for k in ("primal_residual", "dual_residual", "gap")) < 1e-5
        transitions = {(i,j):weight for i,j,weight in r["transitions"]}
        support = set(map(tuple, graph["edges"]))
        if r["edge"]:
            support.add(tuple(r["edge"]))
        sums = [0.0] * n
        for (i,j), weight in transitions.items():
            assert weight >= 0 and weight == transitions[j,i]
            assert i == j or tuple(sorted((i,j))) in support
            sums[i] += weight
        assert max(abs(v-1) for v in sums) < 1e-10
    a = r["analysis"]
    assert a["max_mass_error"] < 1e-9 and a["minimum_probability"] >= 0
    assert min(r["stationary"]) > 0 and isclose(sum(r["stationary"]), 1, abs_tol=1e-12)
    starts = a["starts"]
    assert starts == sorted(set(starts)) and a["diagnostic_start"] in starts
    assert [s["start"] for s in a["per_start"]] == starts
    assert (a["scope"] == "all") == (starts == list(range(r["states"])))
    for snap in a["snapshots"]:
        assert min(snap["distribution"]) >= 0 and abs(sum(snap["distribution"]) - 1) < 1e-9
    curves = {(c["space"], c["step"]): c for c in a["curves"]}
    borderline = []
    for j, c in enumerate(a["crossings"]):
        eps, t, space = c["epsilon"], c["steps"], c["space"]
        first = [s["first_crossings"][j] for s in a["per_start"]]
        if c["status"] == "resolved":
            assert 0 <= t <= c["checked_through"] and curves[space, t]["maximum"] <= eps
            assert min(first) >= 0 and max(first) <= t
            if space == "full" and max(first) != t:
                # A plateau equal to the threshold can straddle it by a few
                # floating-point ulps. Verify the earlier maximum and retain
                # both results rather than silently rounding the measurements.
                early = max(first)
                previous = curves[space, early]["maximum"]
                assert early < t and eps < previous <= eps + 16*ulp(eps)
                borderline.append(dict(run=r["id"], method=r["method"], graph=r["graph"],
                                       epsilon=eps, last_individual_first=early,
                                       earlier_maximum_tv=previous, simultaneous_first=t))
        else:
            assert t is None and curves[space, c["checked_through"]]["maximum"] > eps
            if c["status"] == "periodic":
                assert space == "full" and eps < 1-1/a["period"]
            else:
                assert c["status"] == "unresolved" and c["checked_through"] == config.get("max_steps", 100000)
    if "mcmc" in r:
        m = r["mcmc"]
        errors = [(v-m["true_mean"])**2 for v in m["means"]]
        assert isclose(m["rmse"], (sum(errors)/m["replicates"])**.5, abs_tol=1e-12)
    return borderline


def audit(output=OUT, *, presets=None):
    rows, borderline = [], []
    directory = ROOT / "experiments"
    files = sorted(directory.glob("*.json")) if presets is None else [directory / f"{name}.json" for name in presets]
    for file in files:
        folder, cfg = output / file.stem, read(file)
        m = validate(folder, cfg)
        graphs = {g["id"]: g for g in m["graphs"]}
        runs = [read(folder / f"run-{i}.json") for i in range(m["run_count"])]
        assert [r["id"] for r in runs] == list(range(len(runs)))
        if "wall_budget" in cfg:
            audit_budget(m, runs)
        else:
            audit_design(m, runs)
            for g in graphs.values():
                if g["source"]:
                    src = g["source"]
                    assert src["sha256"] == hashlib.sha256((ROOT/"data"/src["edges"]).read_bytes()).hexdigest()
            for r in runs:
                for entry in audit_run(r, graphs[r["graph"]], cfg) or []:
                    borderline.append(dict(preset=file.stem, **entry))
        status = Counter(r["status"] for r in runs)
        crossings = Counter(c["status"] for r in runs for c in r.get("analysis", {}).get("crossings", []))
        rows.append(dict(preset=file.stem, graphs=len(graphs), planned=len(runs), measured=status["ok"],
                         skipped=status["skipped_sdp_limit"], resolved=crossings["resolved"],
                         periodic=crossings["periodic"], unresolved=crossings["unresolved"]))
    with (output/"dana-parity.csv").open() as f:
        parity = list(csv.DictReader(f))
    assert [int(r["n"]) for r in parity] == [10,20,50,100,200,500]
    assert max(float(r["max_curve_error"]) for r in parity) < 1e-12
    assert read(output/"student-parity.json")["caroline_matched_crossings"] == 960
    write_csv(output / "threshold-roundoff.csv", borderline, ("preset","run","epsilon","last_individual_first","earlier_maximum_tv","simultaneous_first"))
    write_csv(output / "inventory.csv", rows)
    return rows


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("directory", type=Path, nargs="?", default=OUT)
    for row in audit(p.parse_args().directory):
        print(row)
