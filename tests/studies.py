"""Small end-to-end fixtures for the four comparison experiments."""
import csv
import json
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"scripts"))
from study_reports import requested_reports, budget_report


def execute(folder, name, config):
    file, out = folder/(name+".json"), folder/name
    file.write_text(json.dumps(config))
    subprocess.run([sys.argv[1], "run", file, "--output", out, "--jobs", "2"], check=True)
    manifest = json.loads((out/"manifest.json").read_text())
    assert manifest["status"] == "complete"
    runs = [json.loads((out/f"run-{i}.json").read_text()) for i in range(manifest["run_count"])]
    requested_reports(out, manifest, runs, config.get("reports", ["wall-budget"]), figures=False)
    return out, manifest, runs


def table(out, name):
    with (out/name).open() as f:
        return list(csv.DictReader(f))


def edge_ranking(folder):
    policies = [{"policy": p} for p in ("none", "all", "optimal_connectivity")]
    s = dict(name="cycle", graph=dict(family="cycle", sizes=[8]), methods=["clockwise"], policies=policies)
    out, _, _ = execute(folder, "edge", dict(studies=[s], reports=["edge-ranking"], epsilons=[.25,.1],
                                             max_steps=500, snapshots=[]))
    agreements = table(out, "edge-agreement.csv")
    assert len(agreements) == 4
    assert all(int(r["connectivity_best_ties"]) == 20 and
               r["shared_best_edges"] == r["fastest_ties"] for r in agreements)



def path_oracle(folder):
    path = json.loads((ROOT/"experiments/path-oracle.json").read_text())
    path["studies"][0]["graph"]["sizes"] = [4, 8, 16, 32]
    path["max_steps"] = 700
    out, _, _ = execute(folder, "path", path)
    assert all(r["oracle_pass"] == "True" and r["steps"] for r in table(out, "path-controls.csv"))
    assert {r["dropped_smallest"] for r in table(out, "path-sensitivity.csv")} == {"0","1"}



def nested_networks(folder):
    nested = json.loads((ROOT/"experiments/nested-networks.json").read_text())
    nested["max_steps"] = 10
    for study, data in zip(nested["studies"], ["glasgow","power"]):
        study["graph"].update(dataset=str(ROOT/"data"/(data+".json")), sizes=[5,8])
        study.update(replicates=2, methods=["max_degree","metropolis"], policies=[{"policy":"none"}])
        if data == "glasgow":
            study["policies"].append({"policy": "random", "samples": 2})
    out, manifest, runs = execute(folder, "nested", nested)
    prefixes = {}
    for g in manifest["graphs"]:
        parent = g["parent_vertex_ids"]
        edges = {tuple(sorted(map(int, l.split()))) for l in
                 (ROOT/"data"/g["source"]["edges"]).read_text().splitlines() if l.strip() and not l.startswith("#")}
        assert g["edges"] == [[i,j] for i in range(len(parent)) for j in range(i+1,len(parent))
                              if tuple(sorted((parent[i],parent[j]))) in edges]
        key = (g["name"],g["replicate"])
        assert parent[:len(prefixes.get(key, []))] == prefixes.get(key, [])
        prefixes[key] = parent
        paired = [[(r["edge"], r["analysis"]["starts"]) for r in runs
                   if r["graph"] == g["id"] and r["method"] == method]
                  for method in ("max_degree", "metropolis")]
        assert paired[0] and paired[0] == paired[1]
    assert all(r["realizations"] == "2" for r in table(out, "nested-variation.csv"))



def wall_budget(folder):
    s = dict(name="budget", graph=dict(family="barbell",sizes=[8]), methods=["max_degree","metropolis"],
             policies=[dict(policy="random")])
    out, manifest, runs = execute(folder, "budget", dict(studies=[s], wall_budget=dict(seconds=[1e-12,.02], repeats=2)))
    assert manifest["jobs"] == 1
    for trial in (0,1):
        a,b = [r for r in runs if r["trial"] == trial]
        assert a["seed"] == b["seed"] and a["edge"] == b["edge"]
        assert a["checkpoints"][0]["status"] == "setup_exhausted"
        assert a["checkpoints"][0]["draws"] == 0
        for run in (a, b):
            final = run["checkpoints"][-1]
            assert final["status"] == "ok" and final["draws"] > 0
            assert all(0 <= final[k] <= 1 for k in ("rare_estimate", "partition_estimate"))
        assert a["rare_truth"] == .125 and a["partition_truth"] == .5
    for r in table(out,"budget-comparison.csv"):
        if float(r["budget_seconds"]) < 1e-10:
            assert r["setup_exhausted"] == "2" and r["rmse"] == ""
    # One observed repetition cannot establish RMSE for a two-repetition design.
    partial = budget_report(out, manifest, [r for r in runs if r["trial"] == 0])
    assert all(r["rmse"] is None and r["bias"] is None and r["planned_trials"] == 2 for r in partial)


SUITES = {f.__name__: f for f in (edge_ranking, path_oracle, nested_networks, wall_budget)}

if __name__ == "__main__":
    with TemporaryDirectory() as folder:
        for name in sys.argv[2:] or SUITES:
            SUITES[name](Path(folder))
