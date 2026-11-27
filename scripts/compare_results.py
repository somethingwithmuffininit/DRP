"""Combine completed experiment directories into case and policy comparison tables."""
import argparse
from collections import defaultdict
import json
from pathlib import Path
from statistics import median

from comparison_stats import summary_rows, write_csv
from size_scaling import write_scaling
from study_reports import budget_rows


def policy_rows(rows):
    groups = defaultdict(list)
    keys = ("experiment", "study", "vertices", "method", "policy", "dhn_c", "space", "epsilon", "scope")
    for row in rows:
        groups[tuple(row.get(k) for k in keys)].append(row)
    result = []
    for key, group in groups.items():
        stats = dict(zip(keys, key))
        stats.update(cases=len(group), graph_realizations=len({r["graph"] for r in group}),
                     resolved=sum(r["status"] == "resolved" for r in group),
                     unresolved=sum(r["status"] == "unresolved" for r in group),
                     periodic=sum(r["status"] == "periodic" for r in group),
                     skipped=sum(r["run_status"] != "ok" for r in group))
        # Give each case one vote, including each seeded edge draw. A missing
        # observation makes this aggregate unavailable; preserve counts above.
        for metric in ("steps", "start_median", "start_p90", "start_mean", "slem", "asymptotic_rate",
                       "observed_tv_log_decay", "lambda2", "delta_lambda2", "vs_unedited",
                       "vs_unedited_median", "vs_max_degree", "vs_max_degree_median", "setup_seconds"):
            values = [r.get(metric) for r in group]
            stats["median_" + metric] = median(values) if all(v is not None for v in values) else None
        result.append(stats)
    return result


def combine(output, directories, figures=False):
    rows, timed = [], []
    if len({p.name for p in directories}) != len(directories):
        raise ValueError("Experiment directory names must be unique")
    for folder in directories:
        manifest = json.loads((folder / "manifest.json").read_text())
        if manifest["status"] != "complete":
            raise ValueError(f"Incomplete experiment: {folder}")
        runs = [json.loads((folder / f"run-{i}.json").read_text()) for i in range(manifest["run_count"])]
        if manifest.get("kind") == "wall_budget":
            timed.extend(dict(experiment=folder.name, **r) for r in budget_rows(manifest, runs))
            continue
        graphs = {g["id"]: g for g in manifest["graphs"]}
        rows.extend(dict(experiment=folder.name, jobs=manifest["jobs"], **r)
                    for r in summary_rows(runs, graphs))
    if not rows and not timed:
        raise ValueError("No completed experiments found")
    output.mkdir(exist_ok=True, parents=True)
    if timed:
        write_csv(output / "budget-trials.csv", timed)
    if not rows:
        return
    write_csv(output / "cases.csv", rows)
    write_csv(output / "policies.csv", policy_rows(rows))
    write_scaling(output, rows, figures=figures)
    print(f"{len(directories)} experiments, {len(rows)} comparison rows in {output}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("directories", type=Path, nargs="*",
                        help="defaults to the experiment directories directly under OUTPUT")
    parser.add_argument("--figures", action="store_true", help="also render graph-size comparisons")
    args = parser.parse_args()
    combine(args.output, args.directories or sorted(p.parent for p in args.output.glob("*/manifest.json")), figures=args.figures)
