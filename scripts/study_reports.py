"""Reports for exhaustive edge rankings, path controls and timed sampling."""
import argparse
from collections import defaultdict
import json
import math
from pathlib import Path
from statistics import fmean, median

from comparison_stats import summary_rows, write_csv


def ranks(values, *, descending=False, tolerance=0.):
    """Competition ranks; tie groups are measured from the group's leading value."""
    ordered = sorted((key for key in values if values[key] is not None),
                     key=values.get, reverse=descending)
    result = {}
    for index, key in enumerate(ordered, 1):
        if index == 1 or abs(values[key] - anchor) > tolerance:
            rank, anchor = index, values[key]
        result[key] = rank
    return result


def edge_report(directory, runs, graphs):
    groups, controls = defaultdict(list), {}
    for row in summary_rows(runs, graphs):
        if row["run_status"] != "ok":
            continue
        key = (row["graph"], row["method"], row["space"], row["epsilon"])
        if row["edge"] is None:
            controls[key] = row
        else:
            groups[key].append(row)
    ranked, starts, agreements = [], [], []
    by_id = {r["id"]: r for r in runs}
    for (gid, method, space, eps), rows in groups.items():
        candidates = {tuple(r["edge"]): r for r in rows if r["policy"] == "all"}
        if not candidates:
            continue
        graph = graphs[gid]
        n = graph["vertices"]
        if len(candidates) != n*(n-1)//2 - len(graph["edges"]):
            raise ValueError("Edge ranking requires exhaustive coverage of all missing edges")
        meta = dict(graph=gid, study=graph["name"], vertices=n, method=method, space=space, epsilon=eps)
        gain = {e: r["delta_lambda2"] for e, r in candidates.items()}
        gain_rank = ranks(gain, descending=True, tolerance=1e-10)
        scores = {}
        for metric in ("steps", "start_median"):
            vals = {e: r[metric] for e, r in candidates.items()}
            complete = all(v is not None for v in vals.values())
            scores[metric] = (ranks(vals), min(vals.values())) if complete else ({}, None)

        def hits(row):
            a = by_id[row["run"]]["analysis"]
            i = next(i for i, c in enumerate(a["crossings"]) if c["space"] == space and c["epsilon"] == eps)
            return {s["start"]: s["first_crossings"][i] for s in a["per_start"]}

        control = controls.get((gid, method, space, eps))
        baseline = hits(control) if control else {}
        first = {e: hits(r) for e, r in candidates.items()}
        ids = next(iter(first.values())).keys()
        if any(times.keys() != ids for times in first.values()):
            raise ValueError("Edge ranking requires matched starting states")
        best = {}
        for s in ids:
            values = [times[s] for times in first.values()]
            best[s] = min(values) if min(values) >= 0 else None
        best_gain = {e for e,v in gain_rank.items() if v == 1}
        for metric, (ranking, optimum) in scores.items():
            fastest = {e for e,v in ranking.items() if v == 1}
            times = [candidates[e][metric] for e in best_gain]
            agreements.append(dict(meta, metric=metric, candidates=len(candidates), fastest_steps=optimum,
                                   connectivity_best_ties=len(best_gain), fastest_ties=len(fastest) if ranking else None,
                                   shared_best_edges=len(best_gain & fastest) if ranking else None,
                                   best_extra_among_connectivity_ties=min(times)-optimum if ranking else None,
                                   worst_extra_among_connectivity_ties=max(times)-optimum if ranking else None))
        for row in rows:
            edge = tuple(row["edge"])
            record = dict(meta, run=row["run"], policy=row["policy"], edge=edge,
                          candidates=len(candidates), scope=row["scope"], status=row["status"],
                          delta_lambda2=gain[edge], connectivity_rank=gain_rank[edge],
                          connectivity_best_ties=sum(v == 1 for v in gain_rank.values()))
            for metric, (ranking, optimum) in scores.items():
                value = row[metric]
                record.update({metric: value, metric+"_rank": ranking.get(edge),
                               metric+"_optimum": optimum,
                               metric+"_best_ties": sum(v == 1 for v in ranking.values()) if ranking else None,
                               metric+"_extra_steps": value-optimum if optimum is not None and value is not None else None})
            ranked.append(record)
            if row["policy"] != "all":
                times = hits(row)
                if times.keys() != ids:
                    raise ValueError("Heuristic uses different starts")
                for s, t in times.items():
                    starts.append(dict(meta, run=row["run"], policy=row["policy"], edge=edge, start=s,
                                       steps=t if t >= 0 else None, fastest_steps=best[s],
                                       saved_vs_unedited=baseline[s]-t if t >= 0 and baseline.get(s, -1) >= 0 else None,
                                       extra_steps=t-best[s] if t >= 0 and best[s] is not None else None))
    if not ranked:
        raise ValueError("No exhaustive missing-edge runs in this experiment")
    write_csv(directory/"edge-ranks.csv", ranked)
    write_csv(directory/"edge-starts.csv", starts)
    write_csv(directory/"edge-agreement.csv", agreements)
    return ranked


def path_report(directory, runs, graphs, figures=True):
    rows, slopes = [], []
    grouped = defaultdict(dict)
    for r in runs:
        if r["status"] == "ok" and r["edge"] is None:
            grouped[r["graph"]][r["method"]] = r
    for gid, methods in grouped.items():
        if not {"dhn", "max_degree"} <= methods.keys():
            continue
        graph = graphs[gid]
        n = graph["vertices"]
        if graph["generator"].get("family") != "path":
            continue
        expected_gap = 4*math.sin(math.pi/(2*n))**2
        expected_eig = math.cos(math.pi/n)
        reversible = methods["max_degree"]
        eig_error = abs(reversible["transition_lambda2"] - expected_eig)
        for method, space, exponent in (("dhn", "position", 1), ("max_degree", "full", 2)):
            r = methods[method]
            a = r["analysis"]
            start = a["diagnostic_start"]
            if r["state_to_position"][start] != 0 or (method == "dhn" and start != 1):
                raise ValueError("Path control requires (+, first position) and the matching reversible endpoint")
            diagnostic = next(s for s in a["per_start"] if s["start"] == start)
            lap_error = abs(r["lambda2_after"] - expected_gap)
            for i, c in enumerate(a["crossings"]):
                if c["space"] != space:
                    continue
                t = diagnostic["first_crossings"][i]
                rows.append(dict(graph=gid, study=graph["name"], vertices=n, method=method, space="position",
                                 epsilon=c["epsilon"], internal_start=start, checked_through=c["checked_through"],
                                 steps=t if t >= 0 else None, expected_exponent=exponent,
                                 scaled_steps=t/n**exponent if t >= 0 else None,
                                 laplacian_lambda2=r["lambda2_after"], expected_laplacian_lambda2=expected_gap,
                                 laplacian_abs_error=lap_error, reversible_lambda2=reversible["transition_lambda2"],
                                 expected_reversible_lambda2=expected_eig, reversible_abs_error=eig_error,
                                 oracle_pass=max(eig_error, lap_error) <= 1e-10))
    if not rows:
        raise ValueError("No matched path controls")
    series = defaultdict(list)
    for r in rows:
        series[r["study"], r["method"], r["epsilon"]].append(r)
    for (study, method, eps), values in series.items():
        values.sort(key=lambda r: r["vertices"])
        if len({r["vertices"] for r in values}) != len(values):
            raise ValueError("Path oracle requires one deterministic graph per size")
        for drop in range(len(values)-2):
            subset = values[drop:]
            slope = None
            if all(r["steps"] is not None and r["steps"] > 0 for r in subset):
                xs = [math.log(r["vertices"]) for r in subset]
                ys = [math.log(r["steps"]) for r in subset]
                x0, y0 = fmean(xs), fmean(ys)
                slope = sum((x-x0)*(y-y0) for x,y in zip(xs,ys))/sum((x-x0)**2 for x in xs)
            slopes.append(dict(study=study, method=method, epsilon=eps, dropped_smallest=drop,
                               sizes=[r["vertices"] for r in subset], exponent=slope))
    write_csv(directory/"path-controls.csv", rows)
    write_csv(directory/"path-sensitivity.csv", slopes)
    if figures:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(1, 2, figsize=(10, 4))
        for (study, method, eps), values in series.items():
            label = f"{method} · ε={eps:g}"
            axes[0].plot([r["vertices"] for r in values], [r["scaled_steps"] for r in values], "o-", label=label)
            fits = [r for r in slopes if (r["study"], r["method"], r["epsilon"]) == (study, method, eps)]
            axes[1].plot([r["sizes"][0] for r in fits], [r["exponent"] for r in fits], "o-", label=label)
        axes[0].set(xlabel="Path vertices", ylabel="Steps / n (DHN) or steps / n² (reversible)")
        axes[1].set(xlabel="Smallest path included", ylabel="Log-log exponent")
        for ax in axes:
            ax.grid(alpha=.2)
        axes[0].legend(fontsize=8)
        fig.tight_layout()
        fig.savefig(directory/"path-scaling.png", dpi=150)
        plt.close(fig)
    if any(not r["oracle_pass"] for r in rows):
        raise ValueError("Path eigenvalue oracle failed; inspect path-controls.csv")
    return rows


def budget_rows(manifest, runs):
    graphs = {g["id"]: g for g in manifest["graphs"]}
    rows = []
    for r in runs:
        graph = graphs[r["graph"]]
        checkpoints = r.get("checkpoints", [dict(budget_seconds=b, draws=0, status=r["status"])
                                              for b in manifest["config"]["wall_budget"]["seconds"]])
        for check in checkpoints:
            for key in ("partition", "rare"):
                estimate = check.get(key+"_estimate")
                truth = r.get(key+"_truth")
                observable = ("left_clique" if graph["generator"]["family"] == "barbell" else "first_half"
                              ) if key == "partition" else "lowest_degree_vertex"
                rows.append(dict(graph=r["graph"], study=graph["name"], vertices=graph["vertices"], trial=r["trial"],
                                 method=r["method"], policy=r["policy"], edge=r.get("edge"), seed=r.get("seed"),
                                 observable=observable, rare_vertex=r.get("rare_vertex"), jobs=manifest["jobs"],
                                 budget_seconds=check["budget_seconds"], setup_seconds=r.get("setup_seconds"),
                                 status=check["status"], draws=check["draws"], truth=truth, estimate=estimate,
                                 error=estimate-truth if estimate is not None else None))
    return rows


def budget_report(directory, manifest, runs):
    rows = budget_rows(manifest, runs)
    groups = defaultdict(list)
    keys = ("graph", "study", "vertices", "method", "policy", "budget_seconds", "observable")
    for row in rows:
        groups[tuple(row[k] for k in keys)].append(row)
    summary = []
    for key, group in groups.items():
        errors = [r["error"] for r in group if r["error"] is not None]
        setup = [r["setup_seconds"] for r in group if r["setup_seconds"] is not None]
        planned = manifest["config"]["wall_budget"]["repeats"]
        full_trials = {r["trial"] for r in group} == set(range(planned))
        complete = full_trials and len(group) == planned and len(errors) == planned
        summary.append(dict(zip(keys, key), trials=len(group), planned_trials=planned, estimates=len(errors),
                            setup_exhausted=sum(r["status"] == "setup_exhausted" for r in group),
                            skipped=sum(r["status"] == "skipped_sdp_limit" for r in group),
                            median_setup_seconds=median(setup) if setup else None,
                            median_draws=median(r["draws"] for r in group),
                            bias=fmean(errors) if complete else None,
                            rmse=math.sqrt(fmean(e**2 for e in errors)) if complete else None))
    write_csv(directory/"budget-trials.csv", rows)
    write_csv(directory/"budget-comparison.csv", summary)
    return summary


def nested_report(directory, runs, graphs):
    cases = summary_rows(runs, graphs)
    keys = ("study", "vertices", "method", "policy", "space", "epsilon", "scope")
    groups = defaultdict(list)
    for r in cases:
        groups[tuple(r.get(k) for k in keys), r["graph"]].append(r)
    metrics = ("steps", "start_median", "lambda2", "delta_lambda2", "vs_unedited", "vs_max_degree")
    paired = []
    cohorts = defaultdict(list)
    for (key, gid), rows in groups.items():
        graph = graphs[gid]
        row = dict(zip(keys, key), graph=gid, replicate=graph["replicate"],
                   parent_root=graph["parent_vertex_ids"][0], draws=len(rows),
                   run_status=rows[0]["run_status"])
        for metric in metrics:
            values = [r.get(metric) for r in rows]
            row[metric] = median(values) if all(x is not None for x in values) else None
        paired.append(row)
        cohorts[key].append(row)
    variation = []
    for key, rows in cohorts.items():
        for metric in metrics:
            values = [r[metric] for r in rows if r[metric] is not None]
            full = len(values) == len(rows)
            variation.append(dict(zip(keys, key), metric=metric, realizations=len(rows), available=len(values),
                                  minimum=min(values) if full else None, median=median(values) if full else None,
                                  maximum=max(values) if full else None))
    write_csv(directory/"nested-realizations.csv", paired)
    write_csv(directory/"nested-variation.csv", variation)


def link_policy_report(directory, runs, graphs):
    rows = summary_rows(runs, graphs)
    by_crossing = {(r["run"], r["space"], r["epsilon"]): r for r in rows if "space" in r}
    keys = ("study", "vertices", "method", "policy", "space", "epsilon", "scope")
    groups = defaultdict(list)
    for row in rows:
        if row["edge"] is None or "space" not in row:
            continue
        control = by_crossing.get((row["vs_unedited_run"], row["space"], row["epsilon"]))
        saved = (control["steps"] - row["steps"]
                 if control and control["steps"] is not None and row["steps"] is not None else None)
        groups[tuple(row[key] for key in keys)].append((row["normalized_gain"], saved))
    summary = []
    for key, pairs in groups.items():
        gains = [gain for gain, _ in pairs if gain is not None]
        savings = [saved for _, saved in pairs if saved is not None]
        summary.append(dict(zip(keys, key), trials=len(pairs), scored=len(gains),
                            paired_crossings=len(savings),
                            mean_normalized_gain=fmean(gains) if gains else None,
                            mean_steps_saved=fmean(savings) if savings else None,
                            slower_trials=sum(saved < 0 for saved in savings)))
    write_csv(directory / "link-policy-summary.csv", summary)


def requested_reports(directory, manifest, runs, requested, figures=True):
    graphs = {g["id"]: g for g in manifest["graphs"]}
    for name in requested:
        if name == "edge-ranking":
            edge_report(directory, runs, graphs)
        elif name == "path-oracle":
            path_report(directory, runs, graphs, figures)
        elif name == "nested-networks":
            nested_report(directory, runs, graphs)
        elif name == "wall-budget":
            budget_report(directory, manifest, runs)
        elif name == "link-policies":
            link_policy_report(directory, runs, graphs)
        else:
            raise ValueError(f"Unknown report {name}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("--report", action="append",
                        choices=("edge-ranking", "path-oracle", "wall-budget", "nested-networks", "link-policies"))
    args = parser.parse_args()
    manifest = json.loads((args.directory/"manifest.json").read_text())
    if manifest["status"] != "complete":
        parser.error("Experiment is incomplete")
    reports = args.report or manifest["config"].get("reports", [])
    if manifest.get("kind") == "wall_budget" and not args.report:
        reports = ["wall-budget"]
    if not reports:
        parser.error("Choose a report for this experiment")
    runs = [json.loads((args.directory/f"run-{i}.json").read_text()) for i in range(manifest["run_count"])]
    requested_reports(args.directory, manifest, runs, reports)
