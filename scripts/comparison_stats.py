"""Comparable transition, topology and finite-time statistics from MixBench runs."""
import csv
import math
import json
from statistics import fmean, pstdev


def crossing(run, epsilon, space):
    return next(c for c in run["analysis"]["crossings"] if c["epsilon"] == epsilon and c["space"] == space)


def write_csv(path, rows, empty_fields=()):
    fields = list(dict.fromkeys(key for row in rows for key in row)) if rows else empty_fields
    with path.open("w", newline="") as output:
        writer = csv.DictWriter(output, fields)
        writer.writeheader()
        writer.writerows(rows)


def start_statistics(analysis, index):
    hits = [s["first_crossings"][index] for s in analysis["per_start"]]
    known = sorted(t for t in hits if t >= 0)
    count = len(hits)

    def quantile(q):
        rank = q * (count - 1)
        hi = math.ceil(rank)
        if hi >= len(known):
            return None
        lo = math.floor(rank)
        return known[lo] + (rank - lo) * (known[hi] - known[lo])

    diagnostic = hits[analysis["starts"].index(analysis["diagnostic_start"])]
    return dict(starts_count=count, starts_resolved=len(known), starts_unresolved=count-len(known),
                start_min=quantile(0), start_p25=quantile(.25), start_median=quantile(.5),
                start_p75=quantile(.75), start_p90=quantile(.9), start_max=quantile(1),
                start_mean=fmean(known) if len(known) == count else None,
                start_sd=pstdev(known) if len(known) == count else None,
                diagnostic_steps=diagnostic if diagnostic >= 0 else None)


def comparable_space(run, other, space):
    if other is None:
        return None
    a, b = run["analysis"], other["analysis"]
    starts, others = a["starts"], b["starts"]
    mapping, other_mapping = run.get("state_to_position"), other.get("state_to_position")
    if mapping == other_mapping and run["states"] == other["states"]:
        return space if starts == others else None
    if space == "position" and mapping and other_mapping == list(range(run["vertices"])):
        if a["scope"] == b["scope"] == "all" or (len(starts) == 1 and others == [mapping[starts[0]]]):
            return "full"
    return None


def summary_rows(runs, graphs):
    lookup = {(r["graph"], r["method"], tuple(r["edge"] or [])): r
              for r in runs if r["status"] == "ok"}
    by_id = {r["id"]: r for r in runs}
    by_run = {}
    rows = []
    for r in runs:
        data = graphs[r["graph"]]
        generator = {k: v for k, v in data["generator"].items() if k != "sizes"}
        common = dict(family=generator.get("family", "dataset"),
                      graph_definition=json.dumps(generator, sort_keys=True, separators=(",", ":")),
                      run=r["id"], graph=r["graph"], study=data["name"], replicate=data["replicate"],
                      vertices=r["vertices"], support_edges=len(data["edges"]) + int(r["edge"] is not None),
                      method=r["method"], policy=r["policy"], edge=r["edge"], dhn_c=r.get("dhn_c"),
                      states=r.get("states"), positions=r.get("positions"),
                      reversible=r.get("reversible"), run_status=r["status"],
                      lambda2=r["lambda2_after"], lambda2_before=r["lambda2_before"],
                      lambda3_before=r["lambda3_before"], lambda3_after=r.get("lambda3_after"),
                      delta_lambda2=r["delta_lambda2"], normalized_gain=r["normalized_gain"])
        for key in ("transition_lambda2", "transition_lambda_min", "dominant_eigenvalue_real",
                    "dominant_eigenvalue_imag", "slem", "spectral_gap", "absolute_spectral_gap",
                    "asymptotic_rate", "spectral_time", "spectral_status", "start_seed",
                    "edge_selection_seconds", "construction_seconds", "connectivity_seconds",
                    "setup_seconds", "spectrum_seconds", "propagation_seconds"):
            common[key] = r.get(key)
        for key in ("draws", "replicates", "estimate", "true_mean", "rmse", "ess", "ess_per_second", "seconds"):
            common["mcmc_" + key] = r.get("mcmc", {}).get(key)
        if r["status"] != "ok":
            rows.append(dict(common, status=r["status"]))
            continue
        a = r["analysis"]
        points = {(p["space"], p["step"]): p["maximum"] for p in a["curves"]}
        for index, c in enumerate(a["crossings"]):
            space = c["space"]
            end = c["steps"] if c["steps"] is not None else c["checked_through"]
            initial, final = points[space, 0], points[space, end]
            decay = (math.log(initial) - math.log(final))/end if end and initial > 0 and final > 1e-14 else None
            logical_updates = len(a["starts"]) * c["checked_through"]
            row = dict(common, **c, logical_updates=logical_updates,
                       updates_per_second=logical_updates / r["propagation_seconds"]
                       if logical_updates and r.get("propagation_seconds") else None,
                       scope=a["scope"], period=a["period"],
                       diagnostic_start=a["diagnostic_start"], tv_rate_end_step=end,
                       tv_initial_max=initial, tv_end_max=final, observed_tv_log_decay=decay,
                       **start_statistics(a, index))
            by_run[r["id"], space, c["epsilon"]] = row
            rows.append(row)

    for row in rows:
        if row["run_status"] != "ok":
            continue
        r = by_id[row["run"]]
        controls = {
            "vs_unedited": lookup.get((r["graph"], r["method"], ())),
            "vs_max_degree": lookup.get((r["graph"], "max_degree", tuple(r["edge"] or []))),
        }
        for label, control in controls.items():
            space = comparable_space(r, control, row["space"])
            other = by_run.get((control["id"], space, row["epsilon"])) if space else None
            row[label + "_run"] = control["id"] if other else None
            for suffix, metric in (("", "steps"), ("_median", "start_median"), ("_p90", "start_p90")):
                value = row[metric]
                row[label + suffix] = other[metric]/value if other and other[metric] is not None and value else None
    return rows


def tables(directory, runs, graphs):
    rows = summary_rows(runs, graphs)
    write_csv(directory / "summary.csv", rows)
    return rows
