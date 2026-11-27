"""Plot measured chain distances and spectral changes from a completed MixBench run."""
import argparse
import csv
from collections import defaultdict
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter

from comparison_stats import crossing, start_statistics, tables
from size_scaling import write_scaling
from study_reports import requested_reports


def read(path):
    return json.loads(path.read_text())


def save(fig, path):
    fig.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def name(r):
    label = r["method"].replace("_", " ")
    if r["edge"] is not None:
        i, j = r["edge"]
        label += f" · {r['policy']} ({i},{j})"
    return label


def mcmc_table(directory, runs):
    measured = [r for r in runs if "mcmc" in r]
    if not measured:
        return
    records = []
    for r in measured:
        baseline = next((b for b in measured if b["graph"] == r["graph"] and
                         (b["method"], b["edge"]) == ((r["method"], None) if r["edge"]
                                                      else ("max_degree", None))), None)
        m = r["mcmc"]
        setup = r.get("setup_seconds", r["construction_seconds"])
        crossover = None
        if baseline and m["ess_per_second"] and baseline["mcmc"]["ess_per_second"]:
            saved = 1 / baseline["mcmc"]["ess_per_second"] - 1 / m["ess_per_second"]
            extra = setup - baseline.get("setup_seconds", baseline["construction_seconds"])
            if saved > 0:
                crossover = max(0., extra / saved)
        records.append(dict(run=r["id"], graph=r["graph"], method=r["method"], edge=r["edge"],
                            baseline_run=baseline["id"] if baseline else None, estimate=m["estimate"],
                            exact=m["true_mean"], rmse=m["rmse"], ess=m["ess"], ess_per_second=m["ess_per_second"],
                            setup_seconds=setup, break_even_effective_draws=crossover))
    with (directory / "mcmc.csv").open("w", newline="") as output:
        writer = csv.DictWriter(output, records[0].keys())
        writer.writeheader()
        writer.writerows(records)


def method_figure(path, graph, runs, epsilon):
    groups = []
    for run in runs:
        analysis = run["analysis"]
        for i, c in enumerate(analysis["crossings"]):
            if c["epsilon"] == epsilon:
                groups.append((run, c, start_statistics(analysis, i)))
    fig, (times, rates) = plt.subplots(1, 2, sharey=True, figsize=(10, 2.3 + .5*len(groups)))
    labels = []
    for y, (run, c, stats) in enumerate(groups):
        label = name(run)
        if run["positions"] != run["states"]:
            label += " · " + c["space"]
        if run["analysis"]["scope"] == "sampled":
            label += f" · {stats['starts_count']} sampled"
        labels.append(label)
        for value, marker, color, legend in (
            (c["steps"], "o", "C0", "Maximum TV crossing"),
            (stats["start_p90"], "|", "C1", "90th percentile"),
            (stats["start_median"], "x", "black", "Median start"),
        ):
            times.plot(value if value is not None else [], y if value is not None else [],
                       marker=marker, ms=8, color=color, ls="", label=legend if y == 0 else None)
        if c["steps"] is None:
            times.plot(c["checked_through"], y, ">", color="gray")
        rate = run.get("asymptotic_rate")
        if rate is not None:
            rates.plot(rate, y, "o", color="C2")
    times.set_yticks(range(len(labels)), labels)
    times.invert_yaxis()
    times.set_xlabel(f"Steps to TV ≤ {epsilon:g}")
    rates.set_xlabel("Asymptotic log decay per step (full chain)")
    for ax in (times, rates):
        ax.grid(axis="x", alpha=.25)
        ax.set_xlim(left=0)
    times.legend(fontsize=8, loc="upper center", bbox_to_anchor=(.5, -.23), ncol=2)
    if any(c["steps"] is None for _, c, _ in groups):
        times.set_title("▷  Unreached by displayed horizon")
    fig.suptitle(f"{graph['name']} · {graph['vertices']} vertices")
    save(fig, path)


def graph_figures(folder, graph, runs):
    ok = [r for r in runs if r["status"] == "ok"]
    if not ok:
        return
    eps = ok[0]["analysis"]["crossings"][0]["epsilon"]
    prefix = graph["id"]
    baseline = {r["method"]: r for r in ok if r["edge"] is None}
    chosen = []
    for method in dict.fromkeys(r["method"] for r in ok):
        group = [r for r in ok if r["method"] == method]
        chosen.append(baseline.get(method, group[0]))
        edited = [r for r in group if r["edge"] is not None]
        if edited:
            best = min(edited, key=lambda r: crossing(r, eps, "full")["steps"]
                       if crossing(r, eps, "full")["steps"] is not None else math.inf)
            if best != chosen[-1]:
                chosen.append(best)
    method_figure(folder / f"{prefix}-methods.png", graph, chosen, eps)
    fig, axes = plt.subplots(1, 2, figsize=(12, 4), sharey=True)
    scopes = {r["analysis"]["scope"] for r in chosen}
    diagnostics = {r["analysis"]["diagnostic_start"] for r in chosen}
    for r in chosen:
        for space in dict.fromkeys(c["space"] for c in r["analysis"]["curves"]):
            points = [c for c in r["analysis"]["curves"] if c["space"] == space]
            label = name(r) + (" · position" if space == "position" else "")
            for ax, field in zip(axes, ("maximum", "diagnostic")):
                series_label = label
                if field == "maximum" and len(scopes) > 1:
                    series_label += f" · {r['analysis']['scope']} starts"
                if field == "diagnostic" and len(diagnostics) > 1:
                    series_label += f" · start {r['analysis']['diagnostic_start']}"
                ax.plot([c["step"] for c in points], [c[field] for c in points], label=series_label)
    if len(scopes) > 1:
        axes[0].set_title("Maximum over measured starts")
        axes[0].legend(fontsize=7)
    else:
        axes[0].set_title("Worst starting state" if scopes == {"all"} else "Maximum over sampled starts")
    axes[1].set_title(f"Starting state {next(iter(diagnostics))}" if len(diagnostics) == 1 else "Selected starting states")
    for ax in axes:
        ax.axhline(eps, color="0.5", ls="--", lw=1)
        ax.set_xlabel("Markov-chain transitions")
        ax.grid(alpha=0.2)
    axes[0].set_ylabel("Total variation distance")
    axes[1].legend(fontsize=7)
    fig.suptitle(f"{graph['name']} · {graph['vertices']} vertices")
    save(fig, folder / f"{prefix}-tv.png")

    first = chosen[0]
    snapshots = first["analysis"]["snapshots"]
    fig, ax = plt.subplots(figsize=(7, 4))
    for snap in snapshots:
        ax.plot(snap["distribution"], label=f"t = {snap['step']}", alpha=.8)
    ax.set(xlabel="Internal state", ylabel="Probability",
           title=f"{graph['name']} · {name(first)} · start {first['analysis']['diagnostic_start']}")
    ax.legend(fontsize=8)
    save(fig, folder / f"{prefix}-snapshots.png")

    methods = sorted({r["method"] for r in ok})
    fig, axes = plt.subplots(1, 2, figsize=(12, 4))
    for method in methods:
        group = [r for r in ok if r["method"] == method]
        for index, r in enumerate(group):
            x = r["id"]
            c = crossing(r, eps, "full")
            axes[0].scatter(x, r["lambda2_after"], color=f"C{methods.index(method)}",
                            label=method if index == 0 else None, s=15)
            if c["steps"] is not None:
                axes[1].scatter(x, c["steps"], color=f"C{methods.index(method)}", s=15)
            else:
                axes[1].scatter(x, c["checked_through"], facecolors="none", edgecolors="gray", marker="^")
    axes[0].set_ylabel("Algebraic connectivity of undirected support")
    levels = [r["lambda2_after"] for r in ok]
    if max(levels) - min(levels) < 1e-10 * max(1, max(levels)):
        axes[0].set_ylim(0, 1.4 * max(levels))
    axes[0].legend(fontsize=8)
    axes[1].set_ylabel(f"TV crossing at ε = {eps:g}")
    for ax in axes:
        ax.set_xlabel("Run ID (edge and method in summary.csv)")
        ax.ticklabel_format(useOffset=False, axis="y")
    save(fig, folder / f"{prefix}-connectivity.png")

    for method, original in baseline.items():
        t0 = crossing(original, eps, "full")["steps"]
        edited = [r for r in ok if r["method"] == method and r["edge"] is not None
                  and crossing(r, eps, "full")["steps"] is not None]
        if not edited or t0 is None:
            continue
        fig, ax = plt.subplots(figsize=(6.5, 4))
        for policy in dict.fromkeys(r["policy"] for r in edited):
            subset = [r for r in edited if r["policy"] == policy]
            ax.scatter([r["delta_lambda2"] for r in subset],
                       [t0 - crossing(r, eps, "full")["steps"] for r in subset],
                       label=policy, alpha=0.75)
        best = min(edited, key=lambda r: crossing(r, eps, "full")["steps"])
        ax.scatter(best["delta_lambda2"], t0 - crossing(best, eps, "full")["steps"],
                   marker="*", s=110, facecolor="none", edgecolor="black", label="Fastest measured edge")
        ax.set(xlabel="Change in algebraic connectivity", ylabel="Mixing steps saved",
               title=f"{graph['name']} · {method}")
        if max(abs(r["delta_lambda2"]) for r in edited) < 1e-10:
            ax.set_xlim(-.01, .01)
            ax.set_xticks([0])
        ax.legend(fontsize=8)
        save(fig, folder / f"{prefix}-{method}-edges.png")

    placements = [r for r in ok if r["method"] == "clockwise" and r["policy"] == "all"
                  and r["edge"][1] == graph["vertices"] - 1]
    if placements:
        fig, ax = plt.subplots(figsize=(6.5, 4))
        for r in placements:
            c = crossing(r, eps, "full")
            ax.scatter(r["edge"][0]+1, c["steps"] if c["steps"] is not None else c["checked_through"],
                       color="C0", marker="o" if c["steps"] is not None else "^")
        if "clockwise" in baseline:
            time = crossing(baseline["clockwise"], eps, "full")["steps"]
            if time is not None:
                ax.axhline(time, color="0.4", ls="--", label="No shortcut")
                ax.legend()
        ax.set(xlabel="Shortcut endpoint k (paper numbering)", ylabel="Mixing time in steps",
               title=f"{graph['name']} · {graph['vertices']} vertices")
        save(fig, folder / f"{prefix}-placement.png")


def scaling_figures(folder, rows):
    by_study = defaultdict(list)
    for row in rows:
        if row.get("space") in ("full", "position"):
            by_study[row["study"]].append(row)
    for index, (study, group) in enumerate(by_study.items()):
        if len({r["vertices"] for r in group}) < 2:
            continue
        epsilon = group[0]["epsilon"]
        fig, axes = plt.subplots(1, 4, figsize=(16, 4))
        series = defaultdict(list)
        for r in group:
            if r["epsilon"] == epsilon:
                series[r["method"], r["policy"], r["space"]].append(r)
        for index_in_plot, ((method, policy, space), points) in enumerate(series.items()):
            speed_metric = "vs_max_degree" if policy == "none" and method != "clockwise" else "vs_unedited"
            metrics = ("steps", speed_metric)
            if space == "full":
                metrics += ("setup_seconds", "lambda2")
            label = f"{method} · {policy}" + (" · position" if space == "position" else "")
            for ax, metric in zip(axes, metrics):
                pairs = [(p["vertices"], p.get(metric)) for p in points if p.get(metric) is not None]
                if pairs:
                    xs, ys = zip(*pairs)
                    ax.scatter(xs, ys, label=label, alpha=.65, color=f"C{index_in_plot}")
            missing = [p for p in points if p["status"] == "unresolved"]
            axes[0].scatter([p["vertices"] for p in missing], [p["checked_through"] for p in missing],
                            marker="^", facecolors="none", edgecolors="gray")
        for ax, label in zip(axes, ("TV threshold crossing (steps)", "Speedup over baseline", "Chain setup (s)",
                                   "Algebraic connectivity")):
            low, high = ax.dataLim.intervaly
            ax.set(xlabel="Vertices", ylabel=label, xscale="log",
                   yscale="log" if high >= low > 0 else "linear")
            sizes = sorted({r["vertices"] for r in group})
            ticks = [sizes[0]]
            for size in sizes[1:]:
                if size >= 1.5 * ticks[-1]:
                    ticks.append(size)
            ax.set_xticks(ticks, labels=[str(size) for size in ticks])
            ax.xaxis.set_minor_formatter(NullFormatter())
        axes[0].legend(fontsize=7)
        fig.suptitle(f"{study} · ε = {epsilon:g}")
        save(fig, folder / f"scaling-{index}.png")


def main(directory):
    manifest = read(directory / "manifest.json")
    if manifest["status"] != "complete":
        raise ValueError("Results are incomplete")
    graphs = {g["id"]: g for g in manifest["graphs"]}
    runs = [read(directory / f"run-{i}.json") for i in range(manifest["run_count"])]
    if manifest.get("kind") == "wall_budget":
        requested_reports(directory, manifest, runs, ["wall-budget"])
        return
    requested_reports(directory, manifest, runs, manifest["config"].get("reports", []))
    folder = directory / "figures"
    folder.mkdir(exist_ok=True)
    rows = tables(directory, runs, graphs)
    write_scaling(directory, [dict(r, experiment=directory.name, jobs=manifest["jobs"]) for r in rows], figures=True)
    mcmc_table(directory, runs)
    for key, graph in graphs.items():
        graph_figures(folder, graph, [r for r in runs if r["graph"] == key])
    scaling_figures(folder, rows)
    print(f"{len(rows)} summary rows and figures in {directory}")
    return rows


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    main(parser.parse_args().directory)

