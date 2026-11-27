"""Size sweeps with per-size case variation and descriptive log-log slopes."""
from collections import defaultdict
import math
from statistics import fmean, median

from comparison_stats import write_csv

# Column, axis label, whether this belongs to the internal chain only.
METRICS = (
    ("steps", "Maximum TV crossing (steps)", False),
    ("start_median", "Median start (steps)", False),
    ("start_p90", "90th-percentile start (steps)", False),
    ("observed_tv_log_decay", "Measured TV log decay / step", False),
    ("vs_unedited", "Speedup over unedited walk", False),
    ("vs_max_degree", "Speedup over maximum-degree walk", False),
    ("slem", "Nonstationary spectral radius", True),
    ("transition_lambda2", "Second transition eigenvalue", True),
    ("transition_lambda_min", "Smallest transition eigenvalue", True),
    ("spectral_gap", "Reversible spectral gap", True),
    ("absolute_spectral_gap", "Absolute spectral gap", True),
    ("asymptotic_rate", "Asymptotic log decay / step", True),
    ("spectral_time", "Spectral e-folding time (steps)", True),
    ("lambda2", "Algebraic connectivity", True),
    ("delta_lambda2", "Gain in algebraic connectivity", True),
    ("setup_seconds", "Setup (s)", True),
    ("spectrum_seconds", "Transition eigensolve (s)", True),
    ("propagation_seconds", "TV propagation (s)", True),
    ("updates_per_second", "Start-transitions / s", True),
    ("mcmc_ess_per_second", "MCMC effective draws / s", True),
    ("mcmc_rmse", "MCMC observable RMSE", True),
)
KEYS = ("experiment", "jobs", "study", "graph_definition", "cycle_parity", "method", "policy",
        "dhn_c", "space", "epsilon", "scope", "sample_count", "metric")


def quartile(values, q):
    pos = (len(values)-1)*q
    i = int(pos)
    return values[i] + (pos-i)*(values[math.ceil(pos)]-values[i])


def scaling_tables(rows):
    # Only compare sizes within a study and a fixed generator configuration.
    sizes = defaultdict(set)
    skipped = defaultdict(set)
    for r in rows:
        study = (r.get("experiment"), r["study"], r["graph_definition"])
        sizes[study].add(r["vertices"])
        if r["run_status"] != "ok":
            skipped[(*study, r["method"], r["policy"])].add(r["vertices"])
    groups = defaultdict(lambda: defaultdict(list))
    for r in rows:
        if len(sizes[r.get("experiment"), r["study"], r["graph_definition"]]) < 2:
            continue
        if r["run_status"] != "ok":
            continue  # skipped run and reason are retained in the case and policy tables
        header = dict(r, cycle_parity=("odd" if r["vertices"] % 2 else "even") if r["family"] == "cycle" else "",
                      sample_count=r["starts_count"] if r["scope"] == "sampled" else None)
        for metric, _, internal in METRICS:
            if internal and r["space"] != "full":
                continue
            header["metric"] = metric
            groups[tuple(header.get(k) for k in KEYS)][r["vertices"]].append(r)
    points, fits = [], []
    for key, by_size in groups.items():
        series = dict(zip(KEYS, key))
        metric = series["metric"]
        usable, omitted = [], []
        for n, cases in sorted(by_size.items()):
            values = sorted(r[metric] for r in cases if r.get(metric) is not None)
            complete = len(values) == len(cases)
            middle = median(values) if complete else None
            points.append(dict(series, vertices=n, cases=len(cases), available=len(values),
                               graph_realizations=len({r["graph"] for r in cases}),
                               q25=quartile(values, .25) if complete else None, median=middle,
                               q75=quartile(values, .75) if complete else None))
            if complete and middle > 0:
                usable.append((n, middle))
            else:
                omitted.append(n)
        fit = dict(series, sizes_used=" ".join(str(n) for n, _ in usable),
                   unavailable_or_nonpositive_sizes=" ".join(map(str, omitted)),
                   size_count=len(usable), exponent=None, log_r_squared=None)
        unavailable = skipped[series["experiment"], series["study"], series["graph_definition"],
                              series["method"], series["policy"]]
        parity = series["cycle_parity"]
        fit["skipped_sizes"] = " ".join(str(n) for n in sorted(unavailable)
                                        if not parity or (n % 2 == 1) == (parity == "odd"))
        if len(usable) >= 3:
            x, y = zip(*[(math.log(n), math.log(v)) for n, v in usable])
            mx, my = fmean(x), fmean(y)
            xx = sum((t-mx)**2 for t in x)
            yy = sum((t-my)**2 for t in y)
            xy = sum((a-mx)*(b-my) for a, b in zip(x, y))
            fit.update(exponent=xy/xx, log_r_squared=min(1., xy**2/(xx*yy)) if yy else None)
        fits.append(fit)
    return points, fits


def write_scaling(directory, rows, figures=False):
    points, fits = scaling_tables(rows)
    write_csv(directory / "size-statistics.csv", points,
              (*KEYS, "vertices", "cases", "available", "graph_realizations", "q25", "median", "q75"))
    write_csv(directory / "size-fits.csv", fits,
              (*KEYS, "sizes_used", "unavailable_or_nonpositive_sizes", "size_count",
               "exponent", "log_r_squared", "skipped_sizes"))
    if figures and points:
        plot_scaling(directory / "size-figures", points)


def plot_scaling(folder, points):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    folder.mkdir(exist_ok=True)
    panels = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
    for p in points:
        group = (p["experiment"], p["study"], p["graph_definition"], p["epsilon"], p["jobs"])
        series = (p["method"], p["policy"], p["space"], p["scope"], p["sample_count"],
                  p["cycle_parity"], p["dhn_c"])
        panels[group][p["metric"]][series].append(p)
    index = []
    for group, metrics in panels.items():
        for section, fields in (
            ("mixing", ("steps", "start_median", "start_p90", "observed_tv_log_decay")),
            ("spectrum", ("slem", "absolute_spectral_gap", "asymptotic_rate", "spectral_time")),
            ("eigenvalues", ("transition_lambda2", "transition_lambda_min", "lambda2", "delta_lambda2")),
            ("speedup", ("vs_unedited", "vs_max_degree")),
            ("cost", ("setup_seconds", "spectrum_seconds", "propagation_seconds", "updates_per_second")),
            ("sampling", ("mcmc_rmse", "mcmc_ess_per_second")),
        ):
            data = [(metric, label) for metric, label, _ in METRICS if metric in fields
                    and any(p["median"] is not None for line in metrics[metric].values() for p in line)]
            if not data:
                continue
            fig, axes = plt.subplots(1, len(data), squeeze=False, figsize=(4*len(data), 4.6))
            handles = {}
            for ax, (metric, title) in zip(axes.flat, data):
                medians = []
                for series, line in metrics[metric].items():
                    if all(p["median"] is None for p in line):
                        continue
                    method, policy, space, scope, samples, parity, c = series
                    label = " · ".join(str(v) for v in (method, policy, space, parity,
                                       f"{samples} starts" if scope == "sampled" else None,
                                       f"c={c:g}" if c is not None else None) if v)
                    xs = [p["vertices"] for p in line]
                    def column(key):
                        return [p[key] if p[key] is not None else math.nan for p in line]
                    median_line = column("median")
                    medians.extend(v for v in median_line if math.isfinite(v))
                    color = handles[label].get_color() if label in handles else f"C{len(handles)}"
                    artist, = ax.plot(xs, median_line, marker="o", color=color)
                    handles[label] = artist
                    ax.fill_between(xs, column("q25"), column("q75"), color=color, alpha=.12)
                ax.set(xlabel="Graph vertices", ylabel=title, xscale="log")
                if min(medians) > 0 and max(medians)/min(medians) >= 10:
                    ax.set_yscale("log")
                ax.grid(alpha=.2)
            fig.suptitle(f"{group[1]} · ε = {group[3]:g}")
            legend = fig.legend(handles.values(), handles.keys(), loc="upper center",
                                bbox_to_anchor=(.5, .01), ncol=min(4, len(handles)), fontsize=8)
            fig.tight_layout()
            filename = f"{len(index):03d}-{section}.png"
            fig.savefig(folder / filename, dpi=140, bbox_inches="tight", bbox_extra_artists=(legend,))
            plt.close(fig)
            index.append(dict(file=filename, experiment=group[0], study=group[1], graph_definition=group[2],
                              epsilon=group[3], jobs=group[4], panels=section))
    write_csv(folder / "index.csv", index)
