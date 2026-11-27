"""Verify the reported speedups, figure labels, and rendered files on actual CLI output."""
import csv
import json
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from comparison_stats import start_statistics
from compare_results import combine
from size_scaling import scaling_tables
import plot_results as plots


def main(folder):
    # Known quadratic law across three usable sizes, two edge draws per size.
    # The fourth size has a censored draw. The fifth exceeds the solver limit.
    rows = []
    for n in (8, 16, 32, 64):
        for multiplier in (.5, 1.5):
            rows.append(dict(graph="g", graph_definition="path", family="path", vertices=n,
                             study="test", method="fmmc", policy="random", jobs=1,
                             space="full", scope="all", epsilon=.25, starts_count=n, run_status="ok",
                             steps=None if n == 64 and multiplier > 1 else multiplier*n*n))
    rows.append(dict(rows[-1], vertices=128, run_status="skipped_sdp_limit"))
    sizes, fits = scaling_tables(rows)
    fit = next(f for f in fits if f["metric"] == "steps")
    assert abs(fit["exponent"] - 2) < 1e-12 and abs(fit["log_r_squared"] - 1) < 1e-12
    assert fit["unavailable_or_nonpositive_sizes"] == "64" and fit["skipped_sizes"] == "128"
    point = next(p for p in sizes if p["metric"] == "steps" and p["vertices"] == 8)
    assert (point["q25"], point["median"], point["q75"]) == (48, 64, 80)
    a = dict(starts=list(range(5)), diagnostic_start=3,
             per_start=[dict(first_crossings=[v]) for v in (3, 1, 5, -1, -1)])
    stats = start_statistics(a, 0)
    assert (stats["start_min"], stats["start_median"], stats["starts_resolved"]) == (1, 5, 3)
    assert stats["start_p90"] is None and stats["start_mean"] is None and stats["diagnostic_steps"] is None
    config = dict(studies=[dict(name="Paths", graph=dict(family="path", sizes=[4, 7]),
                                methods=["max_degree", "dhn"], policies=[dict(policy="none")])],
                  epsilons=[.20001], max_steps=128, all_starts_limit=8, sample_starts=2)
    path = folder / "config.json"
    for scenario, overrides in (("mixed", {}), ("zero", dict(max_steps=0, epsilons=[.99])),
                                ("matched", dict(all_starts_limit=0, sample_starts=1)),
                                ("unmatched", dict(all_starts_limit=0, sample_starts=1, diagnostic_start=1))):
        path.write_text(json.dumps(dict(config, **overrides)))
        out = folder / scenario
        subprocess.run([sys.argv[1], "run", path, "--jobs", "2", "--output", out], check=True)
        runs = [plots.read(f) for f in sorted(out.glob("run-*.json"))]
        graphs = {g["id"]: g for g in plots.read(out/"manifest.json")["graphs"]}
        if scenario in ("mixed", "zero"):
            figures = {}
            save = plots.save

            def capture(fig, dest):
                figures[dest.name] = fig
                save(fig, dest)

            with patch.object(plots, "save", capture):
                plots.main(out)
            for name in ("g0-tv.png", "g1-tv.png", "g1-connectivity.png", "g1-methods.png", "scaling-0.png"):
                assert (out/"figures"/name).read_bytes().startswith(b"\x89PNG\r\n\x1a\n")
            max_panel = figures["g1-tv.png"].axes[0]
            assert "worst" not in max_panel.get_title().lower(), max_panel.get_title()
            assert any("sampled" in line.get_label() for line in max_panel.lines)
            if scenario == "zero":
                axis = figures["scaling-0.png"].axes[0]
                assert axis.get_yscale() == "linear" and axis.get_ylim()[0] <= 0 <= axis.get_ylim()[1]
                continue
            speeds = figures["scaling-0.png"].axes[1].collections
            projected = [s for s in speeds if "dhn" in s.get_label() and "position" in s.get_label()]
            assert len(projected) == 1 and list(projected[0].get_offsets()[:, 0]) == [4]
        else:
            plots.tables(out, runs, graphs)
        with (out/"summary.csv").open() as f:
            rows = list(csv.DictReader(f))
        for n in ("4", "7"):
            base, full, pos = [row for row in rows if row["vertices"] == n]
            assert base["method"] == "max_degree" and full["space"] == "full" and pos["space"] == "position"
            assert full["vs_max_degree"] == ""
            comparable = scenario == "matched" or (scenario == "mixed" and n == "4")
            if comparable:
                assert float(pos["vs_max_degree"]) == int(base["steps"])/int(pos["steps"])
            else:
                assert pos["vs_max_degree"] == ""

    combine(folder, [folder/"mixed", folder/"zero"])
    with (folder/"cases.csv").open() as f:
        cases = list(csv.DictReader(f))
    assert len(cases) == 12 and len({row["experiment"] for row in cases}) == 2
    assert all(row["observed_tv_log_decay"] == "" for row in cases if row["experiment"] == "zero")


if __name__ == "__main__":
    with TemporaryDirectory() as temp:
        main(Path(temp))
