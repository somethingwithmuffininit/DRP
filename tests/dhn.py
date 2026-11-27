"""Independent dense checks of the public DHN experiment and trajectory outputs."""
import json
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory

import numpy as np

from cli import measured


def physical_walk(n, c):
    # Enumerate each heading along the path. At either wall an advance reflects
    # in place. A flip chooses the opposite heading at the same destination.
    right = list(range(1, n + 1))
    left = [0, *range(2 * n - 1, n, -1)]
    position = np.empty(2 * n, dtype=int)
    p = np.zeros((2 * n, 2 * n))
    for j in range(n):
        position[[right[j], left[j]]] = j
        for start, advance, flip in (
            (right[j], right[j+1] if j+1 < n else left[j], left[j+1] if j+1 < n else right[j]),
            (left[j], left[j-1] if j else right[j], right[j-1] if j else left[j]),
        ):
            p[start, advance] = 1 - c/n
            p[start, flip] = c/n
    return p, position


def study(n, c):
    return dict(name=f"path-{n}-{c}", graph=dict(family="path", sizes=[n]),
                methods=["dhn"], policies=[dict(policy="none")], dhn_c=c)


def run(folder, name, config, jobs=4, failure=None):
    path = folder / f"{name}.json"
    path.write_text(json.dumps(config))
    out = folder / name
    proc = subprocess.run([sys.argv[1], "run", path, "--jobs", str(jobs), "--output", out],
                          capture_output=True, text=True, timeout=120)
    assert (out / "manifest.json").is_file(), proc.stderr
    manifest = json.loads((out / "manifest.json").read_text())
    if failure:
        assert proc.returncode and manifest["status"] == "failed" and failure in proc.stderr, proc.stderr
        assert not list(out.glob("run-*.json"))
        return []
    assert proc.returncode == 0 and manifest["status"] == "complete", proc.stderr
    return [json.loads((out / f"run-{i}.json").read_text()) for i in range(manifest["run_count"])]


def near(actual, expected):
    np.testing.assert_allclose(actual, expected, atol=2e-12, rtol=0)


def first(values, epsilon):
    hits = np.flatnonzero(values <= epsilon)
    return int(hits[0]) if len(hits) else -1


def check_walk(r, cfg, c):
    n = r["vertices"]
    p, position = physical_walk(n, c)
    assert r["dhn_c"] == c and r["states"] == 2*n and r["positions"] == n
    kernel = np.zeros_like(p)
    for i, j, w in r["transitions"]:
        kernel[i, j] = w
    np.testing.assert_array_equal(kernel, p)
    np.testing.assert_array_equal(r["state_to_position"], position)
    near(r["stationary"], np.full(2*n, 1/(2*n)))
    near(r["lambda2_after"], 4*np.sin(np.pi/(2*n))**2)

    eigenvalues = np.linalg.eigvals(p - 1/(2*n))
    mode = complex(r["dominant_eigenvalue_real"], r["dominant_eigenvalue_imag"])
    near(abs(mode), max(abs(eigenvalues)))
    near(min(abs(eigenvalues - mode)), 0)
    near(r["asymptotic_rate"], -np.log(max(abs(eigenvalues))))

    a = r["analysis"]
    starts = a["starts"]
    diagnostic = starts.index(cfg["diagnostic_start"])
    if 2*n <= cfg["all_starts_limit"]:
        assert starts == list(range(2*n)) and a["scope"] == "all"
    else:
        assert len(starts) == len(set(starts)) == cfg["sample_starts"] and a["scope"] == "sampled"
    assert a["period"] == 1 and a["max_mass_error"] < 2e-12 and a["minimum_probability"] >= 0

    times = range(cfg["max_steps"] + 1)
    project = np.eye(n)[position]
    dist = np.eye(2*n)[starts]
    tv = {"full": [], "position": []}
    snaps = {s["step"]: s for s in a["snapshots"]}
    assert set(snaps) == ({t for t in cfg["snapshots"] if t in times} | {times[-1]})
    for t in times:
        for space, mass in (("full", dist), ("position", dist @ project)):
            # Sum the excess above uniform, an equivalent definition of TV.
            tv[space].append(np.maximum(mass - 1/mass.shape[1], 0).sum(axis=1))
        if t in snaps:
            near(snaps[t]["distribution"], dist[diagnostic])
            near(snaps[t]["min_probability_over_starts"], dist.min())
        dist = dist @ p

    required_steps = set(times[::cfg["stride"]]) | {times[-1]}
    per_start = [[] for _ in starts]
    crossings = iter(a["crossings"])
    for space, history in tv.items():
        history = np.array(history)
        maximum = history.max(axis=1)
        for epsilon in cfg["epsilons"]:
            hit = first(maximum, epsilon)
            individual = [first(column, epsilon) for column in history.T]
            required_steps.update(t for t in (hit, individual[diagnostic]) if t >= 0)
            for row, value in zip(per_start, individual):
                row.append(value)
            expected = dict(epsilon=epsilon, space=space, checked_through=times[-1],
                            steps=None if hit < 0 else hit, status="unresolved" if hit < 0 else "resolved")
            assert next(crossings) == expected, (n, c, expected)
        rows = [row for row in a["curves"] if row["space"] == space]
        steps = [row["step"] for row in rows]
        assert steps == sorted(set(steps))
        for key, values in (("maximum", maximum), ("minimum", history.min(axis=1)),
                            ("diagnostic", history[:, diagnostic])):
            near([row[key] for row in rows], values[steps])
    assert next(crossings, None) is None
    for space in tv:
        assert [r["step"] for r in a["curves"] if r["space"] == space] == sorted(required_steps)
    assert a["per_start"] == [dict(start=s, first_crossings=row) for s, row in zip(starts, per_start)]


def check_trajectories(r):
    m = r["mcmc"]
    p, position = physical_walk(r["vertices"], r["dhn_c"])
    f = (position < r["vertices"]//2).astype(float)
    truth = f.mean()
    centered = f - truth
    # Stationary finite-run variance from all autocovariances through draws-1.
    draws, reps = m["draws"], m["replicates"]
    lagged = centered
    variance = draws * np.mean(centered**2)
    for lag in range(1, draws):
        lagged = p @ lagged
        variance += 2 * (draws-lag) * np.mean(centered * lagged)
    variance /= draws**2
    means = np.array(m["means"])
    assert len(means) == reps and np.all((0 <= means) & (means <= 1))
    near(m["true_mean"], truth)
    near(m["estimate"], means.mean())
    near(m["rmse"], np.sqrt(np.mean((means - truth)**2)))
    np.testing.assert_allclose(m["ess"], reps * truth*(1-truth)/means.var(ddof=1), rtol=1e-12)
    assert abs(means.mean() - truth) < 6*np.sqrt(variance/reps)
    assert abs(means.var(ddof=1)/variance - 1) < .15


def main(folder):
    studies = [study(n, c) for n in (3, 4, 5, 7, 8, 11, 17)
               for c in sorted({.01*n, .25*n, 1., .5*n, .99*n})]
    cfg = dict(studies=studies, stride=1, max_steps=64, epsilons=[.10003, .90001, .250123, .01003],
               diagnostic_start=2, snapshots=[0, 1, 31, 32, 33, 63, 64, 65], export_kernel=True,
               all_starts_limit=256, sample_starts=3)
    runs = run(folder, "grid", cfg)
    assert len(runs) == len(studies)
    for r, s in zip(runs, studies):
        check_walk(r, cfg, s["dhn_c"])

    # Exact equality at step zero, periodic output stride and block boundaries.
    for horizon in (0, 1, 31, 32, 33):
        case = dict(cfg, studies=[study(4, 2)], max_steps=horizon, stride=13,
                    epsilons=[.875, .75, .031253], diagnostic_start=7)
        check_walk(run(folder, f"horizon-{horizon}", case)[0], case, 2)

    sampled = dict(cfg, studies=[study(9, 1)], seed=98, stride=11, diagnostic_start=17,
                   all_starts_limit=0, mcmc=dict(replicates=2048, draws=128))
    one = run(folder, "serial", sampled, jobs=1)[0]
    four = run(folder, "parallel", sampled)[0]
    assert measured(one) == measured(four)
    check_walk(one, sampled, 1)
    check_trajectories(one)
    other = run(folder, "other-seed", dict(sampled, seed=99))[0]
    assert one["mcmc"]["means"] != other["mcmc"]["means"]
    check_trajectories(other)

    for i, c in enumerate((-1, 0, 5, 5.01)):
        run(folder, f"invalid-rate-{i}", dict(cfg, studies=[study(5, c)]), failure="0 < c < n")
    for i, override in enumerate((dict(graph=dict(family="cycle", sizes=[5])),
                                  dict(policies=[dict(policy="fixed", edge=[0, 2])]))):
        run(folder, f"invalid-path-{i}", dict(cfg, studies=[dict(study(5, 1), **override)]),
            failure="unedited, ordered undirected path")
    run(folder, "invalid-start", dict(cfg, studies=[study(5, 1)], diagnostic_start=10),
        failure="invalid start options")


if __name__ == "__main__":
    with TemporaryDirectory() as temp:
        main(Path(temp))
