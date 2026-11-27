"""Execution helpers shared by the notebook and focused tests."""
from concurrent.futures import ThreadPoolExecutor
import fcntl
import json
import os
from pathlib import Path
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]
BUILD = Path(os.environ.get("MIXBENCH_BUILD_DIR", ROOT / "build")).resolve()
EXE = BUILD / "mixbench"
OUT = Path(os.environ.get("MIXBENCH_RESULTS", ROOT / "results")).resolve()


def read(path):
    return json.loads(path.read_text())


def single_thread_env():
    return dict(os.environ, OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1", OMP_NUM_THREADS="1")


def execute(config, folder, jobs, cpu=None):
    cmd = [str(EXE), "run", str(config), "--jobs", str(jobs), "--output", str(folder)]
    if cpu is not None:
        cmd = ["taskset", "-c", str(cpu), *cmd]
    subprocess.run(cmd, check=True, env=single_thread_env())


def validate(folder, config):
    manifest = read(folder / "manifest.json")
    if manifest["status"] != "complete" or manifest["config"] != config:
        raise ValueError(f"Incomplete or different experiment in {folder}. Choose a fresh MIXBENCH_RESULTS.")
    if any(not (folder / f"run-{i}.json").is_file() for i in range(manifest["run_count"])):
        raise ValueError(f"Missing measurements in {folder}")
    return manifest


def run_preset(name, jobs=4):
    config = ROOT / "experiments" / (name + ".json")
    folder = OUT / name
    OUT.mkdir(parents=True, exist_ok=True)
    with (OUT / ("." + name + ".lock")).open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        if not folder.exists():
            execute(config, folder, jobs)
        manifest = validate(folder, read(config))
    print(f"{name}: {manifest['run_count']} records, {len(manifest['graphs'])} graphs")
    return folder


def budget_repetitions(config, output, cpus):
    """Run whole repetitions on separate physical cores, one worker per timed trial."""
    cfg = read(config)
    repeats = cfg["wall_budget"]["repeats"]
    if len(cpus) != repeats or len(set(cpus)) != repeats:
        raise ValueError("Supply one distinct allowed physical CPU per repetition")
    physical = set()
    for cpu in cpus:
        topo = Path(f"/sys/devices/system/cpu/cpu{cpu}/topology")
        if cpu not in os.sched_getaffinity(0):
            raise ValueError(f"CPU {cpu} unavailable")
        physical.add(tuple((topo / k).read_text().strip() for k in ("physical_package_id", "core_id")))
    if len(physical) != repeats:
        raise ValueError("Choose CPUs on distinct physical cores")
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=output.parent) as tmp:
        stage = Path(tmp)
        def run_rep(rep):
            scoped = dict(cfg, wall_budget=dict(cfg["wall_budget"], only_repetition=rep))
            file = stage / f"rep-{rep}.json"
            file.write_text(json.dumps(scoped))
            folder = stage / f"rep-{rep}"
            execute(file, folder, 1, cpus[rep])
            return folder, validate(folder, scoped)
        with ThreadPoolExecutor(max_workers=repeats) as pool:
            parts = list(pool.map(run_rep, range(repeats)))
        manifest = dict(parts[0][1], config=cfg, repetition_cpus=cpus)
        merged = stage / "merged"
        merged.mkdir()
        count = 0
        for graph in manifest["graphs"]:
            for rep, (folder, part) in enumerate(parts):
                if part["graphs"] != manifest["graphs"]:
                    raise ValueError("Budget repetitions built different graphs")
                for i in range(part["run_count"]):
                    run = read(folder / f"run-{i}.json")
                    if run["graph"] == graph["id"]:
                        if run["trial"] != rep:
                            raise ValueError("Unexpected repetition in budget shard")
                        run.update(id=count, cpu=cpus[rep])
                        (merged / f"run-{count}.json").write_text(json.dumps(run))
                        count += 1
        manifest["run_count"] = count
        (merged / "manifest.json").write_text(json.dumps(manifest))
        merged.rename(output)


def run_budget():
    name = "wall-budget"
    file = ROOT / "experiments" / (name + ".json")
    folder = OUT / name
    OUT.mkdir(parents=True, exist_ok=True)
    with (OUT / ".wall-budget.lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        if not folder.exists():
            cpus = os.environ.get("MIXBENCH_BUDGET_CPUS")
            if cpus:
                budget_repetitions(file, folder, list(map(int, cpus.split(","))))
            else:
                execute(file, folder, 1)
        manifest = validate(folder, read(file))
    from study_reports import budget_report
    runs = [read(folder / f"run-{i}.json") for i in range(manifest["run_count"])]
    summary = budget_report(folder, manifest, runs)
    print(f"{len(runs)} trials, {len(summary)} observable / deadline summaries")
    return folder
