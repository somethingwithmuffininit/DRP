"""Summarize experiments and save the notebook's tables, figures, and tests."""
from collections import Counter
import json
import os
import shutil
import subprocess
import sys
import xml.etree.ElementTree as ET

from IPython.display import Image, Markdown, display

from notebook_runs import BUILD, OUT, ROOT, read, single_thread_env
from comparison_stats import tables, write_csv
from plot_results import main as plot_all, mcmc_table
from study_reports import requested_reports
from size_scaling import write_scaling


def report(folder):
    if os.environ.get("MIXBENCH_ALL_PLOTS") == "1":
        rows = plot_all(folder)
    else:
        manifest = read(folder / "manifest.json")
        runs = [read(folder / f"run-{i}.json") for i in range(manifest["run_count"])]
        graphs = {g["id"]: g for g in manifest["graphs"]}
        rows = tables(folder, runs, graphs)
        requested_reports(folder, manifest, runs, manifest["config"].get("reports", []))
        write_scaling(folder, [dict(r, experiment=folder.name, jobs=manifest["jobs"]) for r in rows])
        mcmc_table(folder, runs)
    print("Crossings and skip records:", dict(Counter(r["status"] for r in rows)))
    for r in [r for r in rows if r.get("epsilon") == .25 and r["policy"] == "none"][:8]:
        print(f"  {r['study']} (n={r['vertices']}) {r['method']} / {r['space']}: "
              f"{r['steps'] if r['steps'] is not None else r['status']} ({r['scope']} starts)")
    return rows


def show_table(name, rows):
    """Display the same data that is written to the table CSV."""
    directory = OUT / "tables"
    directory.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0])
    write_csv(directory / f"{name}.csv", rows)

    def line(values):
        return "| " + " | ".join(str(v).replace("|", r"\|") for v in values) + " |"

    markdown = [line(fields), line(["---"] * len(fields))]
    markdown.extend(line(row.get(field, "") for field in fields) for row in rows)
    display(Markdown("\n".join(markdown)))


def show_plot(source, edition):
    """Save a vector figure and a notebook preview under its source edition."""
    folder = OUT / "figures" / edition
    folder.mkdir(parents=True, exist_ok=True)
    target = folder / source.name
    if target.resolve() != source.resolve():
        shutil.copyfile(source, target)
    prefix = target.with_suffix("")
    subprocess.run(["pdftoppm", "-singlefile", "-scale-to", "900", "-png", target, prefix],
                   check=True)
    display(Image(filename=str(prefix.with_suffix(".png")), width=800))


def release_tests(descriptions):
    import matplotlib.pyplot as plt

    env = single_thread_env()
    OUT.mkdir(parents=True, exist_ok=True)
    configure = ["cmake", "-S", ROOT, "-B", BUILD, "-DCMAKE_BUILD_TYPE=Release",
                 f"-DPython3_EXECUTABLE={sys.executable}"]
    compile_release = ["cmake", "--build", BUILD, "--target", "mixbench", "mix_tests", "--parallel", "2"]
    with (OUT / "build.log").open("w") as log:
        for command in (configure, compile_release):
            done = subprocess.run(command, env=env, stdout=log, stderr=subprocess.STDOUT)
            if done.returncode:
                raise RuntimeError(f"Build failed. See {OUT / 'build.log'}.")

    ctest = ["ctest", "--test-dir", str(BUILD)]
    registered = json.loads(subprocess.check_output([*ctest, "--show-only=json-v1"], env=env, text=True))
    names = {test["name"] for test in registered["tests"]}
    if missing := descriptions.keys() - names:
        raise RuntimeError(f"Missing CTest suites: {', '.join(sorted(missing))}")

    junit = OUT / "tests.xml"
    junit.unlink(missing_ok=True)
    run = subprocess.run([*ctest, "--parallel", "2", "--output-on-failure", "--output-junit", junit],
                         env=env, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    (OUT / "release-tests.log").write_text(run.stdout)
    if not junit.is_file():
        raise RuntimeError(f"CTest did not write JUnit results. See {OUT / 'release-tests.log'}.")
    rows = []
    for case in ET.parse(junit).iter("testcase"):
        name = case.attrib["name"]
        status = case.attrib["status"]
        rows.append({"Suite": name, "What it checks": descriptions.get(name, "See CTest configuration"),
                     "Result": {"run": "Passed", "fail": "Failed", "notrun": "Skipped"}.get(status, status),
                     "Seconds": round(float(case.attrib["time"]), 3)})
    passed = sum(row["Result"] == "Passed" for row in rows)
    print(f"{passed} of {len(names)} CTest suites passed.")
    show_table("test-results", rows)

    fig, ax = plt.subplots(figsize=(8, 5.3), layout="constrained")
    ordered = sorted(rows, key=lambda row: row["Seconds"])
    ax.barh([row["Suite"] for row in ordered], [row["Seconds"] for row in ordered],
            color=["#087f73" if row["Result"] == "Passed" else "#ab512f" for row in ordered])
    ax.set(xlabel="Elapsed seconds per suite", title="Test runtimes")
    ax.grid(axis="x", alpha=.2)
    folder = OUT / "figures" / "tests"
    folder.mkdir(parents=True, exist_ok=True)
    target = folder / "test-runtimes.pdf"
    fig.savefig(target)
    plt.close(fig)
    show_plot(target, "tests")
    if run.returncode or passed != len(names):
        print(run.stdout)
        raise RuntimeError(f"CTest checks did not all pass. See {OUT / 'release-tests.log'}.")
