"""Small LaTeX tables calculated from the audited notebook outputs."""
import csv
from notebook_runs import OUT, ROOT


def table(path):
    with path.open() as f:
        return list(csv.DictReader(f))


def number(value, places=0):
    return f"{float(value):,.{places}f}" if value not in ("", None) else "--"


def write_values():
    inv = table(OUT/"inventory.csv")
    budget = table(OUT/"wall-budget/budget-comparison.csv")
    timed_runs = table(OUT/"wall-budget/budget-trials.csv")
    nested = table(OUT/"nested-networks/nested-variation.csv")
    selected = [r for r in nested if r["study"].startswith("Western") and r["vertices"]=="96"
                and r["epsilon"]=="0.25" and r["metric"]=="steps" and r["space"]=="full"
                and r["method"]=="max_degree" and r["policy"]=="none"]
    grid = selected[0]
    # A trial appears twice per deadline (one row per observable).
    exhausted = sum(r["status"]=="setup_exhausted" and r["observable"]=="lowest_degree_vertex" for r in timed_runs)
    vals = {"PlannedRuns":sum(int(r["planned"]) for r in inv),
            "MeasuredRuns":sum(int(r["measured"]) for r in inv),
            "SkippedRuns":sum(int(r["skipped"]) for r in inv),
            "PeriodicRows":sum(int(r["periodic"]) for r in inv),
            "CensoredRows":sum(int(r["unresolved"]) for r in inv),
            "Exhausted":exhausted,
            "GridLow":number(grid["minimum"]), "GridHigh":number(grid["maximum"])}
    path=ROOT/"presentation/measurements.tex"
    with path.open("w") as f:
        for key,value in vals.items():
            print(r"\newcommand{\%s}{%s}" % (key, value),file=f)
        def macro(name, lines):
            print(r"\newcommand{\%s}{" % name + "%",file=f)
            for line in lines:
                print(" & ".join(line)+r" \\",file=f)
            print("}",file=f)
        for label, group in (("InventoryA",inv[:8]),("InventoryB",inv[8:])):
            macro(label, [[r["preset"], r["planned"],r["skipped"],r["periodic"],r["unresolved"]] for r in group])
        data=[]
        for study,method,policy,label in [("Lifted path","max_degree","none","Path, reversible"),
                    ("Lifted path","dhn","none","Path, lifted"),
                    ("Cycle shortcut","clockwise","none","Ring, clockwise"),
                    ("Cycle shortcut","clockwise","halves","Ring, diameter")]:
            group={(float(r["budget_seconds"]),r["observable"]):r for r in budget
                   if (r["study"],r["method"],r["policy"])==(study,method,policy)}
            row=[label]
            for secs in (10,60):
                for obs in ("first_half","lowest_degree_vertex"):
                    value=group[secs,obs]["rmse"]
                    row.append(f"{float(value):.2e}" if value not in ("", None) else "--")
            data.append(row)
        macro("BudgetRows",data)
        m=table(OUT/"mcmc/mcmc.csv")
        r=next(r for r in m if r["method"]=="fmmc")
        if r["break_even_effective_draws"]:
            phrase="The estimated karate-club FMMC setup break-even was "+number(float(r["break_even_effective_draws"])/1e6,1)+" million effective draws against maximum degree on this server."
        else:
            phrase="The karate-club FMMC throughput did not yield a finite setup break-even against maximum degree in this run."
        print(r"\newcommand{\BreakEven}{"+phrase+"}",file=f)
