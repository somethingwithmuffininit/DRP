"""Build the slide illustrations from the notebook's measured CSV and JSON files."""
from collections import defaultdict
import csv
import json
import os
from pathlib import Path
from statistics import median

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from notebook_runs import OUT, ROOT, read

FIG = Path(os.environ.get("MIXBENCH_FIGURE_DIR", ROOT / "presentation/figures")).resolve()
NAVY, TEAL, GOLD, PLUM, GRAY = "#173B59", "#008475", "#BC572C", "#795594", "#667786"
PALETTE = [NAVY, TEAL, GOLD, PLUM, GRAY]
plt.rcParams.update({"font.size": 11, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.labelcolor": NAVY, "text.color": NAVY, "axes.titlecolor": NAVY,
                     "axes.prop_cycle": plt.cycler(color=PALETTE), "pdf.fonttype": 42})
METHOD = dict(clockwise="Clockwise", max_degree="Max degree", metropolis="Metropolis",
              fmmc="FMMC", dhn="Lifted")
POLICY = dict(none="Original", random="Random", min_degree="Min degree", fiedler="Fiedler",
              halves="Diameter", optimal_connectivity="Best connectivity")


def rows(name, file="summary.csv"):
    with (OUT/name/file).open() as f:
        return list(csv.DictReader(f))


def save(fig, name):
    FIG.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG / (name + ".pdf"), bbox_inches="tight", pad_inches=.07)
    plt.close(fig)


def pair():
    return plt.subplots(1, 2, figsize=(10.4, 3.05), layout="constrained")


def curve(ax, data, xlabel, ylabel, **style):
    ax.plot(*zip(*data), **style)
    ax.set(xlabel=xlabel, ylabel=ylabel)
    ax.grid(axis="y", alpha=.17)


def riffle():
    fig, (a, b) = pair()
    data = rows("", "riffle.csv")
    for n in (8, 16, 32, 52):
        xy = [(int(r["shuffles"]), float(r["tv"])) for r in data if int(r["cards"]) == n]
        curve(a, xy, "Independent riffles", "TV from uniform order", label=f"{n} cards", marker=".")
    a.legend(fontsize=9)
    values = [float(r["tv"]) for r in data if int(r["cards"]) == 52 and 5 <= int(r["shuffles"]) <= 10]
    bars = b.bar(range(5,11), values, color=[GRAY, GRAY, GOLD, TEAL, TEAL, TEAL])
    b.bar_label(bars, fmt="%.3f", fontsize=9, padding=3)
    b.axhline(.25, color=NAVY, lw=1, ls="--")
    b.set(xlabel="Riffles of a 52-card deck", ylabel="TV from uniform order", ylim=(0,1.04))
    save(fig, "riffle")


def placement():
    data = rows("cycle-placement")
    edited = [r for r in data if r["edge"]]
    ks = [json.loads(r["edge"])[0]+1 for r in edited]
    fig, (a,b) = pair()
    curve(a, zip(ks, (int(r["steps"]) for r in edited)), "Chord endpoint k (other endpoint 32)",
          "Worst-start steps to TV 0.25", marker="o", ms=3, label="Reciprocal chord")
    a.axhline(int(data[0]["steps"]), color=GOLD, ls="--", label="Original ring")
    a.legend(fontsize=9)
    curve(b, zip(ks, (float(r["lambda2"]) for r in edited)), "Chord endpoint k",
          "Support algebraic connectivity", marker="o", ms=3, color=TEAL)
    b.set_ylim(0, .05)
    save(fig, "placement")


def start_effects():
    fig, (a,b) = pair()
    records = rows("real-networks")
    for row in [r for r in records if r["graph"]=="g1" and r["epsilon"]=="0.25" and r["policy"]=="none"]:
        run = read(OUT/"real-networks"/f"run-{row['run']}.json")
        label = METHOD[row["method"]]
        t = [s["first_crossings"][0] for s in run["analysis"]["per_start"]]
        curve(a, enumerate(sorted(t),1), "Starting vertices, ranked within method", "Steps to TV 0.25", label=label)
        tv = [(p["step"],p["maximum"]) for p in run["analysis"]["curves"] if p["step"] <= 75]
        curve(b, tv, "Transitions", "Maximum TV over 34 starts", label=label, marker=".", ms=5)
    b.axhline(.25, color=GRAY, ls="--", lw=1)
    a.legend(fontsize=9)
    save(fig, "starts")


def glasgow_map():
    info = read(ROOT / "data/glasgow.json")
    study = next(s for s in read(ROOT / "experiments/real-networks.json")["studies"]
                 if s["graph"]["dataset"].endswith("glasgow.json"))
    shortcut = next(p["edge"] for p in study["policies"] if p["policy"] == "fixed")
    angles = np.deg2rad(120 - 360 * np.arange(info["vertices"]) / info["vertices"])
    locations = np.column_stack((np.cos(angles), np.sin(angles)))
    fig, ax = plt.subplots(figsize=(6, 3))
    edges = (ROOT / "data" / info["edges"]).read_text().splitlines()
    for edge in edges:
        pair = [int(i) for i in edge.split()]
        ax.plot(*locations[pair].T, color=NAVY, lw=1.4)
    ax.plot(*locations[shortcut].T, color=GOLD, ls="--", label="Hypothetical link")
    ax.scatter(*locations.T, s=40, facecolors="white", edgecolors=NAVY, zorder=3)
    ax.scatter(*locations[shortcut].T, s=40, color=GOLD, zorder=4)
    for index in [0, *shortcut]:
        right = locations[index, 0] > 0
        ax.annotate(info["mapping"][index], locations[index],
                    xytext=(7 if right else -7, 0), textcoords="offset points",
                    ha="left" if right else "right", va="center")
    ax.set(xlim=(-1.8, 1.8), ylim=(-1.15, 1.15), aspect="equal")
    ax.axis("off")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, 0), frameon=False)
    save(fig, "glasgow-map")


def bridges_and_edges():
    fig,(a,b)=pair()
    data=rows("bridge-count")
    for method in ("max_degree","fmmc"):
        xy=[(json.loads(r["graph_definition"])["bridges"], int(r["steps"])) for r in data
            if r["method"]==method and r["epsilon"]=="0.25"]
        curve(a,xy,"Disjoint bridges between two 8-cliques","Worst-start steps to TV 0.25",
              marker="o",label=METHOD[method])
    a.set_xticks([1,2,4,8])
    a.legend(fontsize=9)
    ranked=[r for r in rows("edge-comparison","edge-ranks.csv") if r["study"]=="Karate club" and r["epsilon"]=="0.25"]
    for policy,marker,color in [("all",".",GRAY),("random","x",GOLD),("min_degree","+",PLUM),
                                ("fiedler","s",TEAL),("optimal_connectivity","o",NAVY)]:
        ps=[r for r in ranked if r["policy"]==policy]
        b.scatter([float(r["delta_lambda2"]) for r in ps],[48-int(r["steps"]) for r in ps],
                  marker=marker, s=15 if policy=="all" else 48,
                  label=POLICY.get(policy,"Every candidate"), alpha=.38 if policy=="all" else .9,
                  **(dict(facecolors="none",edgecolors=color) if marker=="o" else dict(color=color)))
    b.axhline(0, color=GRAY, lw=.8)
    b.set(xlabel="Karate club · increase in support λ₂", ylabel="Worst-start steps saved")
    b.legend(fontsize=8)
    save(fig,"edges")


def cycle_scaling():
    fig,(a,b)=pair()
    data=[r for r in rows("scaling") if r["study"]=="Cycle scaling" and int(r["vertices"])%2==0]
    for pol in ("none","halves","random"):
        group=defaultdict(list)
        for r in data:
            if r["epsilon"]=="0.25" and r["policy"]==pol:
                group[int(r["vertices"])].append(float(r["steps"]))
        x=sorted(group)
        curve(a,[(n,median(group[n])) for n in x],"Cycle vertices","Worst-start steps to TV 0.25",
              marker=".",label=POLICY[pol])
        if pol=="random":
            a.fill_between(x,[min(group[n]) for n in x],[max(group[n]) for n in x],alpha=.16,color=GOLD)
    a.set(xscale="log",yscale="log")
    a.legend(fontsize=9)
    for pol in ("none","halves"):
        d=[r for r in data if r["vertices"]=="256" and r["policy"]==pol]
        curve(b,[(float(r["epsilon"]),int(r["steps"])) for r in d],"TV tolerance ε (stricter →)",
              "Worst-start steps at n = 256",label=POLICY[pol],marker="o")
    b.set(xscale="log",yscale="log")
    b.invert_xaxis()
    b.legend(fontsize=9)
    save(fig,"cycle-scaling")


def path_scaling():
    fig,(a,b)=pair()
    data=rows("path-oracle","path-controls.csv")
    for method in ("dhn","max_degree"):
        points=[r for r in data if r["epsilon"]=="0.1" and r["method"]==method]
        curve(a,[(int(r["vertices"]),int(r["steps"])) for r in points],"Path positions",
              "Endpoint position steps to TV 0.1",label=METHOD[method],marker="o")
        curve(b,[(int(r["vertices"]),float(r["scaled_steps"])) for r in points],"Path positions",
              "Normalized endpoint crossing",label="Lifted / n" if method=="dhn" else "Reversible / n²",marker="o")
    a.set(xscale="log",yscale="log")
    b.set_xscale("log")
    a.legend(fontsize=9)
    b.legend(fontsize=9)
    save(fig,"path")


def dana():
    fig,(a,b)=pair()
    run=read(OUT/"dana-theorem1/run-0.json")
    history=run["analysis"]["curves"]
    for space in ("full","position"):
        curve(a,[(p["step"],p["maximum"]) for p in history if p["space"]==space and p["step"]<=500],
              "Transitions at n = 100","Maximum TV across 200 starts",label=space.capitalize())
    a.axhline(.1,color=GRAY,ls="--")
    a.legend(fontsize=9)
    data=rows("","dana-parity.csv")
    curve(b,[(int(r["n"]),float(r["measured_C"])) for r in data],"Path positions",
          "Measured 4n-step minorization constant",color=TEAL,marker="o")
    b.set(xscale="log",ylim=(.8,1.0))
    save(fig,"dana")


def family():
    data=rows("paper-families")
    fig,axes=plt.subplots(1,3,figsize=(10.4,3),layout="constrained")
    for ax,st,title in zip(axes,("er","ba","tree"),("Connected Erdős–Rényi","Preferential attachment","Binary tree")):
        for p in ("none","random","min_degree","fiedler"):
            grouped=defaultdict(list)
            for r in data:
                if r["study"]==st and r["policy"]==p:
                    grouped[int(r["vertices"])].append(int(r["steps"]))
            curve(ax,[(n,median(v)) for n,v in sorted(grouped.items())],"Vertices",
                  "Steps to TV 0.25" if st=="er" else "",marker="o",label=POLICY[p])
        ax.set_title(title,fontsize=11)
        ax.set_xticks(sorted(grouped))
    axes[0].legend(fontsize=8)
    save(fig,"families")


def method_scaling():
    data=rows("method-scaling","size-statistics.csv")
    fig,(a,b)=pair()
    for ax,study in [(a,"Connectivity heuristics on trees"),(b,"Small-graph connectivity optimum")]:
        for policy in ("none","random","min_degree","fiedler","optimal_connectivity"):
            subset=sorted((r for r in data if r["study"]==study and r["epsilon"]=="0.25" and
                           r["metric"]=="steps" and r["policy"]==policy),key=lambda r:int(r["vertices"]))
            if subset:
                x=[int(r["vertices"]) for r in subset]
                curve(ax,zip(x,[float(r["median"]) if r["median"] else np.nan for r in subset]),
                      "Tree vertices" if ax is a else "Two-clique vertices","Worst-start steps to TV 0.25",
                      label=POLICY[policy],marker=".")
                if policy=="random":
                    ax.fill_between(x,[float(r["q25"]) for r in subset],[float(r["q75"]) for r in subset],
                                    alpha=.15,color=TEAL)
        ax.legend(fontsize=8)
    save(fig,"method-scaling")


def nested():
    data=rows("nested-networks","nested-variation.csv")
    fig,(a,b)=pair()
    for ax,graph,policy,methods in [(a,"Glasgow","none",("dhn","max_degree")),
                                    (b,"Western","none",("max_degree","metropolis","fmmc"))]:
        for idx,method in enumerate(methods):
            space="position" if method=="dhn" else "full"
            subset=sorted((r for r in data if r["study"].startswith(graph) and r["epsilon"]=="0.25" and
                    r["metric"]=="steps" and r["method"]==method and r["policy"]==policy and r["space"]==space
                    and r["median"]),key=lambda r:int(r["vertices"]))
            x=[int(r["vertices"]) for r in subset]
            curve(ax,zip(x,[float(r["median"]) for r in subset]),"Stations in section" if ax is a else "Grid vertices",
                  "Steps to TV 0.25",label=METHOD[method],marker="o")
            ax.fill_between(x,[float(r["minimum"]) for r in subset],[float(r["maximum"]) for r in subset],
                            alpha=.13,color=PALETTE[idx])
        ax.legend(fontsize=9)
    save(fig,"nested")


def mcmc():
    records=rows("mcmc","mcmc.csv")
    fig,(a,b)=pair()
    labels=["Maximum\ndegree","Metro-\npolis","FMMC","Directed\nring","Ring +\nchord"]
    for ax,metric,scale,label in [(a,"rmse",1,"RMSE of replicate mean"),(b,"ess_per_second",1e-6,"Effective draws per second (millions)")]:
        ax.bar(labels,[float(r[metric])*scale for r in records],color=[NAVY,NAVY,NAVY,TEAL,TEAL])
        ax.set_ylabel(label)
        ax.grid(axis="y",alpha=.15)
    save(fig,"mcmc")


def budget():
    records=rows("wall-budget","budget-comparison.csv")
    fig,(a,b)=pair()
    for ax,obs,title in [(a,"left_clique","Left-clique membership"),(b,"lowest_degree_vertex","One low-degree vertex")]:
        selection=[r for r in records if r["study"]=="Two cliques" and r["observable"]==obs]
        for i,method in enumerate(("max_degree","metropolis","fmmc")):
            for seconds,marker in [(10,"o"),(60,"s")]:
                subset={r["policy"]:r for r in selection if r["method"]==method and float(r["budget_seconds"])==seconds}
                policies=("none","random","min_degree","fiedler","optimal_connectivity")
                ax.plot(np.arange(5)+.12*(i-1)+.04*(seconds==60),
                        [float(subset[p]["rmse"]) if subset[p]["rmse"] else np.nan for p in policies],
                        ls="none",marker=marker,ms=6,color=PALETTE[i],
                        markerfacecolor="none" if seconds==10 else PALETTE[i],label=f"{METHOD[method]} · {seconds}s")
        ax.set_xticks(range(5),["Original","Random","Min deg.","Fiedler","Optimal"],rotation=15)
        ax.set(yscale="log",ylabel="RMSE over five trials",title=title)
        ax.grid(axis="y",alpha=.2)
    a.legend(fontsize=7,ncol=2)
    save(fig,"budget")


def make_figures(*, include_budget=True):
    for draw in (riffle,placement,start_effects,glasgow_map,bridges_and_edges,cycle_scaling,path_scaling,
                 dana,family,method_scaling,nested,mcmc,connectivity_map,work_scaling,cycle_distributions,nested_map):
        draw()
    if include_budget:
        budget()


def connectivity_map():
    from scipy.linalg import eigh
    data = rows("edge-comparison","edge-ranks.csv")
    r = next(r for r in data if r["study"] == "Karate club" and r["policy"] == "fiedler")
    run = read(OUT/"edge-comparison"/f"run-{r['run']}.json")
    graph = next(g for g in read(OUT/"edge-comparison/manifest.json")["graphs"] if g["id"]==run["graph"])
    n = graph["vertices"]
    adj = np.zeros((n,n))
    for i,j in graph["edges"]:
        adj[i,j]=adj[j,i]=1
    before, vecs = eigh(np.diag(adj.sum(axis=1))-adj)
    # Draw a fixed layout from the parent network on both panels.
    # Kamada–Kawai uses graph distances and separates pendant vertices.
    import networkx as nx
    g = nx.from_numpy_array(adj)
    xy = np.array([v for _,v in sorted(nx.kamada_kawai_layout(g).items())])
    i,j = run["edge"]
    fig,(a,b)=pair()
    for ax,edit in [(a,False),(b,True)]:
        vals = before
        if edit:
            adj[i,j]=adj[j,i]=1
            vals, vecs_after = eigh(np.diag(adj.sum(axis=1))-adj)
            if vecs_after[:,1]@vecs[:,1]<0:
                vecs_after[:,1]*=-1
        colors=(vecs_after if edit else vecs)[:,1]
        for u,v in graph["edges"]:
            ax.plot(*xy[[u,v]].T,color="#B9C4CB",lw=.8,zorder=0)
        if edit:
            ax.plot(*xy[[i,j]].T,color=GOLD,lw=3,zorder=1)
        circles=ax.scatter(*xy.T,c=colors,cmap="BrBG",vmin=-.45,vmax=.45,s=85,edgecolor="white",lw=.5)
        for v in (i,j):
            ax.annotate(str(v),xy[v],xytext=(5,6),textcoords="offset points",color=NAVY,weight="bold")
        ax.set_title(f"{'Fiedler edge added' if edit else 'Original friendships'}   ·   λ₂ = {vals[1]:.3f}")
        ax.set_aspect("equal")
        ax.axis("off")
    fig.colorbar(circles,ax=[a,b],shrink=.7,label="Fiedler coordinate")
    save(fig,"karate-map")


def work_scaling():
    data=rows("method-scaling","size-statistics.csv")
    fig,(a,b)=pair()
    for ax,metric,ylabel in ((a,"setup_seconds","Kernel setup (seconds)"),(b,"asymptotic_rate","Asymptotic log decay per step")):
        for m in ("dhn","max_degree","fmmc"):
            sub=sorted((r for r in data if r["study"]=="Ordered paths" and r["epsilon"]=="0.25"
                        and r["metric"]==metric and r["method"]==m and r["median"]),key=lambda r:int(r["vertices"]))
            curve(ax,[(int(r["vertices"]),float(r["median"])) for r in sub],"Path positions",ylabel,
                  marker="o",label=METHOD[m])
        ax.set(xscale="log",yscale="log")
        ax.legend(fontsize=9)
    save(fig,"work-scaling")


def cycle_distributions():
    runs = [read(f) for f in (OUT/"cycle-placement").glob("run-*.json")]
    fig,(a,b)=pair()
    for ax,r in ((a,next(r for r in runs if r["edge"] is None)),
                 (b,next(r for r in runs if r["edge"]==[15,31]))):
        for s in r["analysis"]["snapshots"]:
            if s["step"] in (10,100):
                curve(ax,enumerate(s["distribution"]),"Cycle vertex", "Probability from vertex 0",
                      label=f"t = {s['step']}",marker=".",ms=4)
        ax.axhline(1/32,color=GOLD,ls="--",label="Uniform")
        ax.set_title("Original ring" if r["edge"] is None else "Ring with diameter shortcut")
        ax.set_ylim(0,.27)
        ax.legend(fontsize=9)
    save(fig,"cycle-mass")


def nested_map():
    import networkx as nx
    m=read(OUT/"nested-networks/manifest.json")
    graphs=sorted((g for g in m["graphs"] if g["name"].startswith("Western") and g["replicate"]==0),
                  key=lambda g:g["vertices"])
    parent=nx.Graph()
    parent.add_edges_from(graphs[-1]["edges"])
    layout=nx.kamada_kawai_layout(parent)
    coords=np.array([layout[i] for i in range(len(layout))])
    fig,axes=plt.subplots(1,4,figsize=(10.6,2.7),layout="constrained")
    old=0
    for ax,g in zip(axes,graphs,strict=True):
        n=g["vertices"]
        for i,j in g["edges"]:
            ax.plot(*coords[[i,j]].T,lw=1,color=TEAL if max(i,j)>=old else "#A8B9C4",zorder=0)
        ax.scatter(*coords[:n].T,s=18,c=[NAVY]*old+[TEAL]*(n-old))
        ax.scatter(*coords[0],s=40,marker="*",color=GOLD)
        ax.set(aspect="equal",title=f"{n} vertices",xlim=(-1.1,1.1),ylim=(-1.1,1.1))
        ax.axis("off")
        old=n
    save(fig,"nested-map")
