"""A repetition filter preserves the serial plan, edit RNG and graph realization."""
import json
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory

def main(temp):
    cfg = dict(seed=19, wall_budget=dict(repeats=3, seconds=[1e-12, .005]), studies=[
        dict(name=family, graph=dict(family=family, sizes=[9,12]), methods=["max_degree","metropolis"],
             policies=[dict(policy="none"), dict(policy="random")]) for family in ("er","tree")])
    def rows(index):
        f = temp / "config.json"
        f.write_text(json.dumps(cfg))
        out = temp / str(index)
        subprocess.run([sys.argv[1], "run", f, "--output", out, "--jobs", "1"], check=True)
        manifest = json.loads((out / "manifest.json").read_text())
        plan = []
        for i in range(manifest["run_count"]):
            r = json.loads((out / f"run-{i}.json").read_text())
            assert r["checkpoints"][0]["draws"] == 0
            plan.append(tuple(r[k] for k in ("trial","graph","seed","method","policy","edge","rare_vertex")))
        return manifest["graphs"], plan
    graphs, full = rows("all")
    for rep in range(3):
        cfg["wall_budget"]["only_repetition"] = rep
        other, part = rows(rep)
        assert other == graphs and part == [v for v in full if v[0] == rep]
    assert len(full) == 48

if __name__ == "__main__":
    with TemporaryDirectory() as tmp:
        main(Path(tmp))
