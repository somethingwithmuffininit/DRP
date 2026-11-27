"""Test the 64-cycle chord that exhausted the initial SCS iteration budget."""
import argparse
import json
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/"scripts"))
from notebook_runs import EXE, read, single_thread_env


def verify(r):
    assert (r["vertices"], r["method"], r["edge"]) == (64, "fmmc", [51,63])
    p = np.zeros((64,64))
    for i,j,v in r["transitions"]:
        assert i == j or (i-j)%64 in (1,63) or {i,j} == {51,63}
        p[i,j] = v
    np.testing.assert_allclose(p.sum(axis=1), 1, atol=1e-12, rtol=0)
    np.testing.assert_array_equal(p,p.T)
    rho = np.max(abs(np.linalg.eigvalsh(p-1/64)))
    assert p.min() >= 0 and abs(rho-r["solver"]["objective"]) < 1e-5
    for key in ("primal_residual","dual_residual","gap"):
        assert abs(r["solver"][key]) < 1e-5
    # The unedited, optimally lazy cycle is feasible on the larger support.
    rho_original = (1+np.cos(2*np.pi/64))/(3-np.cos(2*np.pi/64))
    assert rho <= rho_original + 1e-5
    return r["solver"]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-json", type=Path, help="independently verify a saved run instead of solving again")
    args = parser.parse_args()
    if args.run_json:
        print(verify(read(args.run_json)))
    else:
        with TemporaryDirectory() as tmp:
            folder=Path(tmp)
            file=folder/"case.json"
            file.write_text(json.dumps(dict(max_steps=2,studies=[
                dict(name="Chorded 64-cycle",graph=dict(family="cycle",sizes=[64]),methods=["fmmc"],
                     policies=[dict(policy="fixed",edge=[51,63])])])))
            subprocess.run([EXE,"run",file,"--jobs","1","--output",folder/"run"],check=True,env=single_thread_env())
            print(verify(read(folder/"run/run-0.json")))
