"""Fixed-input cardinality experiment for the planar I, X and H certificates.

Run from any directory with --run to collect new timings, or --check to replay
the grids and audit the archived energies without re-running timed-out searches.
Each (panel, seed) fixes the coordinates, row order, grid origin and accuracy.
"""
import argparse
from collections import Counter, defaultdict
from fractions import Fraction as F
from math import comb
from multiprocessing import get_context
from pathlib import Path
from statistics import median
from time import perf_counter
import json
import platform
import sys
import numpy as np
from paired import make_input, certificates, grid
from epas_refine_v2 import D, energy, fingerprint, solve_reduced

OUT = Path(__file__).resolve().parent.parent / "results"
RECORD = OUT / "fixed-input-runs.json"
POINTS = OUT / "fixed-input-points.json"
EPS = F(1, 2)
SEEDS = (17, 29, 43)
METHODS = ("incumbent", "exchange", "hybrid")
T = 10**7
MAX_M = 160
SECONDS = 5.0
PANELS = {
    "small": (dict(n=2048, c=24, family="irregular", rho="1/4096"),
              (4, 5, 8, 10, 11, 12, 13, 16)),
    "large": (dict(n=1024, c=256, family="lattice", rho="1/256"),
              (4, 8, 12, 13, 16, 32, 64, 96, 128, 192, 256)),
}


def write(path, data):
    path.write_text(json.dumps(data, separators=(",", ":"))+"\n", encoding="utf-8")


def worker(pipe, X, g):
    try:
        pipe.send("ready")
        pipe.recv()
        # Enforce one external wall limit including construction of pair weights.
        # The larger internal timeout cannot fire before the external limit.
        value = solve_reduced(X, g, timeout=10*SECONDS, max_M=g["M"],
                              max_subsets=comb(g["M"], g["k"]))
        pipe.send(("ok", value))
    except Exception as exc:
        pipe.send(("error", repr(exc)))
    finally:
        pipe.close()


def timed(X, grid_data):
    ctx = get_context("spawn")
    left, right = ctx.Pipe()
    args = {key: grid_data[key] for key in ("k","M","denominator","epsilon","representatives")}
    proc = ctx.Process(target=worker, args=(right, X, args))
    proc.start()
    right.close()
    try:
        if not left.poll(60) or left.recv() != "ready":
            raise RuntimeError("Worker startup failed.")
        start = perf_counter()
        left.send("go")
        if not left.poll(SECONDS):
            return dict(status="timeout", seconds=SECONDS)
        tag, data = left.recv()
        elapsed = perf_counter()-start
        if tag != "ok" or data["status"] != "complete":
            raise RuntimeError(str(data))
        return dict(status="complete", seconds=elapsed, energy=data["energy_exact"],
                    selected=data["selected"], lower=data["certificate_lower_exact"],
                    nodes=data["search_stats"]["nodes"])
    finally:
        if proc.is_alive():
            proc.terminate()
        proc.join()
        left.close()


def audited_split():
    counts = defaultdict(Counter)
    solved_together = 0
    for line in (OUT/"paired.jsonl").read_text(encoding="utf-8").splitlines():
        row = json.loads(line)
        m = {x["method"]:x for x in row["methods"]}
        i,x,h = (m[j] for j in METHODS)
        order = "fewer" if x["M"]<i["M"] else "more" if x["M"]>i["M"] else "equal"
        counts["planned" if row["config"]["search"] else "occupancy_only"][order]+=1
        if h["status"]==i["status"]=="complete":
            assert h["representatives_sha256"]==i["representatives_sha256"] and h["energy_exact"]==i["energy_exact"]
            solved_together+=1
    return dict(counts={k:dict(v) for k,v in counts.items()}, identical_completed_HI=solved_together)


def grids(X, k):
    bounds, meta = certificates(X,k)
    gs = {m:grid(X,k,EPS,bounds[m],F(0)) for m in METHODS}
    assert gs["exchange"]["M"]>=2*k and gs["hybrid"]["M"]==min(gs[m]["M"] for m in METHODS[:2])
    if k<=11:
        assert bounds["incumbent"]<=bounds["exchange"] and gs["hybrid"]["representatives"]==gs["incumbent"]["representatives"]
    summary = {m:dict(B=str(bounds[m]), h=g["h"], M=g["M"],
                      subsets=comb(g["M"],k), reps_sha256=fingerprint(np.array(g["representatives"],dtype=np.int64)))
               for m,g in gs.items()}
    return gs, dict(rk2=meta["rk2"],r2k2=meta["r2k2"],methods=summary)


def admit(info):
    return info["M"]<=MAX_M and info["subsets"]<=T


def generate():
    points={}
    result=[]
    for panel,(cfg,ks) in PANELS.items():
        for si,seed in enumerate(SEEDS):
            X = make_input(dict(cfg,seed=seed))
            key = f"{panel}-{seed}"
            points[key] = dict(config=dict(cfg,seed=seed),denominator=D,sha256=fingerprint(X),numerators=X.tolist())
            for ki,k in enumerate(ks):
                gs,row = grids(X,k)
                row.update(panel=panel,seed=seed,k=k,input_sha256=points[key]["sha256"])
                # Rotate the three methods to avoid assigning one method a fixed position.
                shift=(si+ki)%3
                order=METHODS[shift:]+METHODS[:shift]
                row["primary"]={m:(timed(X,gs[m]) if admit(row["methods"][m])
                                   else {"status":"size_excluded"}) for m in order}
                # Preselected high-cardinality stress test: no binomial or occupancy gate.
                if panel=="large" and k==128:
                    row["ungated"]={m:timed(X,gs[m]) for m in order}
                result.append(row)
                print(key,k,[(row["methods"][m]["M"],row["primary"][m]["status"]) for m in METHODS],flush=True)
    write(POINTS,points)
    rec=dict(epsilon=str(EPS),origin=[0,0],max_M=MAX_M,budget=T,seconds=SECONDS,
             environment=dict(python=platform.python_version(),numpy=np.__version__,system=platform.platform()),
             units="Integer coordinate numerators divided by 2^20; seconds exclude worker startup and grid construction.",
             paired_split=audited_split(),rows=result)
    write(RECORD,rec)
    return rec


def verify():
    data=json.loads(RECORD.read_text()); points=json.loads(POINTS.read_text())
    assert data["paired_split"]==audited_split()
    assert data["epsilon"]==str(EPS) and data["origin"]==[0,0]
    assert (data["max_M"],data["budget"],data["seconds"])==(MAX_M,T,SECONDS)
    rows={(r["panel"],r["seed"],r["k"]):r for r in data["rows"]}
    assert len(rows)==len(data["rows"])==sum(len(ks)*len(SEEDS) for _,ks in PANELS.values())
    for panel,(cfg,ks) in PANELS.items():
        for seed in SEEDS:
            ref=points[f"{panel}-{seed}"]; X=make_input(dict(cfg,seed=seed))
            assert ref["denominator"]==D and np.array_equal(X,np.asarray(ref["numerators"])) and ref["sha256"]==fingerprint(X)
            for k in ks:
                r=rows[panel,seed,k]; gs,calc=grids(X,k)
                assert r["input_sha256"]==ref["sha256"] and all(r[key]==val for key,val in calc.items())
                assert ("ungated" in r)==(panel=="large" and k==128)
                for group in ("primary","ungated"):
                    if group not in r: continue
                    for m in METHODS:
                        rec=r[group][m]
                        tested=group=="ungated" or admit(calc["methods"][m])
                        assert (rec["status"]!="size_excluded")==tested
                        if not tested: continue
                        assert rec["status"] in ("timeout","complete")
                        if rec["status"]=="timeout":
                            assert rec["seconds"]==SECONDS and "energy" not in rec
                        else:
                            ids=rec["selected"]
                            assert len(set(ids))==k and set(ids)<=set(gs[m]["representatives"])
                            assert energy(X,ids)==F(rec["energy"])==F(rec["lower"])*(1+EPS)
                for u,v in (("incumbent","hybrid"),("exchange","hybrid")):
                    a,b=r["primary"][u],r["primary"][v]
                    if a["status"]==b["status"]=="complete" and gs[u]["representatives"]==gs[v]["representatives"]:
                        assert a["energy"]==b["energy"]
    return data


def show(data):
    print("Original paired split:",data["paired_split"])
    for panel,(_,ks) in PANELS.items():
        for k in ks:
            rows=[r for r in data["rows"] if r["panel"]==panel and r["k"]==k]
            print(panel,k,[(median(r["methods"][m]["M"] for r in rows),
                  dict(Counter(r["primary"][m]["status"] for r in rows))) for m in METHODS])
    print("Ungated:",[(r["seed"],{m:r["ungated"][m]["status"] for m in METHODS}) for r in data["rows"] if "ungated" in r])


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    modes=p.add_mutually_exclusive_group(required=True)
    for name in ("run","check","summary"):modes.add_argument("--"+name,action="store_true")
    args=p.parse_args()
    if args.run:generate()
    data=verify() if args.run or args.check else json.loads(RECORD.read_text())
    show(data)
