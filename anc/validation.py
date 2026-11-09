"""Reproduce the nine validation runs at epsilon one half."""
import csv
import json
import platform
import sys
from collections import defaultdict
from fractions import Fraction
from itertools import combinations
from pathlib import Path
from statistics import median
import numpy
from epas_refine_v2 import dataset, energy, fingerprint, prepare, refine, solve_reduced

ROOT = Path(__file__).resolve().parent
CAPS = (0, 1, 2, 3, 5, 10)

def rerun():
    runs = []
    for n, kind, k in ((20,"uniform",4),(20,"clusters",4),(1000,"clusters",5)):
        for seed in (17,29,43):
            spec = dict(n=n,k=k,c=2*k,kind=kind,epsilon="1/2",rho="1/4096",seed=seed)
            points = dataset(spec)
            grid = prepare(points,k,spec["epsilon"])
            reduced = solve_reduced(points,grid)
            assert reduced["status"] == "complete"
            initial = tuple(reduced["selected"])
            # Every budget uses precisely the same initial state. A deterministic
            # single trace yields its prefixes including early stopping.
            trajectory = refine(points,grid,initial,max(CAPS))
            optimum = None
            if n == 20:
                optimum = min(energy(points,choice) for choice in combinations(range(n),k))
                assert optimum <= Fraction(reduced["energy_exact"]) <= (1+Fraction(spec["epsilon"]))*optimum
            observations = []
            for cap in CAPS:
                state = trajectory["trace"][min(cap,trajectory["actual_sweeps"])]
                value = Fraction(state["energy_exact"])
                assert value <= Fraction(reduced["energy_exact"])
                observations.append(dict(cap=cap,actual_passes=state["sweep"],energy_exact=str(value),
                    energy=float(value),gain_percent=float(100*(1-value/Fraction(reduced["energy_exact"]))),
                    extra_s=state["seconds"],total_s=grid["prep_seconds"]+reduced["search_seconds"]+state["seconds"]))
            runs.append(dict(config=spec,sha256=fingerprint(points),M=grid["M"],M_over_n=grid["M"]/n,
                        cell_side=grid["h"],initial=initial,representatives=grid["representatives"],
                        prep_s=grid["prep_seconds"],search_s=reduced["search_seconds"],
                        optimum_exact=str(optimum) if optimum else None,observations=observations))
    return runs

if __name__ == "__main__":
    path = ROOT.parent / "measurements" / "validation.json"
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(rerun(),indent=2),encoding="utf-8")
