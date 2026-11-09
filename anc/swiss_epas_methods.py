"""Paired planar inverse-square experiment on fixed perturbed coordinates."""
from fractions import Fraction as F
from collections import defaultdict
from itertools import combinations
from pathlib import Path
import argparse, csv, json, math
import numpy as np
from exact_cities import minimize_cities, energy
from planar_riesz_epas_v2 import farthest_first, dyadic_inverse_sqrt

ROOT=Path(__file__).resolve().parent
K=10

def geometry(x,epsilon=F(1,2)):
    pts=[tuple(F(int(v)) for v in row) for row in x]
    centers,r2,c2,trace=farthest_first(pts,2*K)
    bx=16*sum((F(1,j) for j in range(1,2*K+1)),F())/r2
    inc=energy(x,centers[:K],1,2)
    hx=epsilon*dyadic_inverse_sqrt(bx)/16
    plans={}
    for name,b in [("exchange",bx),("hybrid",min(bx,inc))]:
        h=epsilon*dyadic_inverse_sqrt(b)/16
        buckets=defaultdict(list)
        for i,pt in enumerate(pts):buckets[tuple((v+hx/2)//h for v in pt)].append(i)
        plans[name]=(h,{cell[0]:cell for cell in buckets.values()})
    return (dict(centers=centers,radii_m=[None]+[math.sqrt(t["insertion_radius_squared"]) for t in trace[1:]],
                 last_insertion_m=math.sqrt(r2),actual_cover_m=math.sqrt(c2),exchange_bound=str(bx),incumbent_energy=str(inc),
                 origin_m=str(-hx/2)),plans)

def refine(x,reps,cells):
    ids=list(reps);changes=[]
    for sweep in range(1,11):
        before=ids[:]
        for pos,rep in enumerate(reps):
            companions=ids[:pos]+ids[pos+1:]
            def potential(i):return sum((energy(x,(i,j),1,2) for j in companions),F())
            cost,i=min((potential(i),i) for i in cells[rep])
            if cost<potential(ids[pos]):
                changes.append([ids[pos],i]);ids[pos]=i
        if before==ids:break
    return ids,changes


