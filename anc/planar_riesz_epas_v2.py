"""Exact-arithmetic reference EPAS for planar inverse-square Riesz selection.

Python 3.10+, standard library only. Default: the manuscript's 2k-center mode.
The optional compact mode uses 2k-2 centers (k>=2); the same proof applies
because at least k-1 free centers remain, matching the k-1 other selected sites.

Input JSON: {"points": [["0", "0"], ...], "k": 5, "epsilon": "1/2"}.
Coordinates may be integers, decimal strings, or rational strings. JSON decimal
tokens are parsed as decimal strings, not binary floats. Every algorithmic
decision uses Fraction; display decimals are outside the solver.

Usage: python planar_riesz_epas.py input.json --mode compact --output result.json
The implementation enumerates only representative subsets. It never computes an
original optimum or a local optimum to choose its cover, scale or representatives.
The optional --search branch_bound accelerates only the exact search using
integer lower bounds and exact rational incumbents (exact_subset_search.py).
"""
from __future__ import annotations
import argparse
from fractions import Fraction as F
from itertools import combinations
import json
from math import comb
from pathlib import Path


def rational(x):
    if isinstance(x, float):
        raise TypeError('Use decimal strings, not binary floating-point coordinates.')
    return F(x)


def distance2(p, q):
    return (p[0]-q[0])**2+(p[1]-q[1])**2


def energy(points, ids):
    return sum((1/distance2(points[i], points[j]) for i,j in combinations(ids,2)), F(0))


def enumerate_best(points, ids, k):
    """Lexicographic input-index tie breaking; build weights only for ids."""
    ids=tuple(sorted(ids))
    weights={(i,j):1/distance2(points[i],points[j]) for i,j in combinations(ids,2)}
    best=None
    for selection in combinations(ids,k):
        value=sum((weights[i,j] for i,j in combinations(selection,2)),F(0))
        if best is None or value<best[0]:
            best=(value,selection)
    return best[1],best[0],comb(len(ids),k)


def farthest_first(points, m):
    centers=[0]
    nearest=[distance2(p,points[0]) for p in points]
    trace=[{'center':0,'insertion_radius_squared':None}]
    chosen={0}
    while len(centers)<m:
        nxt=max((i for i in range(len(points)) if i not in chosen),key=lambda i:(nearest[i],-i))
        r2=nearest[nxt]
        trace.append({'center':nxt,'insertion_radius_squared':r2})
        centers.append(nxt)
        chosen.add(nxt)
        nearest=[min(v,distance2(p,points[nxt])) for p,v in zip(points,nearest)]
    return centers,trace[-1]['insertion_radius_squared'],max(nearest),trace


def dyadic_inverse_sqrt(bound):
    """Return the largest power of two z with bound*z*z<=1; no square roots."""
    e=(bound.denominator.bit_length()-bound.numerator.bit_length())//2
    z=F(2**e) if e>=0 else F(1,2**(-e))
    if bound*z*z>1:
        z/=2
    assert bound*z*z<=1<bound*(2*z)**2
    return z


def solve(raw_points,k,epsilon=F(1,2),cover_mode='paper',search='enumerate'):
    points=[tuple(rational(v) for v in p) for p in raw_points]
    n=len(points)
    if any(len(p)!=2 for p in points) or len(set(points))!=n:
        raise ValueError('Candidates must be distinct planar points.')
    if not isinstance(k,int) or not 0<=k<=n:
        raise ValueError('Require integer 0<=k<=n.')
    requested=rational(epsilon)
    if requested<=0 or cover_mode not in ('paper','compact') or search not in ('enumerate','branch_bound'):
        raise ValueError('Require epsilon>0, mode paper or compact, and a valid search method.')
    eps=min(requested,F(1))
    result={'points':points,'n':n,'k':k,'s':2,'dimension':2,'epsilon_requested':requested,
            'epsilon_used':eps,'mode':cover_mode,'indexing':'zero-based input indices'}
    m=2*k if cover_mode=='paper' else max(k,2*k-2)
    if k<=1 or k==n:
        selection=tuple(range(k))
        value=energy(points,selection)
        result.update(branch='trivial',representatives=list(range(n)),M=n,
                      selected=selection,energy=value,certificate=[value,value],subsets_evaluated=1)
        return result
    if n<=m:
        selection,value,count=enumerate_best(points,range(n),k)
        result.update(branch='small_instance_exact',representatives=list(range(n)),M=n,
                      selected=selection,energy=value,certificate=[value,value],subsets_evaluated=count)
        return result
    centers,r2,cover2,trace=farthest_first(points,m)
    harmonic=sum((F(1,j) for j in range(1,m+1)),F(0))
    bound=16*harmonic/r2
    scale=dyadic_inverse_sqrt(bound)
    h=eps*scale/16
    # A translated global grid is equally valid. Integer-coordinate anchors
    # lie at cell centers when h divides 1, avoiding a misleading boundary case.
    origin=(-h/2,-h/2)
    keyed=[]
    for i,p in enumerate(points):
        key=tuple((p[a]-origin[a])//h for a in range(2))
        keyed.append((key,i))
    keyed.sort()  # O(n log n), no expected-time hashing assumption.
    cells=[]
    mapping=[None]*n
    for key,i in keyed:
        if not cells or cells[-1]['key']!=key:
            cells.append({'key':key,'representative':i,'members':[]})
        cells[-1]['members'].append(i)
        mapping[i]=cells[-1]['representative']
    representatives=sorted(c['representative'] for c in cells)
    assert len(representatives)>=k
    search_stats=None
    if search=='branch_bound':
        from exact_subset_search_v2 import minimize
        weights={(i,j):1/distance2(points[i],points[j]) for i,j in combinations(representatives,2)}
        selection,value,search_stats=minimize(weights,representatives,k)
        count=None  # Partial-tree visits are not full subset evaluations.
    else:
        selection,value,count=enumerate_best(points,representatives,k)
    result.update(branch='exchange_grid',center_count=m,centers=centers,traversal=trace,
                  last_insertion_radius_squared=r2,actual_cover_radius_squared=cover2,
                  harmonic=harmonic,potential_bound=bound,dyadic_separation_scale=scale,
                  cell_side=h,grid_origin=origin,cells=cells,replacement_map=mapping,
                  representatives=representatives,M=len(representatives),selected=selection,
                  energy=value,certificate=[value/(1+eps),value],subsets_evaluated=count)
    if search_stats is not None:
        result.update(search_method='branch_bound',search_stats=search_stats)
    return result


def serializable(value):
    if isinstance(value,F):
        return str(value)
    if isinstance(value,dict):
        return {str(k):serializable(v) for k,v in value.items()}
    if isinstance(value,(tuple,list)):
        return [serializable(v) for v in value]
    return value


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input',type=Path)
    parser.add_argument('--mode',choices=['paper','compact'],default='paper')
    parser.add_argument('--output',type=Path)
    parser.add_argument('--search',choices=['enumerate','branch_bound'],default='enumerate')
    args=parser.parse_args()
    data=json.loads(args.input.read_text(encoding='utf-8'),parse_float=str)
    result=solve(data['points'],data['k'],data.get('epsilon','1/2'),args.mode,args.search)
    output=json.dumps(serializable(result),indent=2)+'\n'
    if args.output:
        args.output.write_text(output,encoding='utf-8')
    else:
        print(output,end='')


if __name__=='__main__':
    main()
