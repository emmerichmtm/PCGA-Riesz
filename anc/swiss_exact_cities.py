"""Exact search accelerator for separated city groups (no approximation).

The EPAS construction, objective and fixed-cell refinement are unchanged.
For this input a same-city pair alone costs more than a feasible k-set;
the assertion below proves that at most one selected point belongs to a city.
Enumerate every k-city combination, then bisect Cartesian products of their
candidate sets. Each pair contributes floor(scale*denominator**s / max_squared_distance**(s/2)),
a rigorous lower bound. Prune only above the exact incumbent's ceiling.
Surviving singleton products are compared as Fractions, including exact ties.
"""
from fractions import Fraction as F
from itertools import combinations
from math import comb
from heapq import heappush, heappop
from time import perf_counter
import numpy as np
def energy(X, selection, denominator, exponent):
    half = exponent // 2
    total = F(0)
    for i,j in combinations(selection,2):
        delta = X[i] - X[j]
        total += F(denominator**2,int(delta @ delta))**half
    return total



def minimize_cities(X,ids,k,city_ids,denominator=1000,exponent=2,required_ids=None):
    if exponent not in (2,4): raise ValueError("This runner supports s=2 and s=4")
    half=exponent//2
    scale=1<<(40*exponent)
    start=perf_counter()
    X=np.asarray(X)
    if X.ndim!=2 or X.shape[1]!=2 or X.dtype!=np.int64:
        raise ValueError("Expected an int64 planar array of millimetre coordinates")
    span=[int(X[:,a].max())-int(X[:,a].min()) for a in (0,1)]
    if sum(width**2 for width in span)>=2**63:
        raise OverflowError("Squared distances exceed int64")
    ids=sorted(map(int,ids))
    d2=np.empty((len(X),len(X)),dtype=np.int64)
    for start_row in range(0,len(X),48):
        block=X[start_row:start_row+48]
        dx=block[:,None,0]-X[None,:,0];dy=block[:,None,1]-X[None,:,1]
        d2[start_row:start_row+len(block)]=dx*dx+dy*dy
    required=set(ids) if required_ids is None else set(ids).intersection(required_ids)
    if not required: raise ValueError("No eligible required site")
    # Floating greedy supplies only a feasible upper bound, never a prune.
    w=np.zeros(d2.shape);np.divide(float(denominator**2),d2,out=w,where=d2>0)
    w **= half
    np.fill_diagonal(w,np.inf)
    seed=ids[:k]
    if not required.intersection(seed): seed[-1]=min(required)
    best=tuple(sorted(seed));value=energy(X,best,denominator,exponent)
    id_array=np.asarray(ids)
    for first in ids[::max(1,len(ids)//256)]:
        selected=[first];potential=w[first].copy()
        while len(selected)<k:
            i=ids[int(np.argmin(potential[id_array]))];selected.append(i)
            potential+=w[i]
        trial=tuple(sorted(selected))
        if not required.intersection(trial):continue
        v=energy(X,trial,denominator,exponent)
        if (v,trial)<(value,best):best,value=trial,v
    if k>=8:
        selected=list(best)
        for sweep in range(10):
            old=selected[:]
            for pos in range(k):
                others=selected[:pos]+selected[pos+1:]
                allowed=ids if required.intersection(others) else sorted(required)
                costs=w[np.ix_(allowed,others)].sum(axis=1)
                selected[pos]=allowed[int(np.argmin(costs))]
            if old==selected:break
        trial=tuple(sorted(selected));v=energy(X,trial,denominator,exponent)
        if (v,trial)<(value,best):best,value=trial,v
    del w
    bound=(value.numerator*scale+value.denominator-1)//value.denominator
    cities=sorted(set(int(city_ids[i]) for i in ids))
    groups=[np.array([i for i in ids if city_ids[i]==c],dtype=np.int64) for c in cities]
    same_city=[]
    for a in groups:
        if len(a)>1:
            lower=F(denominator**2,int(d2[np.ix_(a,a)].max()))**half
            assert lower>value, 'Cannot exclude two points in a city on this input'
            same_city.append(lower)
    nodes=[];children={};diameters=[]
    def node(a):
        idx=len(nodes);nodes.append(a)
        diameters.append(int(d2[np.ix_(a,a)].max()))
        return idx
    roots=[node(a) for a in groups]
    eligible_roots={r for r in roots if required.intersection(nodes[r])}
    def split(t):
        if t not in children:
            a=nodes[t];axis=int(np.argmax(np.ptp(X[a],axis=0)))
            a=a[np.lexsort((a,X[a,axis]))];m=len(a)//2
            children[t]=(node(a[:m]),node(a[m:]))
        return children[t]
    cache={}
    def edge(a,b):
        key=(min(a,b),max(a,b))
        if key not in cache:
            far=int(d2[np.ix_(nodes[a],nodes[b])].max())
            cache[key]=(denominator**exponent*scale)//(far**half)
        return cache[key]
    def lower(ts):return sum(edge(a,b) for a,b in combinations(ts,2))
    heap=[];queued=0;pruned=0;leaves=0;visited=0
    # Bound remaining distinct-region choices before expanding their whole subtree.
    # For each eligible region add the (t-1) cheapest edges to other eligible
    # regions. Each internal edge contributes twice in the resulting vertex sum.
    root_count=len(roots)
    edge_rows=[[0 if i==j else edge(i,j) for j in roots] for i in roots]
    internal={}
    for start_id in range(root_count):
        for i in range(start_id,root_count):
            candidates=sorted(edge_rows[i][j] for j in range(start_id,root_count) if j!=i)
            for t in range(1,min(k,root_count-start_id)+1):
                internal[i,start_id,t]=sum(candidates[:t-1])
    rejected_roots=0
    def enqueue(prefix,start_id,t,cost,incident):
        nonlocal rejected_roots, queued
        if t==0:
            if cost<=bound:
                heappush(heap,(cost,prefix));queued+=1
            else:rejected_roots+=1
            return
        remaining=root_count-start_id
        if required_ids is not None and not eligible_roots.intersection(prefix) and not any(r>=start_id for r in eligible_roots):
            rejected_roots+=comb(remaining,t);return
        doubled=sorted(2*incident[i]+internal[i,start_id,t] for i in range(start_id,root_count))
        if 2*cost+sum(doubled[:t])>2*bound:
            rejected_roots+=comb(remaining,t)
            return
        for i in range(start_id,root_count-t+1):
            extended=cost+incident[i]
            if extended>bound:
                rejected_roots+=comb(root_count-i-1,t-1)
                continue
            if t==1:
                if required_ids is not None and not eligible_roots.intersection(prefix+(i,)):
                    rejected_roots+=1;continue
                heappush(heap,(extended,prefix+(i,)));queued+=1
                continue
            new_incident=[incident[j]+edge_rows[i][j] for j in roots]
            enqueue(prefix+(i,),i+1,t-1,extended,new_incident)
    enqueue((),0,k,0,[0]*root_count)
    assert rejected_roots+queued==comb(root_count,k)
    pruned+=rejected_roots
    initial_products=queued
    if k>=8 and exponent==2 and denominator==1:
        from product_bounds import solve_products
        best,value,st=solve_products(X,d2,groups,heap,best,value,energy,required_ids)
        visited=st["visited"];pruned+=st["pruned"];leaves=st["exact_leaves"]
        heap=[]
    while heap:
        b,ts=heappop(heap);visited+=1
        if b>bound:pruned+=1;continue
        if required_ids is not None and not any(required.intersection(nodes[t]) for t in ts):
            pruned+=1;continue
        if all(len(nodes[t])==1 for t in ts):
            leaves+=1;trial=tuple(sorted(int(nodes[t][0]) for t in ts))
            v=energy(X,trial,denominator,exponent)
            if (v,trial)<(value,best):
                best,value=trial,v
                bound=(value.numerator*scale+value.denominator-1)//value.denominator
            continue
        pos=max(range(k),key=lambda p:diameters[ts[p]])
        for child in split(ts[pos]):
            nxt=ts[:pos]+(child,)+ts[pos+1:];lb=lower(nxt)
            if lb<=bound:heappush(heap,(lb,nxt));queued+=1
            else:pruned+=1
    return dict(status='complete',selected=list(best),energy_exact=str(value),energy=float(value),
        seconds=perf_counter()-start,possible_subsets=comb(len(ids),k),
        city_products=comb(len(cities),k),initial_surviving_products=initial_products,
        visited_products=visited,pruned_products=pruned,exact_leaves=leaves,
        minimum_same_city_pair_bound=str(min(same_city)) if same_city else None,
        integer_scale=scale,unary_bound_scale=(1<<60) if k>=8 and exponent==2 and denominator==1 else None,method='exact Cartesian-product branch-and-bound')
