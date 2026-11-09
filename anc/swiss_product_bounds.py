"""Certified product search for s=2, using integer unary lower bounds."""
from fractions import Fraction
from itertools import combinations
from math import prod
import numpy as np

def solve_products(x,d2,groups,root_heap,best,value,energy_fn,required_ids):
    # Each integer pair weight is at most 2**60 times its true contribution.
    scale=1<<60
    weight=np.maximum(d2,1)
    np.floor_divide(scale,weight,out=weight)
    required=None if required_ids is None else set(map(int,required_ids))
    k=len(best)
    pairs=list(combinations(range(k),2))
    stats=dict(visited=0,pruned=0,exact_leaves=0)
    def ceiling():return (value.numerator*scale+value.denominator-1)//value.denominator
    def consider(assignment):
        nonlocal best,value
        trial=tuple(sorted(map(int,assignment)))
        if required is not None and not required.intersection(trial):return
        candidate=energy_fn(x,trial,1,2)
        if (candidate,trial)<(value,best):value,best=candidate,trial
    def split_pair(a,b):
        w=weight[np.ix_(a,b)]
        # w_ab(i,j) >= row_min(i) + min_i(w_ab(i,j)-row_min(i)).
        left=w.min(axis=1)
        right=(w-left[:,None]).min(axis=0)
        return left,right
    root_pair={}
    for i,j in combinations(range(len(groups)),2):
        a,b=groups[i],groups[j]
        # All subsequent vector sums stay in int64 for this geographical input.
        limit=int(weight[np.ix_(a,b)].max())*k*k
        if limit>=2**63:raise OverflowError("Increase group separation or use larger integer accumulators")
        root_pair[i,j]=split_pair(a,b)
    def visit(domains,root_ids=None):
        stats["visited"]+=1
        if required is not None and not any(required.intersection(a) for a in domains):
            stats["pruned"]+=1;return
        unary=[np.zeros(len(a),dtype=np.int64) for a in domains]
        for i,j in pairs:
            if root_ids is None:a,b=split_pair(domains[i],domains[j])
            else:a,b=root_pair[root_ids[i],root_ids[j]]
            unary[i]+=a;unary[j]+=b
        minima=[int(v.min()) for v in unary];lower=sum(minima)
        slack=ceiling()-lower
        if slack<0:stats["pruned"]+=1;return
        # For any fixed value of variable i, its unary excess above its minimum
        # adds to the lower bound. This deletion is therefore safe at equality.
        reduced=[a[u-m<=slack] for a,u,m in zip(domains,unary,minima)]
        if all(len(a)==1 for a in reduced):
            stats["exact_leaves"]+=1;consider([a[0] for a in reduced]);return
        if root_ids is not None:
            # An optional feasible start inside this product. Integer weights
            # here guide the starting choice. Only exact energy updates the bound.
            start=[a[int(u.argmin())] for a,u in zip(domains,unary)]
            for _ in range(2):
                for i,a in enumerate(reduced):
                    others=start[:i]+start[i+1:]
                    start[i]=a[int(weight[np.ix_(a,others)].sum(axis=1).argmin())]
            consider(start)
        if prod(map(len,reduced))*2<prod(map(len,domains)):
            visit(reduced);return
        pos=max(range(k),key=lambda i:np.ptp(x[reduced[i]],axis=0).max())
        a=reduced[pos];axis=int(np.ptp(x[a],axis=0).argmax())
        a=a[np.lexsort((a,x[a,axis]))];mid=len(a)//2
        assert mid>0
        for child in (a[:mid],a[mid:]):
            visit(reduced[:pos]+[child]+reduced[pos+1:])
    for _,root_ids in sorted(root_heap):
        visit([groups[i] for i in root_ids],root_ids)
    return best,value,stats
