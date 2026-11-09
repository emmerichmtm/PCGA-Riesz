"""Exact cardinality search with integer lower bounds (standard library only).

For every positive rational pair weight w, floor(2**48*w) is a lower bound
in scaled units. All pruning is integer arithmetic; surviving leaves are
compared as Fractions. Thus this is an exact accelerator, not a heuristic.
"""
from fractions import Fraction as F
from itertools import combinations
from math import comb

SCALE = 1 << 48

def minimize(weights, ids, k, incumbent=None):
    ids=tuple(sorted(ids))
    if k==0:
        return (),F(0),{'nodes':1,'pruned':0,'exact_leaves':1,'possible_subsets':1}
    exact=lambda s:sum((weights[min(i,j),max(i,j)] for i,j in combinations(s,2)),F(0))
    best=tuple(sorted(incumbent)) if incumbent is not None else ids[:k]
    value=exact(best)
    # A deterministic greedy completion supplies only an upper bound.
    for first in ids:
        trial=[first]
        while len(trial)<k:
            trial.append(min((i for i in ids if i not in trial),
                             key=lambda i:(sum(weights[min(i,j),max(i,j)] for j in trial),i)))
        trial=tuple(sorted(trial)); v=exact(trial)
        if (v,trial)<(value,best):value,best=v,trial
    bound=(value.numerator*SCALE+value.denominator-1)//value.denominator
    w={ij:(x.numerator*SCALE)//x.denominator for ij,x in weights.items()}
    edge=lambda i,j:w[min(i,j),max(i,j)]
    stats={'nodes':0,'pruned':0,'exact_leaves':0,'possible_subsets':comb(len(ids),k),
           'integer_scale':SCALE,'pruning':'integer floors below exact rational weights'}

    def visit(selected,candidates,cost,attached):
        nonlocal best,value,bound
        stats['nodes']+=1
        q=k-len(selected)
        if len(candidates)<q:return
        if q==0:
            stats['exact_leaves']+=1
            val=exact(selected)
            if (val,selected)<(value,best):
                value,best=val,selected
                bound=(value.numerator*SCALE+value.denominator-1)//value.denominator
            return
        # For each possible new vertex, charge its links to the partial set
        # twice and its q-1 cheapest available links once. Any completion's
        # doubled remaining cost dominates the sum of its q vertex charges,
        # hence also the q smallest charges over all candidates.
        if q>=2:
            scores=[]
            for i in candidates:
                incident=sorted(edge(i,j) for j in candidates if j!=i)
                scores.append(2*attached[i]+sum(incident[:q-1]))
            if 2*cost+sum(sorted(scores)[:q])>2*bound:
                stats['pruned']+=1; return
        elif cost+min(attached[i] for i in candidates)>bound:
            stats['pruned']+=1; return
        for pos,i in enumerate(candidates):
            if len(candidates)-pos<q:break
            next_cost=cost+attached[i]
            if next_cost>bound:
                stats['pruned']+=1; continue
            rest=candidates[pos+1:]
            new_att={j:attached[j]+edge(i,j) for j in rest}
            rest=tuple(j for j in rest if next_cost+new_att[j]<=bound)
            visit(selected+(i,),rest,next_cost,new_att)

    visit((),ids,0,{i:0 for i in ids})
    return best,value,stats
