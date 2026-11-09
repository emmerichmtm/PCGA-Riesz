"""Planar s=2 EPAS and fixed-cell refinement on bounded dyadic input.

The mathematical algorithm is the Fraction reference in planar_riesz_epas_v2.
This NumPy implementation uses exact int64 squared distances/grid keys, exact
rational search, and outward binary64 intervals with Fraction fallback for
within-cell minimization. No energy tolerance is used. See README.md.
"""
from fractions import Fraction as F
from itertools import combinations
from math import comb
from time import perf_counter
import hashlib
import numpy as np
from planar_riesz_epas_v2 import dyadic_inverse_sqrt
from exact_subset_search_v2 import minimize

D = 1 << 20

def check_points(X):
    X=np.asarray(X)
    if X.dtype.kind not in 'iu': raise ValueError('Supply integer numerators, not floating coordinates')
    if X.size==0: raise ValueError('Expected nonempty planar input')
    if int(X.min())<-(2**63) or int(X.max())>=2**63: raise ValueError('Coordinates must fit in int64')
    X=X.astype(np.int64,copy=False)
    if X.ndim!=2 or X.shape[1]!=2: raise ValueError('Expected planar integer coordinates')
    span=[int(X[:,a].max())-int(X[:,a].min()) for a in (0,1)]
    if sum(z*z for z in span)>2**53: raise ValueError('Squared distances must fit exactly in binary64')
    if len(np.unique(X,axis=0))!=len(X): raise ValueError('Candidates must be distinct')
    return X

def distances2(X,p):
    z=X-p
    return z[:,0]*z[:,0]+z[:,1]*z[:,1]

def energy(X,selected,denominator=D):
    return sum((F(denominator**2,int(np.dot(X[i]-X[j],X[i]-X[j])))
                for i,j in combinations(selected,2)),F(0))

def prepare(X,k,epsilon,denominator=D,center_count=None):
    """Requires validated distinct points. Returns grouping without n^2 weights."""
    start=perf_counter(); n=len(X); eps=min(F(epsilon),F(1))
    if not 2<=k<n or eps<=0: raise ValueError('Require 2<=k<n and epsilon>0')
    m=min(n,2*k if center_count is None else center_count)
    if m<2*k-2: raise ValueError('Too few centers for the exchange proof')
    centers=[0]; nearest=distances2(X,X[0]); radii=[]
    while len(centers)<m:
        i=int(np.argmax(nearest)); r2=int(nearest[i]); centers.append(i); radii.append(r2)
        nearest=np.minimum(nearest,distances2(X,X[i]))
    H=sum((F(1,j) for j in range(1,m+1)),F(0))
    scale=dyadic_inverse_sqrt(16*H*denominator**2/r2)
    h=eps*scale/16; h_units=h*denominator
    a,b=h_units.numerator,h_units.denominator
    max_abs=max(abs(int(X.min())),abs(int(X.max())))
    if 2*max_abs*b+a>=2**63: raise OverflowError('Grid integer expression exceeds int64')
    keys=(2*X*b+a)//(2*a)
    order=np.lexsort((np.arange(n),keys[:,1],keys[:,0]))
    sorted_keys=keys[order]
    starts=np.flatnonzero(np.r_[True,np.any(sorted_keys[1:]!=sorted_keys[:-1],axis=1)])
    ends=np.r_[starts[1:],n]
    reps=order[starts]; cell_id=np.empty(n,dtype=np.int64)
    cell_id[order]=np.repeat(np.arange(len(starts)),ends-starts)
    assert len(reps)>=k
    return dict(n=n,k=k,epsilon=str(eps),denominator=denominator,centers=centers,
        r2_units=r2,cover2_units=int(nearest.max()),h=str(h),H=str(H),scale=str(scale),
        keys=keys,order=order,starts=starts,ends=ends,cell_id=cell_id,
        representatives=sorted(map(int,reps)),M=len(reps),prep_seconds=perf_counter()-start)

def solve_reduced(X,g,timeout=30,max_M=300,max_subsets=10000000000):
    start=perf_counter(); ids=g['representatives']; k=g['k']; possibilities=comb(len(ids),k)
    out=dict(possible_subsets=possibilities,search_seconds=0.0,status='skipped_search_budget',
             search_budget=dict(max_M=max_M,max_subsets=max_subsets,seconds=timeout))
    if len(ids)>max_M or possibilities>max_subsets:return out
    weights={(i,j):F(g['denominator']**2,int(np.dot(X[i]-X[j],X[i]-X[j])))
             for i,j in combinations(ids,2)}
    try:
        selected,value,stats=minimize(weights,ids,k,timeout=max(.001,timeout-(perf_counter()-start)))
    except TimeoutError:
        out.update(status='timeout',search_seconds=perf_counter()-start);return out
    out.update(status='complete',selected=list(selected),energy_exact=str(value),energy=float(value),
        certificate_lower_exact=str(value/(1+F(g['epsilon']))),search_stats=stats,
        search_seconds=perf_counter()-start)
    return out

def best_candidate(X,candidates,others,incumbent,denominator=D):
    """Certified exact minimization, interval filter only discards provable losers.

    d2 is an exactly representable positive integer. nextafter encloses each
    rounded division and addition. Every true minimizer has lower<=min(upper).
    Evaluate all such survivors as Fractions. Incumbent wins every exact tie.
    """
    candidates=np.asarray(candidates,dtype=np.int64)
    lo=np.zeros(len(candidates));hi=np.zeros(len(candidates))
    for j in others:
        d2=distances2(X[candidates],X[j]); assert np.all(d2>0)
        val=1.0/d2.astype(np.float64)
        lo=np.nextafter(lo+np.nextafter(val,-np.inf),-np.inf)
        hi=np.nextafter(hi+np.nextafter(val,np.inf),np.inf)
    survivors=candidates[lo<=np.min(hi)]
    def potential(i):
        return sum((F(denominator**2,int(np.dot(X[i]-X[j],X[i]-X[j]))) for j in others),F(0))
    old=potential(incumbent);best=incumbent;value=old
    for i in sorted(map(int,survivors)):
        v=potential(i)
        if v<value or (v==value and value<old and i<best):best,value=i,v
    return best,value,old,len(survivors)

def refine(X,g,selected,sweeps=10,global_search=False,timeout=None):
    """Strict sequential minimization. global_search=True is a heuristic baseline.

    Only the default fixed-cell variant is called EPAS-Refine. The global
    baseline can start from any feasible solution and has no standalone bound.
    """
    selected=list(map(int,selected)); fixed=[int(g['cell_id'][i]) for i in selected]
    if not global_search: assert len(set(fixed))==len(selected)
    value=energy(X,selected,g['denominator']);start=perf_counter();updates=[]
    trace=[dict(sweep=0,energy_exact=str(value),energy=float(value),seconds=0.0,
                candidate_evaluations=0,pair_evaluations=0,accepted=0,exact_survivors=0)]
    candidates_total=pair_total=accepted_total=survivor_total=0; first_no_improvement=None
    stop_reason='sweep_cap'
    for t in range(1,sweeps+1):
        if timeout is not None and perf_counter()-start>=timeout:
            stop_reason='time_budget_at_sweep_boundary';break
        sweep_accept=0
        for pos in range(g['k']):
            update_start=perf_counter();old_id=selected[pos];others=selected[:pos]+selected[pos+1:]
            if global_search:
                mask=np.ones(len(X),dtype=bool);mask[others]=False;candidates=np.flatnonzero(mask)
            else:
                cell=fixed[pos];candidates=g['order'][g['starts'][cell]:g['ends'][cell]]
            best,newpot,oldpot,survivors=best_candidate(X,candidates,others,old_id,g['denominator'])
            old_value=value;accepted=(newpot<oldpot)
            if accepted:selected[pos]=best;value+=newpot-oldpot;sweep_accept+=1
            else: assert best==old_id
            assert value<=old_value and len(set(selected))==g['k']
            if not global_search: assert [int(g['cell_id'][i]) for i in selected]==fixed
            candidates_total+=len(candidates);pair_total+=len(candidates)*(g['k']-1)
            accepted_total+=int(accepted);survivor_total+=survivors
            updates.append(dict(sweep=t,position=pos,old=old_id,new=selected[pos],accepted=accepted,
                energy_exact=str(value),energy=float(value),candidate_evaluations=len(candidates),
                exact_survivors=survivors,seconds=perf_counter()-update_start))
        assert value==energy(X,selected,g['denominator'])
        trace.append(dict(sweep=t,energy_exact=str(value),energy=float(value),seconds=perf_counter()-start,
            candidate_evaluations=candidates_total,pair_evaluations=pair_total,accepted=accepted_total,
            exact_survivors=survivor_total))
        if sweep_accept==0:first_no_improvement=t;stop_reason='no_improvement';break
    return dict(selected=selected,trace=trace,updates=updates,requested_sweeps=sweeps,
        actual_sweeps=len(trace)-1,first_no_improvement=first_no_improvement,
        early_stop=first_no_improvement is not None,M=g['M'],fixed_cells=fixed,
        mode='global_heuristic' if global_search else 'fixed_cells',stop_reason=stop_reason)

def dataset(config):
    """Balanced clusters on a unit square lattice; unique dyadic candidates."""
    n=config['n'];rng=np.random.default_rng(config['seed']);kind=config['kind'];c=config.get('c',10)
    width=int(np.ceil(np.sqrt(c)))
    if kind=='uniform':
        # Sampling linear indices without replacement is uniform on the dyadic box.
        side=width*D;idx=rng.choice(side*side,size=n,replace=False)
        X=np.column_stack((idx%side,idx//side)).astype(np.int64)
    else:
        radius=F(config['rho']);R=int(radius*D);side=2*R+1;chunks=[]
        for i in range(c):
            count=n//c+(i<n%c)
            idx=rng.choice(side*side,size=count,replace=False)
            chunks.append(np.column_stack((idx%side-R+(i%width)*D,idx//side-R+(i//width)*D)))
        X=np.concatenate(chunks).astype(np.int64)
        # Overlapping clusters may duplicate lattice points: fill deterministically.
        X=np.unique(X,axis=0)
        while len(X)<n:
            extra=rng.integers(-R,width*D+R+1,size=(n-len(X),2),dtype=np.int64)
            X=np.unique(np.concatenate((X,extra)),axis=0)
        X=X[rng.permutation(n)]
    shift=config.get('translation','0')
    X=X+int(F(shift)*D)
    return check_points(X)

def fingerprint(X):return hashlib.sha256(X.astype('<i8',copy=False).tobytes()).hexdigest()

def public_grid(g):
    return {k:v for k,v in g.items() if k not in ('keys','order','starts','ends','cell_id')}
