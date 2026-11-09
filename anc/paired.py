"""Paired, reproducible planar certificate comparison at epsilon one half.

All methods use the same input order, FFT prefix, floor origin, representative
rule and exact search. The common physical origin is specified in units of the exchange-grid width.
The known-size constants, rather than just their asymptotic orders, are tested.
"""
from pathlib import Path
from fractions import Fraction as F
from itertools import combinations
from math import comb
from time import perf_counter
import json, hashlib, csv, argparse, sys
import numpy as np
from epas_refine_v2 import (D,check_points,distances2,energy,solve_reduced,refine)
from planar_riesz_epas_v2 import dyadic_inverse_sqrt, enumerate_best

OUT=Path(__file__).resolve().parents[1]/'measurements'/'paired'
METHODS=('pair_count','packing','incumbent','exchange','hybrid')

def make_input(cfg):
    rng=np.random.default_rng(cfg['seed']); n=cfg['n']; c=cfg['c']
    if cfg['family']=='uniform':
        X=rng.uniform(0,np.sqrt(c),size=(n,2))
    else:
        w=int(np.ceil(np.sqrt(c)))
        centers=np.array([(i%w,i//w) for i in range(c)],dtype=float)
        if cfg['family']!='lattice':
            centers+=rng.uniform(-.3,.3,size=(c,2))
        counts=np.full(c,n//c);counts[:n%c]+=1
        if cfg['family']=='unequal':
            counts=np.ones(c,dtype=int)+rng.multinomial(n-c,rng.dirichlet(np.ones(c)))
        chunks=[]
        for center,count in zip(centers,counts):
            radius=int(F(cfg['rho'])*D)//4;side=2*radius+1
            sampled=rng.choice(side*side,size=int(count),replace=False)
            offset=np.column_stack((sampled%side-radius,sampled//side-radius))*4/D
            if cfg['family']=='elongated':offset[:,0]*=8
            chunks.append(center+offset)
        X=np.concatenate(chunks)
        if cfg['family'] in ('rotated','elongated'):
            theta=.37; rot=np.array([[np.cos(theta),-np.sin(theta)],[np.sin(theta),np.cos(theta)]])
            X=X@rot.T
    X=np.rint(X*D).astype(np.int64)
    # The four-unit offset lattice prevents within-cluster collisions after
    # rotation and rounding. A final check rejects any accidental collision.
    X=X[rng.permutation(n)]
    return check_points(X)

def certificates(X,k):
    start=perf_counter(); centers=[0];nearest=distances2(X,X[0]);radii=[]
    while len(centers)<2*k:
        idx=int(np.argmax(nearest));radii.append(int(nearest[idx]));centers.append(idx)
        nearest=np.minimum(nearest,distances2(X,X[idx]))
    fft_seconds=perf_counter()-start
    rk2=F(radii[k-2],D*D);r2k2=F(radii[-1],D*D)
    H=lambda m:sum((F(1,j) for j in range(1,m+1)),F(0))
    start=perf_counter();A=1<<48
    u=F(sum((A*D*D+int(distances2(X[[i]],X[j])[0])-1)//int(distances2(X[[i]],X[j])[0]) for i,j in combinations(centers[:k],2)),A)
    energy_seconds=perf_counter()-start
    pc=F(comb(k,2))/rk2;pk=min(pc,4*k*H(k-1)/rk2);ex=16*H(2*k)/r2k2
    # Bounded-denominator upper sum differs from U by < binom(k,2)/2^48.
    u=min(u,pk)
    if k<=8:assert energy(X,centers[:k])<=u
    return dict(pair_count=pc,packing=pk,incumbent=u,exchange=ex,hybrid=min(u,ex)),dict(
        centers=centers,rk2=str(rk2),r2k2=str(r2k2),fft_seconds=fft_seconds,incumbent_seconds=energy_seconds)

def grid(X,k,eps,B,origin='0'):
    start=perf_counter();scale=dyadic_inverse_sqrt(B);h=F(eps)*scale/16
    units=h*D;off=F(origin)*D;a,b=units.numerator,units.denominator
    # floor((X/D-origin)/h), with a shared physical origin for all methods.
    den=a*off.denominator
    mul=b*off.denominator;shift=b*off.numerator
    if max(abs(int(X.min())),abs(int(X.max())))*mul+abs(shift)>=2**63:
        raise OverflowError('Grid intermediate out of checked int64 range')
    keys=(X*mul-shift)//den
    order=np.lexsort((np.arange(len(X)),keys[:,1],keys[:,0]));sk=keys[order]
    starts=np.flatnonzero(np.r_[True,np.any(sk[1:]!=sk[:-1],axis=1)])
    ends=np.r_[starts[1:],len(X)];cell_id=np.empty(len(X),dtype=np.int64)
    cell_id[order]=np.repeat(np.arange(len(starts)),ends-starts)
    reps=sorted(map(int,order[starts]));assert len(reps)>=k
    return dict(n=len(X),k=k,epsilon=str(F(eps)),denominator=D,M=len(reps),
        representatives=reps,keys=keys,order=order,starts=starts,ends=ends,
        cell_id=cell_id,scale=str(scale),h=str(h),grid_seconds=perf_counter()-start)

def configurations():
    configs=[]
    def add(**kw):configs.append(dict(**kw,id=f'paired-v3-half-{len(configs):03d}'))
    for family in ['lattice','irregular','rotated','unequal','elongated','uniform']:
        for rho in ['1/4096','1/1024','1/256']:
            for seed in [17,29,43]:
                add(n=2000,k=5,c=10,epsilon='1/2',family=family,rho=rho,seed=seed,offset='0',group='geometry',search=True)
    for k in [8,16,32,64,128,256,512]:
        for family in ['lattice','irregular']:
            for seed in [17,29,43]:
                add(n=100000,k=k,c=2*k,epsilon='1/2',family=family,rho='1/256',seed=seed,offset='0',group='cardinality',search=False)
    for offset in ['1/4','1/2','3/4']:
        for seed in [17,29,43]:
            add(n=2000,k=5,c=10,epsilon='1/2',family='irregular',rho='1/1024',seed=seed,offset=offset,group='origin',search=True)
    for family in ['irregular','elongated','uniform']:
        for seed in [17,29,43]:
            add(n=24,k=4,c=8,epsilon='1/2',family=family,rho='1/32',seed=seed,offset='0',group='exact',search=True)
    return configs

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--smoke',action='store_true');args=ap.parse_args()
    OUT.mkdir(parents=True,exist_ok=True);cfgs=configurations()
    if args.smoke:cfgs=[cfgs[0],cfgs[-1]]
    (OUT/('smoke-configs-v2.json' if args.smoke else 'configs-v2.json')).write_text(json.dumps(cfgs,indent=2))
    path=OUT/('smoke-results-v2.jsonl' if args.smoke else 'results-v2.jsonl')
    existing={json.loads(l)['id'] for l in path.read_text().splitlines()} if path.exists() else set()
    for cfg in cfgs:
        if cfg['id'] in existing:continue
        X=make_input(cfg);bounds,meta=certificates(X,cfg['k'])
        opt=None
        if cfg['group']=='exact':
            pts=[tuple(F(int(z),D) for z in p) for p in X]
            _,opt,_=enumerate_best(pts,range(len(X)),cfg['k'])
        rows=[]
        # Rotating execution order avoids assigning all cold-start effects to PC.
        shift=cfg['seed']%len(METHODS);order=METHODS[shift:]+METHODS[:shift]
        for method in order:
            origin=F(cfg['offset'])*F(cfg['epsilon'])*dyadic_inverse_sqrt(bounds['exchange'])/16
            g=grid(X,cfg['k'],cfg['epsilon'],bounds[method],origin)
            row=dict(method=method,B=str(bounds[method]),scale=g['scale'],h=g['h'],M=g['M'],
                grid_seconds=g['grid_seconds'],origin=str(origin),representatives=g['representatives'])
            if cfg['search']:
                result=solve_reduced(X,g,timeout=5,max_M=160,max_subsets=10000000)
                row.update(result)
                if result['status']=='complete':
                    tr=refine(X,g,result['selected'],sweeps=2)
                    row.update(refined_energy_exact=tr['trace'][-1]['energy_exact'],refine_seconds=tr['trace'][-1]['seconds'])
                    if opt is not None:
                        assert F(result['energy_exact'])<=(1+F(cfg['epsilon']))*opt
                        row.update(optimum_exact=str(opt),base_opt_ratio=float(F(result['energy_exact'])/opt),refine_opt_ratio=float(F(row['refined_energy_exact'])/opt))
            else:row.update(status='occupancy_only',possible_subsets=comb(g['M'],cfg['k']))
            rows.append(row)
        record=dict(id=cfg['id'],config=cfg,shared=meta,input_sha256=hashlib.sha256(X.astype('<i8').tobytes()).hexdigest(),methods=rows)
        with path.open('a') as f:f.write(json.dumps(record)+'\n')
        print(cfg['id'],cfg['family'],cfg['k'],[(r['method'],r['M'],r['status']) for r in rows],flush=True)
    print('Finished',len(cfgs),'paired inputs',flush=True)

if __name__=='__main__':main()
