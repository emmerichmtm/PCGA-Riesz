"""Reproducible sequential experiment matrix; append-only JSONL checkpoints.

python code/scale.py --smoke
python code/scale.py             # all 3 seeds; resume completed records
python code/scale.py --fresh     # distinct output directory required to rerun
"""
import argparse,json,platform,sys,os,csv
from pathlib import Path
from time import perf_counter
from fractions import Fraction as F
import numpy as np
from epas_refine_v2 import *

OUT=Path(__file__).parent.parent/'measurements'

def configs():
    runs={}
    def add(group,**kw):
        for seed in [17,29,43]:
            cfg=dict(n=100000,k=5,epsilon='1/2',kind='clusters',c=10,rho='1/4096',seed=seed);cfg.update(kw)
            key=json.dumps(cfg,sort_keys=True)
            if key not in runs:runs[key]=dict(config=cfg,groups=[])
            runs[key]['groups'].append(group)
    for n in [1000,10000,100000,1000000]:
        for kind in ['clusters','uniform']:add('input_size',n=n,kind=kind)
    for c in [5,10,25,50]:
        for rho in ['1/4096','1/1024','1/256','1/64','1/16','1/4','1']:
            add('cluster_spread',c=c,rho=rho)
    for epsilon in ['1/2','1']:
        add('accuracy',n=10000,rho='1/1024',epsilon=epsilon)
    for k in [3,4,5,6,8]:add('selected_size',n=10000,k=k,c=2*k)
    for k,n in [(3,15),(4,20),(5,25)]:add('tiny_optimum',n=n,k=k,c=2*k,rho='1/4096')
    for translation in ['1/2048','1/1024','3/2048']:
        add('grid_translation',translation=translation)
    for k in [3,4,5,6,8]:add('fixed_communities',n=10000,k=k,c=25)
    return [dict(id=f'v3-half-{i:03d}',**r) for i,r in enumerate(runs.values())]

def run_one(spec):
    cfg=spec['config'];start=perf_counter();X=dataset(cfg);gen=perf_counter()-start
    g=prepare(X,cfg['k'],cfg['epsilon']);base=solve_reduced(X,g)
    out={**spec,'dataset_sha256':fingerprint(X),'generation_validation_seconds':gen,
         'grid':public_grid(g),'base':base}
    fft=g['centers'][:cfg['k']];out['farthest_first_energy']=float(energy(X,fft))
    if base['status']=='complete':
        trace=refine(X,g,base['selected'],10);out['refine']=trace
        out['ablation']=[dict(requested_sweeps=t,**trace['trace'][min(t,len(trace['trace'])-1)]) for t in [0,1,2,3,5,10]]
        baseline=refine(X,g,fft,2,global_search=True)
        out['global_two_sweep_baseline']=baseline
        assert F(trace['trace'][-1]['energy_exact'])<=F(base['energy_exact'])
    if 'tiny_optimum' in spec['groups']:
        gg={**g,'representatives':list(range(len(X))),'M':len(X)}
        exact=solve_reduced(X,gg,timeout=30,max_M=80,max_subsets=5000000)
        out['original_optimum']=exact
        if exact['status']=='complete' and base['status']=='complete':
            opt=F(exact['energy_exact']);assert F(base['energy_exact'])<=(1+F(cfg['epsilon']))*opt
            out['base_opt_ratio']=base['energy']/float(opt)
            out['refine_opt_ratio']=out['refine']['trace'][-1]['energy']/float(opt)
    out['wall_seconds']=perf_counter()-start
    return out

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--smoke',action='store_true');ap.add_argument('--fresh',action='store_true');a=ap.parse_args()
    OUT.mkdir(exist_ok=True);specs=configs()
    name='smoke-v2.jsonl' if a.smoke else 'results-v2.jsonl';path=OUT/name
    if a.fresh and path.exists():raise SystemExit('Preserve existing measurements; use a separate output directory for fresh runs.')
    if a.smoke:specs=[specs[0],specs[1],next(r for r in specs if 'tiny_optimum' in r['groups'])]
    (OUT/'configs-v2.json').write_text(json.dumps(configs(),indent=2))
    meta=dict(python=sys.version,numpy=np.__version__,platform=platform.platform(),processor=platform.processor(),
        cpu_count=os.cpu_count(),denominator=D,seeds=[17,29,43],clock='time.perf_counter',
        search_budget=dict(max_M=300,max_subsets=10000000000,seconds=30),sequential=True,
        arithmetic='Exact int64 geometry; Fraction search; outward binary64 intervals with exact Fraction fallback; zero tolerance')
    (OUT/'environment-v2.json').write_text(json.dumps(meta,indent=2))
    done={json.loads(line)['id'] for line in path.read_text().splitlines()} if path.exists() else set()
    for j,spec in enumerate(specs):
        if spec['id'] in done:continue
        result=run_one(spec)
        with path.open('a') as f:f.write(json.dumps(result)+'\n');f.flush()
        b=result['base'];tr=result.get('refine',{}).get('trace',[{'energy':b.get('energy',0)}])
        gain=1-tr[-1]['energy']/b['energy'] if b['status']=='complete' else None
        print(json.dumps(dict(index=j+1,total=len(specs),id=spec['id'],n=spec['config']['n'],k=spec['config']['k'],
            c=spec['config']['c'],rho=spec['config']['rho'],M=result['grid']['M'],status=b['status'],gain=gain,
            seconds=round(result['wall_seconds'],3))),flush=True)
    print('Finished',path,flush=True)

if __name__=='__main__':main()
