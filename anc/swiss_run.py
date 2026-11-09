"""Recompute the dense input exactly; verify the two grids and their refinement."""
from pathlib import Path
from fractions import Fraction as F
from itertools import combinations
import csv, json, hashlib
import numpy as np
from exact_cities import minimize_cities, energy
from epas_methods import geometry, refine

root=Path(__file__).resolve().parent
def load(name):return json.loads((root/name).read_text(encoding='utf-8'))
def rows(name):
    with (root/name).open(encoding='utf-8-sig') as f:return list(csv.DictReader(f))

def main():
    data=rows('candidate-sites.csv')
    assert len(data)==5403
    assert hashlib.sha256((root/'candidate-sites.csv').read_bytes()).hexdigest()==load('generation.json')['output_sha256']
    xy=np.array([[int(r['east_lv95_m']),int(r['north_lv95_m'])] for r in data],dtype=np.int64)
    labels=np.unique([r['locality'] for r in data],return_inverse=True)[1]
    mid=np.array(load('central-region.json')['centroid_lv95_m'])
    centre_distance=np.linalg.norm(xy-mid,axis=1)/1000
    def describe(ids):return [dict(id=data[i]['site_id'],locality=data[i]['locality'],centre_km=float(centre_distance[i])) for i in ids]
    # Validate exact search independently by all 10-subsets of 14 sites,
    # including additions and two candidates from the same locality.
    sample=[82,270,378,800,1047,1197,1300,1348,1583,1711]+[1801,1852,2501,5251]
    assert len(set(sample))==14
    sm=xy[sample];smg=labels[sample]
    brute=min((energy(sm,c,1,2),c) for c in combinations(range(14),10))
    chk=minimize_cities(sm,range(14),10,smg,1,2)
    assert (F(chk['energy_exact']),tuple(chk['selected']))==brute
    geo,plans=geometry(xy,F(1))
    print('Counts:',[(name,str(width),len(c)) for name,(width,c) in plans.items()],flush=True)
    full=minimize_cities(xy,range(len(xy)),10,labels,1,2)
    opt=F(full['energy_exact'])
    assert opt==energy(xy,full['selected'],1,2)
    print('Full optimum:',float(opt),flush=True)
    out=dict(n=len(xy),k=10,s=2,epsilon='1',geometry=geo,full=full,selected=describe(full['selected']),methods={})
    for name,(h,cells) in plans.items():
        reps=sorted(cells)
        red=minimize_cities(xy,reps,10,labels,1,2)
        chosen,moves=refine(xy,red['selected'],cells)
        v=energy(xy,red['selected'],1,2);u=energy(xy,chosen,1,2)
        assert v==F(red['energy_exact']) and opt<=u<=v<=2*opt
        out['methods'][name]=dict(h_m=str(h),M=len(reps),representatives=reps,
            reduced=red,refined=chosen,
            refined_sites=describe(chosen),moves=moves,refined_energy=str(u),certificate_lower=str(v/2),
            base_error_percent=float(100*(v/opt-1)),refined_error_percent=float(100*(u/opt-1)))
        print(name,len(reps),out['methods'][name]['base_error_percent'],out['methods'][name]['refined_error_percent'],flush=True)
    outside=np.flatnonzero(centre_distance>30).tolist()
    no_central=minimize_cities(xy,outside,10,labels,1,2)
    out['best_outside_30km']=no_central
    out['central_exclusion_gap_percent']=float(100*(F(no_central['energy_exact'])/opt-1))
    (root/'dense-results.json').write_text(json.dumps(out,indent=2,ensure_ascii=False),encoding='utf-8')
    print('Central exclusion gap:',out['central_exclusion_gap_percent'],flush=True)
if __name__=='__main__':main()
