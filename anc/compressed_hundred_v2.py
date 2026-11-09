"""A compressed 100-point variant of the spread example, with exact audit."""
from pathlib import Path
from fractions import Fraction as F
from itertools import combinations
import json
import numpy as np
from planar_riesz_epas_v2 import solve,serializable,distance2,energy
from exact_subset_search_v2 import minimize
from epas_refine_v2 import D,prepare,refine
ROOT=Path(__file__).resolve().parent

def main():
    old=json.loads((ROOT/'epas_100_results.json').read_text())['compact_run']
    points=[tuple(F(z) for z in p) for p in old['points']]
    base=points[:10];centers=old['centers'];trace=old['traversal']
    # Only old candidates with a positive margin to every traversal boundary
    # become cloud anchors. Five dyadic offsets remain in the covered region.
    offsets=[(F(0),F(0)),(F(1,512),F(0)),(-F(1,512),F(0)),(F(0),F(1,512)),(F(0),-F(1,512))]
    eligible=[]
    for p in points[10:]:
        cloud=[(p[0]+a,p[1]+b) for a,b in offsets]
        if all(min(distance2(q,base[c]) for c in centers[:j])<F(trace[j]['insertion_radius_squared']) for q in cloud for j in range(1,8)):
            eligible.append(p)
    anchors=[];selected=base[:]
    while len(anchors)<18:
        p=min(eligible,key=lambda p:(-min(distance2(p,q) for q in selected),p))
        anchors.append(p);selected.append(p);eligible.remove(p)
    new=base+[tuple(p[i]+off[i] for i in range(2)) for p in anchors for off in offsets]
    assert len(new)==len(set(new))==100
    run=solve(new,5,F(1,2),'paper','branch_bound')
    assert run['centers'][:8]==centers
    weights={(i,j):1/distance2(new[i],new[j]) for i,j in combinations(range(100),2)}
    optset,opt,stats=minimize(weights,range(100),5,run['selected'],timeout=120)
    X=np.array([[int(z*D) for z in p] for p in new],dtype=np.int64)
    g=prepare(X,5,F(1,2));assert g['representatives']==run['representatives']
    tr=refine(X,g,run['selected'],sweeps=10)
    x=optset[0];others=[i for i in optset if i!=x];r2=run['last_insertion_radius_squared']
    free=[c for c in run['centers'] if all(distance2(new[c],new[y])>=r2/4 for y in others)]
    assert len(free)>=6 and run['energy']<=F(3,2)*opt
    audit=dict(original_optimum=str(opt),original_optimal_set=optset,ratio=str(run['energy']/opt),
               search_stats=stats,free_centers=free,exchange_x=x,exchange_others=others,
               preserved_first_eight_centers=True,all_assertions_passed=True)
    out=dict(construction='Original ten sites plus five-point crosses at 18 spread admissible anchors; offsets 1/512.',
             anchors=anchors,offsets=offsets,run=run,audit=audit,refinement=tr)
    (ROOT/'compressed_hundred_v2.json').write_text(json.dumps(serializable(out),indent=2))
    print(json.dumps(dict(M=run['M'],centers=run['centers'],selected=run['selected'],opt=optset,
                         energy=float(run['energy']),ratio=float(run['energy']/opt),h=str(run['cell_side']),free=free,
                         refined=tr['selected'],refined_ratio=float(F(tr['trace'][-1]['energy_exact'])/opt))),flush=True)
if __name__=='__main__':main()
