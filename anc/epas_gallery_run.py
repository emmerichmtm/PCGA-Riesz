"""Regenerate the 5-of-10 trace and independently audit both EPAS modes.

Run alongside planar_riesz_epas.py. Writes epas_gallery_results.json.
Exhaustive original optimization is an audit only, never used by the EPAS.
"""
from fractions import Fraction as F
from itertools import combinations
import json
from pathlib import Path
import random
from planar_riesz_epas import solve, serializable

HERE=Path(__file__).resolve().parent


def squared(p,q):
    return sum((a-b)*(a-b) for a,b in zip(p,q))


def full_table(points,k):
    weights={(i,j):F(1)/squared(points[i],points[j]) for i,j in combinations(range(len(points)),2)}
    return {s:sum((weights[i,j] for i,j in combinations(s,2)),F(0))
            for s in combinations(range(len(points)),k)}


def check_case(raw,k,epsilon,mode):
    result=solve(raw,k,epsilon,mode)
    points=result['points']
    table=full_table(points,k)
    opt=min(table.values())
    assert len(result['selected'])==k and table[tuple(result['selected'])]==result['energy']
    assert opt<=result['energy']<=(1+result['epsilon_used'])*opt
    assert result['certificate'][0]<=opt<=result['certificate'][1]
    if result['branch']=='exchange_grid':
        centers=result['centers']; r2=result['last_insertion_radius_squared']
        assert min(squared(points[i],points[j]) for i,j in combinations(centers,2))==r2
        assert max(min(squared(p,points[c]) for c in centers) for p in points)<=r2
        mapping=result['replacement_map']; h=result['cell_side']; B=result['potential_bound']
        for s,val in table.items():
            if val!=opt:
                continue
            mapped=tuple(sorted(mapping[i] for i in s))
            assert len(set(mapped))==k and table[mapped]<=(1+result['epsilon_used'])*opt
            for x in s:
                others=[y for y in s if y!=x]
                free=[c for c in centers if all(squared(points[c],points[y])>=r2/4 for y in others)]
                assert len(free)>=len(centers)-(k-1)>=k-1
                potential=sum((F(1)/squared(points[x],points[y]) for y in others),F(0))
                alternatives=[sum((F(1)/squared(points[c],points[y]) for y in others),F(0)) for c in free]
                assert potential<=min(alternatives)<=sum(alternatives)/len(free)<=B
                assert squared(points[x],points[mapping[x]])<=2*h*h
    return result,table


def main():
    inp=json.loads((HERE/'epas_gallery_input.json').read_text(),parse_float=str)
    compact,table=check_case(inp['points'],inp['k'],inp['epsilon'],'compact')
    paper,_=check_case(inp['points'],inp['k'],inp['epsilon'],'paper')
    opt=min(table.values()); originals=[s for s,v in table.items() if v==opt]
    assert len(originals)==1
    star=originals[0]
    snap=tuple(sorted(compact['replacement_map'][i] for i in star))
    assert compact['M']==8 and compact['subsets_evaluated']==56 and len(table)==252
    assert compact['energy']>opt and table[snap]==compact['energy']
    points=compact['points']; centers=compact['centers']; r2=compact['last_insertion_radius_squared']
    x=star[0]; others=[i for i in star if i!=x]
    free=[c for c in centers if all(squared(points[c],points[y])>=r2/4 for y in others)]
    local_minima=[]
    for s,value in table.items():
        if all(table[tuple(sorted((set(s)-{a})|{b}))]>=value
               for a in s for b in range(len(points)) if b not in s):
            local_minima.append(s)
            for a in s:
                rest=[y for y in s if y!=a]
                potential=sum((F(1)/squared(points[a],points[y]) for y in rest),F(0))
                assert potential<=compact['potential_bound']
    checks=[]
    rng=random.Random(20261101)
    for t in range(30):
        n=3+t%10
        cloud=set()
        while len(cloud)<n:
            cloud.add((F(rng.randrange(-20,21),4),F(rng.randrange(-20,21),4)))
        pts=sorted(cloud)
        k=2+t%min(5,n-1)
        eps=[F(1,8),F(1,2),F(1),F(2)][t%4]
        for mode in ('paper','compact'):
            run,tt=check_case(pts,k,eps,mode)
            checks.append({'n':n,'k':k,'mode':mode,'branch':run['branch'],'subsets':len(tt)})
    for pts,k in [([],0),([(0,0)],0),([(0,0)],1),([(0,0),(1,0)],2),
                  (inp['points'],0),(inp['points'],1),(inp['points'],10)]:
        for mode in ('paper','compact'):
            run,tt=check_case(pts,k,F(1,2),mode)
            checks.append({'n':len(pts),'k':k,'mode':mode,'branch':run['branch'],'subsets':len(tt)})
    for scale,translation in [(F(1,2**80),(F(2**100),F(-2**100))),
                              (F(2**60),(F(-13,7),F(19,11)))]:
        transformed=[tuple(scale*x+translation[a] for a,x in enumerate(p)) for p in points]
        for mode in ('paper','compact'):
            run,tt=check_case(transformed,5,F(1,2),mode)
            checks.append({'n':10,'k':5,'mode':mode,'branch':run['branch'],'subsets':len(tt),'scaled':True})
    invalid=[([(0,0),(0,0)],1,F(1,2)), ([(0,0)],2,F(1,2)),
             ([(0,0)],1,F(0)), ([(0.1,0)],1,F(1,2))]
    for pts,k,eps in invalid:
        try: solve(pts,k,eps)
        except (ValueError,TypeError): pass
        else: raise AssertionError('Invalid input accepted')
    output={'purpose':'Exact reproducible teaching example and finite verification, not a speed benchmark.',
            'compact_run':compact,'paper_run':paper,
            'audit':{'original_subsets':len(table),'original_optimal_sets':originals,'original_optimum':opt,
                     'snapped_optimum':snap,'snapped_energy':table[snap],
                     'observed_ratio':compact['energy']/opt,
                     'exchange_example_x':x,'exchange_example_others':others,'exchange_example_free_centers':free,
                     'all_one_swap_local_minima':local_minima},
            'tests':{'all_assertions_passed':True,'additional_cases':len(checks),
                     'additional_original_subsets':sum(c['subsets'] for c in checks),
                     'invalid_inputs_rejected':len(invalid),'cases':checks}}
    (HERE/'epas_gallery_results.json').write_text(json.dumps(serializable(output),indent=2)+'\n',encoding='utf-8')
    print(json.dumps(serializable({'centers':centers,'r2':r2,'h':compact['cell_side'],
           'bound':compact['potential_bound'],'scale':compact['dyadic_separation_scale'],
           'opt':opt,'Sstar':star,'returned':compact['selected'],
           'V':compact['energy'],'ratio_display':str(float(compact['energy']/opt)),
           'free':free,'test_cases':len(checks),'all_assertions_passed':True}),indent=2))


if __name__=='__main__': main()
