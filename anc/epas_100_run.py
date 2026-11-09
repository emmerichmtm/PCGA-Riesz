"""Spread 100 candidates across the admissible covered region, preserving centers.

Keeps the first ten input points, k=5, epsilon=1/2, and all eight farthest-first
choices. Exact branch-and-bound searches representatives; an independent run
on all original candidates certifies the optimum. Standard library only.
"""
from fractions import Fraction as F
from itertools import combinations,product
from math import comb
from pathlib import Path
import json
from planar_riesz_epas import solve,serializable,distance2,energy
from exact_subset_search import minimize
HERE=Path(__file__).resolve().parent

def make_input():
    base=json.loads((HERE/'epas_gallery_input.json').read_text())
    points=[tuple(F(v) for v in p) for p in base['points']]
    original=solve(points,5,F(1,2),'compact')
    ip=[tuple(int(x*512) for x in p) for p in points]
    centers=original['centers']
    radii=[None]+[int(t['insertion_radius_squared']*512**2) for t in original['traversal'][1:]]
    d2=lambda p,q:sum((x-y)**2 for x,y in zip(p,q))
    pool=[]
    for a,b in product(range(-28,45),range(-16,49)):
        q=(a*128,b*128)
        if q not in ip and all(min(d2(q,ip[c]) for c in centers[:j])<=radii[j] for j in range(1,8)):
            pool.append(q)
    nearest={q:min(d2(q,p) for p in ip) for q in pool}
    for _ in range(90):
        q=min(nearest,key=lambda q:(-nearest[q],q))
        del nearest[q]; ip.append(q)
        nearest={p:min(v,d2(p,q)) for p,v in nearest.items()}
    return [tuple(F(v,512) for v in p) for p in ip],original,len(pool)

def main():
    points,original,pool_size=make_input()
    run=solve(points,5,F(1,2),'compact','branch_bound')
    assert len(points)==len(set(points))==100 and points[:10]==original['points']
    for key in ['centers','traversal','last_insertion_radius_squared','harmonic',
                'potential_bound','dyadic_separation_scale','cell_side','grid_origin']:
        assert run[key]==original[key],key
    assert run['M']==98 and run['representatives']==list(range(8))+list(range(10,100))
    assert sorted(len(c['members']) for c in run['cells'])==[1]*96+[2]*2
    assert run['energy']==energy(points,run['selected'])
    assert all(distance2(points[i],points[j])>=13 for i,j in combinations(run['centers'],2))
    assert max(min(distance2(p,points[c]) for c in run['centers']) for p in points)<=13
    weights={(i,j):1/distance2(points[i],points[j]) for i,j in combinations(range(100),2)}
    star,opt,stats=minimize(weights,range(100),5,run['selected'])
    assert opt==run['energy'] and star==run['selected']
    x=star[0];others=tuple(i for i in star if i!=x)
    free=[c for c in run['centers'] if all(distance2(points[c],points[y])>=F(13,4) for y in others)]
    assert len(free)>=4
    rho2=max(distance2(points[i],points[run['replacement_map'][i]]) for i in range(100))
    result={'purpose':'Spread-out exact-coordinate illustration; no claimed large compression or measured speedup.',
            'indexing':'zero-based; p_i in the paper is one-based',
            'construction':{'pool_lattice_denominator':4,'a_range':[-28,44],'b_range':[-16,48],
                            'admissible_pool_size':pool_size,'added_points':90,
                            'rule':'maximin distance to existing candidates; lexicographic coordinate ties'},
            'compact_run':run,
            'audit':{'all_assertions_passed':True,'same_center_order_and_grid':True,
                     'original_subsets':comb(100,5),'representative_subsets':comb(98,5),
                     'representative_pair_weights':comb(98,2),'original_optimum_computed':True,
                     'original_optimum':opt,'original_optimal_set':star,
                     'exact_optimum_interval':[opt,opt],'certified_output_ratio':run['energy']/opt,
                     'search_stats':stats,'max_replacement_squared':rho2,
                     'exchange_example_x':x,'exchange_example_others':others,'exchange_example_free_centers':free}}
    for name,data in [('epas_100_input.json',{'points':points,'k':5,'epsilon':F(1,2)}),
                      ('epas_100_results.json',result)]:
        (HERE/name).write_text(json.dumps(serializable(data),indent=2)+'\n',encoding='utf-8')
    print(json.dumps(serializable({'M':run['M'],'centers':[i+1 for i in run['centers']],
          'selected':[i+1 for i in run['selected']],'value':str(float(run['energy'])),
          'exact_value':run['energy'],'ratio':run['energy']/opt,
          'free':[i+1 for i in free],'cover_squared':run['actual_cover_radius_squared'],
          'possible_representative_subsets':comb(98,5),'all_assertions_passed':True}),indent=2))

if __name__=='__main__':main()
