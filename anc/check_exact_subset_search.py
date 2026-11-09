"""Compare exact pruning with independent exhaustive rational minimization."""
from fractions import Fraction as F
from itertools import combinations
import random,json
from exact_subset_search import minimize
from planar_riesz_epas import solve,energy
from pathlib import Path
rng=random.Random(20261103)
cases=0;sets=0
for n in range(3,13):
    for k in range(2,min(n,5)+1):
        for mode in range(3):
            ids=list(range(n))
            w={(i,j):F(rng.randrange(1,1000),rng.randrange(1,1000)) for i,j in combinations(ids,2)}
            if mode==1:w={ij:v/F(2**70) for ij,v in w.items()}
            if mode==2:w={ij:F(1,3) for ij in w}
            exact=min((sum((w[ij] for ij in combinations(s,2)),F(0)),s) for s in combinations(ids,k))
            found,val,stats=minimize(w,ids,k)
            assert (val,found)==exact,(n,k,mode)
            cases+=1;sets+=stats['possible_subsets']
base=json.loads(Path(__file__).with_name('epas_gallery_input.json').read_text())
for mode in ['paper','compact']:
    a=solve(base['points'],5,F(1,2),mode)
    b=solve(base['points'],5,F(1,2),mode,'branch_bound')
    assert a['selected']==b['selected'] and a['energy']==b['energy']
print(json.dumps({'all_assertions_passed':True,'weighted_cases':cases,'independent_exhaustive_subsets':sets,'original_example_modes_checked':2},indent=2))
