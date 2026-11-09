"""Independent exact checks of the paired certificate and grid code."""
import json,sys
from pathlib import Path
from fractions import Fraction as F
from itertools import combinations
import numpy as np
from paired import certificates,grid,make_input,configurations,METHODS,OUT
from epas_refine_v2 import D,energy,solve_reduced

def main():
    rng=np.random.default_rng(20261106);checks=0
    for trial in range(16):
        n=11;k=3 if trial%2==0 else 4
        X=rng.choice(1000000,size=(n,2),replace=False).astype(np.int64)-500000
        bounds,meta=certificates(X,k)
        vals=[(energy(X,S),S) for S in combinations(range(n),k)]
        opt,star=min(vals)
        for method in METHODS:
            B=bounds[method]
            assert all(F(D*D,int(np.dot(X[i]-X[j],X[i]-X[j])))<=B for i,j in combinations(star,2))
            g=grid(X,k,'1/2',B,F(1,777));h=F(g['h'])
            independent=[tuple((F(int(z),D)-F(1,777))//h for z in p) for p in X]
            assert [tuple(map(int,row)) for row in g['keys']]==independent
            result=solve_reduced(X,g,timeout=30)
            assert result['status']=='complete'
            assert F(result['energy_exact'])<=F(3,2)*opt
            checks+=1
        for U,S in vals:
            local=all(energy(X,tuple(y for y in S if y!=x)+(z,))>=U for x in S for z in range(n) if z not in S)
            if local:
                for x in S:
                    potential=sum((F(D*D,int(np.dot(X[x]-X[y],X[x]-X[y]))) for y in S if y!=x),F(0))
                    assert potential<=bounds['exchange']
                    checks+=1
    # Paired datasets, including rotated and unequal clusters, must be stable.
    for cfg in configurations()[::8]:
        assert np.array_equal(make_input(cfg),make_input(cfg));checks+=1
    result=dict(status='passed',exact_checks=checks,random_instances=16,
                tests=['exact optimum pair bounds','independent rational cell floors','original-objective approximation','every one-swap local optimum on tiny inputs','deterministic input generation'])
    OUT.mkdir(exist_ok=True,parents=True);(OUT/'verification-v2.json').write_text(json.dumps(result,indent=2))
    print(json.dumps(result))
if __name__=='__main__':main()
