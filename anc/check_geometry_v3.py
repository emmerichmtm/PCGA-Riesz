"""Small exact checks of polygon witnesses, continuous traversal, and metric rounding."""
from fractions import Fraction as F
from itertools import combinations, product
from pathlib import Path
import json
from planar_riesz_epas_v2 import dyadic_inverse_sqrt

def squared(p, q):
    return sum((a-b)**2 for a,b in zip(p,q))

def clip(vertices, a, b, c):
    """Intersect a closed convex polygon (also a segment/point) with ax+by <= c."""
    result = []
    for p,q in zip(vertices,vertices[1:]+vertices[:1]):
        u,v = a*p[0]+b*p[1]-c, a*q[0]+b*q[1]-c
        if u <= 0:
            result.append(p)
        if (u < 0 < v) or (v < 0 < u):
            t = F(u,u-v)
            result.append(tuple(x+t*(y-x) for x,y in zip(p,q)))
    return list(dict.fromkeys(result))

def witness(triangle, low, high):
    poly = list(triangle)
    for boundary in ((-1,0,-low[0]),(1,0,high[0]),(0,-1,-low[1]),(0,1,high[1])):
        poly = clip(poly,*boundary)
    if not poly:
        return None
    mean = tuple(F(sum(z),len(poly)) for z in zip(*poly))
    return mean if all(low[t] <= mean[t] < high[t] for t in (0,1)) else None

def square_traversal():
    square = [(0,0),(1,0),(1,1),(0,1)]
    centers = [square[0]]
    radii = []
    for _ in range(4):
        candidates = set()
        for p in centers:
            cell = square[:]
            for q in centers:
                if p != q:
                    cell = clip(cell,2*(q[0]-p[0]),2*(q[1]-p[1]),squared(q,(0,0))-squared(p,(0,0)))
            candidates.update(cell)
        radius, point = max((min(squared(v,c) for c in centers),v) for v in candidates)
        radii.append(radius)
        centers.append(point)
    assert set(centers[:4]) == set(square) and centers[4] == (F(1,2),F(1,2))
    assert radii == [2,1,1,F(1,2)] and F(1,max(squared(p,q) for p,q in combinations(square,2))) == F(1,2)
    return [str(v) for v in radii]

def metric_trial(points, k, epsilon):
    # The metric is handed to the algorithm only as a distance matrix. Finite
    # subsets of the line admit a doubling bound lambda=2 with input-centered balls.
    distances = [[abs(p-q) for q in points] for p in points]
    ids = list(range(len(points)))
    for i in ids:
        for radius in set(distances[i]):
            ball = [j for j in ids if distances[i][j] <= radius]
            cover = []
            for j in sorted(ball,key=points.__getitem__):
                if all(2*distances[j][c] > radius for c in cover):
                    cover.append(j)
            assert len(cover)<=4
    centers = [0]
    for _ in range(2*k-1):
        radius,i = max((min(distances[j][c] for c in centers),j) for j in ids)
        centers.append(i)
    assert min(distances[i][j] for i,j in combinations(centers,2))>=radius
    assert max(min(distances[i][j] for j in centers) for i in ids)<=radius
    harmonic=sum(F(1,j) for j in range(1,2*k+1))
    bound=16*harmonic/radius**2
    eta=epsilon*dyadic_inverse_sqrt(bound)/8
    net=[]
    for j in ids:
        if all(distances[j][c]>eta for c in net):
            net.append(j)
    assert all(distances[i][j]>eta for i,j in combinations(net,2))
    assert len(net) <= 2*k*(4*radius/eta)**2
    assignment={i:min(net,key=lambda j:distances[i][j]) for i in ids}
    assert all(distances[i][j]<=eta for i,j in assignment.items())
    def energy(selection):
        return sum(F(1,distances[i][j]**2) for i,j in combinations(selection,2))
    value,optimal=min((energy(c),c) for c in combinations(ids,k))
    assert all(F(1,distances[i][j]**2)<=bound for i,j in combinations(optimal,2))
    snap=[assignment[i] for i in optimal]
    assert len(set(snap))==k
    reduced=min(energy(c) for c in combinations(net,k))
    assert value<=reduced<=energy(snap)<=(1+epsilon)*value
    return dict(n=len(points),k=k,epsilon=str(epsilon),M=len(net),opt=str(value),reduced=str(reduced))

if __name__ == "__main__":
    triangle=[(0,0),(1,0),(0,1)]
    tests=[((-1,-1),(0,0),False),((-1,0),(0,1),False),((1,0),(2,1),True),
           ((0,1),(1,2),True),((0,0),(1,1),True),((F(1,2),F(1,2)),(1,1),True),
           ((1,1),(2,2),False)]
    for low,high,expected in tests:
        p=witness(triangle,low,high)
        assert (p is not None)==expected
        if p is not None:assert min(p)>=0 and sum(p)<=1
    samples=([0,1,3,1000,1001,2000,2010,4000,4010,6000,9000,12000],
             [13*i for i in range(14)])
    trials=[metric_trial(points,k,e) for points,k,e in product(samples,(2,3,4),(F(1,2),F(1)))]
    result=dict(square_radii_squared=square_traversal(),square_optimum="1/2",cell_tests=len(tests),metric_trials=trials)
    assert any(r["M"]<r["n"] for r in trials)
    Path(__file__).with_name("geometry-checks-v3.json").write_text(json.dumps(result,indent=2))
    print("Verified square traversal, seven half-open witnesses and",len(trials),"exact metric trials.")
