"""Add nearby synthetic sites to the frozen Swiss input. No optimizer or EPAS grid is consulted."""
import argparse, collections, csv, hashlib, json, sqlite3
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from shapely import force_2d, from_wkb
from shapely.geometry import Point
from shapely.prepared import prep

ROOT=Path(__file__).resolve().parent
SEED=20261107
def digest(p): return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('address_zip',type=Path); ap.add_argument('boundary_gpkg',type=Path)
    args=ap.parse_args()
    base=ROOT/'baseline-1801.csv'
    assert digest(base)=='141570a8559da8a513a5b2f7404755e33058b31069d60041652d014b47d332aa'
    with base.open(encoding='utf-8-sig') as f: rows=list(csv.DictReader(f))
    with sqlite3.connect(args.boundary_gpkg) as db:
        blob=db.execute("select geom from tlm_landesgebiet where icc='CH'").fetchone()[0]
    border=prep(force_2d(from_wkb(blob[8+{0:0,1:32,2:48,3:48,4:64}[(blob[3]>>1)&7]:])))
    cols=['ADR_EASTING','ADR_NORTHING']
    frames=[a[cols].dropna().to_numpy(dtype=float) for a in pd.read_csv(args.address_zip,sep=';',usecols=cols,chunksize=200000)]
    registered=np.concatenate(frames);del frames
    tree=cKDTree(registered)
    xy=lambda row:(int(row['east_lv95_m']),int(row['north_lv95_m']))
    occupied=set(xy(r) for r in rows)
    assert len(occupied)==1801
    rng=np.random.Generator(np.random.PCG64(SEED))
    ids=[i for i,row in enumerate(rows) if row['group_type']!='eastern']
    assert len(ids)==1800
    # Two extra draws elsewhere keep the isolated border site a singleton.
    schedule=list(ids)*2+list(rng.choice(ids,2,replace=False))
    buckets=collections.defaultdict(list)
    def bucket(p):return p[0]//15,p[1]//15
    for p in occupied:buckets[bucket(p)].append(p)
    radii=[]
    for i in schedule:
        source=rows[i]; origin=np.array(xy(source))
        for _ in range(100):
            theta=rng.uniform(0,2*np.pi,32)
            r=np.sqrt(rng.uniform(75**2,250**2,32))
            offsets=np.rint(np.column_stack((np.cos(theta)*r,np.sin(theta)*r))).astype(np.int64)
            points=origin+offsets
            r2=np.sum(offsets*offsets,axis=1)
            ok=np.flatnonzero((r2>=75**2)&(r2<=250**2)&(tree.query(points)[0]>=15))
            for j in ok:
                p=tuple(map(int,points[j])); a,b=bucket(p)
                if any((x-p[0])**2+(y-p[1])**2<225 for c in range(a-1,a+2) for d in range(b-1,b+2) for x,y in buckets.get((c,d),())):continue
                if border.covers(Point(*p)):break
            else:continue
            break
        else:raise RuntimeError('Failed to sample: '+source['site_id'])
        rows.append(dict(source,site_id=f'P{len(rows)+1:04}',east_lv95_m=p[0],north_lv95_m=p[1]))
        occupied.add(p);buckets[a,b].append(p);radii.append(float(np.sqrt(r2[j])))
    points=np.array([xy(r) for r in rows])
    address_min=float(tree.query(points)[0].min())
    add_min=float(cKDTree(points).query(points[1801:],k=2)[0][:,1].min())
    assert len(occupied)==len(rows)==5403 and address_min>=15 and add_min>=15 and all(border.covers(Point(*p)) for p in points)
    assert sum(r['group_type']=='eastern' for r in rows)==1
    dst=ROOT/'candidate-sites.csv'
    with dst.open('w',newline='',encoding='utf-8-sig') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    info=dict(seed=SEED,n=len(rows),original=1801,added=len(radii),extra_parents=[rows[i]['site_id'] for i in schedule[-2:]],
              min_offset_m=min(radii),max_offset_m=max(radii),min_distance_from_directory_m=address_min,
              min_new_to_other_site_m=add_min,directory_coordinates_checked=len(registered),all_inside_country=True,
              input_sha256=digest(base),output_sha256=digest(dst),address_zip_sha256=digest(args.address_zip),
              boundary_gpkg_sha256=digest(args.boundary_gpkg))
    (ROOT/'generation.json').write_text(json.dumps(info,indent=2))
    print(json.dumps(info,indent=2))
if __name__=='__main__':main()
