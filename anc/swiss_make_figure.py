"""Write the standalone-data TikZ body of Figure 1 (all coordinates projected in LV95)."""
import csv,json,math
from pathlib import Path
import numpy as np

root=Path(__file__).resolve().parent
def load(f):return json.loads((root/f).read_text(encoding="utf8"))
with (root/"candidate-sites.csv").open(encoding="utf-8-sig") as f:rows=list(csv.DictReader(f))
xy=np.array([[int(r["east_lv95_m"]),int(r["north_lv95_m"])] for r in rows])
d=load("dense-results.json");h=d["methods"]["hybrid"];reps=h["representatives"];opt=d["full"]["selected"]
loc=[i for i,r in enumerate(rows) if r["locality"]=="Hasliberg Wasserwendi"]
here=sorted(set(loc)&set(reps));start,win=next((a,b) for a,b in h["moves"] if a in loc)
assert set(opt)==set(h["refined"]) and len(here)==8 and len(loc)==150
geometry=load("boundary-lv95.json")["geometry"]
polys=geometry["coordinates"] if geometry["type"]=="MultiPolygon" else [geometry["coordinates"]]
allb=np.concatenate([np.asarray(ring) for poly in polys for ring in poly])
def transform(lo,hi,box):
    p,q,w,t=box;scale=min(w/(hi[0]-lo[0]),t/(hi[1]-lo[1]))
    shift=np.array([p+w/2,q+t/2])-scale*(lo+hi)/2
    return lambda v: np.asarray(v)*scale+shift
nat=transform(allb.min(0),allb.max(0),(0,0,236,146))
lc=xy[loc];mid=(lc.min(0)+lc.max(0))/2;dw,dh=3440,2360
lo=mid-[dw/2,dh/2];hi=mid+[dw/2,dh/2]
zoom=transform(lo,hi,(249,20,145,145*dh/dw))
def coord(z):return "(%.3f,%.3f)"%tuple(z)
def vpt(z):return "{%.3f}{%.3f}"%tuple(z)
lines=[r"\resizebox{\linewidth}{!}{\begin{tikzpicture}[x=1pt,y=1pt]",
       r"\definecolor{dsLand}{HTML}{F3F4EF}\definecolor{dsBorder}{HTML}{91A39D}",
       r"\definecolor{dsGray}{HTML}{A5AFAD}\definecolor{dsTeal}{HTML}{007574}",
       r"\definecolor{dsRed}{HTML}{A13F20}\definecolor{dsTint}{HTML}{FAEFE3}",
       r"\path[use as bounding box] (-1,-1) rectangle (395,156);",
       r"\def\dsM#1#2{\pgfpathmoveto{\pgfpoint{#1pt}{#2pt}}}",
       r"\def\dsL#1#2{\pgfpathlineto{\pgfpoint{#1pt}{#2pt}}}",
       r"\def\dsDot#1#2{\pgfpathcircle{\pgfpoint{#1pt}{#2pt}}{\dsR pt}}"]
def batches(cmd,items,n=12):
    for i in range(0,len(items),n):lines.append("".join(cmd(t) for t in items[i:i+n]))
for poly in polys:
    lines.append(r"\begin{scope}\pgfsetfillcolor{dsLand}\pgfsetstrokecolor{dsBorder}\pgfsetlinewidth{.4pt}\pgfseteorule")
    for ring in poly:
        pp=nat(ring);lines.append(r"\dsM"+vpt(pp[0]));batches(lambda q:r"\dsL"+vpt(q),pp[1:]);lines.append(r"\pgfpathclose")
    lines.append(r"\pgfusepath{fill,stroke}\end{scope}")
def dots(ids,proj,color,radius):
    lines.append(r"\begin{scope}\pgfsetfillcolor{"+color+r"}\def\dsR{"+str(radius)+"}")
    for i in range(0,len(ids),72):
        batches(lambda q:r"\dsDot"+vpt(q),proj(xy[ids[i:i+72]]))
        lines.append(r"\pgfusepath{fill}")
    lines.append(r"\end{scope}")
def diamond(proj,i,r):
    z=proj(xy[i])
    path="--".join(coord(z+np.array(u)*r) for u in [(0,1),(1,0),(0,-1),(-1,0)])
    lines.append(r"\path[fill=dsRed,draw=white,line width=.55pt]"+path+"--cycle;")
dots(list(range(len(xy))),nat,"dsGray",.43);dots(reps,nat,"dsTeal",.48)
for i in opt:diamond(nat,i,3.4)
# Give the central site a leader to its detail panel identification.
z=nat(xy[win])
lines.append(r"\draw[dsRed,line width=.35pt]"+coord(z+[-1.7,-1.7])+"--"+coord(z+[-12,-14])+r" node[below left,inner sep=.6pt,font=\scriptsize]{(b)};")
lines+=[r"\node[anchor=west,inner sep=0pt,font=\small\bfseries] at (1,150){(a) Switzerland};",
        r"\node[anchor=west,inner sep=0pt,font=\small\bfseries] at (249,150){(b) Hasliberg Wasserwendi};"]
b0=coord(zoom(lo));b1=coord(zoom(hi))
lines+= [r"\begin{scope}\clip "+b0+" rectangle "+b1+";"]
side=int(h["h_m"]);o=int(d["geometry"]["origin_m"])
key=(xy[start]-o)//side;cl=key*side+o
assert np.all((xy[win]-o)//side==key)
lines.append(r"\fill[dsTint]"+coord(zoom(cl))+" rectangle "+coord(zoom(cl+side))+";")
for dim in (0,1):
    for k in range(math.ceil((lo[dim]-o)/side),math.floor((hi[dim]-o)/side)+1):
        aa=lo.copy();bb=hi.copy();aa[dim]=bb[dim]=o+k*side
        lines.append(r"\draw[dsBorder,line width=.34pt]"+coord(zoom(aa))+"--"+coord(zoom(bb))+";")
dots(loc,zoom,"dsGray",.95)
for i in here:lines.append(r"\draw[dsTeal,line width=.9pt]"+coord(zoom(xy[i]))+" circle[radius=2.15pt];")
lines.append(r"\draw[-{Stealth[length=3pt,width=2.6pt]},dsRed,line width=.65pt,shorten <=2.5pt,shorten >=3pt]"+coord(zoom(xy[start]))+"--"+coord(zoom(xy[win]))+";")
diamond(zoom,win,2.95)
lines.append(r"\end{scope}\draw[dsBorder,line width=.34pt]"+b0+" rectangle "+b1+";")
lines.append(r"\end{tikzpicture}}")
(root/"figure-1-inline.tex").write_text("\n".join(lines)+"\n",encoding="utf8")
