from pathlib import Path
from fractions import Fraction as F
import math, json
# City-like variant: the global points are exact, with an explicitly magnified cell.
data=json.loads((Path(__file__).with_name('compressed_hundred_v2.json')).read_text());run=data['run'];audit=data['audit']
pts=[tuple(map(F,p)) for p in run['points']];chosen=audit['original_optimal_set'];centers=run['centers'];free=audit['free_centers']
r=math.sqrt(float(F(run['last_insertion_radius_squared'])));h=F(run['cell_side'])
tex=[r'\begin{figure}[t]\centering',r'\definecolor{cityblue}{RGB}{0,114,178}\definecolor{citygreen}{RGB}{0,130,105}\definecolor{cityorange}{RGB}{213,94,0}\definecolor{citypurple}{RGB}{145,70,160}',r'\tikzset{citypoint/.style={circle,fill=black!55,inner sep=0pt,minimum size=1.7pt},citycenter/.style={circle,draw=cityblue,line width=.7pt,inner sep=0pt,minimum size=6pt},cityfree/.style={regular polygon,regular polygon sides=3,draw=citypurple,line width=.55pt,inner sep=0pt,minimum size=9pt},citychosen/.style={diamond,draw=cityorange,fill=cityorange,inner sep=0pt,minimum size=5pt},cityrep/.style={rectangle,draw=citygreen,line width=.7pt,inner sep=0pt,minimum size=5pt}}',r'\begin{minipage}{.61\linewidth}\centering\begin{tikzpicture}[x=.40cm,y=.40cm]',r'\path[use as bounding box] (-6,-4) rectangle (11,12);',r'\node[anchor=north west,font=\small\bfseries] at (-6,12) {(a)};']
for i in centers:
 x,y=map(float,pts[i]);tex.append(f'\\draw[cityblue!22,line width=.25pt] ({x:.9f},{y:.9f}) circle[radius={r:.9f}];')
for x,y in pts:tex.append(f'\\node[citypoint] at ({float(x):.9f},{float(y):.9f}) {{}};')
for ids,style in [(centers,'citycenter'),(free,'cityfree'),(chosen,'citychosen')]:
 for i in ids:
  x,y=map(float,pts[i]);tex.append(f'\\node[{style}] at ({x:.9f},{y:.9f}) {{}};')
tex += [r'\end{tikzpicture}\end{minipage}\hfill',r'\begin{minipage}{.37\linewidth}\centering\begin{tikzpicture}[x=3.35cm,y=3.35cm]',r'\path[use as bounding box] (-.6,-.7) rectangle (.6,.8);',r'\node[anchor=north west,font=\small\bfseries] at (-.6,.8) {(b)};',r'\draw[black!35] (-.5,-.5) rectangle (.5,.5);']
anchor=pts[15];rep=15;target=17
for i in range(15,20):
 x,y=[float((pts[i][j]-anchor[j])/h) for j in range(2)];tex.append(f'\\node[citypoint,minimum size=3pt] at ({x},{y}) {{}};')
x,y=[float((pts[target][j]-anchor[j])/h) for j in range(2)]
tex += [r'\node[cityrep,minimum size=7pt] at (0,0) {};',f'\\draw[->,citypurple,line width=.8pt] (0,0) -- ({x},{y});',f'\\node[citychosen,minimum size=7pt] at ({x},{y}) {{}};',r'\end{tikzpicture}\end{minipage}',r'\par\smallskip\begin{tikzpicture}[font=\scriptsize]',r'\node[citypoint] at (0,0) {};\node[anchor=west] at (.1,0) {candidate};',r'\node[citycenter] at (1.55,0) {};\node[anchor=west] at (1.7,0) {cover center};',r'\node[cityfree] at (3.65,0) {};\node[anchor=west] at (3.85,0) {free center};',r'\node[cityrep] at (5.7,0) {};\node[anchor=west] at (5.85,0) {representative};',r'\node[citychosen] at (8.05,0) {};\node[anchor=west] at (8.2,0) {selected site};',r'\end{tikzpicture}',r'\caption{A synthetic facility-dispersion map with $n=100$, $k=5$, and $\varepsilon=1/2$. (a) All candidates at their actual coordinates: the original ten sites plus five nearby sites in each of 18 spread-out cities. The ten farthest-first centers and their radius-$r$ cover are blue; purple marks the nine centers free for replacing $p_{10}$ in the independently certified optimum (orange). Small city clouds overlap at this scale. (b) One entire occupied cell, magnified by $1/h$, with all five candidates; its representative $p_{16}$ moves to $p_{18}$ during refinement. The exact grid has $h=1/64$ and $M=26$. Reduced search has relative error $0.05966\%$; refinement reaches the independently computed optimum. These are observed errors, separate from the certified factor $1.5$.}\label{fig:hundred-clusters}',r'\end{figure}']
Path(__file__).with_name('hundred-clusters.tex').write_text('\n'.join(tex))
print('City',dict(r=r,h=str(h),V=run['energy'],opt=audit['original_optimum'],subsets=math.comb(26,5)))
