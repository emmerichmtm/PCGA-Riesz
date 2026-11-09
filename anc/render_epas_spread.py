"""Exact-coordinate TikZ for the ten-point and spread-out hundred-point runs.

Only coordinates converted to decimal text for TikZ rendering; no point jitter.
All classifications (centers, free centers, optimum, grid) read from audits.
"""
from pathlib import Path
from fractions import Fraction as F
from itertools import combinations
import json,math
HERE=Path(__file__).resolve().parent
styles=r'''\definecolor{epasblue}{RGB}{0,114,178}
\definecolor{epasgreen}{RGB}{0,130,105}
\definecolor{epasorange}{RGB}{213,94,0}
\definecolor{epaspurple}{RGB}{145,70,160}
\tikzset{
spinput/.style={circle,draw=black!55,fill=white,line width=.35pt,inner sep=0pt,minimum size=2.8pt},
spcenter/.style={circle,draw=epasblue,line width=.9pt,inner sep=0pt,minimum size=7.4pt},
spfree/.style={regular polygon,regular polygon sides=3,draw=epaspurple,line width=.8pt,inner sep=0pt,minimum size=12pt},
spkeep/.style={rectangle,draw=epasgreen,line width=.6pt,inner sep=0pt,minimum size=4.4pt},
spchosen/.style={diamond,draw=epasorange,fill=epasorange,line width=.5pt,inner sep=0pt,minimum size=5.6pt},
spcover/.style={draw=epasblue!35,line width=.3pt,fill=epasblue,fill opacity=.018},
spexcluded/.style={draw=epasorange!65,densely dashed,line width=.45pt,fill=epasorange,fill opacity=.035}
}
'''
def xy(p):return '('+','.join(f'{float(v):.10f}' for v in p)+')'
def load(name):
    d=json.loads((HERE/name).read_text());r=d['compact_run'];a=d['audit']
    p=[tuple(F(v) for v in q) for q in r['points']]
    return r,a,p
def nodes(P,ids,style):return '\n'.join(rf'\node[{style}] at {xy(P[i])} {{}};' for i in ids)
def picture(letter,body,scale):
    return '\n'.join([rf'\begin{{tikzpicture}}[x={scale}cm,y={scale}cm]',
      r'\path[use as bounding box] (-7,-4.4) rectangle (11,12.0);',
      rf'\node[anchor=north west,font=\small\bfseries] at (-7,12) {{({letter})}};',body,r'\end{tikzpicture}'])
def world(run,a,P,kind):
    r=math.sqrt(float(F(run['last_insertion_radius_squared'])))
    star=a.get('original_optimal_set',a.get('original_optimal_sets',[None])[0])
    body=[]
    if kind=='cover':
        body += [rf'\draw[spcover] {xy(P[i])} circle[radius={r:.10f}];' for i in run['centers']]
    if kind=='free':
        body += [rf'\draw[spexcluded] {xy(P[i])} circle[radius={r/2:.10f}];' for i in a['exchange_example_others']]
    if kind=='output':
        body += [rf'\draw[epasorange!20,line width=.3pt] {xy(P[i])}--{xy(P[j])};' for i,j in combinations(run['selected'],2)]
    body += [nodes(P,range(len(P)),'spinput')]
    if kind in ('output','free') and len(P)==100:body += [nodes(P,run['representatives'],'spkeep')]
    body += [nodes(P,run['centers'],'spcenter')]
    body += [nodes(P,run['selected'] if kind=='output' else star,'spchosen')]
    if kind=='free':
        body += [nodes(P,a['exchange_example_free_centers'],'spfree')]
        x=a['exchange_example_x']
        labelpos=(P[x][0]-F(4,5),P[x][1]-F(1,5)) if len(P)==100 else (P[x][0]+F(6,5),P[x][1])
        anchor='east' if len(P)==100 else 'west'
        body += [rf'\node[anchor={anchor},font=\scriptsize] at {xy(labelpos)} {{$x$}};']
    if len(P)==10 and kind=='cover':
        for i,label in [(0,'A'),(4,'B')]:
            body += [rf'\draw[black!40,densely dotted] {xy(P[i])} circle[radius=.45];',
                     rf'\node[anchor=east,font=\scriptsize] at {xy((P[i][0]-F(3,5),P[i][1]))} {{${label}$}};']
    return '\n'.join(body)
def zoom(run,P,side=1.55):
    h=F(run['cell_side']);out=[]
    # Shift is in outer plot coordinates; local x/y units are physical cm.
    for shift,label,i,j in [((-2.6,4),'A',0,8),((6.0,4),'B',4,9)]:
        rel=tuple((P[j][a]-P[i][a])/h for a in range(2))
        out += [rf'\begin{{scope}}[shift={{{xy(shift)}}},x={side}cm,y={side}cm]',
            r'\fill[epasgreen!4] (-.5,-.5) rectangle (.5,.5);',
            r'\draw[black!40,line width=.4pt] (-.5,-.5) rectangle (.5,.5);',
            rf'\node[font=\scriptsize] at (0,.68) {{${label}$}};',
            rf'\draw[-{{Stealth[length=1.1mm]}},epasorange,line width=.55pt,shorten <=1pt,shorten >=1pt] {xy(rel)} -- (0,0);',
            r'\node[spinput] at (0,0) {};',rf'\node[spinput] at {xy(rel)} {{}};',
            r'\node[spkeep] at (0,0) {};',
            # In A the center is p1, in B it is p10.
            rf'\node[spcenter] at {"(0,0)" if i==0 else xy(rel)} {{}};',
            rf'\node[spchosen] at {xy(rel)} {{}};' if len(P)==10 else '',
            r'\draw[<->,>=Stealth,black!50,line width=.35pt] (-.5,-.72)--(.5,-.72) node[midway,below,font=\scriptsize] {$h$};',
            r'\end{scope}']
    return '\n'.join(out)
legend=r'''\par\smallskip
\begin{tikzpicture}[font=\scriptsize]
\node[spinput] at (0,0) {};\node[anchor=west] at (.14,0) {candidate};
\node[spcenter] at (1.85,0) {};\node[anchor=west] at (2.03,0) {cover center};
\node[spfree] at (4.2,0) {};\node[anchor=west] at (4.42,0) {free center};
\node[spkeep] at (6.5,0) {};\node[anchor=west] at (6.65,0) {representative};
\node[spchosen] at (9.15,0) {};\node[anchor=west] at (9.32,0) {selected point};
\end{tikzpicture}
'''
def panel(pic,width):return rf'\begin{{minipage}}[c]{{{width}\textwidth}}\centering'+'\n'+pic+'\n'+r'\end{minipage}'
def small_main():
    r,a,P=load('epas_gallery_results.json')
    body=r'\hfill'.join(panel(picture(letter,b,.23),'.32') for letter,b in [
        ('a',world(r,a,P,'cover')),('b',world(r,a,P,'free')),('c',zoom(r,P))])
    caption=r'''The ten-point run, $k=5$, $d=s=2$, $\varepsilon=1/2$. All global views plot all ten candidates at their actual coordinates.
(a) The eight farthest-first centers and radius-$\sqrt{13}$ cover; orange diamonds show the independently verified optimum $S^\star=\{p_3,p_4,p_7,p_9,p_{10}\}$, not the returned set. The near pairs at $A,B$ overlap only at this global scale.
(b) With $x=p_3$, the four dashed radius-$r/2$ disks around $S^\star\setminus\{x\}$ leave the purple-triangle centers $p_2,p_3,p_6,p_8$ free; $x$ itself is free. All eight centers remain visible as blue rings.
(c) True magnifications of the two occupied grid cells, $h=1/64$: $p_9\mapsto p_1$ and $p_{10}\mapsto p_5$. No candidate is displaced for legibility. The EPAS searches $56$ representative subsets and returns $\{p_1,p_3,p_4,p_5,p_7\}$, with $V\approx0.267086330$ and certificate $[2V/3,V]$. The optimum and free-center classification illustrate the proof; the solver does not use the original optimum.'''
    return '% BEGIN COMPACT EPAS GALLERY\n'+r'\begin{figure}[!htb]\centering'+'\n'+styles+body+legend+'\n\\caption{'+caption+'}\n\\label{fig:epas-run}\n\\end{figure}\n% END COMPACT EPAS GALLERY\n'
def large(supp=False):
    r,a,P=load('epas_100_results.json')
    scale=.34 if supp else .33
    first=panel(picture('a',world(r,a,P,'cover'),scale),'.49')+r'\hfill'+panel(picture('b',world(r,a,P,'free'),scale),'.49')
    if supp:
        first+=r'\par\medskip'+panel(picture('c',zoom(r,P,2.25),scale),'.49')+r'\hfill'+panel(picture('d',world(r,a,P,'output'),scale),'.49')
    caption=r'''The spread-out $n=100$, $k=5$ run, with the same eight centers, $r=\sqrt{13}$, and $h=1/64$.
(a) All 100 candidates at their computed global positions, the actual covering disks, and the returned five points (orange diamonds), independently certified optimal.
(b) All 98 representatives (green squares) over the original candidates, all eight centers (blue rings), and the seven free centers (purple triangles) for $x=p_{12}$. The four dashed radius-$r/2$ disks are centered at $p_{19},p_{30},p_{33},p_{87}$; only $p_{10}$ is excluded. Thus free centers need not be selected points.'''
    if supp:
        caption+=r''' (c) The only nonsingleton cells, $\{p_1,p_9\}$ and $\{p_5,p_{10}\}$, shown at a common magnification with their actual offsets and replacement arrows. No selected optimum point is discarded. (d) The returned optimum among the representatives, with all ten selected-pair interactions. Every global panel plots all 100 original candidates; only the two original near pairs are unresolved there. All panels use the stored coordinates without jitter; the two local windows use the explicitly different scale $h$. This instance illustrates the conservative grid and the exchange geometry, not a large compression or a running-time improvement.'''
    else:
        caption+=r''' Only the original two near pairs overlap at global scale; their exact cell geometry and replacements are shown in Figure~\ref{fig:epas-run}(c). No added point is jittered, and none of the five selected points is a cover center. This spread-out input retains 98 representatives, so it does not demonstrate large compression.'''
    label='fig:epas-hundred-supp' if supp else 'fig:epas-hundred'
    return r'\begin{figure}['+('p' if supp else '!htb')+r']\centering'+'\n'+styles+first+legend+'\n\\caption[Spread-out hundred-point EPAS run]{'+caption+'}\n\\label{'+label+'}\n\\end{figure}\n'
if __name__=='__main__':
    (HERE/'epas_100_supp_inline.tex').write_text(large(True).replace('[p]','[!htbp]'), encoding='utf-8')
    print('Generated the spread-out hundred-site TikZ figure.')
