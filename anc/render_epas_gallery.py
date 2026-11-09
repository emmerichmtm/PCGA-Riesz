"""Generate the two inline TikZ figures for the ten-site walkthrough."""
from pathlib import Path
from fractions import Fraction as F
from itertools import combinations
import json, math

PKG=Path(__file__).resolve().parent
data=json.loads((PKG/'epas_gallery_results.json').read_text())
run=data['compact_run']; audit=data['audit']
P=[tuple(F(x) for x in p) for p in run['points']]
C=run['centers']; selected=run['selected']; star=audit['original_optimal_sets'][0]
assert run['k']==5 and run['n']==10 and run['epsilon_used']=='1/2'
assert P==[(F(0),F(0)),(F(4),F(0)),(F(7),F(2)),(F(7),F(6)),(F(3),F(8)),(F(-1),F(6)),(F(-3),F(3)),(F(2),F(3)),(F(-1,512),F(-1,512)),(F(3),F(4097,512))], 'This renderer is specific to the documented gallery input.'
r=math.sqrt(float(F(run['last_insertion_radius_squared'])))
h=F(run['cell_side']); V=F(run['energy']); opt=F(audit['original_optimum'])

def coord(p):return '('+','.join(f'{float(x):.8f}' for x in p)+')'
def dot(i,style):return rf'\node[{style}] at {coord(P[i])} {{}};'
def points(style='ginput',ids=range(10)):return '\n'.join(dot(i,style) for i in ids)
def panel(letter,body,large=False):
    # Equal physical panel sizes and equal scales for all global views.
    return '\n'.join([r'\begin{minipage}[t]{.49\textwidth}\centering',
        r'\begin{tikzpicture}[x=.35cm,y=.35cm]',
        r'\path[use as bounding box] (-7,-4.7) rectangle (11,12.3);',
        rf'\node[anchor=north west,font=\small\bfseries] at (-7,12.3) {{({letter})}};',
        body,r'\end{tikzpicture}',r'\end{minipage}'])

legend=r'''\par\smallskip
\begin{tikzpicture}[x=1cm,y=1cm,font=\footnotesize]
\node[ginput] at (0,0) {}; \node[anchor=west] at (.16,0) {candidate};
\node[gcenter] at (2.3,0) {}; \node[anchor=west] at (2.46,0) {cover center};
\node[gkeep] at (5.3,0) {}; \node[anchor=west] at (5.46,0) {representative};
\node[gchosen] at (8.9,0) {}; \node[anchor=west] at (9.06,0) {selected point};
\node[gfree] at (0,-.55) {}; \node[anchor=west] at (.20,-.55) {free center (exchange panel)};
\end{tikzpicture}
'''
styles=r'''\definecolor{epasblue}{RGB}{0,114,178}
\definecolor{epasgreen}{RGB}{0,130,105}
\definecolor{epasorange}{RGB}{213,94,0}
\definecolor{epaspurple}{RGB}{145,70,160}
\tikzset{
ginput/.style={circle,draw=black!60,fill=white,line width=.6pt,inner sep=0pt,minimum size=4.6pt},
gcenter/.style={circle,draw=epasblue,line width=.8pt,inner sep=0pt,minimum size=9pt},
gfree/.style={regular polygon,regular polygon sides=3,draw=epaspurple,line width=.8pt,inner sep=0pt,minimum size=14pt},
gkeep/.style={rectangle,draw=epasgreen,fill=white,line width=.8pt,inner sep=0pt,minimum size=6pt},
gchosen/.style={diamond,draw=epasorange,fill=epasorange,line width=.5pt,inner sep=0pt,minimum size=7pt},
gmove/.style={-{Stealth[length=1.6mm]},epasorange,line width=.8pt},
gcover/.style={draw=epasblue!45,line width=.35pt,fill=epasblue,fill opacity=.035},
gexcluded/.style={draw=epasorange!70,densely dashed,line width=.5pt,fill=epasorange,fill opacity=.07}
}
'''

def zoom_tags():
    return r'''\draw[black!35,densely dashed] (0,0) circle[radius=.52];
\node[anchor=east,font=\footnotesize] at (-.65,0) {$A$};
\draw[black!35,densely dashed] (3,8) circle[radius=.52];
\node[anchor=south,font=\footnotesize] at (3,8.7) {$B$};'''

def numbered(ids):
    offsets={0:(.35,-.6),3:(.4,.5),5:(-.5,.5),9:(.1,.65),1:(.4,-.6),2:(.4,-.6),6:(-.6,-.5),7:(.5,.4)}
    return '\n'.join(dot(i,'gcenter')+rf'\node[text=epasblue,font=\scriptsize] at {coord((P[i][0]+F(str(offsets[i][0])),P[i][1]+F(str(offsets[i][1]))))} {{{j+1}}};' for j,i in enumerate(ids))

firstbody=points()+ '\n'+zoom_tags()
prefix=points()+'\n'+numbered(C[:4])
cover='\n'.join(rf'\draw[gcover] {coord(P[i])} circle[radius={r:.8f}];' for i in C)
cover+='\n'+points()+'\n'+numbered(C)
cover+='\n'+points('gchosen',star)
cover+=rf'\draw[epasblue,line width=.55pt] {coord(P[2])} -- ++(0,-{r:.8f}) node[midway,right,font=\footnotesize] {{$r$}};'
x=audit['exchange_example_x']; others=audit['exchange_example_others']; free=audit['exchange_example_free_centers']
exchange='\n'.join(rf'\draw[gexcluded] {coord(P[i])} circle[radius={r/2:.8f}];' for i in others)
exchange+='\n'+points()+'\n'+points('gcenter',C)
exchange+='\n'+'\n'.join(dot(i,'gchosen') for i in others)
exchange+=rf'\node[diamond,draw=epasorange,fill=white,minimum size=7pt,inner sep=0pt,line width=.9pt] at {coord(P[x])} {{}};'
exchange+=rf'\node[anchor=west,font=\footnotesize] at {coord((P[x][0]+F(1,2),P[x][1]))} {{$x$}};'
exchange+='\n'+points('gfree',free)

fig1=r'''\begin{figure}[!htbp]\centering
'''+styles+panel('a',firstbody)+r'\hfill'+panel('b',prefix)+r'\par\medskip'+panel('c',cover)+r'\hfill'+panel('d',exchange)+legend+r'''
\caption[Five-of-ten EPAS: cover and exchange scale]{Cover and separation for the actual ten-point run, with $k=5$ and $\varepsilon=1/2$.
(a) The ten input candidates; the two pairs at $A$ and $B$ are unresolved at this global scale and are magnified in Figure~\ref{fig:epas-grid-run}. No positions are displaced for visibility.
(b) The first four farthest-first centers; the numerals give traversal order, not input indices.
(c) All ten candidates, all eight centers in the certified $2k-2$ variant (blue rings), and the independently verified optimum $S^\star=\{p_3,p_4,p_7,p_9,p_{10}\}$ (orange diamonds), with true radius-$r$ covering disks, $r=\sqrt{13}$. The center order is $p_1,p_4,p_6,p_{10},p_2,p_3,p_7,p_8$; a cover is distinct from a five-point selection.
(d) The proof mechanism for $x=p_3$. The four orange open-radius-$r/2$ neighborhoods around $S^\star\setminus\{x\}$ exclude four centers, leaving $p_2,p_3,p_6,p_8$ free (purple triangles). All eight centers remain blue rings. The outlined diamond denotes $x$; it is itself free. The algorithm does not know $S^\star$ or perform these exchanges: averaging the feasible replacement potentials proves the separation scale used by the grid. Distinct shapes also identify each role in greyscale.}
\label{fig:epas-cover-run}
\end{figure}
'''

def zoom_panel(arrows):
    out=[]
    for shift,label,base,mate in [(-2.5,'A',0,8),(6.5,'B',4,9)]:
        # One local coordinate unit is exactly one grid side.
        out.append(rf'\begin{{scope}}[shift={{({shift},4.0)}},x=2.7cm,y=2.7cm]')
        out.append(r'\fill[epasgreen!5] (-.5,-.5) rectangle (.5,.5);')
        # Actual cell boundaries (the displayed grid is centered on the base).
        out.append(r'\draw[black!35,line width=.45pt] (-.5,-.58)--(-.5,.58) (.5,-.58)--(.5,.58) (-.58,-.5)--(.58,-.5) (-.58,.5)--(.58,.5);')
        out.append(rf'\node[font=\footnotesize] at (0,.74) {{${label}$}};')
        rel=tuple((P[mate][a]-P[base][a])/h for a in range(2))
        if arrows:
            out.append(rf'\draw[gmove,shorten >=1.8pt,shorten <=1.8pt] {coord(rel)} -- (0,0);')
            out.append(rf'\node[gchosen] at {coord(rel)} {{}};')
            out.append(r'\node[gkeep] at (0,0) {};')
        else:
            out.append(rf'\node[ginput] at {coord(rel)} {{}};')
            out.append(r'\node[ginput] at (0,0) {};')
        out.append(r'\draw[<->,black!60,>=Stealth,line width=.4pt] (-.5,-.72)--(.5,-.72) node[midway,below,font=\footnotesize] {$h$};')
        out.append(r'\end{scope}')
    return '\n'.join(out)

reps=points()+'\n'+points('gkeep',run['representatives'])+'\n'+zoom_tags()
final=points('gkeep',run['representatives'])
final+='\n'+'\n'.join(rf'\draw[epasorange!22,line width=.35pt] {coord(P[i])} -- {coord(P[j])};' for i,j in combinations(selected,2))
final+='\n'+points('gkeep',run['representatives'])+'\n'+points('gchosen',selected)
fig2=r'''\begin{figure}[!htbp]\centering
'''+styles+panel('a',zoom_panel(False))+r'\hfill'+panel('b',zoom_panel(True))+r'\par\medskip'+panel('c',reps)+r'\hfill'+panel('d',final)+legend+rf'''
\caption[Five-of-ten EPAS: grid, replacement, and output]{{Grid compression and the returned five-point selection on the same input.
(a) Two magnified windows of the \emph{{actual}} grid, of side $h=1/64$ and origin $(-1/128,-1/128)$. A central square is one half-open cell; $A$ contains $p_1,p_9$ and $B$ contains $p_5,p_{{10}}$. Magnification is uniform within each window; the remaining six occupied cells are singletons. The grid is not drawn at global scale, where it would be illegible.
(b) The first input candidate in each cell is retained (green square). Arrows show $p_9\mapsto p_1$ and $p_{{10}}\mapsto p_5$ for the independently audited optimum (orange diamonds). The other three optimum points map to themselves. This illustrates the injective replacement used in the proof; the algorithm simply forms $Q$ and searches it.
(c) The eight representatives; a discarded candidate may coincide visually with its representative at this scale.
(d) Exact enumeration of the $\binom85=56$ subsets returns $\{{p_1,p_3,p_4,p_5,p_7\}}$, with full energy $V\approx {float(V):.9f}$. Pale segments show all ten selected-pair interactions. The independent $252$-subset check gives $\OPT\approx {float(opt):.9f}$ and $V/\OPT\approx {float(V/opt):.9f}$; the guaranteed interval is $[2V/3,V]$. This small illustrative run is not a performance benchmark.}}
\label{{fig:epas-grid-run}}
\end{{figure}}
'''

fragment=fig1+fig2
output=Path(__file__).with_name('epas_gallery_inline.tex')
output.write_text(fragment,encoding='utf-8')
print(output)
