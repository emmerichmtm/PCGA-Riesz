# Reproduce *Parameterized Complexity and Geometric Approximation of Riesz Energy Subset Selection*

Michael Emmerich and Soumen Atta, first preprint version, 7 November 2026.

From the directory containing `main_arxiv_v1.tex`, run:

    python anc/reproduce.py
    python anc/reproduce.py --swiss

The second command includes the first; choose one. The driver works on a temporary copy. It requires Python 3 and NumPy (tested with Python 3.12.2 and NumPy 1.26.4). No geodata downloads, commercial solvers or network connections are needed. Allow about 1 GB of free memory for the full Swiss audit. Each test either matches its reference or exits with an error.

| Article material | Directory | What is included |
| --- | --- | --- |
| Geometric rounding, spread and cluster walkthroughs | `walkthroughs/` | Exact rational inputs and selections; independent search audits; generators for the appendix figures. |
| Figure 1 and the Swiss audit | `swiss/` | Frozen synthetic candidates in metres (Swiss LV95), display outline, traversal, both reduced solutions, full optimum, refinement trace, central-exclusion control and TikZ generator. IDs are one-based in the CSV and zero-based in result arrays. See `swiss/README.md` for optional regeneration of the densification. |
| Figure 2, paired table, occupancy and refinement studies | `benchmarks/` | Runners, seeds, configurations, per-run records, exact energies, timings and aggregate statistics. `benchmarks/README.md` explains how to regenerate the reports or rerun every measurement. |
| Compact-domain and doubling-metric checks | `geometry/` | A square traversal, seven half-open cell cases and twelve finite metrics; small exact checks of the implementation examples, not substitutes for the general proofs. |
| Cardinality on fixed inputs | `benchmarks/code/fixed_inputs.py` | Six frozen point arrays; all 171 I/X/H grids; 57 matched comparisons; separate records of size exclusions, completed reduced solves and nine additional attempts without size gates. The default replay verifies counts, grid widths and selected energies without retiming. |

The benchmark replay regenerates summary statistics from the recorded runs. To repeat new timing measurements, run the full `scale.py` and `paired.py` drivers **sequentially**, following the benchmark README. Fresh timings and timeout outcomes can differ from the recorded ones. For a reduced search that times out or exceeds the preset budget, the reported occupancy is still meaningful, but there is no completed-solve approximation certificate.

To keep regenerated figures or inspect an intermediate JSON, copy the relevant directory and run its scripts directly. For instance, in `swiss/`, `python run.py` writes `dense-results.json` and `python make_figure.py` writes `figure-1-inline.tex`. In `walkthroughs/`, the three figure drivers are `render_epas_gallery.py`, `render_epas_spread.py`, and `render_clusters.py`. All TikZ needed to compile the paper is already inline.

To compile the paper, use pdfLaTeX on `main_arxiv_v1.tex` until references and table widths settle, or run `latexmk -pdf main_arxiv_v1.tex`. `PRIMEarxiv.sty` is provided next to the main source. The manuscript needs neither Python nor BibTeX, shell escape or external figure files.
