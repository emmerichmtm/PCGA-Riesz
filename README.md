# PCGA-Riesz

**Parameterized Complexity and Geometric Approximation of Riesz Energy Subset Selection**  
Michael Emmerich and Soumen Atta · Faculty of Information Technology, University of Jyväskylä, Finland

The 35-page PDF is the version arXiv compiled for the update to submission 8351874. The LaTeX source includes the submitted abstract. [Read the paper](main_arxiv_v1.pdf) · [LaTeX source](main_arxiv_v1.tex) · [Submitted PDF and source archive](https://github.com/emmerichmtm/PCGA-Riesz/releases/tag/arxiv-submission-8351874-r1)

The repository also contains the planar Riesz-energy selection algorithm, the candidate sites, the measured runs and the code that checks the reported solutions. [Experimental data release](https://github.com/emmerichmtm/PCGA-Riesz/releases/tag/arxiv-v1)

## Reproduce the results

From the repository root, with Python 3.12 and NumPy, run

```sh
python -m pip install -r requirements.txt
python anc/reproduce.py --swiss
```

This command rebuilds the benchmark summaries and fixed-input grids, checks the geometric examples, and reruns both reduced Swiss searches and the independent full-input search. The hybrid grid should retain **706 of 5,403 sites**. Eight within-cell replacements recover the full optimum. The Swiss search needs about 1 GB of free memory. Use `python anc/reproduce.py` for the shorter replay. Both commands work offline once NumPy is installed.

## Data by experiment

| Experiment | Inputs | Recorded results |
| --- | --- | --- |
| Swiss facilities, Fig. 1 | [5,403 sites in CSV](anc/swiss-candidate-sites.csv) · [1,801 seed sites](anc/swiss-baseline-1801.csv) · [map outline](anc/swiss-boundary-lv95.json) | [Exact energies, selected sites and refinement steps](anc/swiss-dense-results.json) · [data provenance](anc/swiss-guide.md) |
| Scale and refinement | [141 point arrays, up to one million sites](https://github.com/emmerichmtm/PCGA-Riesz/releases/download/arxiv-v1/scale-points.npz) · [configurations](anc/scale-configs.json) | [156 run records](anc/scale.jsonl) · [full logs with all representative indices](https://github.com/emmerichmtm/PCGA-Riesz/releases/download/arxiv-v1/raw-scale.jsonl.gz) |
| Paired certificate comparison, Fig. 2 | [99 point arrays](https://github.com/emmerichmtm/PCGA-Riesz/releases/download/arxiv-v1/paired-points.npz) · [configurations](anc/paired-configs.json) | [114 comparisons of five methods](anc/paired.jsonl) · [full logs with all representative indices](https://github.com/emmerichmtm/PCGA-Riesz/releases/download/arxiv-v1/raw-paired.jsonl.gz) |
| Varying cardinality on fixed sites, Table 8 | [Six frozen arrays in JSON](anc/fixed-input-points.json) | [171 grids, completed solves and timeouts](anc/fixed-input-runs.json) · [all 171 representative lists](https://github.com/emmerichmtm/PCGA-Riesz/releases/download/arxiv-v1/fixed-grid-representatives.jsonl.gz) |
| Validation and occupancy | [Validation and fixed-input arrays](https://github.com/emmerichmtm/PCGA-Riesz/releases/download/arxiv-v1/validation-and-fixed-points.npz) · [input index](anc/points-index.csv) | [Validation](anc/validation.json) · [grid occupancy](anc/occupancy.json) · [geometric checks](anc/geometry-checks-v3.json) |
| Worked examples, Appendix C | [Five of ten](anc/epas_gallery_input.json) · [100 spread sites](anc/epas_100_input.json) · [100 clustered sites](anc/compressed_hundred_v2.json) | [Ten-site walkthrough](anc/epas_gallery_results.json) · [100-site walkthrough](anc/epas_100_results.json) |

The [point index](anc/points-index.csv) maps each benchmark case to its archive, array key, coordinate denominator and SHA-256 hash. Cases that share the same ordered sites can reuse an array. Each NumPy array stores two integer coordinate numerators per row. For example, after downloading `scale-points.npz`, read the first scale input with

```python
import numpy as np

with np.load("scale-points.npz", allow_pickle=False) as inputs:
    xy = inputs["v3-half-000"].astype(np.float64) / 2**20
```

The recorded input hashes use the integer numerators in row order as little-endian signed 64-bit integers. The Switzerland CSV uses integer metres in the Swiss LV95 projection (EPSG:2056). Its candidates are artificial locations displaced from registered addresses. The data-generation record gives the seed, source hashes and minimum distance from registered coordinates. Geodata © Federal Office of Topography swisstopo. [Generation record](anc/swiss-generation.json) · [Official directory of building addresses](https://www.swisstopo.admin.ch/en/official-directory-of-building-addresses) · [January 2026 boundary edition](https://data.geo.admin.ch/api/stac/v1/collections/ch.swisstopo.swissboundaries3d/items/swissboundaries3d_2026-01) · [Geodata terms](https://www.swisstopo.admin.ch/en/terms-of-use-free-geodata-and-geoservices)

## Build the paper or rerun a study

Run `pdflatex main_arxiv_v1.tex` from the root until references settle, or use `latexmk -pdf main_arxiv_v1.tex`. All figures and references are inline. The [PRIME style](PRIMEarxiv.sty) and its [attribution](anc/style-notices.txt) are included.

The single `anc/` directory contains the code and study records. To work with the individual commands and paths printed in Appendix G, export them to a new directory beside the checkout. This also checks the archived file hashes.

```sh
python anc/reproduce.py --export ../riesz-paper
cd ../riesz-paper/anc/benchmarks
python code/fixed_inputs.py --check
```

The [submitted source ZIP](https://github.com/emmerichmtm/PCGA-Riesz/releases/download/arxiv-submission-8351874-r1/arxiv-source-v1.zip) contains the same 62 files as this export. The [benchmark guide](anc/benchmarks-guide.md) gives the commands for fresh timed runs. Wall-clock results can vary by machine. To recreate the downloadable coordinate arrays from the published generators, run `python anc/export_points.py ../riesz-points` in the checkout. [Experimental release checksums](anc/release-SHA256SUMS.txt)
