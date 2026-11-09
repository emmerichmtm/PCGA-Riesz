# Riesz-energy selection: synthetic benchmarks

The paired benchmark runs at epsilon = 1/2; the scale benchmark uses 1/2 and 1.
This is the runnable companion to the scale, refinement and paired-certificate
studies in the article and its appendices. The adjacent swiss directory
supplies the 5,403-site example at epsilon = 1.

## Regenerate tables and plot from the measurements

Use Python 3 with NumPy. The recorded environment was Python 3.12.2 and NumPy
1.26.4 on Windows 11, AMD Ryzen 7 5800H. From the extracted root, run:

    python code/report.py

This reads compact per-run records in results/ and the nine validation runs
in measurements/validation.json. It regenerates results/summary.json, the
TikZ series in results/paired-series.tex and all medians in results/paired-rows.tex.
The cap table is in summary.json, scale.caps; its named time columns use
milliseconds. Main Figure 2 uses only the occupancy records; all range bars
are minima and maxima over three seeds, and all curves use their median.

## Run the experiments again on the original point arrays

Use a fresh extraction separate from the archived copy. The three drivers
below do not require network access. Run the timing studies sequentially:

    python code/scale.py
    python code/paired.py
    python code/validation.py
    python code/checks.py
    python code/report.py

The first two drivers checkpoint full JSON lines for every configuration and
resume by case ID. Report generation prefers fresh raw JSONL when present.
Configuration files, input hashes, integer grid widths (stored as rationals),
certificate parameters, timings and exact reported energies are supplied.
The compact records store ordered representative counts and SHA-256 hashes
in place of bulky lists; full reruns emit the original lists. To verify a
representative hash, encode zero-based row indices as little-endian uint64
in the order emitted by the runner. Repetition may change timings and
wall-clock timeout status. Some internal output basenames contain “v2”
because the original certified solver and raw-record format are reused.

There are 156 scale runs on 52 configurations with seeds 17, 29, 43. The
paired study has 114 cases and 570 method evaluations; 42 of the cases
measure occupancy only. The 20-site validation includes a separate exhaustive
original optimum. File measurements/paired/verification-v2.json records the
157 independent checks. Full input data are generated deterministically from
the configurations; no real addresses or source geodata are used in these
synthetic benchmarks.

## Check occupancy before subset search

Run `python code/occupancy.py` to reproduce the worked five-cluster example.
The script regenerates three of the 100,000-site inputs and verifies their hashes.
It counts occupied cells for exchange and hybrid on each of two common physical
origins, and records the resulting binomial counts in results/occupancy.json.
When the hybrid retains exactly five sites, the sole five-subset supplies its
exact reduced optimum and the corresponding lower bound on the original energy.

## Vary cardinality without changing the data

From this `anc/benchmarks/` directory:

    python code/fixed_inputs.py --summary
    python code/fixed_inputs.py --check
    python code/fixed_inputs.py --run

The first command prints the archived results. The second rebuilds and checks the six actual arrays, all 171 I/X/H grids, the eligibility decisions, the ordered representative hashes, and the exact energy and certificate arithmetic for selections returned by completed searches. It does not repeat or certify timeouts. The third repeats the timed searches and replaces the two fixed-input JSON files; use it in a fresh copy. The main `anc/reproduce.py` driver includes `--check`.

`results/fixed-input-points.json` contains the exact integer coordinate numerators (common denominator 2^20), generator parameters, seeds and ordered input hashes. `results/fixed-input-runs.json` contains one record per (panel, seed, k), with bounds, grid widths, occupancies, binomial counts, feasibility decisions, solution indices and exact energies when available. Its `paired_split` separately audits the search-planned and occupancy-only cases of the original 114-case study. The new study keeps the coordinates, ordering, physical origin (0,0) and epsilon=1/2 fixed when varying k. It uses three seeds (17, 29, 43) for each panel: 24 irregular groups / 2,048 points and 256 lattice groups / 1,024 points. Primary size gates are M <= 160 and binom(M,k) <= 10^7, followed by a five-second wall limit. At k=128 on the lattice panel the three methods are also run without either size gate, with the same five-second wall limit. The clock includes reduced pair-weight construction and excludes worker startup and grid construction. Wall-clock outcomes are machine-dependent.
