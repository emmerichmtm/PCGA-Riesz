# Reproduce Figure 1: ten facilities on 5,403 Swiss settlement sites

All reported optimization results in this archive are for the 5,403-site input, k=10, s=2, epsilon=1. Coordinates are projected integer metres (EPSG:2056). Artificial locations resemble settlement footprints; they do not identify building addresses or imply that a facility could be built there.

## Offline reproduction

Extract the files together. With Python 3.12 and NumPy (tested 3.12.2 and 1.26.4), run:

    python run.py
    python make_figure.py

The first command validates the input hash, checks the search against exhaustive enumeration on a small restriction, and recomputes the exact full optimum, exchange and hybrid grids, reduced optima, within-cell refinement, and the central-exclusion control. It writes dense-results.json; wall times can vary. The second writes the TikZ geometry used in Figure 1 to figure-1-inline.tex. No geodata downloads or optimization packages are required. Allow approximately 1 GB of free memory for the exact reference run.

The expected exchange and hybrid counts are 3,365 and 706. The hybrid reduced energy is 2.8215031282144387e-9, and after eight accepted moves its energy is 2.8123482117390995e-9, exactly equal as a rational number to the independently searched full optimum. The Hasliberg cloud contributes 150 candidates and eight hybrid representatives. CSV site P0001 is result index 0. Full fractions, site IDs, traversal radii, all representatives, replacements and search counts are in dense-results.json.

## Optional regeneration of the dense cloud

The frozen baseline-1801.csv is the displaced seed geometry (not a separate reported optimization). To regenerate candidate-sites.csv using the original PCG64 seed and rejection tests, run:

    python make_dense.py national-address-directory.csv.zip swissBOUNDARIES3D.gpkg

This optional step also uses pandas 2.2.1, SciPy 1.14.1 and Shapely 2.2.0. Use the national LV95 [swisstopo official directory of building addresses](https://www.swisstopo.admin.ch/en/official-directory-of-building-addresses) snapshot of 5 November 2026 and the detailed GeoPackage in the [January 2026 swissBOUNDARIES3D edition](https://data.geo.admin.ch/api/stac/v1/collections/ch.swisstopo.swissboundaries3d/items/swissboundaries3d_2026-01). File hashes, generation seed and clearance checks are in generation.json. The directory and polygon are used for checking offsets and containment; the supplied simplified outline is for display. The first 1,801 rows of the final CSV are the fixed seed positions. Original address identifiers and the displacement correspondence are not released.
