"""Replay the article's examples in a temporary directory; no network required."""
import argparse
import json
from fractions import Fraction
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile


def load(path):
    return json.loads(path.read_text(encoding="utf-8"))


def check(actual, expected, description):
    if actual != expected:
        raise ValueError("Replay disagrees with the record: " + description)


def run(where, *scripts):
    for step in scripts:
        script, *args = (step,) if isinstance(step, str) else step
        print("  " + where.name + "/" + " ".join((script, *args)), flush=True)
        command = [sys.executable, str(where / script), *args]
        proc = subprocess.run(command, cwd=where, stdout=subprocess.PIPE,
                              stderr=subprocess.STDOUT, text=True, errors="replace")
        if proc.returncode:
            raise RuntimeError(script + " failed:\n" + proc.stdout)


def selected_fields(data, *keys):
    return {key: data[key] for key in keys}


def replay(source, dest, swiss):
    b, ref = dest/"benchmarks", source/"benchmarks"
    run(b, "code/report.py", "code/occupancy.py", "code/checks.py")
    for file in ("results/summary.json", "results/occupancy.json",
                 "measurements/paired/verification-v2.json"):
        check(load(b/file), load(ref/file), file)
    print("  Verified the 156-run and 114-case summaries, occupancy and 157 checks.")
    run(b, ("code/fixed_inputs.py", "--check"))
    print("  Verified the six fixed arrays, all 171 grids, search admissions and recorded solution energies.")

    g = dest/"geometry"
    run(g, "check_geometry_v3.py")
    check(load(g/"geometry-checks-v3.json"),
          load(source/"geometry/geometry-checks-v3.json"), "geometry")
    print("  Verified polygon traversal, seven cell tests and twelve metric instances.")

    w = dest/"walkthroughs"
    run(w, "epas_gallery_run.py", "epas_100_run.py", "check_exact_subset_search.py",
        "compressed_hundred_v2.py", "render_epas_gallery.py", "render_epas_spread.py",
        "render_clusters.py")
    for name,mode in (("epas_gallery_results.json","compact_run"),
                      ("epas_100_results.json","compact_run"),
                      ("compressed_hundred_v2.json","run")):
        actual, expected = load(w/name), load(source/"walkthroughs"/name)
        keys=("points","centers","representatives","cell_side","M","selected","energy")
        check(selected_fields(actual[mode],*keys),selected_fields(expected[mode],*keys),name)
        check(actual["audit"]["original_optimum"],expected["audit"]["original_optimum"],name)
    print("  Verified the 8-, 98- and 26-representative examples and their original optima.")
    if not swiss:
        return
    sw = dest/"swiss"
    run(sw, "run.py", "make_figure.py")
    actual = load(sw/"dense-results.json")
    expected = load(source/"swiss/dense-results.json")
    for key in ("full","best_outside_30km"):
        for field in ("selected","energy_exact"):
            check(actual[key][field], expected[key][field], "Swiss "+key+" "+field)
    check(actual["geometry"],expected["geometry"],"Swiss traversal")
    for method in ("exchange","hybrid"):
        a,e=actual["methods"][method],expected["methods"][method]
        keys=("h_m","M","representatives","moves","refined","refined_energy","certificate_lower")
        check(selected_fields(a,*keys),selected_fields(e,*keys),"Swiss "+method)
        check(a["reduced"]["energy_exact"],e["reduced"]["energy_exact"],"Swiss reduced energy")
        check(a["reduced"]["selected"],e["reduced"]["selected"],"Swiss reduced sites")
    h=actual["methods"]["hybrid"]
    check(Fraction(h["refined_energy"]), Fraction(actual["full"]["energy_exact"]), "Swiss recovery")
    print("  Verified Swiss full optimum, both reduced solves, eight hybrid moves and central-exclusion control.")


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--swiss", action="store_true", help="also solve the 5,403-site Swiss instance")
    args=parser.parse_args()
    source=Path(__file__).resolve().parent
    with tempfile.TemporaryDirectory(prefix="riesz-replay-") as scratch:
        dest=Path(scratch)/"anc"
        shutil.copytree(source,dest,ignore=shutil.ignore_patterns("__pycache__","*.pyc"))
        replay(source,dest,args.swiss)
    print("All requested replays agree with the archived results.")
