"""Reproduce the two-origin, five-cluster example using the archived input arrays."""
import json
from fractions import Fraction
from math import comb
from pathlib import Path

from epas_refine_v2 import dataset, energy, fingerprint
from paired import certificates, grid


ROOT = Path(__file__).resolve().parent.parent


def five_cluster_example(record):
    config = record["config"]
    points = dataset(config)
    assert fingerprint(points) == record["dataset_sha256"]
    k = config["k"]
    bounds, _ = certificates(points, k)
    counts = []
    for origin in (Fraction(0), Fraction(-1, 256)):
        for method in ("exchange", "hybrid"):
            cells = grid(points, k, config["epsilon"], bounds[method], origin)
            count = cells["M"]
            result = dict(method=method, origin_on_both_axes=str(origin), h=cells["h"],
                          M=count, subsets=comb(count, k))
            if count == k:
                # There is exactly one reduced k-subset, so its exact energy is V.
                representatives = cells["representatives"]
                value = energy(points, representatives)
                result.update(indices=representatives, reduced_optimum=str(value),
                              original_lower_bound=str(value / (1 + Fraction(config["epsilon"]))))
            counts.append(result)
    return dict(config=config, input_sha256=record["dataset_sha256"], grids=counts)


def main():
    with (ROOT / "results/scale.jsonl").open("r", encoding="utf8") as source:
        required = dict(n=100000, k=5, c=5, epsilon="1/2", kind="clusters", rho="1/4096")
        examples = []
        for line in source:
            record = json.loads(line)
            cfg = record["config"]
            if "translation" not in cfg and all(cfg[key] == value for key, value in required.items()):
                examples.append(five_cluster_example(record))
    assert sorted(example["config"]["seed"] for example in examples) == [17, 29, 43]
    result_path = ROOT / "results/occupancy.json"
    result_path.write_text(json.dumps(examples, indent=2) + "\n", encoding="utf8")
    for example in examples:
        print(example["config"]["seed"], [(r["method"], r["origin_on_both_axes"], r["h"], r["M"])
                                          for r in example["grids"]])


if __name__ == "__main__":
    main()
