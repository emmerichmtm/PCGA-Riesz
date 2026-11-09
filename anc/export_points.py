"""Export the ordered planar inputs, case index and fixed-input representative lists."""
import argparse
import csv
import gzip
import json
from pathlib import Path
from tempfile import TemporaryDirectory
from zipfile import ZIP_DEFLATED, ZipFile, ZipInfo

import numpy as np
from epas_refine_v2 import D, dataset, fingerprint
from fixed_inputs import METHODS, grids
from paired import make_input

ANC = Path(__file__).resolve().parent
FIELDS = ("study", "case_id", "n", "k", "epsilon", "seed", "archive", "array", "denominator", "sha256_i64le")


def read_json(name):
    return json.loads((ANC / name).read_text(encoding="utf-8"))


def records(study):
    if study in ("scale", "paired"):
        hash_key = "dataset_sha256" if study == "scale" else "input_sha256"
        generator = dataset if study == "scale" else make_input
        with (ANC / (study + ".jsonl")).open(encoding="utf-8") as source:
            for line in source:
                item = json.loads(line)
                yield item["id"], item["config"], item[hash_key], generator
    elif study == "validation":
        for i, item in enumerate(read_json("validation.json")):
            yield f"validation-{i:02}", item["config"], item["sha256"], dataset
    elif study == "fixed":
        points = read_json("fixed-input-points.json")
        for row in read_json("fixed-input-runs.json")["rows"]:
            name = f"{row['panel']}-{row['seed']}"
            item = points[name]
            config = dict(item["config"], k=row["k"], epsilon="1/2")
            if item["sha256"] != row["input_sha256"]:
                raise ValueError(f"Fixed-input record disagrees with the archived coordinates: {name}")
            yield f"{name}-k{row['k']}", config, item["sha256"], make_input


def write_array(archive, name, points):
    """Store exact numerators as little-endian int32 after a checked conversion."""
    numbers = points.astype("<i4")
    if not np.array_equal(numbers, points):
        raise OverflowError(f"Coordinates exceed int32 for {name}")
    member = ZipInfo(name + ".npy", (1980, 1, 1, 0, 0, 0))
    member.compress_type = ZIP_DEFLATED
    with archive.open(member, "w") as stream:
        np.lib.format.write_array(stream, numbers, allow_pickle=False)


def write_fixed_grids(staging):
    points = read_json("fixed-input-points.json")
    arrays = {key: np.asarray(item["numerators"], dtype=np.int64) for key, item in points.items()}
    with (staging / "fixed-grid-representatives.jsonl.gz").open("wb") as destination:
        with gzip.GzipFile(filename="", fileobj=destination, mode="wb", mtime=0) as compressed:
            for row in read_json("fixed-input-runs.json")["rows"]:
                key = f"{row['panel']}-{row['seed']}"
                all_grids, summary = grids(arrays[key], row["k"])
                if summary["methods"] != row["methods"]:
                    raise ValueError(f"Grid differs from the archived run for {key}, k={row['k']}")
                for method in METHODS:
                    record = dict(panel=row["panel"], seed=row["seed"], k=row["k"], method=method,
                                  input_sha256=row["input_sha256"],
                                  representatives=all_grids[method]["representatives"],
                                  representatives_sha256=row["methods"][method]["reps_sha256"])
                    compressed.write((json.dumps(record, separators=(",", ":")) + "\n").encode("utf-8"))


def export(staging):
    index = []
    files_by_hash = {}
    for filename, studies in (("scale-points.npz", ("scale",)),
                              ("paired-points.npz", ("paired",)),
                              ("validation-and-fixed-points.npz", ("validation", "fixed"))):
        saved = {}
        with ZipFile(staging / filename, "w", compression=ZIP_DEFLATED) as archive:
            for study in studies:
                for case, config, digest, generator in records(study):
                    if digest not in saved:
                        points = generator(config)
                        if points.shape != (config["n"], 2) or fingerprint(points) != digest:
                            raise ValueError(f"Input hash or shape mismatch for {case}")
                        write_array(archive, case, points)
                        saved[digest] = case
                    location = (filename, saved[digest])
                    files_by_hash.setdefault(digest, location)
                    index.append((study, case, config["n"], config["k"], config["epsilon"],
                                  config["seed"], *location, D, digest))
        print(f"{filename}: {len(saved)} distinct ordered arrays", flush=True)
    for item in read_json("occupancy.json"):
        cfg = item["config"]
        digest = item["input_sha256"]
        index.append(("occupancy", f"occupancy-{cfg['seed']}", cfg["n"], cfg["k"], cfg["epsilon"],
                      cfg["seed"], *files_by_hash[digest], D, digest))
    with (staging / "points-index.csv").open("w", encoding="utf-8", newline="") as output:
        writer = csv.writer(output, lineterminator="\n")
        writer.writerow(FIELDS)
        writer.writerows(index)
    write_fixed_grids(staging)
    print(f"Indexed {len(index)} experimental cases.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path, help="directory for the coordinate archives, index and fixed-grid lists")
    destination = parser.parse_args().directory.resolve()
    destination.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(dir=destination) as temp:
        staging = Path(temp)
        export(staging)
        for path in staging.iterdir():
            path.replace(destination / path.name)


if __name__ == "__main__":
    main()
