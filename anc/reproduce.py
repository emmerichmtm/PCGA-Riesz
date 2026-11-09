"""Run the arXiv replay, or export the traditional appendix layout for individual studies."""
import argparse
import csv
from hashlib import sha256
from pathlib import Path
import shutil
import subprocess
import sys
from tempfile import TemporaryDirectory

ANC = Path(__file__).resolve().parent


def export_paper(destination):
    """Copy the three manuscript files and restore each archived ancillary path."""
    with (ANC / "manifest.csv").open(encoding="utf-8", newline="") as manifest:
        for entry in csv.DictReader(manifest):
            source = ANC / entry["file"]
            contents = source.read_bytes()
            if sha256(contents).hexdigest() != entry["sha256"]:
                raise ValueError(f"Archived file differs from its checksum: {source.name}")
            target = destination / entry["appendix_path"]
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes(contents)
    for name in ("main_arxiv_v1.tex", "PRIMEarxiv.sty", "00README.json"):
        shutil.copyfile(ANC.parent / name, destination / name)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--swiss", action="store_true", help="include the exact Swiss search and Figure 1")
    mode.add_argument("--export", type=Path, metavar="DIRECTORY",
                      help="write the paper's original runnable layout to a new directory")
    args = parser.parse_args()
    if args.export is not None:
        destination = args.export.resolve()
        if destination.exists():
            parser.error("--export requires a new directory")
        destination.parent.mkdir(exist_ok=True, parents=True)
        with TemporaryDirectory(dir=destination.parent) as temp:
            staged = Path(temp) / "paper"
            export_paper(staged)
            staged.rename(destination)
        print(f"Exported manuscript, code and data to {destination}")
        return
    with TemporaryDirectory(prefix="riesz-paper-") as temp:
        root = Path(temp)
        export_paper(root)
        command = [sys.executable, str(root / "anc" / "reproduce.py")]
        if args.swiss:
            command.append("--swiss")
        subprocess.run(command, cwd=root, check=True)


if __name__ == "__main__":
    main()
