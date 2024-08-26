import argparse
import os
import sys
from pathlib import Path



def parse_args(args):
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "-i",
        "--tol-path",
        type=Path,
        required=True,
        help="Darwin Tree of Life Root Directory"
    )
    args = args[1:]
    return parser.parse_args(args)


def find_dtol_released_assemblies(tol_path: Path):
    counter = 0 
    with open("dtol_assemblies.txt", "w") as outfile:
        for i in os.listdir(tol_path):
            taxon_path = os.path.join(tol_path, i)
            for sample in os.listdir(taxon_path):
                sample_path = os.path.join(taxon_path, sample)
                sample_assembly_path = os.path.join(sample_path, "assembly", "release")
                if not os.path.exists(sample_assembly_path):
                    continue
                if len(os.listdir(sample_assembly_path)) == 0:
                    continue
                for sample_assembly_dir in os.listdir(sample_assembly_path):
                    fields = sample_assembly_dir.split(".")
                    sample = fields[0]
                    if len(fields) == 1:
                        continue
                    if not fields[1].isdigit():
                        continue
                    sample_assembly_version_path = os.path.join(sample_assembly_path, sample_assembly_dir, "insdc")
                    for file in os.listdir(sample_assembly_version_path):
                        if not file.startswith("GCA"):
                            continue
                        if file.endswith(".fasta.gz"):
                            outfile.write("{}\t{}\t{}\n".format(sample, sample_assembly_dir, os.path.join(sample_assembly_version_path, file)))
                            break


def main():
    options = parse_args(sys.argv)
    find_dtol_released_assemblies(options.tol_path)
    return 0


if __name__ == "__main__":
    exit_code = main()
    sys.exit(exit_code)
