#!/usr/bin/env python3

import subprocess
import argparse
import pandas as pd
import helper
from pathlib import Path
import re


def run_jasmine(input_vcfs, output_vcf, jasmine_path="jasmine",threads=1):
    """
    Merge SV VCFs using Jasmine.

    Args:
        input_vcfs (list[str]): Input VCF files.
        output_vcf (str): Output merged VCF.
        jasmine_path (str): Jasmine executable.
    """

    with open("tmp.txt", "w") as f:
        f.write("\n".join(input_vcfs) + "\n")


    cmd = [
        jasmine_path,
        "file_list=tmp.txt",
        f"out_file={output_vcf}",
        "--output_genotypes",
        f"--threads {threads}"
    ]

    print("Running Jasmine:")
    print(" ".join(cmd))

    result = subprocess.run(
        cmd,
        check=True,
        capture_output=True,
        text=True,
    )

    print(result.stdout)
    print(result.stderr)

    if Path("tmp.txt").exists():
        Path("tmp.txt").unlink()



if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Run Jasmine SV merging and extract variants present in affected samples only."
    )

    parser.add_argument("--input_vcfs",nargs="+",required=True)
    parser.add_argument("--output_vcf",required=True)
    parser.add_argument("--jasmine_path",default="jasmine")
    parser.add_argument("--threads",default=1)

    args = parser.parse_args()

    Path(args.output_vcf).parent.mkdir(parents=True, exist_ok=True)


    run_jasmine(
        input_vcfs=args.input_vcfs,
        output_vcf=args.output_vcf,
        jasmine_path=args.jasmine_path,
        threads=args.threads
    )