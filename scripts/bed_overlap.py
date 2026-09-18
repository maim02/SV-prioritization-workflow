import argparse
import os

from helper import perform_bed_overlap


### perform overlap of two bed files based on minimal overlap


if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--bed1",
        required=True
    )

    parser.add_argument(
        "--bed2",
        required=True
    )

    parser.add_argument(
        "--output",
        required=True
    )

    parser.add_argument(
        "--min_overlap",
        required=True
    )

    args = parser.parse_args()


    out_dir = os.path.dirname(args.output)
    os.makedirs(out_dir, exist_ok=True)


    perform_bed_overlap(
        bed1=args.bed1,
        bed2=args.bed2,
        output_file_name=args.output,
        min_overlap=float(args.min_overlap)
    )