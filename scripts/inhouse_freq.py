import os
import argparse
import pandas as pd
import gzip
import re
import helper
from pathlib import Path


def read_vcf(vcf:str):
    header_lines = []
    if vcf.endswith(".gz"):
        with gzip.open(vcf, "rt") as f:
            for line in f:
                if line.startswith("#CHROM"):
                    header = line.strip().lstrip("#").split("\t")
                    break
                elif line.startswith("#"):
                    header_lines.append(line)
    else:
        with open(vcf) as f:
            for line in f:
                if line.startswith("#CHROM"):
                    header = line.strip().lstrip("#").split("\t")
                    break
                elif line.startswith("#"):
                    header_lines.append(line)

    df = pd.read_csv(
            vcf,
            sep="\t",
            comment="#",
            names=header
        )

    return df, header_lines


def parse_info(info_str:str):
    """
    Parse the INFO field of a VCF file into a dictionary.
    # info col: AN;AC;Het;Hom;Hemi;AF;TEN;calling_methodology;SVLEN;SVTYPE;END
    """
    info_dict = {}
    for item in info_str.split(";"):
        if "=" in item:
            key, value = item.split("=")
            info_dict[key] = value
        else:
            info_dict[item] = True
    return info_dict


def sample_match(row):
    if row["GT"] == "1/1":
        samples = [] if pd.isna(row["SAMPLES_HOM"]) else re.split(r"[;,]", str(row["SAMPLES_HOM"]))
    else:
        samples = [] if pd.isna(row["SAMPLES_HET"]) else re.split(r"[;,]", str(row["SAMPLES_HET"]))

    ten_samples = [] if pd.isna(row["TEN_sample"]) else re.split(r"[;,]", str(row["TEN_sample"]))

    samples = {s.strip() for s in samples if s.strip()}
    ten_samples = {s.strip() for s in ten_samples if s.strip()}

    return ten_samples == samples


def check_inhouse_freq(inhouse_file:str, variants_file:str, out_file:str):
    """
    Check the frequency of inhouse variants in the given variants file.
    """

    # check if dir exists, if not create it
    Path(out_file).parent.mkdir(parents=True, exist_ok=True)

    # reciprocal overlap between inhouse and variants
    helper.perform_bed_overlap(bed1=variants_file, bed2=inhouse_file, output_file_name=out_file, min_overlap=0.7)

    # read the overlap file and add a new column to the variants_df indicating whether the variant is inhouse or not
    overlap_df = pd.read_csv(out_file, sep="\t", header=None)
    overlap_df.columns = [
        "CHROM","AVG_START","AVG_END","GT","ID","SVTYPE","AVG_LEN","SAMPLES_HET","SAMPLES_HOM",
        "CHROM_i","START_i","END_i","SVTYPE_i",
        "Het","Hom","TEN","AF", 
    ]

    # rm all lines where SVTYPE_i != SVTYPE (mismatched svtype)
    overlap_df = overlap_df[overlap_df["SVTYPE_i"] == overlap_df["SVTYPE"]].copy()

    # rm all lines where GT == 1/1 and Hom == 0 or where GT == 0/1 and Het == 0 (mismatched genotype)
    overlap_df = overlap_df[
        ~((overlap_df["GT"] == "1/1") & (overlap_df["Hom"] == 0)) &
        ~((overlap_df["GT"] == "0/1") & (overlap_df["Het"] == 0))
    ]

    # rm all lines samples are the same
    overlap_df["TEN_sample"] = overlap_df["TEN"].str.replace(r"HOM|HET", "", regex=True)
    overlap_df["SAMPLE_MATCH"] = overlap_df.apply(sample_match, axis=1)
    overlap_df = overlap_df[~overlap_df["SAMPLE_MATCH"]].copy()

    #  all ten_samples do not have pa in the sample name -> healthy only -> inhouse variants have healthy containing this variant -> should be removed
    overlap_df = overlap_df[
        overlap_df["TEN_sample"].apply(
            lambda x: not all(s.strip() == "pa" for s in x.split(",") if s.strip())
        )
    ]

    overlap_df_out = overlap_df[["CHROM","AVG_START","AVG_END","GT","ID","SVTYPE","AVG_LEN","SAMPLES_HET","SAMPLES_HOM","CHROM_i","START_i","END_i","SVTYPE_i","Het","Hom","TEN","AF"]].copy()
    overlap_df_out.to_csv(
        os.path.join(out_file),
        sep="\t",
        header=False,
        index=False
    ) # == all variants to be removed


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
            description="Perform inhouse freq filtering"
        )
    
    parser.add_argument("--inhouse_file",required=True)
    parser.add_argument("--variants_file",required=True)
    parser.add_argument("--out_file",required=True)
    args = parser.parse_args()
    check_inhouse_freq(args.inhouse_file, args.variants_file, args.out_file)


    