#!/usr/bin/env python3

import subprocess
import argparse
import pandas as pd
import helper
from pathlib import Path
import re



def parse_jasmine_info(info):
    """
    This function parses the INFO field of a Jasmine VCF entry and extracts relevant information into a pandas Series.
    The INFO field is expected to be a semicolon-separated string of key-value pairs, where each key and value are separated by an equals sign (=). 
    The function extracts specific keys such as AVG_START, AVG_END, AVG_LEN, SVTYPE, converting them to appropriate data types (e.g., integers or floats) as needed. 
    If a key is not present in the INFO field, the corresponding value in the Series will be None.
    Args:
        info (str): The INFO field from a Jasmine VCF entry, formatted as a semicolon-separated string of key-value pairs.
    Returns:
        pd.Series: A pandas Series containing the extracted information.
    """
    d = {}

    for item in info.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            d[key] = value

    return pd.Series({
        "AVG_START": int(float(d["AVG_START"])) if "AVG_START" in d else None,
        "AVG_END": int(float(d["AVG_END"])) if "AVG_END" in d else None,
        "AVG_LEN": abs(int(float(d["AVG_LEN"]))) if "AVG_LEN" in d else None,
        "SVTYPE": d.get("SVTYPE"),
        #"ID_LIST": d.get("IDLIST"),
        #"SUPP_VEC":d.get("SUPP_VEC")
    })


def _get_gt(cell, gt_index=0):
    if pd.isna(cell):
        return "./."
    return cell.split(":")[gt_index]


def analyze_variant(row, affected, all_samples, het="0/1", hom="1/1"):
    """
    This function analyzes a variant row to determine the genotype (GT) of the variant for all samples.
    It checks if the variant is present in affected samples and absent in unaffected samples, 
    and returns the genotype (0/1 or 1/1 or None) along with lists of samples that are heterozygous or homozygous for the variant.
    GT=None indicates that the variant is not present in affected samples or is present in unaffected samples and needs to be filtered out.
    args:
        row (pd.Series): a row from the jasmine output dataframe
        affected (list): list of affected sample names
        all_samples (list): list of all sample names
        het (str): genotype string for heterozygous (default "0/1")
        hom (str): genotype string for homozygous (default "1/1")
    return: pd.Series with (GT (str), SAMPLES_HET (list), SAMPLES_HOM (list))
    """
    het_samples = []
    hom_samples = []

    has_het_aff = False
    has_hom_aff = False
    has_het_unaff = False
    has_hom_unaff = False

    # for each sample, check if it has the variant and if it's affected or unaffected
    for sample in all_samples:
        # get the genotype for the sample
        gt = _get_gt(row[sample])

        # check if the genotype is heterozygous or homozygous
        if gt == het:
            het_samples.append(sample)

            # check if the sample is affected or unaffected
            if sample in affected:
                has_het_aff = True
            else:
                has_het_unaff = True

        elif gt == hom:
            hom_samples.append(sample)

            if sample in affected:
                has_hom_aff = True
            else:
                has_hom_unaff = True

    # if there is at least one affected sample with the variant and no unaffected samples with the variant, add genotype info else none
    if has_het_aff and not has_het_unaff:
        variant_gt = het
    elif has_hom_aff and not has_hom_unaff:
        variant_gt = hom
    else:
        variant_gt = None

    return pd.Series({
        "GT": variant_gt,
        "SAMPLES_HET": ",".join(het_samples) if het_samples else "None",
        "SAMPLES_HOM": ",".join(hom_samples) if hom_samples else "None",
    })


def samples_are_only_pa(value):
    """
    Return True if all listed samples contain 'pa'.

    Empty/NA values are considered True -> so not removed.
    """
    if pd.isna(value):
        return True
    elif value == "None":
            return True
    
    samples = {
        sample.strip()
        for sample in re.split(r"[;,]", str(value))
        if sample.strip()
    }

    return all("pa" in sample for sample in samples)


def remove_heterozygous_when_homozygous(df):
    """
    Remove 0/1 rows when the same SV has a 1/1 call.

    Exception:
        Keep the 0/1 variant if all homozygous samples are 'pa'.

    """
    before = len(df)

    # Determine which to drop
    df["remove"] = (
        (df["GT"] == "0/1") & # is heterozygous
        (~df["SAMPLES_HOM"].apply(samples_are_only_pa)) # & homozygous call has not only pa
    )

    df = df[~df["remove"]]

    df = df.drop(
        columns=["remove"]
    )

    print(
        f"0/1 variants superseded by 1/1: "
        f"{before:,} -> {len(df):,}"
    )

    return df


def extract_affected_sv(jasmine_vcf, out_file, samples):
    """
    Extracts structural variants (SVs) from a Jasmine VCF file that are present in affected samples and absent in unaffected samples.
    Also perform GT based filtering and remove all variants that are called heterozygously but found homozygously on unaffected.
    Args:
        jasmine_vcf (str): Path to the Jasmine VCF file.
        out_file (str): Path to the output BED file where the filtered SVs will be saved.
        samples (list): List of sample names to consider for filtering.
    """

    # samples should contain your sample names
    affected = [
        s for s in samples
        if "pa" in s
    ]

    print("Extract relevant variants in affected-only samples and perform GT based filtering...")

    jasmine_out_df, _ = helper.read_vcf(jasmine_vcf)

    info_columns = [
        "CHROM",
        "POS",
        "ID",
        "REF",
        "ALT",
        "QUAL",
        "FILTER",
        "INFO",
        "FORMAT",
    ]

    sample_columns = helper.get_jasmine_sample_order(
        info_columns,
        jasmine_out_df
    )

    jasmine_out_df.columns = info_columns + sample_columns

    # get gt info 
    gt_analysis = jasmine_out_df.apply(
        analyze_variant,
        axis=1,
        affected=affected,
        all_samples=sample_columns
    )


    jasmine_out_df = pd.concat(
        [jasmine_out_df, gt_analysis],
        axis=1
    )

    # keep only variants present in affected samples
    merged_jasmine_out_df = (
        jasmine_out_df
        .dropna(subset=["GT"])
        .copy()
    )

    # filter based on GT
    merged_jasmine_out_df = remove_heterozygous_when_homozygous(merged_jasmine_out_df)

    # extract relevant info from merged output and write to bed file
    jasmine_info_df = (merged_jasmine_out_df["INFO"].apply(parse_jasmine_info))
    jasmine_info_df["ID"] = merged_jasmine_out_df["ID"]
    jasmine_info_df["CHROM"] = merged_jasmine_out_df["CHROM"]
    jasmine_info_df["AVG_START"] = (jasmine_info_df["AVG_START"] - 1) # VCF -> BED coordinate conversion
    jasmine_info_df["GT"] = merged_jasmine_out_df["GT"]
    jasmine_info_df["SAMPLES_HET"] = (merged_jasmine_out_df["SAMPLES_HET"])
    jasmine_info_df["SAMPLES_HOM"] = (merged_jasmine_out_df["SAMPLES_HOM"])

    jasmine_info_df = jasmine_info_df[
        [
            "CHROM",
            "AVG_START",
            "AVG_END",
            "GT",
            "ID",
            "SVTYPE",
            "AVG_LEN",
            "SAMPLES_HET",
            "SAMPLES_HOM",
        ]
    ] # reorder columns

    jasmine_info_df.to_csv(
        out_file,
        sep="\t",
        header=False,
        index=False
    )


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Extract variants present in affected samples only and perform GT based filtering"
    )

    parser.add_argument("--input_vcf",required=True)
    parser.add_argument("--out_bed",required=True)
    parser.add_argument("--samples",nargs="+",required=True)

    args = parser.parse_args()

    Path(args.out_bed).parent.mkdir(parents=True, exist_ok=True)

    extract_affected_sv(
        args.input_vcf,
        args.out_bed,
        args.samples
    )