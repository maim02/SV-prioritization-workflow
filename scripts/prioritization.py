#!/usr/bin/env python3

"""
Prioritize disease/SCN structural variants.

Expected disease/SCN BED columns
---------------------------------
CHROM
AVG_START
AVG_END
GT
ID
SVTYPE
AVG_LEN
SAMPLES_HET
SAMPLES_HOM
CNV_SUPP
cCRE_CHROM
cCRE_START
cCRE_END
cCRE_ID
cCRE_SCORE
cCRE_STRAND
cCRE_THICK_START
cCRE_THICK_END
cCRE_RGB
cCRE_CLASS
cCRE_BIOSAMPLE
gene_id
gene_name


Usage
-----
python prioritize_variants.py \
    --variants_file filtered_sv.bed \
    --output_scn prioritized_scn_variants.tsv  \
    --output_all prioritized_variants.tsv
"""

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd

# todo:
# filter out: low DNAse
# sort based on cCREs class
# per class:
# sort based on SCN 
# sort based on number of patient per variant gt
# 


# -------------------------------------------------------------------------
# Configuration
# -------------------------------------------------------------------------

CCRE_CLASS_ORDER = [
    "PLS",
    "pELS",
    "dELS",
    "CA-H3K4me3",
    "CA-CTCF",
    "CA-TF",
    "CA-only",
    "Low-DNase",
]


SUBSET_COL = [
    "CHROM",
    "AVG_START",
    "AVG_END",
    "GT",
    "ID",
    "SVTYPE",
    "AVG_LEN",
    "SAMPLES_HET",
    "SAMPLES_HOM",
    "CNV_SUPP",
    "cCRE_ID",
    "cCRE_CLASS",
    "GENE_ID",
    "GENE_NAME",
    #"CNV_SUPP",
    "eQTL_SUPP",
    "PANEL_SUPP",
#    "NEAREST_GENE_STRAND",
#    "NEAREST_GENE_ID",
]


KEY_COLUMNS = [
    "CHROM",
    "AVG_START",
    "AVG_END",
    "GT",
]


# -------------------------------------------------------------------------
# General helpers
# -------------------------------------------------------------------------

def read_bed(path, header=None):
    """Read a tab-delimited file."""
    return pd.read_csv(
        path,
        sep="\t",
        header=header,
    )


def normalize_sample_string(value):
    """
    Convert a sample string into a set.

    Supports comma or semicolon separated sample names.
    """
    if pd.isna(value):
        return set()

    return {
        sample.strip()
        for sample in re.split(r"[;,]", str(value))
        if sample.strip()
    }


def samples_are_only_pa(value):
    """
    Return True if all listed samples contain 'pa'.

    Empty/NA values are considered True -> so not removed.
    """
    samples = normalize_sample_string(value)

    if not samples:
        return True

    return all("pa" in sample for sample in samples)


def count_samples(value):
    """Count unique samples in a comma/semicolon separated string."""
    return len(normalize_sample_string(value))


# -------------------------------------------------------------------------
# Variant filtering
# -------------------------------------------------------------------------

def remove_low_dnase(df):
    """
    Remove rows where the cCRE annotation is Low-DNase.
    """
    before = len(df)

    df = df[
        df["cCRE_CLASS"].fillna("") != "Low-DNase"
    ].copy()

    print(f"Low-DNase filtering: {before:,} -> {len(df):,}")

    return df


def merge_gene_lists(values):
    genes = set()

    for value in values.dropna():
        genes.update(
            gene.strip()
            for gene in str(value).split(",")
            if gene.strip()
        )

    return ",".join(sorted(genes))


def collapse_ccre_annotations(df):
    """
    Collapse duplicate SV rows caused by multiple cCRE annotations.

    cCRE_CLASS values are combined into a comma-separated list.
    """
    before = len(df)

    df_subset = df[SUBSET_COL]

    # Columns defining the same SV/annotation except cCRE class.
    group_cols = [col for col in df_subset.columns if col != "cCRE_ID" and col != "cCRE_CLASS" and col != "GENE_ID" and col != "GENE_NAME"]
    df_agg = (
        df_subset.groupby(group_cols, dropna=False, as_index=False)
        .agg({
            "cCRE_ID": lambda x: ",".join(sorted(set(x.dropna()))),
            "cCRE_CLASS": lambda x: ",".join(sorted(set(x.dropna()))),
            "GENE_NAME": merge_gene_lists,
            "GENE_ID": lambda x: ",".join(sorted(set(x.dropna())))
        })
    ).copy()

    df_agg = df_agg.drop_duplicates()

    print(f"Collapsed annotations: {before:,} -> {len(df_agg):,}")

    return df_agg


# -------------------------------------------------------------------------
# cCRE prioritization
# -------------------------------------------------------------------------

def highest_priority_class(value):
    """
    Return the highest-priority cCRE class in a combined annotation.
    """
    if pd.isna(value):
        return np.nan

    classes = {
        str(c).strip()
        for c in str(value).split(",")
        if str(c).strip()
    }

    for cCRE_class in CCRE_CLASS_ORDER:
        if cCRE_class in classes:
            return cCRE_class

    return np.nan


def add_ccre_priority(df):
    """
    Add sortable cCRE priority.
    """
    df["cCRE_CLASS_SORT"] = (
        df["cCRE_CLASS"]
        .apply(highest_priority_class)
    )

    df["cCRE_CLASS_SORT"] = pd.Categorical(
        df["cCRE_CLASS_SORT"],
        categories=CCRE_CLASS_ORDER,
        ordered=True,
    )

    return df


# -------------------------------------------------------------------------
# Patient prioritization
# -------------------------------------------------------------------------

def count_patients(row):
    """
    Count affected patients according to genotype.

    0/1 -> SAMPLES_HET
    1/1 -> SAMPLES_HOM
    """
    if row["GT"] == "0/1":
        return count_samples(row["SAMPLES_HET"])

    if row["GT"] == "1/1":
        return count_samples(row["SAMPLES_HOM"])

    return 0


def add_patient_count(df):
    """
    Add number of affected patients.
    """
    df["N_PATIENTS"] = (
        df.apply(
            count_patients,
            axis=1,
        )
    )

    return df

# -------------------------------------------------------------------------
# Panel gene prioritization
# -------------------------------------------------------------------------


def count_panel_genes(row):
    """
    Count panel genes
    """

    return len(row["PANEL_SUPP"])


def add_panel_genes_count(df):
    """
    Add number of involved panel genes.
    """
    df["N_PANEL_GENES"] = (
        df.apply(
            count_panel_genes,
            axis=1,
        )
    )

    return df

# -------------------------------------------------------------------------
# Final prioritization
# -------------------------------------------------------------------------

def prioritize(variants_file:str):
    """
    Main prioritization workflow.
    """

    print("=" * 70)
    print("Loading disease variants")
    print("=" * 70)

    df = pd.read_csv(variants_file, sep="\t")

    print(f"Initial variants/annotations: {len(df):,}")

    # -------------------------------------------------------------
    # cCRE filtering
    # -------------------------------------------------------------

    #df = remove_low_dnase(df)

    df = collapse_ccre_annotations(df)


    # -------------------------------------------------------------
    # cCRE priority
    # -------------------------------------------------------------

    df = add_ccre_priority(df)

    # -------------------------------------------------------------
    # Patient count
    # -------------------------------------------------------------

    df = add_patient_count(df)

    # -------------------------------------------------------------
    # Panel gene count
    # -------------------------------------------------------------

    df = add_panel_genes_count(df)

    # -------------------------------------------------------------
    # Final sorting
    # -------------------------------------------------------------

    df = df.sort_values(
        by=[
            "cCRE_CLASS_SORT",
            "N_PATIENTS",
            "N_PANEL_GENES",
            "eQTL_SUPP"
            
        ],
        ascending=[
            True,
            False,
            False,
            False
        ],
        na_position="last",
    )

    return df


# -------------------------------------------------------------------------
# Command-line interface
# -------------------------------------------------------------------------

def main():

    parser = argparse.ArgumentParser(
        description=(
            "Prioritize disease/SCN structural variants "
            "using cCRE class, patient count, involved panel genes and eQTL."
        )
    )

    parser.add_argument(
        "--variants_file",
        required=True,
        help=(
            "Disease/SCN cCRE-gene BED file. "
            "Must contain the 22 expected columns."
        ),
    )

    parser.add_argument(
        "--output_file",
        required=True,
        help="Output prioritized TSV.",
    )

    # parser.add_argument(
    #     "--output_panel",
    #     required=True,
    #     help="Output prioritized panel only variants TSV.",
    # )

    args = parser.parse_args()

    result = prioritize(variants_file=args.variants_file)

    result = result.drop(columns=["cCRE_CLASS_SORT", "N_PATIENTS", "N_PANEL_GENES"])

    result.to_csv(
        args.output_file,
        sep="\t",
        index=False,
    )

    # result[result["PANEL_SUPP"].apply(lambda x: len(x) > 0)].to_csv(
    #     args.output_panel,
    #     sep="\t",
    #     index=False,
    # )

    print("=" * 70)
    print("Prioritization complete")
    print("=" * 70)
    print(f"Final variants: {len(result):,}")
    print(f"Output: {args.output_file}")


if __name__ == "__main__":
    main()
