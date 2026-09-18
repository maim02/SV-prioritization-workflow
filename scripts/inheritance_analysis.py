import argparse
import os
import re
from collections import defaultdict

import pandas as pd

import helper


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

INFO_COLUMNS = [
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

VARIANT_GTS = {"0/1", "1/1"}


TRIO_INHERITANCE = {
    # -----------------------------------------------------------------------
    # Child only
    # -----------------------------------------------------------------------
    ("./.", "./.", "0/1"):
        "De novo candidate (het) (trio analysis)",

    ("./.", "./.", "1/1"):
        "De novo candidate (hom) (trio analysis)",

    # -----------------------------------------------------------------------
    # Maternal inheritance
    # -----------------------------------------------------------------------
    ("0/1", "./.", "0/1"):
        "Maternal inheritance (same het) (trio analysis)",

    ("1/1", "./.", "0/1"):
        "Maternal inheritance (hom→het) (trio analysis)",

    ("0/1", "./.", "1/1"):
        "Maternal inheritance (het→hom) (trio analysis)",

    ("1/1", "./.", "1/1"):
        "Maternal inheritance (same hom) (trio analysis)",

    # -----------------------------------------------------------------------
    # Paternal inheritance
    # -----------------------------------------------------------------------
    ("./.", "0/1", "0/1"):
        "Paternal inheritance (same het) (trio analysis)",

    ("./.", "1/1", "0/1"):
        "Paternal inheritance (hom→het) (trio analysis)",

    ("./.", "0/1", "1/1"):
        "Paternal inheritance (het→hom) (trio analysis)",

    ("./.", "1/1", "1/1"):
        "Paternal inheritance (same hom) (trio analysis)",

    # -----------------------------------------------------------------------
    # Both parents carry
    # -----------------------------------------------------------------------
    ("0/1", "0/1", "0/1"):
        "Familial variant (all het) (trio analysis)",

    ("0/1", "0/1", "1/1"):
        "Recessive candidate (trio analysis)",

    ("0/1", "1/1", "0/1"):
        "Familial variant (trio analysis)",

    ("1/1", "0/1", "0/1"):
        "Familial variant (trio analysis)",

    ("1/1", "1/1", "0/1"):
        "Unexpected child het (trio analysis)",

    ("1/1", "1/1", "1/1"):
        "Familial homozygous variant (trio analysis)",

    ("0/1", "1/1", "1/1"):
        "Inherited from both parents (trio analysis)",

    ("1/1", "0/1", "1/1"):
        "Inherited from both parents (trio analysis)",

    # -----------------------------------------------------------------------
    # Parent only
    # -----------------------------------------------------------------------
    ("0/1", "./.", "./."):
        "Mother only (het) (trio analysis)",

    ("1/1", "./.", "./."):
        "Mother only (hom) (trio analysis)",

    ("./.", "0/1", "./."):
        "Father only (het) (trio analysis)",

    ("./.", "1/1", "./."):
        "Father only (hom) (trio analysis)",

    ("0/1", "0/1", "./."):
        "Parents only (trio analysis)",

    ("0/1", "1/1", "./."):
        "Parents only (trio analysis)",

    ("1/1", "0/1", "./."):
        "Parents only (trio analysis)",

    ("1/1", "1/1", "./."):
        "Parents only (trio analysis)",

    # -----------------------------------------------------------------------
    # Absent
    # -----------------------------------------------------------------------
    ("./.", "./.", "./."):
        "Absent",
}


ONE_PARENT_INHERITANCE = {
    # -----------------------------------------------------------------------
    # Patient only
    # -----------------------------------------------------------------------
    ("./.", "0/1"):
        "Present only in patient (het) (one-parent family)",

    ("./.", "1/1"):
        "Present only in patient (hom) (one-parent family)",

    # -----------------------------------------------------------------------
    # Inherited from available parent
    # -----------------------------------------------------------------------
    ("0/1", "0/1"):
        "Possibly inherited from available parent (one-parent family)",

    ("1/1", "0/1"):
        "Inherited from available parent (one-parent family)",

    ("0/1", "1/1"):
        "Possible recessive inheritance (one-parent family)",

    ("1/1", "1/1"):
        "Likely inherited from both parents (one-parent family)",

    # -----------------------------------------------------------------------
    # Parent only
    # -----------------------------------------------------------------------
    ("0/1", "./."):
        "Parents only (one-parent family)",

    ("1/1", "./."):
        "Parents only (one-parent family)",

    # -----------------------------------------------------------------------
    # Absent
    # -----------------------------------------------------------------------
    ("./.", "./."):
        "Absent",
}


# ---------------------------------------------------------------------------
# Genotype helpers
# ---------------------------------------------------------------------------

def get_gt(cell):
    """
    Extract GT from a VCF sample field.

    Examples:
        '0/1:.:INS:0:.' -> '0/1'
        '1/1:.:INS:0:.' -> '1/1'
        './.:NA:NA:NA:NA' -> './.'
        NaN -> './.'
    """
    if pd.isna(cell):
        return "./."

    return str(cell).split(":", 1)[0]


def extract_gt_columns(df, samples):
    """
    Extract GT for the requested samples.

    Returns a DataFrame containing one GT column per sample.
    """
    gt_df = pd.DataFrame(index=df.index)

    for sample in samples:
        gt_df[sample] = df[sample].map(get_gt)

    return gt_df


# ---------------------------------------------------------------------------
# Family identification
# ---------------------------------------------------------------------------

def get_family_trios(samples):
    """
    Group samples by family and identify available parents and children.

    Returns columns:
        SAMPLE
        MOTHER
        FATHER
        ANALYSIS
    """

    family_pattern = re.compile(r"^(SCN\d+|U\d+|DPLD\d+)")

    families = defaultdict(list)

    for sample in samples:
        match = family_pattern.match(sample)

        if match:
            family_id = match.group(1)
            families[family_id].append(sample)

    results = []

    for members in families.values():

        mother = next(
            (sample for sample in members if "mo" in sample),
            None,
        )

        father = next(
            (sample for sample in members if "fa" in sample),
            None,
        )

        children = [
            sample
            for sample in members
            if "pa" in sample
            and "gf" not in sample
            and "gm" not in sample
            and "mo" not in sample
            and "fa" not in sample
        ]

        if mother and father:
            analysis = "Trio analysis"

        elif mother:
            analysis = "Mother and child only"

        elif father:
            analysis = "Father and child only"

        else:
            analysis = "No parents"

        for child in children:
            results.append(
                {
                    "SAMPLE": child,
                    "MOTHER": mother,
                    "FATHER": father,
                    "ANALYSIS": analysis,
                }
            )

    return pd.DataFrame(
        results,
        columns=[
            "SAMPLE",
            "MOTHER",
            "FATHER",
            "ANALYSIS",
        ],
    )


# ---------------------------------------------------------------------------
# Jasmine sample selection
# ---------------------------------------------------------------------------

def get_jasmine_selection(jasmine_vcf_df, samples):
    """
    Select INFO columns and requested sample columns from a Jasmine VCF.

    Jasmine sample columns may have an index prefix such as:
        0_SCN344pa
        1_SCN344mo
        2_SCN344fa

    The returned DataFrame contains the requested samples using their
    original sample names.
    """

    jasmine_sample_columns = helper.get_jasmine_sample_order(
        INFO_COLUMNS,
        jasmine_vcf_df,
    )

    selected_columns = INFO_COLUMNS + [
        f"{jasmine_sample_columns.index(sample)}_{sample}"
        for sample in samples
    ]

    selection = jasmine_vcf_df[selected_columns].copy()

    selection.columns = INFO_COLUMNS + samples

    return selection


# ---------------------------------------------------------------------------
# Variant filtering
# ---------------------------------------------------------------------------

def keep_patient_variants(df, patient):
    """
    Keep variants where the patient carries the variant.

    A variant is considered present when GT is:
        0/1
        1/1

    Parents are deliberately NOT used for filtering here.

    This is important because inheritance should be determined after
    retaining all variants present in the patient.
    """

    patient_gt = df[patient].map(get_gt)

    mask = patient_gt.isin(VARIANT_GTS)

    return df.loc[mask].copy()


# ---------------------------------------------------------------------------
# Inheritance classification
# ---------------------------------------------------------------------------

def classify_trio_inheritance(df, patient, mother, father):
    """
    Classify inheritance for a trio.
    """

    mother_gt = df[mother].map(get_gt)
    father_gt = df[father].map(get_gt)
    patient_gt = df[patient].map(get_gt)

    keys = zip(
        mother_gt,
        father_gt,
        patient_gt,
    )

    return pd.Series(
        [
            TRIO_INHERITANCE.get(
                key,
                "Complex/Review",
            )
            for key in keys
        ],
        index=df.index,
    )


def classify_one_parent_inheritance(df, patient, parent):
    """
    Classify inheritance when only one parent is available.
    """

    parent_gt = df[parent].map(get_gt)
    patient_gt = df[patient].map(get_gt)

    keys = zip(
        parent_gt,
        patient_gt,
    )

    return pd.Series(
        [
            ONE_PARENT_INHERITANCE.get(
                key,
                "Complex/Review",
            )
            for key in keys
        ],
        index=df.index,
    )


# ---------------------------------------------------------------------------
# Output construction
# ---------------------------------------------------------------------------

def build_output_dataframe(
    sv_df,
    patient,
    mother,
    father,
):
    """
    Create the final inheritance-analysis output.
    """

    info_df = sv_df["INFO"].apply(
        helper.parse_jasmine_info
    )

    output = pd.DataFrame(
        {
            "CHR": sv_df["CHROM"].values,
            "AVG_START": info_df["AVG_START"].values,
            "AVG_END": info_df["AVG_END"].values,
            "VARIANTid": sv_df["ID"].values,
            "AVG_SVLEN": info_df["AVG_LEN"].values,
            "SVTYPE": info_df["SVTYPE"].values,
            "PATIENT": patient,
            "MOTHER": mother,
            "FATHER": father,
            "GT_pa_mo_fa": sv_df["GT_pa_mo_fa"].values,
            "INHERITANCE": sv_df["INHERITANCE"].values,
            "ID_LIST": info_df["ID_LIST"].values,
        }
    )

    return output


# ---------------------------------------------------------------------------
# Family analysis
# ---------------------------------------------------------------------------

def analyze_family(
    jasmine_vcf_df,
    patient,
    mother,
    father,
    analysis,
):
    """
    Perform inheritance analysis for one family.
    """

    family_samples = [patient]

    if mother is not None:
        family_samples.append(mother)

    if father is not None:
        family_samples.append(father)

    # Select only the samples needed for this family.
    sv_df = get_jasmine_selection(
        jasmine_vcf_df,
        family_samples,
    )

    # Keep all variants present in the patient.
    # Do NOT filter based on parental genotypes here.
    sv_df = keep_patient_variants(
        sv_df,
        patient,
    )

    if sv_df.empty:
        return None

    # -----------------------------------------------------------------------
    # Trio
    # -----------------------------------------------------------------------

    if analysis == "Trio analysis":

        sv_df["INHERITANCE"] = classify_trio_inheritance(
            sv_df,
            patient,
            mother,
            father,
        )

        patient_gt = sv_df[patient].map(get_gt)
        mother_gt = sv_df[mother].map(get_gt)
        father_gt = sv_df[father].map(get_gt)

        sv_df["GT_pa_mo_fa"] = (
            patient_gt
            + ":"
            + mother_gt
            + ":"
            + father_gt
        )

    # -----------------------------------------------------------------------
    # Mother + child
    # -----------------------------------------------------------------------

    elif analysis == "Mother and child only":

        sv_df["INHERITANCE"] = classify_one_parent_inheritance(
            sv_df,
            patient,
            mother,
        )

        patient_gt = sv_df[patient].map(get_gt)
        mother_gt = sv_df[mother].map(get_gt)

        sv_df["GT_pa_mo_fa"] = (
            patient_gt
            + ":"
            + mother_gt
            + ":."
        )

    # -----------------------------------------------------------------------
    # Father + child
    # -----------------------------------------------------------------------

    elif analysis == "Father and child only":

        sv_df["INHERITANCE"] = classify_one_parent_inheritance(
            sv_df,
            patient,
            father,
        )

        patient_gt = sv_df[patient].map(get_gt)
        father_gt = sv_df[father].map(get_gt)

        sv_df["GT_pa_mo_fa"] = (
            patient_gt
            + ":."
            + ":"
            + father_gt
        )

    # -----------------------------------------------------------------------
    # Patient only
    # -----------------------------------------------------------------------

    elif analysis == "No parents":

        sv_df["INHERITANCE"] = "No parents, in patient"

        patient_gt = sv_df[patient].map(get_gt)

        sv_df["GT_pa_mo_fa"] = (
            patient_gt
            + ":.:."
        )

    else:
        raise ValueError(
            f"Unknown analysis type: {analysis}"
        )

    return build_output_dataframe(
        sv_df,
        patient,
        mother,
        father,
    )


# ---------------------------------------------------------------------------
# Main inheritance analysis
# ---------------------------------------------------------------------------

def inheritance_analysis(
    samples,
    jasmine_vcf_df,
):
    """
    Run inheritance analysis for all families.
    """

    family_df = get_family_trios(samples)

    if family_df.empty:
        return pd.DataFrame()

    all_results = []

    for _, family in family_df.iterrows():

        result = analyze_family(
            jasmine_vcf_df=jasmine_vcf_df,
            patient=family["SAMPLE"],
            mother=family["MOTHER"],
            father=family["FATHER"],
            analysis=family["ANALYSIS"],
        )

        if result is not None and not result.empty:
            all_results.append(result)

    if not all_results:
        return pd.DataFrame()

    return pd.concat(
        all_results,
        ignore_index=True,
    )


# ---------------------------------------------------------------------------
# Command-line interface
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Perform inheritance analysis on a Jasmine VCF."
    )

    parser.add_argument(
        "--jasmine_vcf",
        required=True,
        help="Jasmine VCF containing all samples.",
    )

    parser.add_argument(
        "--samples",
        nargs="+",
        required=True,
        help="Samples to include in the analysis.",
    )

    parser.add_argument(
        "--output",
        required=True,
        help="Output TSV file.",
    )

    args = parser.parse_args()

    # -----------------------------------------------------------------------
    # Read VCF
    # -----------------------------------------------------------------------

    jasmine_df, _ = helper.read_vcf(
        args.jasmine_vcf
    )

    # -----------------------------------------------------------------------
    # Run analysis
    # -----------------------------------------------------------------------

    result = inheritance_analysis(
        args.samples,
        jasmine_df,
    )

    # -----------------------------------------------------------------------
    # Sort output
    # -----------------------------------------------------------------------

    if not result.empty:

        result = result.sort_values(
            by=[
                "PATIENT",
                "CHR",
                "AVG_START",
                "AVG_END",
            ]
        ).reset_index(drop=True)

    # -----------------------------------------------------------------------
    # Write output
    # -----------------------------------------------------------------------

    output_dir = os.path.dirname(args.output)

    if output_dir:
        os.makedirs(
            output_dir,
            exist_ok=True,
        )

    result.to_csv(
        args.output,
        sep="\t",
        index=False,
    )


if __name__ == "__main__":
    main()
