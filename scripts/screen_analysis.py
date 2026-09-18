import argparse
import pandas as pd
import helper

# -------------------------------------------------------------------------
# Configuration
# -------------------------------------------------------------------------
KEY_COL = [
    "CHROM",
    "AVG_START",
    "AVG_END",
    "GT"
]

SV_COL=[
    "CHROM", 
    "AVG_START",
    "AVG_END", 
    "GT", 
    "ID",
    "SVTYPE",
    "AVG_LEN", 
    "SAMPLES_HET",
    "SAMPLES_HOM", 
    "CNV_SUPP"
]

INTERACTION_COL = [
    "cCRE_ID",
    "GENE_ID",
    "GENE_NAME",
    "gene_type",
    "assay_type",
    "experiment_id",
    "biosample",
    "score",
    "p_value"
]

EQTL_COL = ["cCRE_ID","GENE_ID","GENE_NAME","gene_type","variant_id","Source","Tissue","Slope","P-value"]

GENES_COL = [
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
    "CHROM_gene",
    "START_gene",
    "END_gene",
    "GENE_NAME",
    "cCRE_CLASS"
]

# NEAREST_GENES_COL = [
#     "CHROM",
#     "START",
#     "END",
#     "cCRE_ID_2",
#     "cCRE_ID",
#     "cCRE_TYPE",
#     "GENE_CHROM",
#     "GENE_START",
#     "GENE_END",
#     "TRANSCRIPT_ID",
#     "GENE_SYMBOL",
#     "STRAND",
#     "GENE_ID"
# ]


# -------------------------------------------------------------------------
# eqtl 
# -------------------------------------------------------------------------

def load_eqtl(eqtl_file:str):
    """
    Load SCREEN eQTL table.
    """
    
    eqtl_df = pd.read_csv(eqtl_file, header=None, sep="\t")
    eqtl_df.columns = EQTL_COL
    eqtl_df['GENE_NAME'] = eqtl_df['GENE_NAME'].astype(str).str.strip()

    return eqtl_df


# -------------------------------------------------------------------------
# cCREs 
# -------------------------------------------------------------------------

def load_relevant_cCREs_overlap(output_file:str, variant_bed:str, cCres_file:str):
    """
    Get the relevant cCREs that overlap with the variants in the variant_bed file.
    """
    helper.perform_bed_overlap(variant_bed, cCres_file, output_file, min_overlap=1E-9) # min. overlap of 1bp to consider a variant overlapping with a cCRE

    # read the overlap to a pandas_df
    sv_ccres_df = pd.read_csv(output_file, sep="\t", header=None, dtype={8: str})
    sv_ccres_df.columns = ["CHROM", "AVG_START","AVG_END", "GT", "ID","SVTYPE","AVG_LEN", "SAMPLES_HET","SAMPLES_HOM", "CNV_SUPP", "cCRE_CHROM", "cCRE_START", "cCRE_END", "cCRE_ID", "cCRE_SCORE", "cCRE_STRAND", "cCRE_THICK_START", "cCRE_THICK_END", "cCRE_RGB", "cCRE_CLASS", "cCRE_BIOSAMPLE"]

    return sv_ccres_df[["CHROM", "AVG_START","AVG_END", "GT", "ID","SVTYPE","AVG_LEN", "SAMPLES_HET","SAMPLES_HOM", "CNV_SUPP", "cCRE_CHROM", "cCRE_START", "cCRE_END", "cCRE_ID", "cCRE_STRAND", "cCRE_CLASS", "cCRE_BIOSAMPLE"]]


def load_cCREs_gene_interaction(interaction_file:str):
    """
    Load relevant interaction info to df.
    """
    interaction_df = pd.read_csv(
            interaction_file,
            sep="\t",
            header=None
        )
    
    interaction_df.columns = INTERACTION_COL
    
    # Filter to only keep relevant columns and drop duplicates
    return (interaction_df[["cCRE_ID", "GENE_ID", "GENE_NAME"]].drop_duplicates())


# -------------------------------------------------------------------------
# panel genes
# -------------------------------------------------------------------------


def add_panel_genes(panel_file: str, df: pd.DataFrame):
    """
    Load panel genes data and add the affected panel gene to df
    if a variant is within a panel gene or its associated cCRE region.
    """
    genes_df = pd.read_csv(panel_file, header=None, sep="\t")

    if len(genes_df.columns) < len(GENES_COL):
        raise ValueError(f"Panel file must contain at least columns: {GENES_COL}")

    genes_df.columns = GENES_COL

    # Clean gene names
    genes_df["GENE_NAME"] = genes_df["GENE_NAME"].str.strip()

    # Map each panel region to its unique genes
    genes_df_keys = (
        genes_df.groupby(KEY_COL)["GENE_NAME"]
        .apply(lambda x: list(dict.fromkeys(x)))
        .to_dict()
    )

    # Add affected panel genes
    df["PANEL_SUPP"] = [
        genes_df_keys.get(tuple(row), [])
        for row in df[KEY_COL].to_numpy()
    ]

    return df


# -------------------------------------------------------------------------
# nearest gene
# -------------------------------------------------------------------------

def add_nearest_gene_info(nearest_gene_file, interaction_merge):
    # load nearest_gene
    nearest_gene_df = pd.read_csv(
                nearest_gene_file,
                sep="\t",
                header=None
            )
        
    nearest_gene_df.columns = NEAREST_GENES_COL

    # strip version of gene_id
    nearest_gene_df["NEAREST_GENE_ID"] = nearest_gene_df["GENE_ID"].str.split(".").str[0]

    # merge to interaction_merge and add strand and gene id of nearest gene
    return interaction_merge.merge(
        nearest_gene_df[["cCRE_ID", "STRAND", "NEAREST_GENE_ID"]].drop_duplicates(),
        on="cCRE_ID",
        how="left"
    ).rename(columns={
        "STRAND": "NEAREST_GENE_STRAND" #,
        #"GENE_ID_y": "NEAREST_GENE_ID"
    })


# -------------------------------------------------------------------------
# RUN
# -------------------------------------------------------------------------

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("--variant_bed", required=True)
    parser.add_argument("--ccres", required=True)
    parser.add_argument("--interactions", required=True)
    parser.add_argument("--eqtl", required=True)
    #parser.add_argument("--nearest_genes", required=True)
    parser.add_argument("--panel_file", required=True)
    parser.add_argument("--output", required=True) # all disease variants and their linked cCREs, eqtl annotated, and nearest gene

    args = parser.parse_args()

    # SV overlap with cCREs
    sv_ccres_df = load_relevant_cCREs_overlap(
        output_file=args.output,
        variant_bed=args.variant_bed,
        cCres_file=args.ccres
    )

    # Load cCRE-gene links
    interaction_df = load_cCREs_gene_interaction(args.interactions)

    interaction_merge = sv_ccres_df.merge(
        interaction_df,
        on="cCRE_ID",
        how="inner"
    )

    # merge gene names and gene id into comma sep entry per variant
    interaction_merge["GENE_NAME"] = (interaction_merge["GENE_NAME"].astype(str).str.strip())

    group_cols = [
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
        "cCRE_CHROM",
        "cCRE_START",
        "cCRE_END",
        "cCRE_ID",
        "cCRE_STRAND",
        "cCRE_CLASS",
        "cCRE_BIOSAMPLE"
    ]

    interaction_merge = (
        interaction_merge
        .groupby(group_cols, dropna=False)
        .agg(
            {
                "GENE_ID": lambda x: ",".join(x.astype(str)),
                "GENE_NAME": lambda x: ",".join(x.astype(str))
            }
        )
        .reset_index()
    )

    # add eqtl info
    eqtl_df = load_eqtl(args.eqtl)
    keys = ["cCRE_ID", "GENE_ID", "GENE_NAME"]
    interaction_merge["eQTL_SUPP"] = (
        interaction_merge[keys].astype(str).agg("_".join, axis=1)
        .isin(
            eqtl_df[keys].astype(str).agg("_".join, axis=1)
        )
    )

    # add if variant associated with panel gene
    interaction_merge = add_panel_genes(args.panel_file, interaction_merge.copy())

    # add nearest gene info
    #interaction_merge = add_nearest_gene_info(args.nearest_genes, interaction_merge.copy())

    interaction_merge.to_csv(
        args.output,
        sep="\t",
        index=False,
        #header=False
    )

    