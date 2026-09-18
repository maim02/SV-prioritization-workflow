import os
import pandas as pd
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--eqtl",
        required=True
    )

    parser.add_argument(
        "--interaction",
        required=True
    )

    parser.add_argument(
        "--output",
        required=True
    )

    args = parser.parse_args()

    output_dir = os.path.dirname(args.output)
    os.makedirs(output_dir, exist_ok=True)

    eqtl_df = pd.read_csv(args.eqtl, header=None, sep="\t")
    eqtl_df.columns = ["cCRE_id","gene_id","gene_name","gene_type","variant_id","Source","Tissue","Slope","P-value"]
    eqtl_df['gene_name'] = eqtl_df['gene_name'].astype(str).str.strip()

    interaction_merge = pd.read_csv(args.interaction, sep="\t", header=None, dtype={8: str})
    interaction_merge.columns = ['CHROM', 'AVG_START', 'AVG_END', 'GT', 'ID', 'SVTYPE', 'AVG_LEN', 'SAMPLES_HET', 'SAMPLES_HOM',
        'cCRE_CHROM', 'cCRE_START', 'cCRE_END', 'cCRE_ID', 'cCRE_SCORE',
        'cCRE_STRAND', 'cCRE_THICK_START', 'cCRE_THICK_END', 'cCRE_RGB',
        'cCRE_CLASS', 'cCRE_BIOSAMPLE', 'gene_id', 'gene_name']

    filtered_eqtl_sv = pd.merge(interaction_merge, eqtl_df, left_on=["cCRE_ID","gene_id", "gene_name"], right_on=["cCRE_id","gene_id", "gene_name"],how="inner")
    filtered_eqtl_sv.to_csv(args.output, sep="\t", index=False)
