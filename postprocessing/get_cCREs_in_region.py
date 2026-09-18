import pandas as pd
import argparse
import os
import psycopg2

import warnings

warnings.filterwarnings("ignore")

#######
# Get for an (enformer) region all its ccres within this region

# call:
# python postprocessing/get_cCREs_in_region.py --output /dir/file --chr <chromosome> --start <start> --end <end> 
#######

def get_cCREs(chr, start, end, cCREs_file, interaction_file):
    ccres_df = pd.read_csv(cCREs_file, sep='\t', header=None)
    ccres_df = ccres_df[[0, 1, 2, 3, 9]]
    ccres_df.columns = ['chrom', 'start', 'end', 'cCRE_id', 'cCRE_class']

    # Subset to same chromosome and any overlap with [start, end)
    overlapping = ccres_df[
        (ccres_df['chrom'] == chr) &
        (ccres_df['start'] < end) &
        (ccres_df['end'] > start)
    ].copy()

    # Clip start/end to the query region bounds
    overlapping['start'] = overlapping['start'].clip(lower=start)
    overlapping['end'] = overlapping['end'].clip(upper=end)

    # get gene interaction
    interaction_df = pd.read_csv(interaction_file, sep='\t', header=None)
    interaction_df.columns = ['cCRE_id', 'gene_id', 'gene_name', 'gene_type', 'assay_type', 'experiment_id', 'biosample', 'score', 'p_value']
    interaction_df["gene_name"] = interaction_df["gene_name"].str.strip()

    merged_df = overlapping.merge(interaction_df[['cCRE_id', 'gene_name']], on='cCRE_id', how='left')
    # drop duplicate rows
    merged_df = merged_df.drop_duplicates()

    return merged_df


if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--chr",
        required=True
    )

    parser.add_argument(
        "--start",
        required=True,
        type=int
    )

    parser.add_argument(
        "--end",
        required=True,
        type=int
    )

    parser.add_argument(
            "--interaction",
            default='/archive/emedgene_data/resources/screen_files/V4-hg38.Gene-Links.3D-Chromatin.txt'
        )
    
    parser.add_argument(
        "--ccres",
        default='/archive/emedgene_data/resources/screen_files/blood.all.cCREs.bed'
    )

    parser.add_argument(
        "--output",
        required=True
    )

    args = parser.parse_args()

    out_dir = os.path.dirname(args.output)
    os.makedirs(out_dir, exist_ok=True)

    cCres_df = get_cCREs(args.chr, args.start, args.end, args.ccres, args.interaction)

    cCres_df.to_csv(args.output, sep="\t", index=False)



