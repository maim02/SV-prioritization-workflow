import pandas as pd
import argparse
import numpy as np

def filter_gnomad(variants_df, gnomad_df):
    # subset gnomad_to include only variants that have AF > 0.001 -> common variants
    keys = ['CHROM', 'AVG_START', 'AVG_END', 'GT',"ID"]

    gnomad_common = (
        gnomad_df.loc[gnomad_df['AF'] > 0.001, keys]
        .drop_duplicates()
    )

    selected_variants = variants_df.merge(
        gnomad_common.assign(IN_GNOMAD_COMMON=True),
        on=keys,
        how='left'
    )

    return selected_variants[selected_variants['IN_GNOMAD_COMMON'].isna()].drop(columns='IN_GNOMAD_COMMON')


def filter_inhouse(variants_df, inhouse_df):
    keys = ["CHROM", "AVG_START", "AVG_END", "GT", "ID"]

    return (
        variants_df.merge(inhouse_df[keys].drop_duplicates(),
                on=keys,
                how="left",
                indicator=True)
        .query('_merge == "left_only"')
        .drop(columns="_merge")
    )


if __name__=="__main__":
    parser = argparse.ArgumentParser(
            description="Population Frequency Filtering"
        )
    
    parser.add_argument("--variants_file",required=True)
    parser.add_argument("--inhouse_overlap",required=True)
    parser.add_argument("--gnomad_overlap",required=True)
    parser.add_argument("--out_file",required=True)

    args = parser.parse_args()

    variants_df = pd.read_csv(args.variants_file, sep="\t", header=None, names=["CHROM","AVG_START","AVG_END","GT","ID","SVTYPE","AVG_LEN","SAMPLES_HET","SAMPLES_HOM"])

    print(len(variants_df["ID"].unique()))

    gnomad_df = pd.read_csv(args.gnomad_overlap, sep="\t", header=None, names=["CHROM","AVG_START","AVG_END","GT","ID","SVTYPE","AVG_LEN","SAMPLES_HET","SAMPLES_HOM", "CHROM_g", "START_g","END_g","ID_g","SVTYPE_g","AF", "FREQ_HOMALT", "FREQ_HET"])

    inhouse_df = pd.read_csv(args.inhouse_overlap, sep="\t", header=None, names=["CHROM","AVG_START","AVG_END","GT","ID","SVTYPE","AVG_LEN","SAMPLES_HET","SAMPLES_HOM", "CHROM_i","START_i","END_i","SVTYPE_i","Het","Hom","TEN","AF"])

    selected_variants = filter_inhouse(variants_df, inhouse_df)

    print(len(selected_variants["ID"].unique()))

    selected_variants = filter_gnomad(selected_variants, gnomad_df)

    print(len(selected_variants["ID"].unique()))

    selected_variants = selected_variants.replace("", np.nan)

    selected_variants.to_csv(args.out_file, sep="\t", header=False, index=False, na_rep="None")

