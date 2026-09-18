import pandas as pd
import helper
import os
import argparse


if __name__=="__main__":
    parser = argparse.ArgumentParser()
    
    parser.add_argument(
        "--sv_bed",
        required=True
    )

    parser.add_argument(
        "--cnv_bed",
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


    helper. perform_bed_overlap(
        bed1=args.sv_bed,
        bed2=args.cnv_bed,
        output_file_name=args.output,
        min_overlap=float(args.min_overlap)
    )

    overlap_df = pd.read_csv(args.output, sep="\t", header=None, names=["CHROM","AVG_START","AVG_END","GT","ID","SVTYPE","AVG_LEN","SAMPLES_HET","SAMPLES_HOM", "CHROM_cnv","AVG_START_cnv","AVG_END_cnv","GT_cnv","ID_cnv","SVTYPE_cnv","AVG_LEN_cnv","SAMPLES_HET_cnv","SAMPLES_HOM_cnv"])

    sv_df = pd.read_csv(args.sv_bed, sep="\t", header=None, names=["CHROM","AVG_START","AVG_END","GT","ID","SVTYPE","AVG_LEN","SAMPLES_HET","SAMPLES_HOM"])

    keys= ["CHROM","AVG_START","AVG_END","GT","ID"]
    cnv_supp_annotated = sv_df.merge(
            overlap_df[keys].assign(CNV_SUPP=True),
            on=keys,
            how='left'
        )

    cnv_supp_annotated["CNV_SUPP"] = (
        cnv_supp_annotated["CNV_SUPP"]
        .fillna(False)
        .astype(bool)
    )
    

    cnv_supp_annotated.to_csv(args.output, sep="\t", header=False, index=False)