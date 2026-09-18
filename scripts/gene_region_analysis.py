import pandas as pd
import argparse
from pathlib import Path

# get
# sv_df            # SVs
# exons_df         # extracted exons
# cds_df           # extracted CDS
# transcripts_df   # extracted transcripts
# ccre_df          # relevant cCREs
# interaction_df   # cCRE -> gene interactions

parser = argparse.ArgumentParser()

parser.add_argument(
    "--variants",
    required=True
)


parser.add_argument(
    "--exons",
    required=True
)

parser.add_argument(
    "--transcripts",
    required=True
)


parser.add_argument(
    "--cds",
    required=True
)

# parser.add_argument(
#     "--map",
#     required=True
# )

parser.add_argument(
    "--out",
    required=True
)

args = parser.parse_args()

# read files 
sv_df = pd.read_csv(args.variants, sep="\t", header = None, names=["CHROM","AVG_START","AVG_END","GT","ID","SVTYPE","AVG_LEN","SAMPLES_HET","SAMPLES_HOM", "CNV_SUPP", "CHROM_g", "START_g","END_g", "GENE_NAME", "LABEL"])
exons_df = pd.read_csv(args.exons, sep="\t")
exons_df["gene_id"] = exons_df["gene_id"].str.split(".").str[0] # drop version number
cds_df = pd.read_csv(args.cds, sep="\t")
cds_df["gene_id"] = cds_df["gene_id"].str.split(".").str[0] # drop version number
transcripts_df = pd.read_csv(args.transcripts, sep="\t")
gene_map_df = pd.read_csv(f"{Path(args.transcripts).parent}/gene_ensembl_map.tsv", sep="\t")


results = []

for _, sv in sv_df.iterrows(): # per variant

    chrom = sv.CHROM # var chrom
    start = sv.AVG_START # var start
    end = sv.AVG_END # var end
    gene = sv.GENE_NAME # affected gene name
    SPLICE_WINDOW = 20

    ensg = gene_map_df[gene_map_df["gene_name"] == gene]["ensg"].iloc[0]
    ############################################################
    # transcripts that overlap with variant
    ############################################################

    tx = transcripts_df[
        (transcripts_df.ensg == ensg) &
        (transcripts_df.chr == chrom) &
        (transcripts_df.start_pos <= end) &
        (transcripts_df.end_pos >= start)
    ]

    ############################################################
    # MANE / canonical
    ############################################################

    mane = tx.is_mane_select.any()
    
    mane_clinical = tx.is_mane_plus_clinical.any()

    canonical = tx.is_ensembl_canonical.any()

    ############################################################
    # exons
    ############################################################

    exon_hits = exons_df[
        (exons_df.gene_id == ensg) &
        (exons_df.chr == chrom) &
        (exons_df.start_pos <= end) &
        (exons_df.end_pos >= start)
    ]

    exon_numbers = sorted(exon_hits.exon_number.unique())

    ############################################################
    # is the variant within +-20 bp of the exon
    ############################################################
    gene_exons = exons_df[
        (exons_df["gene_id"] == ensg) &
        (exons_df["chr"] == chrom)
    ]

    splice_hits_left = []
    splice_hits_right = []

    for _, exon in gene_exons.iterrows():

        # Intronic region left of exon
        left_start = exon["start_pos"] - SPLICE_WINDOW
        left_end = exon["start_pos"] - 1

        # Intronic region right of exon
        right_start = exon["end_pos"] + 1
        right_end = exon["end_pos"] + SPLICE_WINDOW

        left_overlap = (
            start <= left_end and
            end >= left_start
        )

        right_overlap = (
            start <= right_end and
            end >= right_start
        )

        if left_overlap:
            splice_hits_left.append(exon["exon_number"])
        if right_overlap:
            splice_hits_right.append(exon["exon_number"])
    

    ############################################################
    # CDS
    ############################################################

    cds_hits = cds_df[
        (cds_df.gene_id == ensg) &
        (cds_df.chr == chrom) &
        (cds_df.start_pos <= end) &
        (cds_df.end_pos >= start)
    ]

    cds_bp = 0

    for _, cds in cds_hits.iterrows():

        ov_start = max(start, cds.start_pos)
        ov_end = min(end, cds.end_pos)

        cds_bp += ov_end - ov_start + 1



    results.append({
        "SV_CHROM": sv.CHROM,
        "AVG_START": sv.AVG_START,
        "AVG_END": sv.AVG_END,
        "GT": sv.GT,
        "ID": sv.ID,
        "SVTYPE": sv.SVTYPE,
        "AVG_LEN": sv.AVG_LEN,
        "SAMPLES_HET": sv.SAMPLES_HET,
        "SAMPLES_HOM": sv.SAMPLES_HOM,
        "GENE_NAME": sv.GENE_NAME,
        "GENE_ID": ensg,
        "CNV_SUPP": sv.CNV_SUPP,
        "NUMBER_AFFECTED_TRANSCRIPTS": len(tx),
        "MANE": mane,
        "MANE_CLINICAL": mane_clinical,
        "CANONICAL": canonical,
        #"Exons_affected": ",".join(map(str, exon_numbers)),
        "NUMBER_AFFECTED_EXONS": len(exon_hits),
        "AFFECTED_EXON_LEFT": splice_hits_left,
        "AFFECTED_EXON_RIGHT": splice_hits_right,
        "AFFECTED_CDS_bp": cds_bp,
        #"Overlaps_cCRE": len(ccre_hits) > 0,
    })

summary = pd.DataFrame(results)

summary["SPLICE_RELEVANT"] = (
    len(summary["AFFECTED_EXON_LEFT"]) != 0 |
    len(summary["AFFECTED_EXON_RIGHT"]) != 0
)

summary = summary.sort_values(
    by=[
        "NUMBER_AFFECTED_EXONS",
        "SPLICE_RELEVANT",
        "AFFECTED_CDS_bp",
        "NUMBER_AFFECTED_TRANSCRIPTS",
    ],
    ascending=[False, False, False, False]
)

Path(args.out).parent.mkdir(parents=True, exist_ok=True)
summary.to_csv(args.out, sep="\t", index=False)
