#!/usr/bin/env python3

import pandas as pd
import psycopg2
import argparse
from pathlib import Path

# ----------------------------
# Configuration
# ----------------------------

parser = argparse.ArgumentParser()
    
parser.add_argument(
    "--host",
    required=True
)

parser.add_argument(
    "--input",
    required=True
)

parser.add_argument(
    "--output_exons",
    required=True
)

parser.add_argument(
    "--output_transcripts",
    required=True
)

parser.add_argument(
    "--output_cds",
    required=True
)

args = parser.parse_args()

PG_DSN, input, output_exons, output_transcripts, output_cds = args.host, args.input, args.output_exons, args.output_transcripts, args.output_cds

Path(output_exons).parent.mkdir(parents=True, exist_ok=True)
Path(output_transcripts).parent.mkdir(parents=True, exist_ok=True)
Path(output_cds).parent.mkdir(parents=True, exist_ok=True)

# ----------------------------
# Read gene IDs
# ----------------------------

input_df = pd.read_csv(input, sep="\t", header=None, names=["CHROM","AVG_START","AVG_END","GT","ID","SVTYPE","AVG_LEN","SAMPLES_HET","SAMPLES_HOM", "CNV_SUPP", "CHROM_g", "START_g","END_g", "GENE_NAME", "LABEL"])

genes = input_df["GENE_NAME"].unique().tolist()

print(f"Loaded {len(genes)} gene IDs.")

# ----------------------------
# Connect to PostgreSQL
# ----------------------------

conn = psycopg2.connect(PG_DSN)
cur = conn.cursor()

# ----------------------------
# Map gene names -> Ensembl IDs
# ----------------------------

print("Mapping gene names to Ensembl gene IDs...")

query_gene = """
SELECT DISTINCT
    ensg,
    gene_name
FROM know_gencode_gene
WHERE gene_name = ANY(%s);
"""

cur.execute(query_gene, (genes,))

rows = cur.fetchall()

gene_map = pd.DataFrame(
    rows,
    columns=["ensg", "gene_name"]
)

print(f"Found {len(gene_map)} gene name -> Ensembl ID mappings.")

# Remove version from Ensembl IDs if present
gene_map["ensg"] = (
    gene_map["ensg"]
    .astype(str)
    .str.split(".", n=1)
    .str[0]
)

gene_map.to_csv(f"{Path(output_transcripts).parent}/gene_ensembl_map.tsv", sep="\t", index=False)

# ----------------------------
# Check for missing gene names
# ----------------------------

found_gene_names = set(gene_map["gene_name"].dropna())

missing_gene_names = set(genes) - found_gene_names

if missing_gene_names:
    print(
        f"Warning: {len(missing_gene_names)} input gene names "
        f"were not found in know_gencode_gene:"
    )
    print(", ".join(sorted(missing_gene_names)))

# Only use successfully mapped Ensembl IDs
ensg_ids = gene_map["ensg"].dropna().unique().tolist()

print(f"Using {len(ensg_ids)} unique Ensembl gene IDs for downstream extraction.")

# ----------------------------
# Query exons & save to file
# ----------------------------
TABLE_NAME = "know_gencode_exon"

query = f"""
SELECT
    enst,
    exon_number,
    exon_id,
    chr,
    start_pos,
    end_pos,
    strand,
    raw->>'gene_id' AS gene_id,
    raw->>'gene_name' AS gene_name
FROM {TABLE_NAME}
WHERE split_part(raw->>'gene_id', '.', 1) = ANY(%s)
ORDER BY raw->>'gene_id', enst, exon_number;
"""

cur.execute(query, (ensg_ids,))

rows = cur.fetchall()

columns = [desc[0] for desc in cur.description]

df = pd.DataFrame(rows, columns=columns)

print(f"Found {len(df)} exons.")

df.to_csv(output_exons, sep="\t", index=False)

# check if df has all ensembl IDs from the input file
missing_genes = set(ensg_ids) - set(df['gene_id'].str.split('.').str[0])

if missing_genes:
    print(f"Warning: Missing IDs in exon extraction: {', '.join(sorted(missing_genes))}")

print(f"Saved to {output_exons}")

# ----------------------------
# Query cds & save to file
# ----------------------------
TABLE_NAME = "know_gencode_cds"

query_cds = f"""
SELECT
    enst,
    genome_build,
    cds_rank,
    chr,
    start_pos,
    end_pos,
    strand,
    phase,
    release_id,
    raw->>'gene_id' AS gene_id,
    raw->>'gene_name' AS gene_name,
    raw->>'transcript_id' AS transcript_id,
    raw->>'transcript_name' AS transcript_name,
    raw->>'transcript_type' AS transcript_type,
    raw->>'protein_id' AS protein_id
FROM {TABLE_NAME}
WHERE split_part(raw->>'gene_id', '.', 1) = ANY(%s)
ORDER BY
    split_part(raw->>'gene_id', '.', 1),
    enst,
    cds_rank;
"""

cur.execute(query_cds, (ensg_ids,))

rows = cur.fetchall()

columns = [desc[0] for desc in cur.description]

df = pd.DataFrame(rows, columns=columns)

print(f"Found {len(df)} cds.")

df.to_csv(output_cds, sep="\t", index=False)

# check if df has all ensembl IDs from the input file
missing_genes = set(ensg_ids) - set(df['gene_id'].str.split('.').str[0])

if missing_genes:
    print(f"Warning: Missing IDs in cds extraction: {', '.join(sorted(missing_genes))}")

print(f"Saved to {output_cds}")

# ----------------------------
# Query transcript & save to file
# ----------------------------
TABLE_NAME = "know_gencode_transcript"

query_transcript = f"""
SELECT
    enst,
    ensg,
    genome_build,
    chr,
    start_pos,
    end_pos,
    strand,
    transcript_name,
    transcript_type,
    tags,
    tsl,
    is_mane_select,
    is_mane_plus_clinical,
    is_ensembl_canonical,
    is_gencode_primary,
    release_id
FROM {TABLE_NAME}
WHERE split_part(ensg, '.', 1) = ANY(%s)
ORDER BY
    split_part(ensg, '.', 1),
    is_mane_select DESC,
    is_mane_plus_clinical DESC,
    is_ensembl_canonical DESC,
    transcript_name;
"""
cur.execute(query_transcript, (ensg_ids,))

rows = cur.fetchall()

columns = [desc[0] for desc in cur.description]

df = pd.DataFrame(rows, columns=columns)

print(f"Found {len(df)} transcripts.")

df.to_csv(output_transcripts, sep="\t", index=False)

# check if df has all ensembl gene IDs from the input file
missing_genes = set(ensg_ids) - set(df['ensg'].str.split('.').str[0])

if missing_genes:
    print(f"Warning: Missing gene IDs in transcript extraction: {', '.join(sorted(missing_genes))}")

print(f"Saved to {output_transcripts}")