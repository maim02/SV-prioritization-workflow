import pandas as pd
import argparse
import os
import psycopg2

import warnings

warnings.filterwarnings("ignore")

#### This script extracts for the diseasemodulediscovery gene selection the gene region and the location of its associated cCREs

# call:
# python get_dmd_panel.py --genes <file> --host <host> --output /path/to/file

def get_diseasemodulediscovery_genes(genes_file):
    print("=" * 70)
    print("Load Genes...")
    print("=" * 70)
    df = pd.read_csv(genes_file, sep="\t")

    # remove seed genes
    df = df[df["is_seed"] == 0].copy()

    #df["name"].to_csv(outfile, sep="\t", index=False, header=False)
    return df["name"].to_list()


def get_gene_regions(host, genes, outfile): 
    print("=" * 70)
    print("Get Gene Region...")
    print("=" * 70)

    conn = psycopg2.connect(host)

    placeholders = ",".join(["%s"] * len(genes))

    query = f"""
        SELECT chr, start_pos, end_pos, gene_name, ensg, strand
        FROM public.know_gencode_gene
        WHERE gene_name IN ({placeholders})
    """

    gene_region_df = pd.read_sql_query(query, conn, params=genes)

    conn.close()

    # Report genes that weren't found
    found = set(gene_region_df["gene_name"])
    missing = set(genes) - found

    for gene in sorted(missing):
        print(f"Gene not found: {gene}")

    
    gene_region_df = gene_region_df[["chr", "start_pos", "end_pos", "gene_name", "ensg","strand"]]

    print("Adding window to start and end...")
    
    gene_region_df.loc[
        gene_region_df["strand"] == 1,
        "start_pos"
    ] -= 1000

    gene_region_df.loc[
        gene_region_df["strand"] == -1,
        "end_pos"
    ] += 1000

    gene_region_df["cCRE_class"] = "gene_region"

    print("Write gene region data to file...")
    out_df = gene_region_df[["chr", "start_pos", "end_pos", "gene_name", "cCRE_class"]]
    out_df.to_csv(outfile, sep="\t", header=False, index=False)

    return gene_region_df


def get_cCREs_regions(interaction_file, cCres_file, genes_df, outfile):
    print("=" * 70)
    print("Get Gene Associated cCREs Region...")
    print("=" * 70)

    interaction_df = pd.read_csv(interaction_file, sep='\t', header=None)
    interaction_df.columns = ['cCRE_id', 'gene_id', 'gene_name', 'gene_type', 'assay_type', 'experiment_id', 'biosample', 'score', 'p_value']
    interaction_df["gene_name"] = interaction_df["gene_name"].str.strip()

    # get all gene_ids in gene_df which are not in interaction_df
    genes_in_interaction = set(interaction_df['gene_id'])
    genes_in_genes_df = set(genes_df['ensg'])
    genes_not_in_interaction = genes_in_genes_df - genes_in_interaction
    print(f'Number of panel genes: {len(genes_in_genes_df)}')
    print(f'Number of genes without interaction data: {len(genes_not_in_interaction)}')
    print(f"Missing genes: {genes_df.loc[genes_df['ensg'].isin(genes_not_in_interaction), 'gene_name'].tolist()}")

    filtered_interaction_df = interaction_df[interaction_df['gene_id'].isin(genes_in_genes_df)]
    print(f'Number of interactions with panel genes: {len(filtered_interaction_df)}')

    relevant_ccres = set(filtered_interaction_df['cCRE_id'])

    ccres_gene = filtered_interaction_df[['cCRE_id', 'gene_name']].drop_duplicates().groupby('cCRE_id')['gene_name'].apply(",".join).reset_index()

    ccres_df = pd.read_csv(cCres_file, sep='\t', header=None)
    ccres_df = ccres_df[[0,1,2,3,9]]
    ccres_df.columns = ['chrom', 'start', 'end', 'cCRE_id', 'cCRE_class']
    filtered_blood_df = ccres_df[ccres_df['cCRE_id'].isin(relevant_ccres)]

    # add gene_name to filtered_blood_df by merging with filtered_interaction_df on cCRE_id
    filtered_blood_df = filtered_blood_df.merge(ccres_gene[['cCRE_id', 'gene_name']], on='cCRE_id', how='left')
    # drop duplicate rows
    filtered_blood_df = filtered_blood_df.drop_duplicates()

    # shuffle to match gene_df
    filtered_blood_df = filtered_blood_df[['chrom', 'start', 'end', 'gene_name', 'cCRE_class']]
    print(f'Number of cCREs regions: {len(filtered_interaction_df)}')

    # add gene_region as label for cCREs_class
    subset_genes_df = genes_df[["chr", "start_pos", "end_pos", "gene_name","cCRE_class"]]
    subset_genes_df.columns = ['chrom', 'start', 'end', 'gene_name', 'cCRE_class'] 


    merged_df = pd.concat([filtered_blood_df, subset_genes_df], ignore_index=True)
    
    merged_df = merged_df.sort_values(["chrom", "start", "end"])

    print("=" * 70)
    print("Save to file...")
    print("=" * 70)

    merged_df.to_csv(outfile, sep="\t", header=False, index=False)


if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--genes",
        required=True
    )

    parser.add_argument(
        "--host", 
        required=True
    )
    
    parser.add_argument(
        "--interaction",
        required=True
    )

    parser.add_argument(
        "--ccres",
       required=True
    )
    
    parser.add_argument(
        "--output_region",
        required=True
    )

    parser.add_argument(
        "--output",
        required=True
    )

    args = parser.parse_args()

    out_dir = os.path.dirname(args.output)
    os.makedirs(out_dir, exist_ok=True)

    genes = get_diseasemodulediscovery_genes(args.genes)

    gene_region_df = get_gene_regions(args.host, genes, args.output_region)
    get_cCREs_regions(args.interaction, args.ccres, gene_region_df, args.output)



