import os
import psycopg2
import pandas as pd
import argparse

from get_alternative_panel import get_cCREs_regions

### get hpo based genes in internal database


QUERY = r"""
WITH RECURSIVE hpo_scope(term_id) AS
(
    SELECT 'HP:0001875'::text

    UNION

    SELECT relation.child_term_id
    FROM public.know_term_relation relation
    JOIN hpo_scope parent
      ON parent.term_id = relation.parent_term_id
    WHERE relation.rel_type = 'is_a'
),
current_hpoa_release AS
(
    SELECT release.id
    FROM public.know_release release
    WHERE release.source_name = 'HPOA'
    ORDER BY
        release.release_date DESC NULLS LAST,
        release.imported_at DESC,
        release.id DESC
    LIMIT 1
),
matched_disease_hpo AS
(
    SELECT DISTINCT
        disease_hpo.disease_term_id,
        disease_hpo.hpo_term_id
    FROM public.know_disease_hpo disease_hpo
    JOIN current_hpoa_release
      ON current_hpoa_release.id = disease_hpo.release_id
    JOIN hpo_scope
      ON hpo_scope.term_id = disease_hpo.hpo_term_id
)
SELECT
    evidence.gene_id,
    evidence.approved_gene_symbol AS gene_symbol,
    evidence.ensg,
    hgnc.curie AS hgnc_id,
    evidence.disease_term_id,
    evidence.disease_name,
    jsonb_agg(
        DISTINCT jsonb_build_object(
            'term_id', matched_disease_hpo.hpo_term_id,
            'label', matched_hpo.label
        )
    ) AS matching_hpo_terms,
    COALESCE(
        array_agg(
            DISTINCT evidence.assertion_mode_of_inheritance
        ) FILTER (
            WHERE evidence.assertion_mode_of_inheritance IS NOT NULL
        ),
        ARRAY[]::text[]
    ) AS modes_of_inheritance,
    COALESCE(
        array_agg(
            DISTINCT evidence.mode_of_inheritance_curie
        ) FILTER (
            WHERE evidence.mode_of_inheritance_curie IS NOT NULL
        ),
        ARRAY[]::text[]
    ) AS mode_of_inheritance_curies,
    COALESCE(
        array_agg(
            DISTINCT evidence.classification
        ) FILTER (
            WHERE evidence.classification IS NOT NULL
        ),
        ARRAY[]::text[]
    ) AS classifications,
    array_agg(
        DISTINCT evidence.evidence_source_name
    ) AS evidence_sources,
    jsonb_agg(
        DISTINCT jsonb_strip_nulls(
            jsonb_build_object(
                'source', evidence.evidence_source_name,
                'classification', evidence.classification,
                'classification_scheme', evidence.classification_scheme,
                'moi', evidence.assertion_mode_of_inheritance,
                'moi_curie', evidence.mode_of_inheritance_curie
            )
        )
    ) AS source_evidence,
    count(DISTINCT evidence.assertion_id)::integer AS assertion_count
FROM public.report_gene_disease_evidence_current evidence
JOIN matched_disease_hpo
  ON matched_disease_hpo.disease_term_id = evidence.disease_term_id
LEFT JOIN public.know_term matched_hpo
  ON matched_hpo.term_id = matched_disease_hpo.hpo_term_id
LEFT JOIN public.know_gene gene
  ON gene.id = evidence.gene_id
LEFT JOIN public.know_xref hgnc
  ON hgnc.id = gene.hgnc_xref_id
WHERE evidence.gene_id IS NOT NULL
  AND evidence.disease_term_id IS NOT NULL
  AND evidence.count_as_independent_evidence
GROUP BY
    evidence.gene_id,
    evidence.approved_gene_symbol,
    evidence.ensg,
    hgnc.curie,
    evidence.disease_term_id,
    evidence.disease_name
ORDER BY
    evidence.approved_gene_symbol,
    evidence.disease_name,
    evidence.disease_term_id;
"""

def load_panel_genes(panel_genes_file):
    return 


def extract_panel_genes(panel_genes_file, genes):
    # load panel data
    panel_df = pd.read_csv(panel_genes_file, sep="\t", header=None, names=["chr", "start_pos", "end_pos", "gene_name", "cCRE_class"])

    # get panel genes
    panel_genes = panel_df["gene_name"].dropna().unique().tolist()

    # remove panel genes from all genes
    return [gene for gene in genes if gene not in panel_genes]
     


def get_gene_region(host, panel_genes_file, db_output, data_output):
    
    try:
        # Execute query directly into a DataFrame
        conn = psycopg2.connect(host)
        df = pd.read_sql_query(QUERY, conn)

        print(f"Returned {len(df)} rows")
        print(f"Columns: {list(df.columns)}")
        #print()
        #print(df.to_string(index=False))

        # save the DataFrame
        df.to_csv(db_output, sep="\t", index=False)

        genes = df["gene_symbol"].dropna().unique().tolist()

        genes = extract_panel_genes(panel_genes_file, genes)

        placeholders = ",".join(["%s"] * len(genes))

        query = f"""
            SELECT chr, start_pos, end_pos, gene_name, ensg, strand
            FROM public.know_gencode_gene
            WHERE gene_name IN ({placeholders})
        """

        genes_bed = pd.read_sql_query(query, conn, params=genes)

        conn.close()

        # Report genes that weren't found
        found = set(genes_bed["gene_name"])
        missing = set(genes) - found

        for gene in sorted(missing):
            print(f"Gene not found: {gene}")

        print("Adding window to start and end...")
        
        genes_bed.loc[
            genes_bed["strand"] == 1,
            "start_pos"
        ] -= 1000

        genes_bed.loc[
            genes_bed["strand"] == -1,
            "end_pos"
        ] += 1000
    
        genes_bed["cCRE_class"] = "gene_region"

        # Write BED
        genes_bed[["chr", "start_pos", "end_pos", "gene_name", "cCRE_class"]].drop_duplicates().to_csv(
            data_output,
            sep="\t",
            header=False,
            index=False
        )
        print(f"Wrote {len(genes_bed)} entries to genes.bed")
        return genes_bed
    
    except psycopg2.Error as e:
        print(f"Database error: {e}")
        raise


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
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
        "--genes",
         required=True
    )

    parser.add_argument(
        "--db_output",
        required=True
    )

    parser.add_argument(
        "--output",
        required=True
    )

    parser.add_argument(
        "--region_output",
        required=True
    )

    args = parser.parse_args()

    # create out dirs
    out_dir = os.path.dirname(args.db_output)
    os.makedirs(out_dir, exist_ok=True)

    out_dir = os.path.dirname(args.output)
    os.makedirs(out_dir, exist_ok=True)

    out_dir = os.path.dirname(args.region_output)
    os.makedirs(out_dir, exist_ok=True)

    # get hpo related genes and their gene info
    gene_region_df = get_gene_region(args.host, args.genes, args.db_output, args.output)

    # get for genes cCREs
    get_cCREs_regions(args.interaction, args.ccres, gene_region_df, args.region_output)