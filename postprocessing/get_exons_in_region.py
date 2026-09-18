import pandas as pd
import argparse
import os
import psycopg2

import warnings

warnings.filterwarnings("ignore")

#######
# Get for an (enformer) region all its exons within the region

# call:
# python postprocessing/get_exons_in_region.py --output /dir/file --chr <chromosome> --start <start> --end <end> --host <host>
#######


def get_exons(chr, start, end, host):
    print("=" * 70)
    print("Get Exons in Region...")
    print("=" * 70)

    conn = psycopg2.connect(host)

    query = """
        SELECT
            chr,
            start_pos,
            end_pos,
            split_part(raw->>'gene_id', '.', 1) AS gene_id,
            raw->>'gene_name' AS gene_name,
            pos_range
        FROM public.know_gencode_exon
        WHERE chr = %s
          AND pos_range && int4range(%s, %s, '[)')
        ORDER BY chr, start_pos;
    """

    exons_df = pd.read_sql_query(
        query,
        conn,
        params=(chr, start, end + 1)
    )

    conn.close()

    # Clip exon coordinates to requested region
    exons_df["start_pos"] = exons_df["start_pos"].clip(lower=start)
    exons_df["end_pos"] = exons_df["end_pos"].clip(upper=end)

    return exons_df


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
        "--host", 
        required=True
    )

    parser.add_argument(
        "--output",
        required=True
    )

    args = parser.parse_args()

    out_dir = os.path.dirname(args.output)
    os.makedirs(out_dir, exist_ok=True)

    exons_df = get_exons(args.chr, args.start, args.end, args.host)

    exons_df.to_csv(args.output, sep="\t", index=False)



