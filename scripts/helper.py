import gzip
import pandas as pd
import subprocess
import os


def read_vcf(vcf:str):
    header_lines = []
    if vcf.endswith(".gz"):
        with gzip.open(vcf, "rt") as f:
            for line in f:
                if line.startswith("#CHROM"):
                    header = line.strip().lstrip("#").split("\t")
                    break
                elif line.startswith("#"):
                    header_lines.append(line)
    else:
        with open(vcf) as f:
            for line in f:
                if line.startswith("#CHROM"):
                    header = line.strip().lstrip("#").split("\t")
                    break
                elif line.startswith("#"):
                    header_lines.append(line)

    df = pd.read_csv(
            vcf,
            sep="\t",
            comment="#",
            names=header
        )

    return df, header_lines


def parse_jasmine_info(info):
    d = {}

    for item in info.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            d[key] = value

    return pd.Series({
        "AVG_START": int(float(d["AVG_START"])) if "AVG_START" in d else None,
        "AVG_END": int(float(d["AVG_END"])) if "AVG_END" in d else None,
        "AVG_LEN": abs(int(float(d["AVG_LEN"]))) if "AVG_LEN" in d else None,
        "SVTYPE": d.get("SVTYPE"),
        "ID_LIST": d.get("IDLIST"),
        "SUPP_VEC":d.get("SUPP_VEC")
    })


def get_jasmine_sample_order(info_columns:list, jasmine_vcf_df:pd.DataFrame):
    """
    Jasmine returns sample with a numbered prefix, remove the leading number and return a list of the ordered sample names
    """
    sample_nr_columns = [col for col in jasmine_vcf_df.columns if col not in info_columns]
    return ["_".join(sample.split("_")[1:]) for sample in sample_nr_columns]


def perform_bed_overlap(bed1:str, bed2:str, output_file_name:str, min_overlap=1E-9):
    """
    Perform bedtools intersect between two bed files and write the output to a file.
    """
    try:
        output_dir = os.path.dirname(output_file_name)
        os.makedirs(output_dir, exist_ok=True)
        with open(output_file_name, "w") as out:
            subprocess.run(
                [
                    "bedtools",
                    "intersect",
                    "-a",
                    bed1,
                    "-b",
                    bed2,
                    "-f",
                    str(min_overlap),
                    "-r",
                    "-wa",
                    "-wb",
                ],
                stdout=out,
                text=True,
                check=True,
            )
    except subprocess.CalledProcessError as e:
        #print("ERROR:")
        print(f"Occured in function perform_bed_overlap(bed1={bed1}, bed2={bed2}, output_file_name={output_file_name}")
        print(f"Exit code {e.returncode}")
        raise