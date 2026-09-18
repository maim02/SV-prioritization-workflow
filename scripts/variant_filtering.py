import pandas as pd

def get_affected_variants(samples:list, jasmine_vcf_df:pd.DataFrame):
    """
    Filter out all variants that are only in unaffected. Keeping variants where at least one affected sample carries the variant and all unaffected samples do not.
    args:
        samples (list): of all samples to analyze
        jasmine_vcf_df (pd.DataFrame): df with jasmine merged vcf data (columns: CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT [n_SAMPLES])
    return: subsetted jasmine_vcf_df containing only variants (rows) which are in diseased and not unaffected.
    """
    # remove all variants that are in unaffected only
    unaffected = []
    affected = []
    for s in samples:
        if "pa" in s:
            affected.append(s)
        else:
            unaffected.append(s)

    # keep variants where at least one affected sample carries the variant and all unaffected samples do not
    mask = (
        jasmine_vcf_df[affected].map(_has_variant).any(axis=1)
        &
        ~jasmine_vcf_df[unaffected].map(_has_variant).any(axis=1)
    )

    return jasmine_vcf_df[mask]

def _has_variant(sample_data):
    """
    Check if vcf output for sample variant shows that it is there.
    args:
        sample_data (str): of format "GT:IS:OT:DV:DR"
    return: boolean -> True if sample has variant else False
    """
    return pd.notna(sample_data) and not sample_data.startswith(("0/0", "./."))


