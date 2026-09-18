from pathlib import Path

configfile: "config/config.yaml"


###############################################################################
# Get samples
###############################################################################

samples = []
cnv_samples = []

for vcf in Path(config["data_dir"]).glob("*/*/*.vcf.gz"):
    if vcf.name.endswith(config["sv_suffix"]):
        samples.append(
            {
                "family": vcf.parent.parent.name,
                "sample": vcf.parent.name,
                "vcf": str(vcf),
            }
        )
    elif vcf.name.endswith(config["cnv_suffix"]): # _sv
        cnv_samples.append(
            {
                "family": vcf.parent.parent.name,
                "sample": vcf.parent.name,
                "vcf": str(vcf),
            }
        )

VCFS = {
    (s["family"], s["sample"]): s["vcf"]
    for s in samples
}

CNV_VCFS = {
    (s["family"], s["sample"]): s["vcf"]
    for s in cnv_samples
}


###############################################################################
# Parsing Helpers
###############################################################################

def input_sv(wc):
    return VCFS[(wc.family, wc.sample)]


def input_cnv(wc):
    return CNV_VCFS[(wc.family, wc.sample)]



###############################################################################
# Final target
###############################################################################

rule all:
    input:
        # Merged SV/CNV VCFs
        "output/sv_jasmine/all_samples_jasmine.vcf",
        "output/cnv_jasmine/all_samples_jasmine.vcf",

        # Inheritance Output
        "output/inheritance/inheritance_affected.tsv",

        # Final SV prioritization
        "output/prioritized/prioritized_filtered_sv.tsv",

        # Final gene-region result
        "output/gene_region/panel/gene_region_summary.bed",


###############################################################################
# Filter SVs & CNVs
###############################################################################

rule filter_sv:
    input:
        vcf=input_sv

    output:
        "output/sv_filtered/{family}/{sample}_filtered.vcf"

    params:
        filterexp=config["filter_expression"]

    threads: 4
    conda: 
        "sv_env"

    shell:
        """
        mkdir -p $(dirname {output})

        bcftools view \
            -i '{params.filterexp}' \
            {input.vcf} \
            -Ov \
            -o {output}
        """


rule filter_cnv:
    input:
        vcf=input_cnv

    output:
        "output/cnv_filtered/{family}/{sample}_filtered.vcf"

    params:
        filterexp=config["filter_expression"]

    threads: 4
    conda: 
        "sv_env"

    shell:
        """
        mkdir -p $(dirname {output})

        bcftools view \
            -i '{params.filterexp}' \
            {input.vcf} \
            -Ov \
            -o {output}
        """ 

###############################################################################
# Perform Jasmine merge of all samples
###############################################################################

rule jasmine_sv:
    input:
        [
            f"output/sv_filtered/{s['family']}/{s['sample']}_filtered.vcf"
            for s in samples
        ]

    output:
        "output/sv_jasmine/all_samples_jasmine.vcf",

    params:
        jasmine="jasmine",
    conda: 
        "jasmine"

    shell:
        """
        python scripts/jasmine_analysis.py \
            --input_vcfs {input} \
            --output_vcf {output} \
            --jasmine_path {params.jasmine} \
        """

rule jasmine_cnv:
    input:
        [
            f"output/cnv_filtered/{s['family']}/{s['sample']}_filtered.vcf"
            for s in cnv_samples
        ]

    output:
        "output/cnv_jasmine/all_samples_jasmine.vcf",

    params:
        jasmine="jasmine",

    conda: 
        "jasmine"

    shell:
        """
        python scripts/jasmine_analysis.py \
            --input_vcfs {input} \
            --output_vcf {output} \
            --jasmine_path {params.jasmine} \
        """

###############################################################################
# Inheritance analysis
###############################################################################

rule inheritance_analysis:
    input:
        "output/sv_jasmine/all_samples_jasmine.vcf"

    output:
        "output/inheritance/inheritance_affected.tsv"

    params:
        samples=[
            s["sample"] for s in samples
        ],
        script="scripts/inheritance_analysis.py"

    shell:
        """
        python {params.script} \
            --jasmine_vcf {input} \
            --samples {params.samples} \
            --output {output}
        """

###############################################################################
# Cohort and GT based filtering
###############################################################################

rule cohort_gt_filtering_sv:
    input:
        "output/sv_jasmine/all_samples_jasmine.vcf",

    output:
        "output/sv_jasmine/affected_sv_gtFiltered.bed"

    params:
        samples=[
            s["sample"] for s in samples
        ]
    conda: 
        "sv_env"

    shell:
        """
        python scripts/cohort_gt_filtering.py \
            --input_vcf {input} \
            --out_bed {output} \
            --samples {params.samples}
        """

rule cohort_gt_filtering_cnv:
    input:
        "output/cnv_jasmine/all_samples_jasmine.vcf"

    output:
        "output/cnv_jasmine/affected_cnv_gtFiltered.bed"

    params:
        samples=[
            s["sample"] for s in cnv_samples
        ]
    conda: 
        "sv_env"

    shell:
        """
        python scripts/cohort_gt_filtering.py \
            --input_vcf {input} \
            --out_bed {output} \
            --samples {params.samples}
        """

###############################################################################
# GnomAD overlap analysis
###############################################################################

rule gnomad_overlap:
    input:
        jasmine_bed="output/sv_jasmine/affected_sv_gtFiltered.bed",
        
    output:
        "output/population_freq_filtered/jasmine_gnomad_overlap_70perc.bed"

    params:
        script="scripts/bed_overlap.py",
        gnomad_bed=config["gnomad"]
    conda: 
        "sv_env"
    shell:
        """
        python {params.script} \
            --bed1 {input.jasmine_bed} \
            --bed2 {params.gnomad_bed}\
            --output {output} \
            --min_overlap 0.7
        """

###############################################################################
# inhouse overlap analysis
###############################################################################

rule inhouse_freq:
    input:
        "output/sv_jasmine/affected_sv_gtFiltered.bed"

    output:
        "output/population_freq_filtered/jasmine_inhouse_overlap_70perc.bed"

    params:
        script="scripts/inhouse_freq.py",
        inhouse_file=config["inhouse"]

    conda: 
        "sv_env"

    shell:
        """
        python {params.script} \
            --variants_file {input} \
            --inhouse_file {params.inhouse_file} \
            --out_file {output}
        """

###############################################################################
# Population-Freq Filtering
###############################################################################

rule variant_filtering:
    input:
        variants_file="output/sv_jasmine/affected_sv_gtFiltered.bed",
        inhouse_overlap="output/population_freq_filtered/jasmine_inhouse_overlap_70perc.bed",
        gnomad_overlap="output/population_freq_filtered/jasmine_gnomad_overlap_70perc.bed",

    output:
        "output/population_freq_filtered/pop_freq_filtered_sv.bed"

    params:
        script="scripts/population_freq_filter.py"

    shell:
        """
        python {params.script} \
            --variants_file {input.variants_file} \
            --inhouse_overlap {input.inhouse_overlap} \
            --gnomad_overlap {input.gnomad_overlap} \
            --out_file {output}
        """

###############################################################################
# annotate CNV-SV support
###############################################################################

rule cnv_sv_support:
    input:
        sv="output/population_freq_filtered/pop_freq_filtered_sv.bed",
        cnv="output/cnv_jasmine/affected_cnv_gtFiltered.bed"

    output:
        "output/cnv_support/filtered_sv.bed"

    params:
        script="scripts/cnv_support.py"
    conda: 
        "sv_env"

    shell:
        """
        python {params.script} \
            --sv_bed {input.sv} \
            --cnv_bed {input.cnv} \
            --min_overlap 0.7 \
            --output {output}
        """


###############################################################################
# get SCN related SV
###############################################################################
rule scn_extraction:
    input:
        "output/cnv_support/filtered_sv.bed"

    output:
        "output/panel_associated/SV_SCN_associated.bed"

    params:
        script="scripts/bed_overlap.py",
        panel_genes=config["regulatory_region_gene"]
        min_overlap=1E-9
    conda: 
        "sv_env"

    shell:
        """
        python {params.script} \
            --bed1 {input} \
            --bed2 {params.panel_genes} \
            --output {output} \
            --min_overlap {params.min_overlap} \
            
        """

###############################################################################
# SCREEN analysis: cCREs, eQTL, nearest gene data on SV
###############################################################################

rule screen_analysis:
    input:
        sv_bed="output/cnv_support/filtered_sv.bed",
        panel_associated="output/panel_associated/SV_SCN_associated.bed"

    output:
        "output/screen_analysis/screen_merge_sv.tsv"

    params:
        ccre_file="resources/screen_files/blood.all.cCREs.bed",
        interaction_file="resources/screen_files/V4-hg38.Gene-Links.3D-Chromatin.txt",
        eqtl_file="resources/screen_files/V4-hg38.Gene-Links.eQTLs.txt",
        #nearest_genes_file="resources/screen_files/GRCh38-Closest-Genes-PC.tsv",
        script="scripts/screen_analysis.py"
    conda: 
        "sv_env"
    shell:
        """
        python {params.script} \
            --variant_bed {input.sv_bed} \
            --ccres {params.ccre_file} \
            --interactions {params.interaction_file} \
            --eqtl {params.eqtl_file} \
            --panel_file {input.panel_associated} \
            --output {output} \
        """##--nearest_genes {params.nearest_genes_file} \

###############################################################################
# variants in gene region extraction
###############################################################################

rule gene_region_extraction:
    input:
        "output/cnv_support/filtered_sv.bed"
        

    output:
        "output/gene_region/panel/sv_in_gene_region.bed"

    params:
        genes_file=config["gene_region"],
        script="scripts/bed_overlap.py",
        min_overlap=1E-9
    conda: 
        "sv_env"

    shell:
        """
        python {params.script} \
            --bed1 {input} \
            --bed2 {params.genes_file}\
            --output {output} \
            --min_overlap {params.min_overlap}
        """

###############################################################################
# extract exon, cds, transcripts
###############################################################################

rule exon_cds_transcript_extraction:
    input:
        "output/gene_region/panel/sv_in_gene_region.bed"
        

    output:
        exons="output/gene_region/panel/extracted_exons.bed",
        cds="output/gene_region/panel/extracted_cds.bed",
        transcripts="output/gene_region/panel/extracted_transcripts.bed",

    params:
        host=config["db_host"],
        script="/archive/emedgene_data/variant_analysis/scripts/cds_exon_transcript_extraction.py",
    conda: 
        "sv_env"

    shell:
        """
        python {params.script} \
            --host {params.host} \
            --input {input}\
            --output_exons {output.exons} \
            --output_transcripts {output.transcripts} \
            --output_cds {output.cds}
        """


###############################################################################
# gene region analyisis
###############################################################################

rule gene_region_analysis:
    input:
        variants="output/gene_region/panel/sv_in_gene_region.bed",
        exons="output/gene_region/panel/extracted_exons.bed",
        cds="output/gene_region/panel/extracted_cds.bed",
        transcripts="output/gene_region/panel/extracted_transcripts.bed",
        

    output:
        "output/gene_region/panel/gene_region_summary.bed"

    params:
        script="/archive/emedgene_data/variant_analysis/scripts/gene_region_analysis.py",
    conda: 
        "sv_env"

    shell:
        """
        python {params.script} \
            --variants {input.variants} \
            --exons {input.exons} \
            --transcripts {input.transcripts} \
            --cds {input.cds} \
            --out {output}
        """

###############################################################################
# prioritization
###############################################################################
rule screen_prioritization:
    input:
        "output/screen_analysis/screen_merge_sv.tsv"

    output:
        out_file="output/prioritized/prioritized_filtered_sv.tsv",

    params:
        script="scripts/prioritization.py"
    conda: 
        "sv_env"
    shell:
        """
        python {params.script} \
            --variants_file {input} \
            --output_file {output.out_file} \
        """
