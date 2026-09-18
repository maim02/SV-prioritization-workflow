# SV-prioritization-workflow

A Snakemake workflow for the prioritization and analysis of structural variants (SVs) identified from short-read whole-genome sequencing (WGS) data.

## Overview

Structural variants can affect coding as well as non-coding regions of the genome and may therefore contribute to disease phenotypes even when standard variant analysis does not identify a genetic cause.

This repository contains the workflow developed to process and prioritize structural variants from short-read WGS data. The workflow integrates SV data with genomic annotations and regulatory information to identify variants that may be relevant to the phenotype under investigation.

The workflow was developed in the context of the analysis of patients with severe congenital neutropenia (SCN), with a particular focus on structural variants overlapping non-coding regulatory regions.

Note: Currently, this pipeline was developed based on the internal database structure of the Dr. von Hauner Children's Hospital. Manual adaptions might be required for other setups

## Workflow

The workflow is implemented using [Snakemake](https://snakemake.readthedocs.io/) and is organized into several processing stages:

```text
Input data
    │
    ▼
Preprocessing
    │
    ▼
SV annotation and filtering
    │
    ▼
Regulatory annotation
    │
    ▼
SV prioritization
    │
    ▼
Postprocessing
    │
    ▼
Evaluation and results
```

The exact workflow configuration and individual processing steps are defined in the `Snakefile` and the corresponding subdirectories.

## Repository structure

```text
SV-prioritization-workflow/
│
├── config/
│   └── Configuration files and workflow parameters
│
├── preprocessing/
│   └── Preparation and preprocessing of input data
│
├── postprocessing/
│   └── Processing and formatting of workflow results
│
├── evaluation/
│   └── Evaluation of prioritization results
│
├── resources/
│   └── Reference data and external resources
│
├── scripts/
│   └── Supporting scripts used by the workflow
│
├── Snakefile
│   └── Main Snakemake workflow
│
├── environment_jasmine.yml
│   └── Conda environment definition
│
├── environment_sv_env.yml
│   └── Conda environment definition
│
└── LICENSE
```

## Requirements

The workflow requires:

* [Snakemake](https://snakemake.readthedocs.io/)
* [Conda](https://docs.conda.io/) or a compatible environment manager
* The software and reference resources specified by the workflow
* Short-read WGS structural variant data
* Appropriate reference genome and annotation resources

Two Conda environment definitions are provided:

* `environment_jasmine.yml`
* `environment_sv_env.yml`


## Installation

Clone the repository:

```bash
git clone https://github.com/maim02/SV-prioritization-workflow.git
cd SV-prioritization-workflow
```

Create the required Conda environments using the provided environment files:

```bash
conda env create -f environment_jasmine.yml
conda env create -f environment_sv_env.yml
```

Depending on the local setup, additional dependencies or reference resources may need to be configured before running the workflow.

## Configuration

Example workflow parameters and paths are defined in the `config/` directory.

Before running the workflow, make sure that:

1. Input files are available at the expected locations.
2. Reference genome resources are correctly configured.
3. Required annotation and regulatory resources are available.
4. Output directories and paths are writable.
5. Configuration files correspond to the genome assembly used for the analysis.

**Important:** Patient-level data and other hospital-specific data are not included in this repository and must not be added to the public repository.

## Running the workflow

The workflow is executed using Snakemake from the repository root.

A typical execution can be started with:

```bash
snakemake --cores <NUMBER_OF_CORES>
```

For a dry run, use:

```bash
snakemake -n
```

This allows the workflow and its dependencies to be inspected without executing the individual rules.

The exact command and configuration required for a specific analysis depend on the available input data, reference resources, and computational environment.

## Input data

The workflow is designed for structural variant data derived from short-read WGS.

Input data used for the analysis may include:

* Structural variant calls
* Reference genome information
* Genomic annotation data
* Regulatory annotations
* cCRE annotations
* cCRE–gene associations
* eQTL information
* Additional resources required by individual workflow steps

The precise input format and file locations are defined by the workflow configuration and scripts.
Example of possibly relevant SCN genes obtained through HPO- or nf-core/DiseaseDiscoveryModule-based approaches, as well as their location based on GENECODE v49 grch38 human reference genome assembly can be found in resources.

## Regulatory annotation

A major component of the workflow is the prioritization of structural variants based on their potential regulatory relevance.

Variants can be assessed with respect to regulatory genomic elements, including candidate cis-regulatory elements (cCREs). Additional information such as cCRE–gene associations and expression quantitative trait loci (eQTLs) can be used to connect regulatory regions to potentially affected genes.

The workflow therefore enables the investigation of structural variants outside protein-coding regions, which may otherwise be missed by analyses focused primarily on coding variants.
Regulatory element annotaion can be downloaded from the ENCODE project via SCREEN: https://screen.wenglab.org/downloads.

## Output

The workflow generates intermediate and final files corresponding to the different processing stages.

These include results related to:

* Processed structural variants
* Genomic and regulatory annotations
* Candidate variant prioritization
* Postprocessed results
* Evaluation of prioritization results

The exact output files depend on the selected workflow configuration.

## Data availability

This repository contains the workflow code and analysis-related resources that can be publicly shared.

Patient-level data are **not** included in the repository. Such data remain within the hospital and are subject to applicable data protection and privacy requirements.

The exact contents of hospital databases are likewise not publicly available.

## Reproducibility

The workflow is intended to support reproducible analysis of structural variants from short-read WGS data.

Software dependencies are documented through the provided Conda environment files and requirements file. Configuration files are separated from the workflow logic to facilitate adaptation to different datasets and computational environments.

Because patient-level data and hospital-specific databases cannot be publicly released, complete reproduction of patient-specific results is not possible from the public repository alone.

## Reference resources

External genomic resources used by the workflow include resources from the ENCODE project and the Genome Reference Consortium (GRC), where applicable.

cCRE annotations, cCRE–gene annotations, and eQTL data can be obtained through the [SCREEN platform](https://screen.wenglab.org/downloads).

Reference genome assemblies are available from the [Genome Reference Consortium](https://www.ncbi.nlm.nih.gov/grc/).

## License

This project is distributed under the MIT License. See [`LICENSE`](LICENSE) for the full license text.

## Citation

If you use this workflow in your research, please cite the associated study and the relevant external resources used by the workflow.
yet to come...


This ReadMe contains AI-generated content.

