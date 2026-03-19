#!/usr/bin/env bash

set -euo pipefail

# ================================================
#   Constants and Variables
# ================================================

CONDA_NEXTFLOW_ENV="nextflow"
CONDA_R_ENV="P06-r-renv"

rnaseq=""
analysis=""

# ================================================
#   Functions
# ================================================

tar_make() {
    conda run -n $CONDA_R_ENV Rscript -e 'targets::tar_make()'
}

nextflow_bin() {
    conda run -n $CONDA_NEXTFLOW_ENV nextflow "$@"
}

help() {
    echo "Usage: ${0} [subcomand]"
    echo "subcomands: rnaseq, analysis, help"
    echo "subcomands can be stacked"
}

# ================================================
#   Parse arguments
# ================================================

if [[ $# == 0 ]]; then
    help
    exit 1
fi

while [[ $# -gt 0 ]]; do
    case $1 in
    rnaseq)
        rnaseq="true"
        shift
        ;;
    analysis)
        analysis="true"
        shift
        ;;
    help)
        help
        exit 0
        ;;
    *)
        help
        exit 1
        ;;
    esac
done

# ================================================
#   RNA-seq processing pipeline
# ================================================

# Implemented as Nextflow pipeline

if [[ -n $rnaseq ]]; then
    nextflow run . -resume -params-file ./params.yml --rnaseq
fi

# ================================================
#   Analysis pipeline
# ================================================

# Mix of R targets workflows and Nextflow pipelines

if [[ -n $analysis ]]; then
    TAR_PROJECT=preprocessing tar_make
    TAR_PROJECT=transformation tar_make
    TAR_PROJECT=loess_residuals tar_make
    TAR_PROJECT=shannon_entropy tar_make
    TAR_PROJECT=measures_comparison tar_make
    TAR_PROJECT=principal_components tar_make
    TAR_PROJECT=gene_sets tar_make
    TAR_PROJECT=functional_enrichment tar_make
    nextflow run . -resume -profile local --chromatin
    TAR_PROJECT=chromatin_states tar_make
    TAR_PROJECT=gene_features tar_make
    nextflow_bin run . -resume -profile local --motifs
    TAR_PROJECT=motifs tar_make
fi
