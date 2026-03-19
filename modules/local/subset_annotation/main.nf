process SUBSET_ANNOTATION {

    label 'small'
    container "https://depot.galaxyproject.org/singularity/pygtftk:1.6.0--py39h2add14b_0"

    input:
    tuple val(meta), path(genelist)
    path(annotation)

    output:
    tuple val(meta), path("*_input.txt"), emit: inputs
    tuple val(meta), path("*_genes.bed"), emit: genes
    tuple val(meta), path("*_exons.bed"), emit: exons
    tuple val(meta), path("*_tss.bed"), emit: tss
    tuple val(meta), path("*_tes.bed"), emit: tes
    tuple val(meta), path("*_tss2kb.bed"), emit: tss2kb

    script:
    """
    subset-annotation.py \\
        --annotation $annotation \\
        --genelist $genelist

    echo "${meta.decile}_genes.bed" > ${meta.decile}_input.txt
    echo "${meta.decile}_exons.bed" >> ${meta.decile}_input.txt
    echo "${meta.decile}_tss2kb.bed" >> ${meta.decile}_input.txt
    """
}
