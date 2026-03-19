process COMBINE_COUNTS {

    label 'small'

    container "https://depot.galaxyproject.org/singularity/bioconductor-tximport:1.34.0--r44hdfd78af_0"
    publishDir "${params.outdir}/counts",
        mode: "copy",
        saveAs: { filename -> filename.equals("versions.yml") ? null : filename }

    input:
    path "counts/*"

    output:
    path("all_gene_counts.rds"), emit: counts
    path("versions.yml"),        emit: versions

    script:
    template "concatenate_counts.R"
}
