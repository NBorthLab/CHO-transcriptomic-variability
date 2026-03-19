process FASTQC {

    tag "$meta.id"
    label 'medium'

    container "https://depot.galaxyproject.org/singularity/trim-galore:0.6.9--hdfd78af_0"
    publishDir "${params.outdir}/fastqc/${meta.dataset}",
        mode: "copy",
        saveAs: { filename -> filename.equals("versions.yml") ? null : filename }

    input:
    tuple val(meta), path(reads)

    output:
    // reads.simpleName is the same as .getSimpleName()
    tuple val(meta), path("*_fastqc.zip"),  emit: zip
    tuple val(meta), path("*_fastqc.html"), emit: html
    path "versions.yml",                    emit: versions

    script:
    """
    fastqc $reads

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqc: \$( fastqc --version | sed '/FastQC v/!d; s/.*v//' )
    END_VERSIONS
    """
}
