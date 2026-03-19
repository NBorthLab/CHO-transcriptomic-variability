process FEATURECOUNTS {

    tag "$meta.id"
    label 'big'

    container "https://depot.galaxyproject.org/singularity/subread:2.0.8--h577a1d6_0"
    publishDir "${params.outdir}/counts/${meta.dataset}",
        mode: "copy",
        saveAs: { filename -> filename.equals("versions.yml") ? null : filename }

    input:
    tuple val(meta), path(bam)
    each path(annotation)

    output:
    tuple val(meta), path("*counts.txt"),         emit: counts
    tuple val(meta), path("*counts.txt.summary"), emit: summary
    path "versions.yml",                          emit: versions

    script:
    def strandedness = 0
    if (meta.strandedness == "forward") {
        strandedness = 1
    } else if (meta.strandedness == "reverse") {
        strandedness = 2
    }
    def p = meta.paired_end ? "-p --countReadPairs" : ""

    """
    featureCounts \\
        -T $task.cpus \\
        $p \\
        -s $strandedness \\
        -a $annotation \\
        -o ${meta.id}.counts.txt \\
        $bam

    # Rename the column with gene counts to the sample's id
    sed -i -e 's/${bam}/${meta.id}/g' ${meta.id}.counts.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        subread: \$( echo \$(featureCounts -v 2>&1) | sed -e "s/featureCounts v//g")
    END_VERSIONS
    """
}
