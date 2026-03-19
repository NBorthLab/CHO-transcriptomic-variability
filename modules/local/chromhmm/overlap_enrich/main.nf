process OVERLAP_ENRICH {

    label 'small'

    container "https://depot.galaxyproject.org/singularity/chromhmm:1.25--hdfd78af_0"
    publishDir "${params.outdir}/chromatin_states/overlap/${meta2.tp}",
        mode: "copy",
        saveAs: { filename -> filename.equals("versions.yml") ? null : filename }

    input:
    tuple val(meta), path(inputfile), path(annotations), val(meta2), path(segments)

    output:
    path ("${meta.decile}.txt")
    path "versions.yml", emit: versions, optional: true

    script:
    """
    ChromHMM.sh OverlapEnrichment \\
        -f $inputfile \\
        $segments \\
        ./ \\
        ${meta.decile}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        chromhmm: \$( ChromHMM.sh Version | sed -E "s/.*Version.([0-9]+).*)/\\1/")
    END_VERSIONS
    """
}
