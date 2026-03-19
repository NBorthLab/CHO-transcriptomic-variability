process INFER_STRAND {

    tag "$meta.id"
    label 'medium'

    container "https://depot.galaxyproject.org/singularity/rseqc:5.0.4--pyhdfd78af_0"
    publishDir "${params.outdir}/aligned_reads/${meta.dataset}",
        mode: "copy",
        saveAs: { filename -> filename.equals("versions.yml") ? null : filename }

    input:
    tuple val(meta), path(alignment)
    each path(annotation)

    output:
    tuple val(meta), path("*infer_experiment.out.txt"), path(alignment), emit: inferred
    path "versions.yml",                                                 emit: versions

    script:
    """
    infer_experiment.py -r $annotation -i $alignment \\
        > ${meta.id}.infer_experiment.out.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rseqc: \$(infer_experiment.py --version | sed -e "s/infer_experiment.py //g")
    END_VERSIONS
    """
}
