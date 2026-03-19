process MULTIQC {

    tag "$dataset"
    label 'medium'

    clusterOptions "-w rkn07"

    container "https://depot.galaxyproject.org/singularity/multiqc:1.30--pyhdfd78af_0"
    publishDir "${params.outdir}/multiqc/${dataset}",
        mode: "copy",
        saveAs: { filename -> filename.equals("versions.yml") ? null : filename }


    input:
    tuple val(dataset), path(multiqc_files)

    output:
    path "*multiqc_report.html", emit: report
    path "*_data",               emit: data

    script:
    """
    multiqc .
    """
}
