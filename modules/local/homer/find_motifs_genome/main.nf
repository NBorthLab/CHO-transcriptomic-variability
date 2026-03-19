process FIND_MOTIFS_GENOME {

    cpus 8

    container "community.wave.seqera.io/library/homer:5.1--851c883b9d1473f9"
    publishDir "${params.outdir}/motifs/homer/${meta.id}"

    input:
    tuple val(meta), path(tss_bed)
    path all_tss_bed
    path genome

    output:
    path 'motifs_out'

    script:
    // WARN:
    // -useNewBg -> Positional dependent background correction
    //     ---> Segmentation Fault when trying to run that!
    """
    findMotifsGenome.pl \\
        $tss_bed \\
        $genome \\
        motifs_out \\
        -bg $all_tss_bed \\
        -size -300,100 \\
        -len 6,8,10 \\
        -mask \\
        -p ${task.cpus}
    """
}
