process STAR_ALIGN {

    tag "$meta.id"
    label 'big'

    container "https://depot.galaxyproject.org/singularity/star:2.7.11b--h5ca1c30_4"
    publishDir "${params.outdir}/aligned_reads/${meta.dataset}",
        mode: "copy",
        saveAs: { filename -> filename.equals("versions.yml") ? null : filename }

    input:
    tuple val(meta), path(reads)
    path index

    output:
    tuple val(meta), path("*Aligned.out.bam"),                 emit: genome_bam
    tuple val(meta), path("*Aligned.toTranscriptome.out.bam"), emit: transcriptome_bam
    tuple val(meta), path("*Log.final.out"),                   emit: log_final
    tuple val(meta), path("*Log.out"),                         emit: log
    path "versions.yml",                                       emit: versions

    script:
    def reads1 = []
    def reads2 = []
    meta.paired_end
        ? reads.eachWithIndex { v, i -> (i & 1 ? reads2 : reads1) << v }
        : [reads].flatten().each { reads1 << it }

    // STAR options derived from the RSEM calculate expression  execution perl
    // script:
    // https://github.com/deweylab/RSEM/blob/master/rsem-calculate-expression
    """
    STAR \\
        --runThreadN ${task.cpus} \\
        --genomeDir ${index} \\
        --genomeLoad LoadAndRemove \\
        --readFilesCommand gunzip -c \\
        --readFilesIn ${reads1.join(",")} ${reads2.join(",")} \\
        --outFileNamePrefix ${meta.id}. \\
        --outSAMtype BAM Unsorted \\
        --outSAMheaderHD \\@HD VN:1.4 SO:unsorted \\
        --outSAMunmapped Within \\
        --outSAMattributes NH HI AS NM MD \\
        --quantMode TranscriptomeSAM \\
        --outFilterType BySJout \\
        --outFilterMultimapNmax 20 \\
        --outFilterMismatchNmax 999 \\
        --outFilterMismatchNoverLmax 0.04 \\
        --alignIntronMin 20 \\
        --alignIntronMax 1000000 \\
        --alignMatesGapMax 1000000 \\
        --alignSJoverhangMin 8 \\
        --alignSJDBoverhangMin 1 \\
        --sjdbScore 1



    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        star: \$(STAR --version | sed -e "s/STAR_//g")
    END_VERSIONS
    """
}
