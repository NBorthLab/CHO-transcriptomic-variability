process RSEM_CALCULATEEXPRESSION {

    tag "$meta.id"
    label 'big'

    container "https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/23/23651ffd6a171ef3ba867cb97ef615f6dd6be39158df9466fe92b5e844cd7d59/data"
    publishDir "${params.outdir}/counts/${meta.dataset}",
        mode: "copy",
        saveAs: { filename -> filename.equals("versions.yml") ? null : filename }

    input:
    tuple val(meta), path(transcriptome_bam)
    path genome_index

    output:
    tuple val(meta), path("*.genes.results"),           emit: counts_gene
    tuple val(meta), path("*.isoforms.results"),        emit: counts_transcript
    tuple val(meta), path("*.stat"),                    emit: stat
    path  "versions.yml",                               emit: versions


    script:

    def paired_end = meta.paired_end ? "--paired-end" : ""

    """
    rsem-calculate-expression \\
        --num-threads $task.cpus \\
        --strandedness $meta.strandedness \\
        $paired_end \\
        --estimate-rspd \\
        --seed 1 \\
        --temporary-folder ./tmp/ \\
        --alignments \\
        $transcriptome_bam \\
        $genome_index/genome \\
        $meta.id

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rsem: \$(rsem-calculate-expression --version | sed -e "s/Current version: RSEM v//g")
        star: \$(STAR --version | sed -e "s/STAR_//g")
    END_VERSIONS
    """
}
