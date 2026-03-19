process RSEM_PREPAREREFERENCE {

    label 'huge'

    container "https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/23/23651ffd6a171ef3ba867cb97ef615f6dd6be39158df9466fe92b5e844cd7d59/data"
    publishDir "${params.outdir}/index",
        mode: "copy",
        saveAs: { filename -> filename.equals("versions.yml") ? null : filename }

    input:
    path genome_fasta
    path gtf

    output:
    path "genome_index", emit: index
    path "versions.yml", emit: versions

    script:
    """
    mkdir genome_index

    rsem-prepare-reference \\
        -p 38 \\
        --star \\
        --gtf $gtf \\
        $genome_fasta \\
        genome_index/genome

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rsem: \$(rsem-calculate-expression --version | sed -e "s/Current version: RSEM v//g")
        star: \$(STAR --version | sed -e "s/STAR_//g")
    END_VERSIONS
    """
}
