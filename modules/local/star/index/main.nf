process STAR_INDEX {

    label 'huge'

    container "https://depot.galaxyproject.org/singularity/star:2.7.11b--h5ca1c30_4"

    input:
    path genome
    path annotation

    output:
    path "genome_index", emit: index
    path "versions.yml", emit: versions

    script:
    """
    mkdir genome_index

    STAR \\
        --runMode genomeGenerate \\
        --genomeFastaFiles $genome \\
        --sjdbGTFfile $annotation \\
        --genomeDir genome_index \\
        --runThreadN $task.cpus

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        star: \$(STAR --version | sed -e "s/STAR_//g")
    END_VERSIONS
    """
}
