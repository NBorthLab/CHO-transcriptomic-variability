process MEME_FIMO {

    label 'small'

    container "https://depot.galaxyproject.org/singularity/meme%3A5.5.9--pl5321h1ca524f_0"
    publishDir "${params.outdir}/motifs/fimo"

    input:
    tuple val(meta), path(fasta)
    each path(motif)

    output:
    path "*best_site.narrowPeak"
    path "*fimo.tsv"
    path "*fimo.gff"
    path "*fimo.html"
    path "versions.yml", emit: versions, optional: true

    script:
    def m = motif.name =~ /(\w+)_.*\.meme/
    def motif_id = m[0][1]
    """
    fimo $motif $fasta

    mv fimo_out/best_site.narrowPeak ${meta.id}.${motif_id}.best_site.narrowPeak
    mv fimo_out/fimo.tsv ${meta.id}.${motif_id}.fimo.tsv
    mv fimo_out/fimo.gff ${meta.id}.${motif_id}.fimo.gff
    mv fimo_out/fimo.html ${meta.id}.${motif_id}.fimo.html
    """
}
