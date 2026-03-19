include { MEME_FIMO          } from "../modules/local/meme/fimo"
include { ALIAS_FASTA        } from "../modules/local/util"
include { FIND_MOTIFS_GENOME } from "../modules/local/homer/find_motifs_genome"
include { GUNZIP             } from "../modules/local/util"

workflow MOTIFS {

    // ================================================
    //   Find known motifs
    // ================================================

    ch_tss_flanks = channel
        .fromPath("results/gene_features/*_tss_flanks.fa")
        .map { file ->
            def m = file.name =~ /(\w+)_tss_flanks\.fa/
            def meta = [id: m[0][1]]
            [meta, file]
        }

    motifs = channel.value(file("resources/motifs/*.meme"))

    MEME_FIMO(
        ch_tss_flanks,
        motifs
    )


    // ================================================
    //   De novo Motif discovery
    // ================================================

    genome_fasta_gz = file(
        "resources/genome/GCF_003668045.3_CriGri-PICRH-1.0_genomic.fna.gz"
    )

    aliases = file("resources/chromatin_states/PICRH_aliases_rev.txt")

    GUNZIP(
        genome_fasta_gz
    )

    ALIAS_FASTA(
        GUNZIP.out,
        aliases
    )

    ch_tss = channel
        .fromPath("results/gene_features/{1,10}_tss.bed")
        .map { file ->
            def m = file.name =~ /(\w+)_tss\.bed/
            def meta = [id: m[0][1]]
            [meta, file]
        }

    all_tss_bed = file(
        "results/gene_features/all_tss.bed"
    )

    FIND_MOTIFS_GENOME(
        ch_tss,
        all_tss_bed,
        ALIAS_FASTA.out
    )

// oras://community.wave.seqera.io/library/homer:5.1--586e112615c26ca1
}
