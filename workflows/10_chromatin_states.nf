include { OVERLAP_ENRICH } from "../modules/local/chromhmm/overlap_enrich"
include { SUBSET_ANNOTATION } from "../modules/local/subset_annotation"

nextflow.enable.moduleBinaries = true

workflow CHROMATIN_STATES {

    ch_genelist = channel
        .fromPath("results/chromatin_states/decile??.txt")
        .map { file ->
            def m = file.name =~ /(decile..)\.txt/
            def meta = [decile: m[0][1]]
            [meta, file]
        }

    annotation = file("resources/chromatin_states/annotation_aliased.gtf")

    ch_segments = channel
        .fromPath("resources/chromatin_states/*_11_all_segments.bed")
        .map { file ->
            def m = file.name =~ /(.*)_11_all_segments\.bed/
            def meta = [tp: m[0][1]]
            [meta, file]
        }

    SUBSET_ANNOTATION(
        ch_genelist,
        annotation
    )

    ch_input = SUBSET_ANNOTATION.out.inputs
        .join(SUBSET_ANNOTATION.out.genes)
        .join(SUBSET_ANNOTATION.out.exons)
        .join(SUBSET_ANNOTATION.out.tss2kb)
        .map { meta, inputs, genes, exons, tss2kb ->
            [meta, inputs, [genes, exons, tss2kb]]
        }
        .combine(ch_segments)

    OVERLAP_ENRICH(
        ch_input
    )
}
