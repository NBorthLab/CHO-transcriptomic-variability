#!/usr/bin/env nextflow

nextflow.enable.moduleBinaries = true


include { RNASEQ           } from "./workflows/01_rnaseq"
include { CHROMATIN_STATES } from "./workflows/chromatin_states.nf"
include { MOTIFS           } from "./workflows/motifs.nf"

include { processVersionsFromYaml } from "./modules/local/util"
include { formatNextflowVersion   } from "./modules/local/util"
include { sendNotification        } from "./modules/local/util"

workflow {

    // ================================================
    //   RNA-seq processing workflow
    // ================================================

    if (params.rnaseq) {
        log.info "Running RNA-seq workflow"

        // Genome files
        ch_gtf         = channel.fromPath(params.annotation_gtf)
        ch_genome      = channel.fromPath(params.genome)
        ch_samplesheet = channel.fromPath(params.input)

        RNASEQ(
            ch_gtf,
            ch_genome,
            ch_samplesheet
        )

        // Process and save versions
        RNASEQ.out.versions
            .map { version -> processVersionsFromYaml(version) }
            .mix(channel.of(formatNextflowVersion()))
            .unique()
            .collectFile(
                storeDir: "results/pipeline_info",
                name:     "all_versions.yml",
                sort:     true,
                newLine:  true
            )
    }


    // ================================================
    //   Analysis
    // ================================================

    if (params.chromatin) {
        CHROMATIN_STATES()
    }

    if (params.motifs) {
        MOTIFS()
    }


    workflow.onComplete = {
        log.info "Workflow completed at $workflow.complete"
        log.info "Duration: $workflow.duration"
        log.info "Execution status: ${workflow.success ? 'OK' : 'Failed'}"

        if (secrets.NTFY_URL != null) {
            sendNotification("${workflow.success ? 'OK': 'FAIL'}", secrets.NTFY_URL)
        }
    }
}



