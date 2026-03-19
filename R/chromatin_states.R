extract_genes <- function(pc, dec) {
  decile_str <- paste0("decile", str_pad(dec, width = 2, pad = "0"))
  genes <- pc %>%
    filter(decile == dec) %>%
    pull(gene)
  genelist_file <- str_glue("results/chromatin_states/{decile_str}.txt")
  writeLines(genes, genelist_file)
  return(genelist_file)
}

subset_annotation <- function(dec, genelist_file) {
  decile_str <- paste0("decile", str_pad(dec, width = 2, pad = "0"))
  system2(
    "conda",
    str_glue(
      " run -n P06-python python scripts/subset-annotation.py",
      " --annotation resources/chromatin_states/annotation_aliased.gtf",
      " --genelist {genelist_file}",
      " --log tmp/subset_annotation_{decile_str}_log.txt"
    )
  )
  inputs_file <- str_glue("results/chromatin_states/{decile_str}_input.txt")
  writeLines(
    str_glue("{decile_str}_{c('genes', 'tss2kb', 'exons')}.bed"),
    con = inputs_file
  )
  return(inputs_file)
}

overlap_enrich <- function(dec, tp, inputs_file, segment = NULL) {
  decile_str <- paste0("decile", str_pad(dec, width = 2, pad = "0"))
  inputs_dir <- str_glue("results/chromatin_states")

  # FIXME: Remove hard coded segmentation file. Allow computation for all
  # timepoints.
  segment <- str_glue("resources/chromatin_states/{tp}_11_all_segments.bed")

  output_prefix <- str_glue(
    "results/chromatin_states/overlap/{tp}/{decile_str}"
  )

  if (!dir.exists(dirname(output_prefix))) {
    dir.create(dirname(output_prefix))
  }

  command_args <- str_glue(
    " run -n P06-chromhmm ChromHMM.sh OverlapEnrichment",
    " -f {inputs_file}",
    " {segment}",
    " {inputs_dir}",
    " {output_prefix}"
  )

  system2("conda", command_args)

  return(str_glue("{output_prefix}.txt"))
}

segment_chromstates <- function(chromstates) {
  chromstate_segments_list <- GenomicRanges::tile(chromstates, width = 200)
  chromstate_segments_lengths <- S4Vectors::elementNROWS(
    chromstate_segments_list
  )
  chromstate_segments <- unlist(chromstate_segments_list)
  segments_labels <- rep(chromstates$name, chromstate_segments_lengths)
  chromstate_segments$state <- segments_labels

  return(chromstate_segments)
}

calculate_coverage <- function(dec, tp) {
  decile_str <- paste0("decile", str_pad(dec, width = 2, pad = "0"))
  decile_bed <- rtracklayer::import(
    str_glue("results/chromatin_states/{decile_str}_genes.bed"),
    format = "bed"
  )

  # TODO:
  # New approach: Use the coverage of every single gene to calculate mean and
  # CI95 of the coverage of each state in deciles.
  # Then I can do statistics.
  # However, how am I going to deal with the chromatin states timepoints. I
  # should probably average the coverages over a few timepoints, or pick one
  # representative timepoint.

  # Add +- 10kbp around gene bodies
  decile_bed <- decile_bed %>%
    plyranges::anchor_center() %>%
    plyranges::stretch(2 * 10000)

  chromstate_tp <- rtracklayer::import(
    str_glue("resources/chromatin_states/{tp}_11_all_segments.bed")
  ) %>%
    segment_chromstates()

  intersect <- plyranges::join_overlap_intersect(decile_bed, chromstate_tp)

  coverages <- tibble(
    gene = intersect$name,
    state = intersect$state
  ) %>%
    group_by(gene, state) %>%
    summarise(count = n()) %>%
    ungroup() %>%
    group_by(gene) %>%
    mutate(
      coverage = count / sum(count) * 100,
      n_segments = sum(count)
    ) %>%
    mutate(
      tp = tp,
      decile = dec,
      state = str_remove(state, "E"),
    )
  # ungroup() %>%
  # group_by(state) %>%
  # summarise(
  #   avg_cov = mean(coverage),
  #   std_cov = sd(coverage),
  #   sem_cov = sd(coverage) / sqrt(length(coverage))
  # )

  # bind_rows(
  #   var_coverages %>%
  #     mutate(
  #       order = str_rank(state, numeric = TRUE),
  #       type = "var"
  #     ) %>%
  #     arrange(order),
  #   stable_coverages %>%
  #     mutate(
  #       order = str_rank(state, numeric = TRUE),
  #       type = "stable"
  #     ) %>%
  #     arrange(order)
  # ) %>%
  #   tidyplot(x = state, y = avg_cov, color = type) %>%
  #   add_mean_bar() %>%
  #   reorder_x_axis_labels(str_sort(var_coverages$state, numeric = TRUE))

  # overlap <- IRanges::subsetByOverlaps(chromstate_tp, decile_bed)
  # n_segments <- length(overlap)
  # n_state <- c(table(overlap$state)[paste0("E", 1:11)])
  # coverage <- (n_state / n_segments) * 100

  # return(
  #   tibble(
  #     coverage = coverage,
  #     decile = dec,
  #     tp = tp
  #   )
  # )
  return(coverages)
}

load_overlap_results <- function(decile, tp) {
  dec_str <- paste0("decile", str_pad(decile, width = 2, pad = "0"))
  read_tsv(
    str_glue("results/chromatin_states/overlap/{tp}/{dec_str}.txt"),
    show_col_types = FALSE
  ) %>%
    rename(
      state = "State (Emission order)",
      genome = "Genome %",
      genes = str_glue("{dec_str}_genes.bed"),
      tss2kb = str_glue("{dec_str}_tss2kb.bed"),
      exons = str_glue("{dec_str}_exons.bed")
    ) %>%
    mutate(decile = decile, tp = tp) %>%
    pivot_longer(
      c(genome, genes, tss2kb, exons),
      names_to = "region",
      values_to = "enrichment"
    ) %>%
    filter(state != "Base") %>%
    mutate(state_label = case_when(
      state == 1 ~ "PcRepr",
      state == 2 ~ "Quie",
      state == 3 ~ "ReprHet",
      state == 4 ~ "Tx",
      state == 5 ~ "WkGEnh",
      state == 6 ~ "WkEnh",
      state == 7 ~ "AEnh1",
      state == 8 ~ "AEnh2",
      state == 9 ~ "APr",
      state == 10 ~ "TSSUp",
      state == 11 ~ "TSSDn"
    )) %>%
    mutate(
      type = case_when(
        state %in% c(1, 2, 3) ~ "Repressive",
        state %in% c(4) ~ "Transcription",
        state %in% c(5, 6, 7, 8) ~ "Enhancer",
        state %in% c(9, 10, 11) ~ "Promoter"
      )
    ) %>%
    mutate(
      state = factor(state, level = 1:10),
      decile = factor(decile, level = rev(1:10)),
      tp = factor(tp),
      region = factor(region)
    ) %>%
    mutate(
      state_tp = interaction(state, tp)
    ) %>%
    mutate(
      log_enrichment = log2(enrichment)
    )
}

load_coverages <- function(...) {
  bind_rows(...) %>%
    mutate(state_label = case_when(
      state == 1 ~ "PcRepr",
      state == 2 ~ "Quie",
      state == 3 ~ "ReprHet",
      state == 4 ~ "Tx",
      state == 5 ~ "WkGEnh",
      state == 6 ~ "WkEnh",
      state == 7 ~ "AEnh1",
      state == 8 ~ "AEnh2",
      state == 9 ~ "APr",
      state == 10 ~ "TSSUp",
      state == 11 ~ "TSSDn"
    )) %>%
    mutate(
      type = case_when(
        state %in% c(1, 2, 3) ~ "Repressive",
        state %in% c(4) ~ "Transcription",
        state %in% c(5, 6, 7, 8) ~ "Enhancer",
        state %in% c(9, 10, 11) ~ "Promoter"
      )
    ) %>%
    mutate(
      tp = factor(tp, levels = paste0("Tp", 0:17)),
      state = factor(state, levels = 1:11),
      decile = factor(decile, levels = rev(1:10)),
      state_label = factor(state_label, levels = c(
        "PcRepr", "Quie", "ReprHet", "Tx", "WkGEnh", "WkEnh", "AEnh1", "AEnh2",
        "APr", "TSSUp", "TSSDn"
      ))
    ) %>%
    filter(!tp %in% c("Tp7", "Tp10")) %>%
    filter(!is.na(state))
}

load_emissions <- function() {
  read_tsv("resources/chromatin_states/emissions_11_all.txt",
    col_types = "cdddddd",
  ) %>%
    rename(`state` = "State (Emission order)") %>%
    pivot_longer(-state, names_to = "mark", values_to = "emission")
}
