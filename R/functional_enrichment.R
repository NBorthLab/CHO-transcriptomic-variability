get_hamster_annotation <- function() {
  annotation_hub <- AnnotationHub::AnnotationHub()

  orgdb <- AnnotationHub::query(
    annotation_hub,
    c("Cricetulus griseus", "OrgDb"),
  )[[1]]

  annot_db_Cg <<- orgdb

  return(NULL)
}

get_all_orthologs <- function(genes) {
  orthologs <- orthogene::convert_orthologs(
    genes,
    gene_input     = "gene",
    gene_output    = "columns",
    input_species  = "Cricetulus griseus",
    output_species = "Mus musculus",
    method         = "gprofiler"
  )

  return(orthologs)
}

convert_orthologs <- function(genes, orthologs) {
  # data.frame(gene = genes) %>%
  genes %>%
    inner_join(orthologs, by = join_by(gene == input_gene)) %>%
    pull(ortholog_gene)
}

get_mouse_universe <- function(orthologs) {
  # Unvierse is all ortholog genes
  universe <- orthologs$ortholog_gene

  return(universe)
}

# AnnotationDbi::keys(org.Mm.eg.db, "UNIPROT")
#
#
#
# https://www.ebi.ac.uk/QuickGO/annotations?
# goUsage=slim&
# goUsageRelationships=is_a,part_of,occurs_in&
# goId=GO:1901135,GO:0140053,GO:0140014,GO:0140013,GO:0098754,GO:0098542,GO:0072659,GO:0071941,GO:0071554,GO:0065003,GO:0061024,GO:0061007,GO:0055086,GO:0055085,GO:0051604,GO:0050886,GO:0050877,GO:0048870,GO:0048856,GO:0044782,GO:0042254,GO:0042060,GO:0034330,GO:0032200,GO:0031047,GO:0030198,GO:0030163,GO:0030154,GO:0023052,GO:0022600,GO:0022414,GO:0016192,GO:0016073,GO:0016071,GO:0015979,GO:0012501,GO:0009101,GO:0007163,GO:0007155,GO:0007059,GO:0007040,GO:0007031,GO:0007018,GO:0007010,GO:0007005,GO:0006954,GO:0006914,GO:0006913,GO:0006886,GO:0006790,GO:0006766,GO:0006629,GO:0006575,GO:0006520,GO:0006457,GO:0006399,GO:0006355,GO:0006351,GO:0006325,GO:0006310,GO:0006281,GO:0006260,GO:0006091,GO:0005975,GO:0003016,GO:0003014,GO:0003013,GO:0003012,GO:0002376,GO:0002181,GO:0000910,GO:0000278,GO:0043226,GO:0031410,GO:0031012,GO:0030312,GO:0009579,GO:0009536,GO:0005929,GO:0005886,GO:0005856,GO:0005840,GO:0005829,GO:0005815,GO:0005811,GO:0005794,GO:0005783,GO:0005777,GO:0005773,GO:0005768,GO:0005764,GO:0005739,GO:0005730,GO:0005694,GO:0005654,GO:0005635,GO:0005634,GO:0005618,GO:0005615,GO:0005576,GO:0000228,GO:0140657,GO:0140313,GO:0140299,GO:0140223,GO:0140110,GO:0140104,GO:0140098,GO:0140097,GO:0140096,GO:0120274,GO:0098772,GO:0098631,GO:0090729,GO:0060090,GO:0060089,GO:0048018,GO:0045735,GO:0045182,GO:0044183,GO:0042393,GO:0038024,GO:0031386,GO:0016874,GO:0016853,GO:0016829,GO:0016787,GO:0016740,GO:0016491,GO:0016209,GO:0009975,GO:0008289,GO:0008092,GO:0005215,GO:0005198,GO:0003924,GO:0003824,GO:0003774,GO:0003723,GO:0003677,GO:0001618&
# taxonId=10090
# TODO: Continue with the API call to get the slim annotation for mouse:
# https://www.ebi.ac.uk/QuickGO/annotations?goUsage=slim&goUsageRelationships=is_a,part_of,occurs_in&goId=GO:1901135,GO:0140053,GO:0140014,GO:0140013,GO:0098754,GO:0098542,GO:0072659,GO:0071941,GO:0071554,GO:0065003,GO:0061024,GO:0061007,GO:0055086,GO:0055085,GO:0051604,GO:0050886,GO:0050877,GO:0048870,GO:0048856,GO:0044782,GO:0042254,GO:0042060,GO:0034330,GO:0032200,GO:0031047,GO:0030198,GO:0030163,GO:0030154,GO:0023052,GO:0022600,GO:0022414,GO:0016192,GO:0016073,GO:0016071,GO:0015979,GO:0012501,GO:0009101,GO:0007163,GO:0007155,GO:0007059,GO:0007040,GO:0007031,GO:0007018,GO:0007010,GO:0007005,GO:0006954,GO:0006914,GO:0006913,GO:0006886,GO:0006790,GO:0006766,GO:0006629,GO:0006575,GO:0006520,GO:0006457,GO:0006399,GO:0006355,GO:0006351,GO:0006325,GO:0006310,GO:0006281,GO:0006260,GO:0006091,GO:0005975,GO:0003016,GO:0003014,GO:0003013,GO:0003012,GO:0002376,GO:0002181,GO:0000910,GO:0000278,GO:0043226,GO:0031410,GO:0031012,GO:0030312,GO:0009579,GO:0009536,GO:0005929,GO:0005886,GO:0005856,GO:0005840,GO:0005829,GO:0005815,GO:0005811,GO:0005794,GO:0005783,GO:0005777,GO:0005773,GO:0005768,GO:0005764,GO:0005739,GO:0005730,GO:0005694,GO:0005654,GO:0005635,GO:0005634,GO:0005618,GO:0005615,GO:0005576,GO:0000228,GO:0140657,GO:0140313,GO:0140299,GO:0140223,GO:0140110,GO:0140104,GO:0140098,GO:0140097,GO:0140096,GO:0120274,GO:0098772,GO:0098631,GO:0090729,GO:0060090,GO:0060089,GO:0048018,GO:0045735,GO:0045182,GO:0044183,GO:0042393,GO:0038024,GO:0031386,GO:0016874,GO:0016853,GO:0016829,GO:0016787,GO:0016740,GO:0016491,GO:0016209,GO:0009975,GO:0008289,GO:0008092,GO:0005215,GO:0005198,GO:0003924,GO:0003824,GO:0003774,GO:0003723,GO:0003677,GO:0001618&taxonId=10090
#
#
#
# url <- paste0(
# "https://www.ebi.ac.uk/QuickGO/services/annotation/downloadSearch?goId=GO%3A1901135%2CGO%3A0140053%2CGO%3A0140014%2CGO%3A0140013%2CGO%3A0098754%2CGO%3A0098542%2CGO%3A0072659%2CGO%3A0071941%2CGO%3A0071554%2CGO%3A0065003%2CGO%3A0061024%2CGO%3A0061007%2CGO%3A0055086%2CGO%3A0055085%2CGO%3A0051604%2CGO%3A0050886%2CGO%3A0050877%2CGO%3A0048870%2CGO%3A0048856%2CGO%3A0044782%2CGO%3A0042254%2CGO%3A0042060%2CGO%3A0034330%2CGO%3A0032200%2CGO%3A0031047%2CGO%3A0030198%2CGO%3A0030163%2CGO%3A0030154%2CGO%3A0023052%2CGO%3A0022600%2CGO%3A0022414%2CGO%3A0016192%2CGO%3A0016073%2CGO%3A0016071%2CGO%3A0015979%2CGO%3A0012501%2CGO%3A0009101%2CGO%3A0007163%2CGO%3A0007155%2CGO%3A0007059%2CGO%3A0007040%2CGO%3A0007031%2CGO%3A0007018%2CGO%3A0007010%2CGO%3A0007005%2CGO%3A0006954%2CGO%3A0006914%2CGO%3A0006913%2CGO%3A0006886%2CGO%3A0006790%2CGO%3A0006766%2CGO%3A0006629%2CGO%3A0006575%2CGO%3A0006520%2CGO%3A0006457%2CGO%3A0006399%2CGO%3A0006355%2CGO%3A0006351%2CGO%3A0006325%2CGO%3A0006310%2CGO%3A0006281%2CGO%3A0006260%2CGO%3A0006091%2CGO%3A0005975%2CGO%3A0003016%2CGO%3A0003014%2CGO%3A0003013%2CGO%3A0003012%2CGO%3A0002376%2CGO%3A0002181%2CGO%3A0000910%2CGO%3A0000278%2CGO%3A0043226%2CGO%3A0031410%2CGO%3A0031012%2CGO%3A0030312%2CGO%3A0009579%2CGO%3A0009536%2CGO%3A0005929%2CGO%3A0005886%2CGO%3A0005856%2CGO%3A0005840%2CGO%3A0005829%2CGO%3A0005815%2CGO%3A0005811%2CGO%3A0005794%2CGO%3A0005783%2CGO%3A0005777%2CGO%3A0005773%2CGO%3A0005768%2CGO%3A0005764%2CGO%3A0005739%2CGO%3A0005730%2CGO%3A0005694%2CGO%3A0005654%2CGO%3A0005635%2CGO%3A0005634%2CGO%3A0005618%2CGO%3A0005615%2CGO%3A0005576%2CGO%3A0000228%2CGO%3A0140657%2CGO%3A0140313%2CGO%3A0140299%2CGO%3A0140223%2CGO%3A0140110%2CGO%3A0140104%2CGO%3A0140098%2CGO%3A0140097%2CGO%3A0140096%2CGO%3A0120274%2CGO%3A0098772%2CGO%3A0098631%2CGO%3A0090729%2CGO%3A0060090%2CGO%3A0060089%2CGO%3A0048018%2CGO%3A0045735%2CGO%3A0045182%2CGO%3A0044183%2CGO%3A0042393%2CGO%3A0038024%2CGO%3A0031386%2CGO%3A0016874%2CGO%3A0016853%2CGO%3A0016829%2CGO%3A0016787%2CGO%3A0016740%2CGO%3A0016491%2CGO%3A0016209%2CGO%3A0009975%2CGO%3A0008289%2CGO%3A0008092%2CGO%3A0005215%2CGO%3A0005198%2CGO%3A0003924%2CGO%3A0003824%2CGO%3A0003774%2CGO%3A0003723%2CGO%3A0003677%2CGO%3A0001618&goUsage=slim&taxonId=10090"
# )

# httr::GET(url, httr::accept("text/tsv"))

enrich_GO_terms <- function(genes, universe, ont = "ALL") {
  get_hamster_annotation()

  terms_raw <- clusterProfiler::enrichGO(
    gene          = genes,
    OrgDb         = annot_db_Cg,
    keyType       = "SYMBOL",
    ont           = ont,
    pvalueCutoff  = 0.05,
    pAdjustMethod = "BH",
    minGSSize     = 4,
    universe      = universe
  )

  # Simplify the GO terms to reduce redundancy of resulting terms
  terms_returned <- tryCatch(
    {
      # clusterProfiler::simplify(terms_raw)
      terms_raw
    },
    error = function(e) {
      cat("WARNING: Terms not simplified!")
      return(terms_raw)
    }
  )

  terms_returned <- terms_returned %>%
    clusterProfiler::mutate(
      log_p_adjust = -log10(p.adjust),
      gene_pct = map_dbl(GeneRatio, (\(x) eval(parse(text = x)) * 100))
    )

  return(terms_returned)
}

enrich_KEGG_pathways <- function(genes, universe, organism = "cge") {
  db <- get_CH_annotation()

  genes_convert <- clusterProfiler::bitr(genes, "SYMBOL", "ENTREZID", db)
  universe_convert <- clusterProfiler::bitr(universe, "SYMBOL", "ENTREZID", db)

  pathways_raw <- clusterProfiler::enrichKEGG(
    gene = genes_convert$ENTREZID,
    organism = organism,
    keyType = "ncbi-geneid",
    pvalueCutoff = 1,
    pAdjustMethod = "BH",
    universe = universe_convert$ENTREZID,
    minGSSize = 10,
    maxGSSize = 500
  )

  pathways_raw@result <- pathways_raw@result %>%
    mutate(log_p_adjust = -log(p.adjust, 10))

  return(pathways_raw)
}

get_secRecon_ontology <- function(file) {
  secRecon_xlsx <- readxl::read_xlsx(file)

  secRecon <- secRecon_xlsx %>%
    dplyr::select(
      `CHO GENE SYMBOL`,
      `CHO ENTREZID`,
      `Process 1`,
      `Process 2`,
      `Process 3`,
      `Process 4`,
      `Process 5`,
      `Process 6`,
      `Process 7`,
      `Process 8`,
      `Process 9`
    ) %>%
    tidyr::pivot_longer(
      -c("CHO GENE SYMBOL", "CHO ENTREZID"),
      names_to = "Process",
      values_to = "Ontology"
    ) %>%
    dplyr::select(-Process) %>%
    dplyr::rename(
      Entrezid = `CHO ENTREZID`,
      Symbol = `CHO GENE SYMBOL`
    ) %>%
    dplyr::filter(!is.na(Ontology) & !is.na(Entrezid))

  secRecon_ontology <- secRecon %>%
    dplyr::select(Ontology, Symbol)

  return(secRecon_ontology)
}

enrich_secRecon <- function(geneset, universe, secRecon) {
  enriched_terms <- clusterProfiler::enricher(
    geneset,
    TERM2GENE = secRecon,
    universe = universe,
    minGSSize = 1
  )

  enriched_terms@result <- enriched_terms@result %>%
    mutate(log_p_adjust = -log(p.adjust, 10))

  return(enriched_terms)
}

get_slim_annotation <- function(tsv, genes) {
  get_hamster_annotation()

  generic_slim_terms <- read_tsv(tsv) %>%
    mutate(
      term = str_replace(`?x`, ".*GO_([0-9]+).", "GO:\\1")
    ) %>%
    dplyr::rename(label = `?label`) %>%
    dplyr::select(-`?x`)

  term2name <- dplyr::select(generic_slim_terms, term, label)

  genes_full_annotation <- AnnotationDbi::select(
    annot_db_Cg,
    keys = genes,
    columns = c("SYMBOL", "GOALL"),
    keytype = "SYMBOL"
  )

  genes_slim_annotation <- genes_full_annotation %>%
    right_join(generic_slim_terms, by = join_by(GOALL == term))

  term2gene <- dplyr::select(genes_slim_annotation, GOALL, SYMBOL)

  list(term2gene = term2gene, term2name = term2name)
}
# get_slim_annotation("resources/gene_ontology/goslim_generic.tsv", first_pc$gene)
# tsv <- "resources/gene_ontology/goslim_generic.tsv"
run_gsea <- function(gene_list, term2gene, term2name = NULL) {
  clusterProfiler::GSEA(gene_list, TERM2GENE = term2gene, TERM2NAME = term2name)
}

run_gseGO <- function(gene_list) {
  get_hamster_annotation()
  go_gsea <- clusterProfiler::gseGO(
    gene_list,
    OrgDb = annot_db_Cg,
    keyType = "SYMBOL"
  )

  go_gsea
}

run_gseKEGG <- function(gene_list) {
  get_hamster_annotation()

  # Convert gene symbols to NCBI gene IDs
  gene_list_id <- data.frame(
    gene = names(gene_list),
    pc = gene_list
  ) %>%
    mutate(
      gene_id = clusterProfiler::bitr(
        gene, "SYMBOL", "ENTREZID", annot_db_Cg,
        drop = FALSE
      )$ENTREZID
    ) %>%
    filter(!is.na(gene_id)) %>%
    pull(pc, name = gene_id)

  kegg_gsea <- clusterProfiler::gseKEGG(
    gene_list_id,
    organism = "cge",
    keyType = "ncbi-geneid"
  )

  kegg_gsea
}
