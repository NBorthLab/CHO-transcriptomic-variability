library(tibble)
library(targets)
library(tarchetypes)
library(crew)

params <- yaml::read_yaml("targets_config.yaml")

controller <- crew_controller_local(
  name = "my_controller",
  workers = 8
)

global_packages <- c(
  "tidyverse",
  "conflicted",
  "SummarizedExperiment"
)

tar_option_set(
  packages = global_packages,
  controller = controller
)

if (params$use_s3) {
  tar_option_set(
    repository = "aws",
    repository_meta = "aws",
    resources = tar_resources(
      aws = tar_resources_aws(
        bucket = "p06-genevariability",
        prefix = "results/R/preprocessing",
        region = "eu2",
        endpoint = "https://eu2.contabostorage.com"
      )
    )
  )
}

tar_source()

# ================================================
#   Inputs
# ================================================

inputs <- list(
  tar_target(
    counts_rds,
    command = "results/counts/all_gene_counts.rds",
    format = "file"
  ),
  tar_target(
    annotation_gtf,
    command =
      "resources/genome/GCF_003668045.3_CriGri-PICRH-1.0_genomic.gtf.gz",
    format = "file"
  ),
  tar_target(
    metadata_csv,
    command = "resources/metadata.csv",
    format = "file"
  )
)


# ================================================
#   Preprocessing
# ================================================

list(
  inputs,
  tar_target(
    counts_tables,
    command = read_counts_tables(counts_rds)
  ),
  tar_target(
    annotation,
    command = read_annotation(annotation_gtf)
  ),
  tar_target(
    metadata,
    command = get_metadata(metadata_csv)
  ),
  tar_target(
    raw_counts,
    command = make_summarized_experiment(counts_tables, metadata)
  ),
  tar_target(
    filtered_samples,
    command = filter_samples(raw_counts)
  ),
  tar_target(
    coding_counts,
    command = filter_protein_coding_genes(filtered_samples, annotation)
  )
)
