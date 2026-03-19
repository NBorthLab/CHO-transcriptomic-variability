library(targets)
library(tidyverse)

params <- yaml::read_yaml("targets_config.yaml")

if (FALSE) {
  tar_manifest()
  options(browser = "brave")
  tar_manifest() %>% DT::datatable()
}

Sys.setenv(TAR_PROJECT = "preprocessing")
tar_make()

Sys.setenv(TAR_PROJECT = "transformation")
tar_make()

Sys.setenv(TAR_PROJECT = "principal_components")
tar_make()

Sys.setenv(TAR_PROJECT = "loess_residuals")
tar_make()

Sys.setenv(TAR_PROJECT = "shannon_entropy")
tar_make()

Sys.setenv(TAR_PROJECT = "measures_comparison")
tar_make()

Sys.setenv(TAR_PROJECT = "gene_sets")
tar_make()

Sys.setenv(TAR_PROJECT = "functional_enrichment")
tar_make()
# tar_make(callr_function = NULL, use_crew = FALSE, as_job = FALSE)

# NOTE: Nextflow pipeline chromatin
Sys.setenv(TAR_PROJECT = "chromatin_states")
tar_make()

# NOTE: Nextflow pipeline motifs
Sys.setenv(TAR_PROJECT = "motifs")
tar_make()

Sys.setenv(TAR_PROJECT = "gene_features")
tar_make()

Sys.setenv(TAR_PROJECT = "gene_characteristics")
tar_make()

Sys.setenv(TAR_PROJECT = "metadata_loadings")
tar_make()

Sys.setenv(TAR_PROJECT = "versions")
tar_make()

if (params$entropy) {
}

# Sys.setenv(TAR_PROJECT = "exploration")
# tar_make(starts_with("mean_sd_plot"))
# tar_make()

# Sys.setenv(TAR_PROJECT = "variance")
# tar_make()
