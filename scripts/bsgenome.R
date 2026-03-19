setwd("resources/genome/bsgenome")

options(timout = 3600)
Sys.setenv("R_DEFAULT_INTERNET_TIMEOUT" = "3600")

BSgenomeForge::forgeBSgenomeDataPkgFromNCBI(
  "GCF_003668045.3",
  "Markus Riedl <markus.riedl@boku.ac.at>"
)

# NOTE: Moved the package dir to
# resources/genome/BSgenome.Cgriseus.NCBI.CriGriPICRH1.0

devtools::build("BSgenome.Cgriseus.NCBI.CriGriPICRH1.0")

devtools::check_built(
  "BSgenome.Cgriseus.NCBI.CriGriPICRH1.0_1.0.0.tar.gz"
)
devtools::install_local(
  "BSgenome.Cgriseus.NCBI.CriGriPICRH1.0_1.0.0.tar.gz"
)
