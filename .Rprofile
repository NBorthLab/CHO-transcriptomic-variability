# if (Sys.info()["nodename"] == "dev.yawin") {
#   Sys.setenv("RENV_PATHS_ROOT" = "/run/host/data/mriedl/R_renv_cache")
# } else {
#   Sys.setenv("RENV_PATHS_ROOT" = "/data/mriedl/R_renv_cache")
# }
# Sys.setenv("RENV_PATHS_LIBRARY" ="renv/library")

source("renv/activate.R")

# DISPLAY envvar set to xpra session
if (stringr::str_detect(Sys.info()["nodename"], "rkn")) {
  message("DISPLAY environment variable set to :103")
  Sys.setenv("DISPLAY" = ":103")
}

params <- yaml::read_yaml("targets_config.yaml")
