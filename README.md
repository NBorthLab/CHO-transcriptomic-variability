# CHO transcriptomic variability

## File structure

```plain
.
├── analysis/                   # Analysis script in Quarto markdown
├── envs/                       # Conda environment YAML + txt pins
├── modules/                    # Nextflow modules
│   └── local/
├── plots/                      # Plots produced in the analysis
├── R/                          # R function used in the targets workflow
├── renv/                       # renv package stuff
├── resources/                  # Data resources (GO, secRecon, ...)
├── results/                    # Results of preprocessing and analysis
├── scripts/                    # Various utility scripts
├── workflows/                  # Definitions of Nextflow and targets pipelines
│   ├── 01_rnaseq/
│   ├── 02_preprocessing.R
│   ├── 03_transformation.R
│   └── ...
├── main.nf                     # Entrypoint for Nextflow pipelines
├── nextflow.config             # Configuration for Nextflow pipelines
├── params.yml                  # Parameters for Nextflow pipelines
├── README.md                   # You are here
├── renv.lock                   # Lockfile for renv
├── run.sh                      # Entrypoint for the full analysis
├── targets_config.yaml         # targets configurations
├── targets_main.R              # R script for running targets pipelines
└── _targets.yaml               # targets workflows manifest
```


## Workflow

The full processing and analysis workflow can be started by invoking the
`run.sh` script in this repository.

Importantly, Nextflow and R with the renv package need to be installed in conda
environments whose names must be specified on top of `run.sh` for the pipeline
to work. Please see the paper's Materials and methods for recommended versions
of the software. For the R environment, find the conda environment YAML
definition (and pinned conda package) in `envs/r-renv.yaml`
(`*.linux-64.pin.txt`).

Alternatively, if conda is not an option for whatever reason, change the
definition of the `tar_make()` and `nextflow_bin()` function in the script
accordingly.

Finally, run

```bash
bash run.sh [rnaseq|analysis|help]
```

e.g. for running the RNA-seq processing pipeline `bash run.sh rnaseq`, for
running the analysis pipeline `bash run.sh analysis`, or `bash run.sh help` to
print a short help message on usage of the script.


### Raw data download

Downloading of the raw data from the SRA was done using
[nf-core/fetchngs v1.12.0](https://nf-co.re/fetchngs/1.12.0/).


### RNA-seq data processing

**Requirements:**

- Nextflow (>= v25.04.04)
- Singularity

The preprocessing pipeline for RNA-seq raw data was adapted from the publicly
available Nextflow pipeline of the Borth group at BOKU:
[github.com/NBorthLab/nf-rnaseq](https://github.com/NBorthLab/nf-rnaseq)
(doi:[10.5281/zenodo.18433643](https://doi.org/10.5281/zenodo.18433643)).
For details on the preparation of the samplesheet CSV file and configuration
options for Nextflow, please find the information in the `nf-rnaseq` repository.
Adaptations to the this workflow were done to accommodate the meta-analysis
setting of this study, i.e. saving results separately per dataset.

The pipeline requires Singularity (or Apptainer, not tested) for running the
software tools and will by default run using SLURM (which is also recommended).
Other execution and containerization options require adaptation of the
`nextflow.config` and/or pipelines modules.

### Analysis

**Requirements:**

- R v4.4.0
- renv
- quarto

The analysis pipeline is written using the `targets` R package and small
Nextflow pipelines for tools outside the R ecosystem (e.g. HOMER, FIMO,
ChromHMM).



