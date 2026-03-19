process SALMON_QUANT {

    tag "$meta.id"
    label 'big'

    container "https://depot.galaxyproject.org/singularity/salmon:1.10.3--h6dccd9a_2"
}
