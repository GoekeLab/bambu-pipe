process AGGREGATE_BINS_VISIUM_HD {
    publishDir "$params.output_dir/unique_counts", mode: 'copy', pattern: 'unique_counts_*'
    publishDir "$params.output_dir/gene_counts", mode: 'copy', pattern: 'gene_counts_*'
    label "r"
    label "low_cpu"
    label "medium_mem"
    label "short"

    input:
    tuple val(resolution), path(tissue_positions), path(spot_mappings)
    path(se_unique_002um)

    output:
    tuple val(resolution), path("unique_counts_${resolution}"), emit: unique_counts
    tuple val(resolution), path("gene_counts_${resolution}"),   emit: gene_counts
    tuple val(resolution), path("gene_counts_${resolution}/se_gene_counts_${resolution}.rds"), emit: se_gene_counts
    path "versions.yml", topic: 'versions'

    script:
    """
    #!/usr/bin/env Rscript
    if ("$params.bambu_path" == "null") { library("bambu") } else { library("devtools"); load_all("$params.bambu_path") }
    library(DropletUtils)
    library(Matrix)
    library(dplyr)
    source(Sys.which("save_counts.R"))
    source(Sys.which("visium_hd_aggregate_resolution.R"))

    unique002um     <- readRDS("$se_unique_002um")
    tissuePositions <- read.csv("$tissue_positions", stringsAsFactors = FALSE)
    spotMappings    <- read.csv("$spot_mappings", stringsAsFactors = FALSE)

    # aggregate transcript and gene counts
    uniqueAgg <- aggregateResolution(unique002um, tissuePositions, spotMappings)
    geneAgg   <- transcriptToGeneExpression(uniqueAgg)

    saveCounts(uniqueAgg, "unique_counts_$resolution", "Transcript Expression")
    saveCounts(geneAgg,   "gene_counts_$resolution")

    writeLines(c('"${task.process}":', paste0('    R: ', R.Version()\$version.string), paste0('    bambu: ', as.character(packageVersion("bambu")))), "versions.yml")
    """
}
