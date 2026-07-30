process BAMBU_SPOT_LEVEL_QUANTIFICATION {
    publishDir "$params.output_dir", mode: 'copy', pattern: 'transcript_counts_*'
    label "r"
    label "low_cpu"
    label "high_mem"
    label "long"

    input:
    tuple val(resolution), path(tissue_positions), path(spot_mappings)
    path(quant_data)
    path(extended_annotation)
    path(genome)

    output:
    path ("transcript_counts_${resolution}")
    path "versions.yml", topic: 'versions'

    script:
    """
    #!/usr/bin/env Rscript
    if ("$params.bambu_path" == "null") { library("bambu") } else { library("devtools"); load_all("$params.bambu_path") }
    library(DropletUtils)
    library(dplyr)
    source(Sys.which("save_counts.R"))

    extendedAnno <- readRDS("$extended_annotation")
    quantData    <- readRDS("$quant_data")

    # For bin-level quantification, we need to use the clusters argument to collapse the quantData object
    # to the desired resolution, since quantData is generated at 2um resolution.
    if ("$resolution" == "002um") {
        clusters <- NULL
    } else {
        # bambu prefixes the sample name onto each cluster label, so pass the bare bin
        spotMappings <- read.csv("$spot_mappings", stringsAsFactors = FALSE)
        clusters     <- setNames(spotMappings\$bin, spotMappings\$sample_barcode)
    }

    se <- bambu.singlecell(
        reads = quantData,
        output = if (is.null(clusters)) "EM" else "clusteredEM",
        annotations = extendedAnno,
        genome = "$genome",
        ncore = $task.cpus,
        verbose = FALSE,
        opt.em = list(degradationBias = FALSE),
        clusters = clusters
    )

    # At the bin-level resolutions bambu builds a fresh colData, so we need to reattach the bin metadata.
    if (!is.null(clusters)) {
        tissuePositions <- read.csv("$tissue_positions", stringsAsFactors = FALSE)
        binColData      <- as.data.frame(colData(se)) %>%
            left_join(tissuePositions, by = c("cluster" = "barcode")) %>%
            rename(barcode = cluster)
        colData(se) <- DataFrame(binColData, row.names = colnames(se))
    }

    saveCounts(se, "transcript_counts_$resolution", "Transcript Expression")

    writeLines(c('"${task.process}":', paste0('    R: ', R.Version()\$version.string), paste0('    bambu: ', as.character(packageVersion("bambu")))), "versions.yml")
    """
}
