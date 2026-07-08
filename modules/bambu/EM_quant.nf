process BAMBU_EM{
    publishDir "$params.output_dir", mode: 'copy', pattern: 'transcript_counts_singlecell'
    publishDir "$params.output_dir", mode: 'copy', pattern: 'transcript_counts_clusters'
    publishDir "$params.output_dir", mode: 'copy', pattern: 'gene_counts_clusters'
    label "r"
    label "low_cpu"
    label "high_mem"
    label "long"

    input:
    path(quant_data)
    path(extended_annotation)
    tuple val(has_clusters), path(clusters)
    path(genome)

    output:
    path ('transcript_counts_singlecell'), optional: true
    path ('transcript_counts_clusters'), optional: true
    path ('gene_counts_clusters'), optional: true
    path "versions.yml", topic: 'versions'

    script:
    """
    #!/usr/bin/env Rscript
    if ("$params.bambu_path" == "null") { library("bambu") } else { library("devtools"); load_all("$params.bambu_path") }
    library(DropletUtils)
    source(Sys.which("save_counts.R"))

    extendedAnno <- readRDS("$extended_annotation")
    quantData    <- readRDS("$quant_data")
    clusters     <- if ("$has_clusters" == "true") readRDS("$clusters") else NULL
    degBias      <- !is.null(clusters)

    se <- bambu.singlecell(
        reads = NULL,
        annotations = extendedAnno,
        genome = "$genome",
        quantData = quantData,
        assignDist = FALSE,
        ncore = $task.cpus,
        discovery = FALSE,
        quant = TRUE,
        verbose = FALSE,
        opt.em = list(degradationBias = degBias),
        clusters = clusters
    )

    if (is.null(clusters)) {
        save_counts(se, "transcript_counts_singlecell", "Transcript Expression")
    } else {
        save_counts(se, "transcript_counts_clusters", "Transcript Expression")

        se_gene <- transcriptToGeneExpression(se)
        save_counts(se_gene, "gene_counts_clusters")
    }

    writeLines(c('"${task.process}":', paste0('    R: ', R.Version()\$version.string), paste0('    bambu: ', as.character(packageVersion("bambu")))), "versions.yml")
    """
}