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
        write10xCounts("transcript_counts_singlecell", assays(se)\$counts, version = "3", gene.type = "Transcript Expression")
        saveRDS(se, "transcript_counts_singlecell/se_transcript_counts_singlecell.rds")
    } else {
        write10xCounts("transcript_counts_clusters", assays(se)\$counts, version = "3", gene.type = "Transcript Expression")
        saveRDS(se, "transcript_counts_clusters/se_transcript_counts_clusters.rds")

        se_gene <- transcriptToGeneExpression(se)
        write10xCounts("gene_counts_clusters", assays(se_gene)\$counts, version = "3")
        saveRDS(se_gene, "gene_counts_clusters/se_gene_counts_clusters.rds")
    }

    writeLines(c('"${task.process}":', paste0('    R: ', R.Version()\$version.string), paste0('    bambu: ', as.character(packageVersion("bambu")))), "versions.yml")
    """
}