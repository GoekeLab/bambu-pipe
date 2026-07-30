process BAMBU_CLUSTER_LEVEL_QUANTIFICATION {
    publishDir "$params.output_dir", mode: 'copy', pattern: 'transcript_counts_clusters'
    publishDir "$params.output_dir", mode: 'copy', pattern: 'gene_counts_clusters'
    label "r"
    label "low_cpu"
    label "high_mem"
    label "long"

    input:
    path(clusters)
    path(quant_data)
    path(extended_annotation)
    path(genome)

    output:
    path ('transcript_counts_clusters')
    path ('gene_counts_clusters')
    path "versions.yml", topic: 'versions'

    script:
    """
    #!/usr/bin/env Rscript
    if ("$params.bambu_path" == "null") { library("bambu") } else { library("devtools"); load_all("$params.bambu_path") }
    library(DropletUtils)
    source(Sys.which("save_counts.R"))

    extendedAnno <- readRDS("$extended_annotation")
    quantData    <- readRDS("$quant_data")
    clusters <- readRDS("$clusters")

    se <- bambu.singlecell(
        reads = quantData,
        output = "clusteredEM",
        annotations = extendedAnno,
        genome = "$genome",
        ncore = $task.cpus,
        verbose = FALSE,
        opt.em = list(degradationBias = TRUE),
        clusters = clusters
    )

    saveCounts(se, "transcript_counts_clusters", "Transcript Expression")

    seGene <- transcriptToGeneExpression(se)
    saveCounts(seGene, "gene_counts_clusters")

    writeLines(c('"${task.process}":', paste0('    R: ', R.Version()\$version.string), paste0('    bambu: ', as.character(packageVersion("bambu")))), "versions.yml")
    """
}
