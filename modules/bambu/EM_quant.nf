process BAMBU_EM{
	publishDir "$params.output_dir", mode: 'copy', pattern: '*.rds'
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
    path ('se_transcript_counts_*.rds')
    path ('se_gene_counts_clusters.rds'), optional: true
    path "versions.yml", topic: 'versions'

    script:
    """
    #!/usr/bin/env Rscript
    if ("$params.bambu_path" == "null") { library("bambu") } else { library("devtools"); load_all("$params.bambu_path") }
    
    extendedAnno <- readRDS("$extended_annotation")
    quantData = readRDS("$quant_data")
    clusters = if ("$has_clusters" == "true") readRDS("$clusters") else NULL
    degBias <- !is.null(clusters) 

    se = bambu.singlecell(
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

    # Save transcript counts; for clustered EM, also save gene counts
    if (is.null(clusters)) {
        saveRDS(se, "se_transcript_counts_singlecell.rds")
    } else {
        saveRDS(se, "se_transcript_counts_clusters.rds")
        saveRDS(transcriptToGeneExpression(se), "se_gene_counts_clusters.rds")
    }
    
    writeLines(c('"${task.process}":', paste0('    R: ', R.Version()\$version.string), paste0('    bambu: ', as.character(packageVersion("bambu")))), "versions.yml")
    """
}