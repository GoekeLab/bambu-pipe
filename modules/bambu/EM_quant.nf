process BAMBU_EM{
	publishDir "$params.output_dir", mode: 'copy'
    label "r"
    label "low_cpu"
    label "high_mem"
    label "long"

    input:
    path(quant_data)
    path(extended_annotation)
    path(clusters)
    path(genome)

    output:
    path ('se.rds')

    script:
    """
    #!/usr/bin/env Rscript
    if ("$params.bambu_path" == "null") { library("bambu") } else { library("devtools"); load_all("$params.bambu_path") }
    
    extendedAnno <- readRDS("$extended_annotation")
    quantData = readRDS("$quant_data")
    clusters = readRDS("$clusters")
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
    saveRDS(se, "se.rds")
    """
}