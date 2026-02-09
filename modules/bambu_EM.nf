process BAMBU_EM{
	publishDir "$params.output_dir", mode: 'copy'
	container "ghcr.io/ch99l/bambu-pipe:latest"
    label "medium_cpu"
    label "high_mem"
    label "medium"

    input:
    tuple val(sample), path(rds_files), val(meta) 
    path(quant_data)
    path(extended_annotation)
    path(clusters)
    path(genome)

    output:
    path ('*_se.rds')

    script:
    """
    #!/usr/bin/env Rscript
    #.libPaths("/usr/local/lib/R/site-library")
    library(devtools)
    load_all("$params.bambu_path")
    
    idNames <- strsplit("${sample.join(',')}", ",")[[1]]
	runName = if (length(idNames) == 1) idNames[1] else "combined"
    readClassFile <- strsplit("${rds_files.join(',')}", ",")[[1]]
    extendedAnno <- readRDS("$extended_annotation")
    quantDatas = readRDS("$quant_data")
    clusters = readRDS("$clusters")
	print(clusters)
    degBias = TRUE
    if(is.null(clusters)){degBias = FALSE}

    se = bambu(reads = readClassFile, annotations = extendedAnno, genome = "$genome", quantData = quantDatas, 
    assignDist = FALSE, ncore = $task.cpus, discovery = FALSE, quant = TRUE, demultiplexed = TRUE, 
    verbose = FALSE, opt.em = list(degradationBias = degBias), clusters = clusters)
    saveRDS(se, paste0(runName, "_se.rds"))
    # writeBambuOutput(se, path = ".", prefix = paste0(runName, "_EM_"),outputExtendedAnno = FALSE, outputAll = FALSE, outputBambuModels = FALSE, outputNovelOnly = FALSE)
    """
}