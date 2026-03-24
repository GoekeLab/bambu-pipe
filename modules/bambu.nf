process BAMBU{ 
    publishDir "$params.output_dir", mode: 'copy', pattern: '*extended_annotations.rds'
    container "ghcr.io/ch99l/bambu-pipe-r:latest"
    label "medium_cpu"
    label "high_mem"
    label "medium"

	input:
    tuple val(sample), path(rds_files), val(meta), path(spatial_metadata_files)
	path(genome)
	path(bambu_annotation)
    val(ndr)
    val(run_clustering)
    
	output: 
    path ('*quantData.rds'), emit: quant_data
	path ('*extended_annotations.rds'), emit: extended_annotations
    path ('*_clusters.rds'), emit: clusters

	script:
	""" 
	#!/usr/bin/env Rscript
    #.libPaths("/usr/local/lib/R/site-library")
    library("bambu")

	## Transcript discovery and quantification without EM
    idNames <- strsplit("${sample.join(',')}", ",")[[1]]
    runName = if (length(idNames) == 1) idNames[1] else "combined"
    annotation <- readRDS("$bambu_annotation")
    readClassFile <- strsplit("${rds_files.join(',')}", ",")[[1]]
    sampleData <- strsplit("${spatial_metadata_files.join(',')}", ",")[[1]]
    
    # Set sampleData to NA/NULL for non-spatial samples
    if (!as.logical("$params.visium_hd")) {
        contains_visium_standard <- grepl("visium-v", sampleData)
        sampleData[!contains_visium_standard] <- NA
        sampleData <- if (all(is.na(sampleData))) NULL else sampleData
    }

    # Transcript discovery
    extendedAnno = bambu(reads = readClassFile, annotations = annotation, genome = "$genome", ncore = $task.cpus, 
    discovery = TRUE, quant = FALSE, demultiplexed = TRUE, verbose = FALSE, assignDist = FALSE, NDR = $ndr)
    saveRDS(extendedAnno, paste0(runName, "_extended_annotations.rds"))

    # Quantification without EM
    se = bambu(reads = readClassFile, annotations = extendedAnno, genome = "$genome", ncore = $task.cpus, 
    discovery = FALSE, quant = FALSE, demultiplexed = TRUE, verbose = FALSE, 
    opt.em = list(degradationBias = FALSE), assignDist = TRUE, sampleData = sampleData)
    saveRDS(se, paste0(runName, "_quantData.rds"))

    # Seurat Clustering (if no clustering provided, automatically cluster)
	path <- Sys.getenv("PATH") |> strsplit(":")
    bin_path <- tail(path[[1]], n=1)
    clusters = NULL
    if(as.logical("$run_clustering")){
        clusters = list()
        cellMixs = list()
        source(file.path(bin_path,"/utilityFunctions.R"))
        for(quantData in se){
            quantData.gene = transcriptToGeneExpression(quantData)
            for(sample in unique(colData(quantData)\$sampleName)){
                i = which(colData(quantData)\$sampleName == sample)
                counts = assays(quantData.gene)\$counts[,i]
                cellMix = clusterCells(counts, resolution = $params.resolution)
                x = setNames(names(cellMix@active.ident), cellMix@active.ident)
                names(x) = paste0(sample,"_",names(x))
                clusters = c(clusters, splitAsList(unname(x), names(x)))
                cellMixs = c(cellMixs, cellMix)
            }
        }
        saveRDS(cellMixs, paste0(runName, "_cellMixs.rds"))
    }
    
    saveRDS(clusters, paste0(runName, "_clusters.rds"))
	"""
}