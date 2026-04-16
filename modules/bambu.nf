process BAMBU{ 
    publishDir "$params.output_dir", mode: 'copy', pattern: '*extended_annotations.gtf'
    container "ghcr.io/ch99l/bambu-pipe-r:latest"
    label "medium_cpu"
    label "high_mem"
    label "medium"

	input:
    tuple val(sample), path(rds_files), val(meta), path(spatial_metadata_files, stageAs: '?/*')
	path(genome)
	path(bambu_annotation)
    val(ndr)
    
	output:
    path ('*quantData.rds'), emit: quant_data
	path ('*extended_annotations.rds'), emit: extended_annotations
    path ('*extended_annotations.gtf'), emit: extended_annotations_gtf

	script:
	""" 
	#!/usr/bin/env Rscript
    #.libPaths("/usr/local/lib/R/site-library")
    if ("$params.bambu_path" == "null") { library("bambu") } else { library("devtools"); load_all("$params.bambu_path") }

	## Transcript discovery and quantification without EM
    idNames <- strsplit("${sample.join(',')}", ",")[[1]]
    runName = if (length(idNames) == 1) idNames[1] else "combined"
    annotation <- readRDS("$bambu_annotation")
    readClassFile <- strsplit("${rds_files.join(',')}", ",")[[1]]
    sampleData <- strsplit("${spatial_metadata_files.join(',')}", ",")[[1]]
    
    # Set sampleData to NA/NULL for non-spatial samples
    contains_visium_standard <- grepl("visium-v", sampleData)
    sampleData[!contains_visium_standard] <- NA
    sampleData <- if (all(is.na(sampleData))) NULL else sampleData

    # Transcript discovery
    extendedAnno = bambu(reads = readClassFile, annotations = annotation, genome = "$genome", ncore = $task.cpus, 
    discovery = TRUE, quant = FALSE, demultiplexed = TRUE, verbose = FALSE, assignDist = FALSE, NDR = $ndr)
    saveRDS(extendedAnno, paste0(runName, "_extended_annotations.rds"))
    writeToGTF(extendedAnno, paste0(runName, "_extended_annotations.gtf"))

    # Quantification without EM
    se = bambu(reads = readClassFile, annotations = extendedAnno, genome = "$genome", ncore = $task.cpus, 
    discovery = FALSE, quant = FALSE, demultiplexed = TRUE, verbose = FALSE, 
    opt.em = list(degradationBias = FALSE), assignDist = TRUE, sampleData = sampleData)
    saveRDS(se, paste0(runName, "_quantData.rds"))
	"""
}