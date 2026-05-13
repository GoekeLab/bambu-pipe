process BAMBU_TRANSCRIPT_DISCOVERY{ 
    publishDir "$params.output_dir", mode: 'copy', pattern: 'extended_annotations.gtf'
    publishDir "$params.output_dir", mode: 'copy', pattern: '{se_unique_counts,se_gene_counts}.rds'
    publishDir "$params.output_dir/intermediate_R", mode: 'copy', pattern: '*.rds', enabled: params.save_intermediates
    label "r"
    label "medium_cpu"
    label "high_mem"
    label "medium"

	input:
    tuple val(sample), path(rds_files), val(meta), path(spatial_metadata_files, stageAs: '?/*') // stageAs prevents filename collisions when multiple samples share the same metadata filename
	path(genome)
	path(bambu_annotation)
    val(ndr)
    
	output:
    path ('quant_data.rds'), emit: quant_data
    path ('se_unique_counts.rds'), emit: se_unique_counts
    path ('se_gene_counts.rds'), emit: se_gene_counts
	path ('extended_annotations.rds'), emit: extended_annotations
    path ('extended_annotations.gtf'), emit: extended_annotations_gtf
    path ('sample_names.rds'), emit: sample_names

	script:
	""" 
	#!/usr/bin/env Rscript
    if ("$params.bambu_path" == "null") { library("bambu") } else { library("devtools"); load_all("$params.bambu_path") }

    annotation <- readRDS("$bambu_annotation")
    readClassFile <- strsplit("${rds_files.join(',')}", ",")[[1]]
    sampleNames <- strsplit("${sample.join(',')}", ",")[[1]]
    sampleData <- strsplit("${spatial_metadata_files.join(',')}", ",")[[1]]
    chemistry  <- setNames(strsplit("${meta.collect { m -> m.chemistry }.join(',')}", ",")[[1]], sampleNames)
    technology <- setNames(strsplit("${meta.collect { m -> m.technology }.join(',')}", ",")[[1]], sampleNames)

    # Save sampleNames (required for multi-sample Seurat clustering)
    saveRDS(sampleNames, "sample_names.rds")
    
    # Set sampleData to NA/NULL for non-spatial samples
    containsVisiumStandard <- grepl("visium-v", sampleData)
    sampleData[!containsVisiumStandard] <- NA
    sampleData <- if (all(is.na(sampleData))) NULL else sampleData

    # Transcript discovery
    extendedAnno <- bambu.singlecell(reads = readClassFile, annotations = annotation, genome = "$genome", ncore = $task.cpus, 
    discovery = TRUE, quant = FALSE, verbose = FALSE, assignDist = FALSE, NDR = $ndr)
    saveRDS(extendedAnno, "extended_annotations.rds")
    writeToGTF(extendedAnno, "extended_annotations.gtf")

    # Quantification without EM
    quantData <- bambu.singlecell(reads = readClassFile, annotations = extendedAnno, genome = "$genome", ncore = $task.cpus, 
    discovery = FALSE, quant = FALSE, verbose = FALSE, opt.em = list(degradationBias = FALSE), assignDist = TRUE, sampleData = sampleData)
    saveRDS(quantData, "quant_data.rds")

    # Generate unique counts SE from quantData
    seDiscovery <- generateUniqueCountsSEFromQuantData(quantData, extendedAnno)
    colData(seDiscovery)\$chemistry  <- chemistry[colData(seDiscovery)\$sampleName] # Add chemistry into colData (for subsequent batch correction)
    colData(seDiscovery)\$technology <- technology[colData(seDiscovery)\$sampleName] # Add technology into colData (for subsequent batch correction)
    saveRDS(seDiscovery, "se_unique_counts.rds")

    # Generate gene counts SE from unique counts SE
    seDiscovery.gene <- transcriptToGeneExpression(seDiscovery)
    saveRDS(seDiscovery.gene, "se_gene_counts.rds")
	"""
}