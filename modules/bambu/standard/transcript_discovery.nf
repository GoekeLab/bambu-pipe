process BAMBU_TRANSCRIPT_DISCOVERY{
    publishDir "$params.output_dir", mode: 'copy', pattern: 'extended_annotations.gtf'
    publishDir "$params.output_dir", mode: 'copy', pattern: 'unique_counts'
    publishDir "$params.output_dir", mode: 'copy', pattern: 'gene_counts'
    publishDir "$params.output_dir/intermediate_R", mode: 'copy', pattern: '*.rds', enabled: params.save_intermediates
    label "r"
    label "medium_cpu"
    label "high_mem"
    label "medium"

	input:
    tuple val(sample), path(rds_files), val(meta), path(spatial_metadata_files)
	path(genome)
	path(bambu_annotation)
    val(ndr)

	output:
    path ('quant_data.rds'), emit: quant_data
    path ('unique_counts/se_unique_counts.rds')
    path ('gene_counts/se_gene_counts.rds'), emit: se_gene_counts
	path ('extended_annotations.rds'), emit: extended_annotations
    path ('extended_annotations.gtf')
    path ('unique_counts')
    path ('gene_counts')
    path "versions.yml", topic: 'versions'

	script:
	"""
	#!/usr/bin/env Rscript
    if ("$params.bambu_path" == "null") { library("bambu") } else { library("devtools"); load_all("$params.bambu_path") }
    library(DropletUtils)
    source(Sys.which("bambu_discovery.R"))
    source(Sys.which("save_counts.R"))

    annotation <- readRDS("$bambu_annotation")
    sampleNames <- strsplit("${sample.join(',')}", ",")[[1]]
    readClassFile <- setNames(strsplit("${rds_files.join(',')}", ",")[[1]], sampleNames)
    sampleData <- strsplit("${spatial_metadata_files.join(',')}", ",")[[1]]
    chemistry  <- setNames(strsplit("${meta.collect { m -> m.chemistry }.join(',')}", ",")[[1]], sampleNames)
    technology <- setNames(strsplit("${meta.collect { m -> m.technology }.join(',')}", ",")[[1]], sampleNames)
    sampleData <- if (any(startsWith(chemistry, "visium-v"))) sampleData else NULL # Add spatial metadata for visium samples

    result <- bambuDiscovery(reads = readClassFile, annotation = annotation, genome = "$genome",
        ncore = $task.cpus, ndr = $ndr, sampleData = sampleData)
    saveRDS(result\$extendedAnno, "extended_annotations.rds")
    writeToGTF(result\$extendedAnno, "extended_annotations.gtf")
    saveRDS(result\$quantData, "quant_data.rds")

    seDiscovery <- result\$seDiscovery
    colData(seDiscovery)\$chemistry  <- unname(chemistry[colData(seDiscovery)\$sampleName]) # Add chemistry into colData (for subsequent batch correction)
    colData(seDiscovery)\$technology <- unname(technology[colData(seDiscovery)\$sampleName]) # Add technology into colData (for subsequent batch correction)
    saveCounts(seDiscovery, "unique_counts", "Transcript Expression")

    # Generate gene counts SE from unique counts SE
    seDiscoveryGene <- transcriptToGeneExpression(seDiscovery)
    saveCounts(seDiscoveryGene, "gene_counts")

    writeLines(c('"${task.process}":', paste0('    R: ', R.Version()\$version.string), paste0('    bambu: ', as.character(packageVersion("bambu")))), "versions.yml")
	"""
}