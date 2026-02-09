process BAMBU_CONSTRUCT_READ_CLASS{
    container "ghcr.io/ch99l/bambu-pipe:latest"
    label "low_cpu"
    label "high_mem"
    label "medium"

    input:
    tuple val(sample), path(bam), val(meta)
	path(genome)
	path(bambu_annotation)
    
    output:
    tuple val(sample), path("${sample}_readClassFile.rds"), val(meta)

    script:
    """
    #!/usr/bin/env Rscript
    
    library("devtools")
    load_all("$params.bambu_path")

    if ("$meta.barcode_map" == "true") {
    demultiplexed = TRUE} else {
    demultiplexed = "$meta.barcode_map"}

	annotation <- readRDS("$bambu_annotation")
    readClassFile = bambu(reads = "$bam", annotations = annotation, genome = "$genome", 
        ncore = $task.cpus, discovery = FALSE, quant = FALSE, demultiplexed = demultiplexed, 
        verbose = FALSE, assignDist = FALSE, processByChromosome = as.logical("$params.process_by_chromosome"), 
        yieldSize = 10000000, dedupUMI = as.logical("$params.deduplicate_umis"))
    saveRDS(readClassFile[[1]], "${sample}_readClassFile.rds") 
    """
}