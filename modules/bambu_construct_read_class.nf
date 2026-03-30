process BAMBU_CONSTRUCT_READ_CLASS{
    container "ghcr.io/ch99l/bambu-pipe-r:latest"
    publishDir "${params.output_dir}/read_class", mode: 'copy'
    label "low_cpu"
    label "high_mem"
    label "medium"

    input:
    tuple val(sample), path(bam), val(meta)
	path(genome)
	path(bambu_annotation)
    
    output:
    tuple val(sample), path("${sample}_readClassFile.rds"), val(meta), emit: rds

    script:
    """
    #!/usr/bin/env Rscript
    if ("$params.bambu_path" == "null") { library("bambu") } else { library("devtools"); load_all("$params.bambu_path") }

	annotation <- readRDS("$bambu_annotation")
    readClassFile = bambu(reads = "$bam", annotations = annotation, genome = "$genome", 
        ncore = $task.cpus, discovery = FALSE, quant = FALSE, demultiplexed = TRUE, 
        verbose = FALSE, assignDist = FALSE, processByChromosome = as.logical("$params.process_by_chromosome"), 
        yieldSize = 10000000, dedupUMI = as.logical("$params.deduplicate_umis"))
    saveRDS(readClassFile[[1]], "${sample}_readClassFile.rds") 
    """
}