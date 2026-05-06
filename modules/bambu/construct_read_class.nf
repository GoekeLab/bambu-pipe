process BAMBU_CONSTRUCT_READ_CLASS{
    publishDir "${params.output_dir}/read_class", mode: 'copy'
    label "r"
    label "low_cpu"
    label "high_mem"
    label "medium"

    input:
    tuple val(sample), path(bam), val(meta)
	path(genome)
	path(bambu_annotation)
    
    output:
    tuple val(sample), path("${sample}_read_class.rds"), val(meta), emit: rds

    script:
    """
    #!/usr/bin/env Rscript
    if ("$params.bambu_path" == "null") { library("bambu") } else { library("devtools"); load_all("$params.bambu_path") }

	annotation <- readRDS("$bambu_annotation")
    readClassFile <- bambu.singlecell(reads = "$bam", annotations = annotation, genome = "$genome", 
        ncore = $task.cpus, discovery = FALSE, quant = FALSE, verbose = FALSE, assignDist = FALSE, 
        processByChromosome = as.logical("$params.process_by_chromosome"), yieldSize = 10000000)
    saveRDS(readClassFile[[1]], "${sample}_read_class.rds") 
    """
}