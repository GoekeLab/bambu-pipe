process BAMBU_CONSTRUCT_READ_CLASS{
    publishDir "$params.output_dir/intermediate_R/read_class", mode: 'copy', pattern: '*_read_class.rds', enabled: params.save_intermediates
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
    # Rename BAM to val(sample) so bambu derives sampleName from val(sample), regardless of the original filename
    file.symlink("$bam", "${sample}.bam")
    readClassFile <- bambu.singlecell(reads = "${sample}.bam", annotations = annotation, genome = "$genome",
        ncore = $task.cpus, discovery = FALSE, quant = FALSE, verbose = FALSE, assignDist = FALSE, 
        processByChromosome = as.logical("$params.process_by_chromosome"), yieldSize = 10000000)
    saveRDS(readClassFile[[1]], "${sample}_read_class.rds") 
    """
}