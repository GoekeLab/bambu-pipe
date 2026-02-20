process MINIMAP_BUILD_INDEX{
    container "ghcr.io/ch99l/bambu-pipe-alignment:latest"
    label "low_cpu"
    label "medium_mem"
    label "short"

    input: 
    path(genome)
    val(fastq_count)

    when: fastq_count > 0 // only build index if there are fastq samples to process

    output: 
    path('ref.mmi')

    script:
    """ 
    minimap2 -k15 -w5 -d ref.mmi $genome # -k and -w flags are used for both splice:hq and splice presets
    """
}

process PAFTOOLS_GFF2BED {
    container "ghcr.io/ch99l/bambu-pipe-alignment:latest"
    label "low_cpu"
    label "low_mem"
    label "short"

    input: 
    path(gtf)
    val(fastq_count)

    when: fastq_count > 0 // only convert annotation if there are fastq samples to process

    output: 
    path('anno.bed')

    script:
    """ 
    paftools.js gff2bed $gtf > anno.bed
    """
}

process MINIMAP_ALIGNMENT{
    publishDir "$params.output_dir", mode: 'copy' 
    container "ghcr.io/ch99l/bambu-pipe-alignment:latest"
    label "high_cpu"
    label "high_mem"
    label "long"
	
	input: 
	tuple val(sample), path(newfastq), val(meta)
    path(ref_mmi)
    path(bed)

	output: 
	tuple val(sample), path("${sample}_demultiplexed.bam"), val(meta)

	script:
	""" 
	if [[ $meta.technology == "PacBio" ]]; then 
		preset="splice:hq"
	else
		preset="splice"  
	fi

    minimap2 -ax \$preset -uf --junc-bed $bed -t $task.cpus $ref_mmi $newfastq > demultiplexed.sam
	samtools sort -@ $task.cpus demultiplexed.sam -o ${sample}_demultiplexed.bam 
	samtools index -@ $task.cpus ${sample}_demultiplexed.bam 

    rm demultiplexed.sam
	"""
}

workflow ALIGNMENT {
    take: 
    ch_unaligned_fastq
    ch_genome
    ch_annotation
    ch_input_fastq_count // fastq count is used to ensure paftools and minimap build index are skipped when there are no fastq samples

    main:
    // Build minimap2 index based on reference genome
    MINIMAP_BUILD_INDEX(ch_genome, ch_input_fastq_count)

    // Convert gtf/gff annotation to bed format
    PAFTOOLS_GFF2BED(ch_annotation, ch_input_fastq_count)

    // Minimap alignment
    MINIMAP_ALIGNMENT(ch_unaligned_fastq, MINIMAP_BUILD_INDEX.out, PAFTOOLS_GFF2BED.out)
    
    emit:
    bam = MINIMAP_ALIGNMENT.out
}