process MINIMAP_BUILD_INDEX{
    container "community.wave.seqera.io/library/minimap2_samtools:b09096fc890429ce"
    label "low_cpu"
    label "medium_mem"
    label "short"

    input:
    path(genome)

    output:
    path('ref.mmi')

    script:
    """ 
    minimap2 -k15 -w5 -d ref.mmi $genome # -k and -w flags are used for both splice:hq and splice presets
    """
}

process PAFTOOLS_GFF2BED {
    container "community.wave.seqera.io/library/minimap2_samtools:b09096fc890429ce"
    label "low_cpu"
    label "low_mem"
    label "short"

    input:
    path(gtf)

    output:
    path('anno.bed')

    script:
    """ 
    paftools.js gff2bed $gtf > anno.bed
    """
}

process MINIMAP_ALIGNMENT{
    publishDir "${params.output_dir}/bam", mode: 'copy' 
    container "community.wave.seqera.io/library/minimap2_samtools:b09096fc890429ce"
    label "high_cpu"
    label "high_mem"
    label "long"
	
	input: 
	tuple val(sample), path(newfastq), val(meta)
    path(ref_mmi)
    path(bed)

	output: 
	tuple val(sample), path("${sample}_demultiplexed.bam"), val(meta), emit: bam
    path("${sample}_demultiplexed.bam.bai")

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

    main:
    // ch_gate emits one item if fastq channel is non-empty, else emits nothing. Used to prevent MINIMAP_BUILD_INDEX and PAFTOOLS_GFF2BED 
    // from running when there are no fastq samples to process
    ch_gate = ch_unaligned_fastq.first()

    // Build minimap2 index based on reference genome
    MINIMAP_BUILD_INDEX(ch_genome.combine(ch_gate).map { g, _ -> g })

    // Convert gtf/gff annotation to bed format
    PAFTOOLS_GFF2BED(ch_annotation.combine(ch_gate).map { a, _ -> a })

    // Minimap alignment
    MINIMAP_ALIGNMENT(ch_unaligned_fastq, MINIMAP_BUILD_INDEX.out, PAFTOOLS_GFF2BED.out)
    
    emit:
    bam = MINIMAP_ALIGNMENT.out.bam
}