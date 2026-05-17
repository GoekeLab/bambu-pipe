process MINIMAP_ALIGNMENT{
    publishDir "${params.output_dir}/bam", mode: 'copy', pattern: '*.bam*'
    label "minimap2_samtools"
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
    path "versions.yml", topic: 'versions'

    script:
    """
    if [[ $meta.technology == "PacBio" ]]; then
        preset="splice:hq"
    else
        preset="splice"
    fi

    minimap2 -ax \$preset -uf --junc-bed $bed -t $task.cpus $ref_mmi $newfastq | \
    samtools sort -@ $task.cpus -o ${sample}_demultiplexed.bam
    samtools index -@ $task.cpus ${sample}_demultiplexed.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version 2>&1)
        samtools: \$(samtools --version 2>&1 | head -1)
    END_VERSIONS
    """
}
