process FILTER_BARCODED_BAM {
    publishDir "$params.output_dir/intermediate_bam", mode: 'copy', pattern: '*_filtered.bam*', enabled: params.save_intermediates
    label "minimap2_samtools"
    label "medium_cpu"
    label "medium_mem"
    label "medium"

    input:
    tuple val(sample), path(bam), val(meta)
    path(barcodes)

    output:
    tuple val(sample), path("${sample}_filtered.bam"), val(meta), emit: bam
    path("${sample}_filtered.bam.bai")
    path "versions.yml", topic: 'versions'

    script:
    """
    samtools view -@ $task.cpus -D CB:$barcodes -o ${sample}_filtered.bam $bam
    samtools index -@ $task.cpus ${sample}_filtered.bam

    if [[ \$(samtools view -c -@ $task.cpus ${sample}_filtered.bam) -eq 0 ]]; then
        echo "No reads remain after barcode filtering -- check that the BAM carries CB tags matching $barcodes" >&2
        exit 1
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version 2>&1 | head -1)
    END_VERSIONS
    """
}
