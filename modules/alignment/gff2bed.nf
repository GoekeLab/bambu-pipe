process PAFTOOLS_GFF2BED {
    label "minimap2_samtools"
    label "low_cpu"
    label "low_mem"
    label "short"

    input:
    path(gtf)

    output:
    path('anno.bed'), emit: bed
    path "versions.yml", topic: 'versions'

    script:
    """
    paftools.js gff2bed $gtf > anno.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version 2>&1)
    END_VERSIONS
    """
}
