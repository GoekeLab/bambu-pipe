process PAFTOOLS_GFF2BED {
    label "minimap2_samtools"
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
