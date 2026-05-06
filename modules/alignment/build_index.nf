process MINIMAP_BUILD_INDEX{
    label "minimap2_samtools"
    label "low_cpu"
    label "medium_mem"
    label "short"

    input:
    path(genome)

    output:
    path('ref.mmi'), emit: index
    path "versions.yml", emit: versions

    script:
    """
    minimap2 -k15 -w5 -d ref.mmi $genome # -k and -w flags are used for both splice:hq and splice presets
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version 2>&1)
    END_VERSIONS
    """
}
