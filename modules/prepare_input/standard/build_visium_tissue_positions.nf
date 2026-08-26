process BUILD_VISIUM_TISSUE_POSITIONS {
    publishDir "$params.output_dir/intermediate_visium", mode: 'copy', pattern: 'tissue_positions.csv', enabled: params.save_intermediates
    publishDir "$params.output_dir/intermediate_visium", mode: 'copy', pattern: 'tissue_barcodes.txt', enabled: params.save_intermediates
    label "pyarrow_pandas"
    label "low_cpu"
    label "low_mem"
    label "short"

    input:
    path(spatial_coordinates)
    path(loupe_alignment)

    output:
    path('tissue_positions.csv'), emit: tissue_positions
    path('tissue_barcodes.txt'),  emit: tissue_barcodes
    path "versions.yml", topic: 'versions'

    script:
    """
    build_visium_tissue_positions.py $loupe_alignment $spatial_coordinates tissue_positions.csv tissue_barcodes.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1)
    END_VERSIONS
    """
}
