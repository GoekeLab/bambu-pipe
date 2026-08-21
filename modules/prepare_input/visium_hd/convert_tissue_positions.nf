process CONVERT_TISSUE_POSITIONS {
    publishDir "$params.output_dir/intermediate_visium_hd", mode: 'copy', pattern: '*.csv', enabled: params.save_intermediates
    publishDir "$params.output_dir/intermediate_visium_hd", mode: 'copy', pattern: 'barcodes_in_tissue.txt', enabled: params.save_intermediates
    label "pyarrow_pandas"
    label "low_cpu"
    label "low_mem"
    label "short"

    input:
    tuple val(resolution), path(tissue_positions)

    output:
    tuple val(resolution), path("tissue_positions_${resolution}.csv"), emit: csv
    path("barcodes_in_tissue.txt"), optional: true, emit: barcodes
    path "versions.yml", topic: 'versions'

    script:
    """
    visium_hd_convert_tissue_positions.py $tissue_positions ${resolution} tissue_positions_${resolution}.csv barcodes_in_tissue.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1)
        pandas: \$(python3 -c 'import pandas; print(pandas.__version__)')
    END_VERSIONS
    """
}
