process CONVERT_BARCODE_MAPPINGS {
    publishDir "$params.output_dir/intermediate_visium_hd", mode: 'copy', pattern: '*.csv', enabled: params.save_intermediates
    label "pyarrow_pandas"
    label "low_cpu"
    label "low_mem"
    label "short"

    input:
    val(resolution)
    path(barcode_mappings)

    output:
    tuple val(resolution), path("barcode_mappings_${resolution}.csv"), emit: csv
    path "versions.yml", topic: 'versions'

    script:
    """
    visium_hd_convert_barcode_mappings.py $barcode_mappings ${resolution} barcode_mappings_${resolution}.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1)
        pandas: \$(python3 -c 'import pandas; print(pandas.__version__)')
    END_VERSIONS
    """
}
