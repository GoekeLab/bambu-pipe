process SPOT_BIN_MAPPINGS {
    publishDir "$params.output_dir/intermediate_visium_hd", mode: 'copy', pattern: '*.csv.gz', enabled: params.save_intermediates
    label "pyarrow_pandas"
    label "low_cpu"
    label "low_mem"
    label "short"

    input:
    tuple val(resolution), path(barcode_mappings)
    path(barcodes_002um)
    val(sample)

    output:
    tuple val(resolution), path("spot_bin_mappings_${resolution}.csv.gz"), emit: csv
    path "versions.yml", topic: 'versions'

    script:
    """
    visium_hd_spot_bin_mappings.py $barcode_mappings $barcodes_002um $sample spot_bin_mappings_${resolution}.csv.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1)
        pandas: \$(python3 -c 'import pandas; print(pandas.__version__)')
    END_VERSIONS
    """
}
