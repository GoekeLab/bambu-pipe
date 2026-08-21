process EXTRACT_10X_BARCODES {
    label "spaceranger"
    executor 'local'

    input:
    val(chemistry)
    path(barcode_coordinate_config)

    output:
    tuple val(chemistry), path("${chemistry}_barcode.txt"), emit: barcodes

    script:
    """
    # extract 10x barcode file path from config csv
    IFS=',' read -r _ bc_filename _ < <(awk -F',' -v chem=$chemistry '\$1 == chem' $barcode_coordinate_config)

    # extract 10x barcode file from spaceranger container (used for subsequent processes)
    if [[ \$bc_filename == *.gz ]]; then
        gunzip -c $params.cellranger_dir/\$bc_filename > ./${chemistry}_barcode.txt
    else
        cp $params.cellranger_dir/\$bc_filename ./${chemistry}_barcode.txt
    fi
    """
}
