process EXTRACT_10X_SPATIAL_COORDINATES {
    label "spaceranger"
    executor 'local'

    input:
    val(chemistry)
    path(barcode_coordinate_config)

    output:
    tuple val(chemistry), path("${chemistry}_spatial_coordinates.txt"), emit: spatial_coordinates

    script:
    """
    # extract spatial coordinate file path from config csv
    IFS=',' read -r _ _ sc_filename < <(awk -F',' -v chem=$chemistry '\$1 == chem' $barcode_coordinate_config)

    cp $params.cellranger_dir/\$sc_filename ./${chemistry}_spatial_coordinates.txt
    sed -i '1ibarcode\tx_coordinate\ty_coordinate' ./${chemistry}_spatial_coordinates.txt
    """
}
