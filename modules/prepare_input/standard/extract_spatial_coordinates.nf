process EXTRACT_10X_SPATIAL_COORDINATES {
    label "spaceranger"
    executor 'local'

    input:
    val(chemistry)
    path(barcode_coordinate_config)

    output:
    path("spatial_coordinates.txt"), emit: spatial_coordinates

    script:
    """
    # extract spatial coordinate file path from config csv
    IFS=',' read -r _ _ sc_filename < <(awk -F',' -v chem=$chemistry '\$1 == chem' $barcode_coordinate_config)

    cp $params.cellranger_dir/\$sc_filename ./spatial_coordinates.txt
    """
}
