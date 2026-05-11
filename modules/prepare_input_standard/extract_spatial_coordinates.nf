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

    # extract spatial coordinate file from spaceranger container (used for subsequent processes)
    # TODO: All visium coordinates file are not gzipped. In the future if there are gzipped files, include decompression)
    if [[ $chemistry == visium* ]]; then
        cp $params.cellranger_dir/\$sc_filename ./${chemistry}_spatial_coordinates.txt
        sed -i '1ibarcode\tx_coordinate\ty_coordinate' ./${chemistry}_spatial_coordinates.txt # adding header to spatial coordinates file as it is required for bambu
    else
        touch ./${chemistry}_spatial_coordinates.txt # create empty file for non-visium chemistries
    fi
    """
}
