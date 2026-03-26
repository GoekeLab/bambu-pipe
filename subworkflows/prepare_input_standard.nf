process EXTRACT_10X_BARCODES {
    executor 'local'
    container "quay.io/nf-core/spaceranger:9c5e7dc93c32448e"

    input:
    val(chemistry)
    path(barcode_coordinate_config)

    output:
    tuple val(chemistry), path("${chemistry}_barcode.txt")

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

process EXTRACT_10X_SPATIAL_COORDINATES {
    executor 'local'
    container "quay.io/nf-core/spaceranger:9c5e7dc93c32448e"

    input:
    val(chemistry)
    path(barcode_coordinate_config)

    output:
    tuple val(chemistry), path("${chemistry}_spatial_coordinates.txt")

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

workflow PREPARE_INPUT_STANDARD {
    take:
    ch_input
    ch_barcode_coordinate_config

    main:
    // read samplesheet CSV into channel of tuples (sample, path, metadata)
    ch_samples = ch_input.splitCsv(header:true, sep:',')
    .map { row ->
        // validate samplesheet format and required columns
        if (!row.containsKey("sample"))
            error "Samplesheet is missing a required 'sample' column"
        if (!row.sample)
            error "A row in the samplesheet has an empty sample name"

        if (!row.containsKey("path"))
            error "Samplesheet is missing a required 'path' column"
        if (!row.path)
            error "Sample '${row.sample}' has an empty 'path' value"

        def sample_path = file(row.path, checkIfExists: true)  // check if file exists at path specified
        
        // validate chemistry and technology
        def chemistry = row.containsKey("chemistry") ? row.chemistry : params.chemistry
        def technology = row.containsKey("technology") ? row.technology : params.technology

        if (!chemistry)
            error "Sample '${row.sample}' is missing a chemistry — set it in the samplesheet or via params.chemistry"
        if (!params.valid_chemistries.contains(chemistry))
            error "Sample '${row.sample}' has invalid chemistry '${chemistry}' — must be one of: ${params.valid_chemistries.join(', ')}"

        if (!technology)
            error "Sample '${row.sample}' is missing a technology — set it in the samplesheet or via params.technology"
        if (!params.valid_technologies.contains(technology))
            error "Sample '${row.sample}' has invalid technology '${technology}' — must be one of: ${params.valid_technologies.join(', ')}"

        def meta = [chemistry: chemistry, technology: technology]

        [row.sample, sample_path, meta]
    }

    // extract distinct chemistries from metadata
    ch_unique_chem = ch_samples.map { sample, path, meta -> meta.chemistry }.unique()

    // extract 10x barcodes and spatial coordinates for each chemistry (copy files from spaceranger container)
    EXTRACT_10X_BARCODES(ch_unique_chem, ch_barcode_coordinate_config)
    EXTRACT_10X_SPATIAL_COORDINATES(ch_unique_chem, ch_barcode_coordinate_config)

    // update metadata with barcode and spatial coordinate paths
    ch_updated_samples = ch_samples.map { sample, path, meta -> [meta.chemistry, sample, path, meta] }
        .combine(EXTRACT_10X_BARCODES.out, by: 0)
        .combine(EXTRACT_10X_SPATIAL_COORDINATES.out, by: 0)
        .map { chem, sample, path, meta, bc, sc ->
            def updated_meta = meta + [
                barcode: bc, 
                spatial_metadata: sc
            ]
            return [sample, path, updated_meta]
    }

    // split samples by file type
    def fastq_exts = ['.fastq', '.fq', '.fastq.gz', '.fq.gz']
    ch_fastq = ch_updated_samples.filter { sample, path, meta -> fastq_exts.any { ext -> path.name.endsWith(ext) } }
    ch_bam = ch_updated_samples.filter { sample, path, meta -> path.name.endsWith('.bam') }
    ch_rds = ch_updated_samples.filter { sample, path, meta -> path.name.endsWith('.rds') }

    emit:
    fastq = ch_fastq
    bam = ch_bam
    rds = ch_rds
}