include { EXTRACT_10X_BARCODES } from '../modules/prepare_input/standard/extract_barcodes.nf'
include { EXTRACT_10X_SPATIAL_COORDINATES } from '../modules/prepare_input/standard/extract_spatial_coordinates.nf'
include { BUILD_VISIUM_TISSUE_POSITIONS } from '../modules/prepare_input/standard/build_visium_tissue_positions.nf'

workflow PREPARE_INPUT_STANDARD {
    take:
    ch_rows
    ch_barcode_coordinate_config
    loupe_alignment

    main:
    // parse samplesheet rows into channel of tuples (sample, path, metadata)
    ch_samples = ch_rows.map { row ->
        Validation.validateRow(row, params, log) // validate samplesheet structure and input values
        def sample_path = file(row.path, checkIfExists: true)
        def meta = [chemistry: row.chemistry, technology: row.technology]
        [row.sample, sample_path, meta]
    }

    // validate: for visium runs, ensure only one sample is processed at a time
    ch_samples.collect(flat: false).map { samples -> Validation.validateVisiumSampleCount(samples, loupe_alignment) }

    // extract distinct chemistries from metadata
    ch_unique_chem = ch_samples.map { _sample, _path, meta -> meta.chemistry }.unique()

    // extract 10x barcodes for standard chemistries only (copy files from spaceranger container)
    ch_unique_standard = ch_unique_chem.filter { chem -> params.valid_chemistries.contains(chem) }
    ch_unique_custom   = ch_unique_chem.filter { chem -> !params.valid_chemistries.contains(chem) }
    EXTRACT_10X_BARCODES(ch_unique_standard, ch_barcode_coordinate_config)
    ch_barcodes = EXTRACT_10X_BARCODES.out.barcodes
        .mix(ch_unique_custom.map { chem -> [chem, null] })

    // extract spatial coordinates for visium chemistries only
    ch_unique_visium = ch_unique_chem.filter { chem -> chem.startsWith('visium') }
    EXTRACT_10X_SPATIAL_COORDINATES(ch_unique_visium, ch_barcode_coordinate_config)

    // update sample tuples with the chemistry's barcode whitelist
    ch_keyed_samples = ch_samples.map { sample, path, meta -> [meta.chemistry, sample, path, meta] }
        .combine(ch_barcodes, by: 0)

    // build the tissue positions (spatial metadata) from the Loupe manual alignment file and the spatial
    // coordinate whitelist file. Also generate the in-tissue barcode list which will be used to filter
    // the BAM file.
    if (loupe_alignment != null) {
        BUILD_VISIUM_TISSUE_POSITIONS(EXTRACT_10X_SPATIAL_COORDINATES.out.spatial_coordinates, loupe_alignment)
        ch_tissue_barcodes = BUILD_VISIUM_TISSUE_POSITIONS.out.tissue_barcodes.first()
        ch_updated_samples = ch_keyed_samples
            .combine(BUILD_VISIUM_TISSUE_POSITIONS.out.tissue_positions.first())
            .map { _chem, sample, path, meta, bc, sc ->
                [sample, path, meta + [barcode: bc, spatial_metadata: sc]]
            }
    } else {
        ch_tissue_barcodes = channel.empty()
        ch_updated_samples = ch_keyed_samples
            .map { _chem, sample, path, meta, bc ->
                [sample, path, meta + [barcode: bc, spatial_metadata: null]]
            }
    }

    // split samples by file type
    def fastq_exts = ['.fastq', '.fq', '.fastq.gz', '.fq.gz']
    ch_fastq = ch_updated_samples.filter { _sample, path, _meta -> fastq_exts.any { ext -> path.name.endsWith(ext) } }
    ch_bam = ch_updated_samples.filter { _sample, path, _meta -> path.name.endsWith('.bam') }

    emit:
    fastq = ch_fastq
    bam = ch_bam
    tissue_barcodes = ch_tissue_barcodes
}