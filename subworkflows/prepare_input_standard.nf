include { EXTRACT_10X_BARCODES            } from '../modules/prepare_input_standard/extract_barcodes.nf'
include { EXTRACT_10X_SPATIAL_COORDINATES } from '../modules/prepare_input_standard/extract_spatial_coordinates.nf'

workflow PREPARE_INPUT_STANDARD {
    take:
    ch_rows
    ch_barcode_coordinate_config

    main:
    // parse samplesheet rows into channel of tuples (sample, path, metadata)
    ch_samples = ch_rows.map { row ->
        // validate required columns exist and are non-empty
        ["sample", "path"].each { col ->
            if (!row.containsKey(col))
                error "Samplesheet is missing a required '${col}' column"
            if (!row[col])
                error "A row in the samplesheet has an empty '${col}' value"
        }

        def sample_path = file(row.path, checkIfExists: true)

        // resolve chemistry and technology from row or params, then validate
        def valid_options = [chemistry: params.valid_chemistries, technology: params.valid_technologies]
        def meta = valid_options.collectEntries { col, valid_list ->
            def val = row.containsKey(col) ? row[col] : params[col]
            if (!val)
                error "Sample '${row.sample}' is missing a ${col} — set it in the samplesheet or via params.${col}"
            if (!valid_list.contains(val))
                error "Sample '${row.sample}' has invalid ${col} '${val}' — must be one of: ${valid_list.join(', ')}"
            [col, val]
        }

        [row.sample, sample_path, meta]
    }

    // extract distinct chemistries from metadata
    ch_unique_chem = ch_samples.map { sample, path, meta -> meta.chemistry }.unique()

    // extract 10x barcodes and spatial coordinates for each chemistry (copy files from spaceranger container)
    EXTRACT_10X_BARCODES(ch_unique_chem, ch_barcode_coordinate_config)
    EXTRACT_10X_SPATIAL_COORDINATES(ch_unique_chem, ch_barcode_coordinate_config)

    // update metadata with barcode and spatial coordinate paths
    ch_updated_samples = ch_samples.map { sample, path, meta -> [meta.chemistry, sample, path, meta] }
        .combine(EXTRACT_10X_BARCODES.out.barcodes, by: 0)
        .combine(EXTRACT_10X_SPATIAL_COORDINATES.out.spatial_coordinates, by: 0)
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
    versions = EXTRACT_10X_BARCODES.out.versions.first()
}