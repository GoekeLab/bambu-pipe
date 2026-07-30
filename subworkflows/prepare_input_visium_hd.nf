include { FILTER_BARCODED_BAM }      from '../modules/prepare_input/visium_hd/filter_barcoded_bam.nf'
include { CONVERT_BARCODE_MAPPINGS } from '../modules/prepare_input/visium_hd/convert_barcode_mappings.nf'
include { CONVERT_TISSUE_POSITIONS } from '../modules/prepare_input/visium_hd/convert_tissue_positions.nf'

workflow PREPARE_INPUT_VISIUM_HD {
    take:
    ch_rows  // raw samplesheet rows

    main:
    // Visium HD: single sample, starting from a pre-aligned, barcode-tagged BAM file
    ch_rows.collect(flat: false).map { rows -> Validation.validateVisiumHDRows(rows) }

    ch_sample = ch_rows.map { row ->
        def sample_path = file(row.path, checkIfExists: true)
        def meta = [chemistry: 'visium-hd', technology: 'NA']
        [row.sample, sample_path, meta]
    }

    // parse the bins samplesheet into one [resolution, tissue_positions] tuple per row (includes mandatory 2um base)
    ch_resolutions = channel.of(params.bins)
        .splitCsv(header: true, sep: ',')
        .map { row ->
            // Convert integer values in the resolution column of the CSV into Spaceranger resolution, e.g. 8 -> "008um"
            def resolution = String.format('%03dum', row.resolution as Integer)
            def tissue_positions = file(row.tissue_positions, checkIfExists: true)
            [resolution, tissue_positions]
        }
    ch_barcode_mappings = channel.value(params.barcode_mappings)

    // convert every resolution's tissue positions parquet to CSV; the 2um task also extracts the in-tissue barcode list
    CONVERT_TISSUE_POSITIONS(ch_resolutions)

    // the 2um in-tissue barcodes are used to filter out-of-tissue reads from the BAM file
    FILTER_BARCODED_BAM(ch_sample, CONVERT_TISSUE_POSITIONS.out.barcodes)

    // convert the barcode mappings parquet to one 'barcode,bin' CSV per bin resolution level (e.g., 8um/16um)
    ch_bin_resolutions = ch_resolutions
        .map { resolution, _tissue_positions -> resolution }
        .filter { resolution -> resolution != '002um' } 
    CONVERT_BARCODE_MAPPINGS(ch_bin_resolutions, ch_barcode_mappings)

    // extract the 2um tissue position separately since transcript discovery is performed on 2um resolution
    // first, before aggregating the SE object to the other lower resolutions
    ch_tissue_positions_002um = CONVERT_TISSUE_POSITIONS.out.csv.filter { resolution, _csv -> resolution == '002um' }.map { _resolution, csv -> csv }
    ch_tissue_positions_bins  = CONVERT_TISSUE_POSITIONS.out.csv.filter { resolution, _csv -> resolution != '002um' }

    emit:
    bam                    = FILTER_BARCODED_BAM.out.bam
    sample_name            = ch_sample.map { sample, _path, _meta -> sample }.first()
    tissue_positions_002um = ch_tissue_positions_002um
    tissue_positions_bins  = ch_tissue_positions_bins
    barcode_mappings       = CONVERT_BARCODE_MAPPINGS.out.csv
}
