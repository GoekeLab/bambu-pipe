#! /usr/bin/env nextflow

nextflow.enable.types = true

// subworkflows
include { PREPARE_INPUT_STANDARD }               from './subworkflows/prepare_input_standard.nf'
include { PREPARE_INPUT_VISIUM_HD }              from './subworkflows/prepare_input_visium_hd.nf'
include { ALIGNMENT }                            from './subworkflows/alignment.nf'
include { CLUSTERING }                           from './subworkflows/clustering_standard.nf'

// modules
include { DECOMPRESS as DECOMPRESS_GENOME }      from './modules/decompress.nf'
include { DECOMPRESS as DECOMPRESS_ANNOTATION }  from './modules/decompress.nf'
include { PREPROCESS_FASTQ }                     from './modules/preprocess_fastq.nf'
include { BAMBU_PREPARE_ANNOTATION }             from './modules/bambu/shared/prepare_annotation.nf'
include { BAMBU_CONSTRUCT_READ_CLASS }           from './modules/bambu/shared/construct_read_class.nf'
include { BAMBU_CLUSTER_LEVEL_QUANTIFICATION }   from './modules/bambu/shared/cluster_level_quantification.nf'
include { BAMBU_TRANSCRIPT_DISCOVERY }           from './modules/bambu/standard/transcript_discovery.nf'
include { BAMBU_SINGLE_CELL_QUANTIFICATION }     from './modules/bambu/standard/single_cell_quantification.nf'
include { BAMBU_TRANSCRIPT_DISCOVERY_VISIUM_HD } from './modules/bambu/visium_hd/transcript_discovery.nf'
include { BAMBU_SPOT_LEVEL_QUANTIFICATION }      from './modules/bambu/visium_hd/spot_level_quantification.nf'
include { AGGREGATE_BINS_VISIUM_HD }             from './modules/bambu/visium_hd/aggregate_bins.nf'
include { EXTRACT_SPOT_BARCODES }                from './modules/bambu/visium_hd/extract_spot_barcodes.nf'
include { MAP_CLUSTERS_TO_2UM_SPOTS }            from './modules/bambu/visium_hd/map_clusters_to_2um_spots.nf'
include { CONVERT_BARCODE_MAPPINGS }             from './modules/prepare_input/visium_hd/convert_barcode_mappings.nf'
include { SPOT_BIN_MAPPINGS }                    from './modules/prepare_input/visium_hd/spot_bin_mappings.nf'
include { SEURAT_VISIUM_HD }                     from './modules/seurat/visium_hd/clustering.nf'

params {
    input: Path
    genome: Path
    annotation: Path
    output_dir: Path
    chemistry: String?
    technology: String?
    bam_only: Boolean
    qscore_filtering: Boolean
    ndr: Float?
    deduplicate_umis: Boolean
    quantification_mode: String
    seurat_resolution: Float
    visium_hd: Boolean
    barcode_mappings: Path?
    bins: Path?
    clustering_bin: Integer
    banksy: Boolean
    banksy_lambda: Float
    banksy_k_geom: Integer
    manual_clustering: Boolean
}

workflow STANDARD {
    take:
    ch_rows: Channel<Map>
    ch_genome: Path
    ch_annotation: Path
    ndr: Float?
    manual_clustering: Boolean

    main:
    if (!manual_clustering) {
        def ndrArg = ndr != null ? ndr : 'NULL'

        // load config files
        ch_barcode_coordinate_config = file("${projectDir}/assets/10x_config/barcode_coordinate_config.csv", checkIfExists: true)
        ch_adapter_seq_config = file("${projectDir}/assets/10x_config/adapter_seq_config.csv", checkIfExists: true)
        ch_flank_seq_config = file("${projectDir}/assets/10x_config/flank_seq_config.csv", checkIfExists: true)

        ch_n_samples = ch_rows.count()

        PREPARE_INPUT_STANDARD(ch_rows, ch_barcode_coordinate_config)

        // input files are split by type (fastq, bam)
        ch_input_fastq = PREPARE_INPUT_STANDARD.out.fastq
        ch_input_bam = PREPARE_INPUT_STANDARD.out.bam

        // process fastq samples
        ch_preprocess_fastq_in = ch_input_fastq.map { sample, path, meta -> [sample, path, meta, meta.barcode] } // add whitelist path to fastq input tuple
        PREPROCESS_FASTQ(ch_preprocess_fastq_in, ch_flank_seq_config, ch_adapter_seq_config)
        ALIGNMENT(PREPROCESS_FASTQ.out.fastq, ch_genome, ch_annotation)

        if (!params.bam_only) {
            // process bam samples
            ch_bam_files = ALIGNMENT.out.bam.mix(ch_input_bam)
            BAMBU_PREPARE_ANNOTATION(ch_annotation)
            BAMBU_CONSTRUCT_READ_CLASS(ch_bam_files, ch_genome, BAMBU_PREPARE_ANNOTATION.out.annotation)

            // collect the read class files for joint clustering
            ch_rds_files_collect = BAMBU_CONSTRUCT_READ_CLASS.out.rds
                .map { sample, path, meta -> [sample, path, meta, meta.spatial_metadata] }
                .collect(flat:false)
                .map { collected_tup ->
                    def (samples, paths, metas, spatial_metadatas) = collected_tup.transpose()
                    def has_spatial = metas.any { meta -> meta.chemistry.startsWith('visium') } // for non-visium samples set the spatial metadata to an empty list (for staging)
                    [samples, paths, metas, has_spatial ? spatial_metadatas : []]
                }
            BAMBU_TRANSCRIPT_DISCOVERY(ch_rds_files_collect, ch_genome, BAMBU_PREPARE_ANNOTATION.out.annotation, ndrArg)

            // cluster the cells first, then pool each cluster's cells for the EM
            if (params.quantification_mode == 'clusteredEM') {
                CLUSTERING(BAMBU_TRANSCRIPT_DISCOVERY.out.gene_counts, BAMBU_TRANSCRIPT_DISCOVERY.out.col_data, ch_n_samples)
                BAMBU_CLUSTER_LEVEL_QUANTIFICATION(CLUSTERING.out.clusters, BAMBU_TRANSCRIPT_DISCOVERY.out.quant_data, BAMBU_TRANSCRIPT_DISCOVERY.out.extended_annotations, ch_genome)
            } else if (params.quantification_mode == 'EM') {
                BAMBU_SINGLE_CELL_QUANTIFICATION(BAMBU_TRANSCRIPT_DISCOVERY.out.quant_data, BAMBU_TRANSCRIPT_DISCOVERY.out.extended_annotations, ch_genome)
            }
        }

    } else {
        // restarting pipeline from quantData and manual cluster map
        ch_clusters   = ch_rows.map { row -> file(row.clusters_path, checkIfExists: true) }
        ch_quant_data = ch_rows.map { row -> file(row.quant_data_path, checkIfExists: true) }

        BAMBU_CLUSTER_LEVEL_QUANTIFICATION(ch_clusters, ch_quant_data, ch_annotation, ch_genome)
    }
}

workflow VISIUM_HD {
    take:
    ch_rows: Channel<Map>
    ch_genome: Path
    ch_annotation: Path
    ndr: Float?
    manual_clustering: Boolean

    main:
    def requested_bin = String.format('%03dum', params.clustering_bin) // convert clustering_bin specified as an integer into Spaceranger format

    if (!manual_clustering) {
        def ndrArg = ndr != null ? ndr : 'NULL'

        PREPARE_INPUT_VISIUM_HD(ch_rows)
        BAMBU_PREPARE_ANNOTATION(ch_annotation)
        BAMBU_CONSTRUCT_READ_CLASS(PREPARE_INPUT_VISIUM_HD.out.bam, ch_genome, BAMBU_PREPARE_ANNOTATION.out.annotation)
        BAMBU_TRANSCRIPT_DISCOVERY_VISIUM_HD(BAMBU_CONSTRUCT_READ_CLASS.out.rds, ch_genome, BAMBU_PREPARE_ANNOTATION.out.annotation, ndrArg, PREPARE_INPUT_VISIUM_HD.out.tissue_positions_002um)

        // perform transcript discovery at 2um first
        ch_quant_data     = BAMBU_TRANSCRIPT_DISCOVERY_VISIUM_HD.out.quant_data.first()
        ch_extended_anno  = BAMBU_TRANSCRIPT_DISCOVERY_VISIUM_HD.out.extended_annotations.first()
        ch_unique_002um   = BAMBU_TRANSCRIPT_DISCOVERY_VISIUM_HD.out.se_unique_002um.first()   // unique counts at 2um resolution
        ch_barcodes_002um = BAMBU_TRANSCRIPT_DISCOVERY_VISIUM_HD.out.barcodes_002um.first()    // list of 2um barcodes (same as the column names)

        // map every 2um spot in the SE to its bin, once per resolution; every module below reads this file
        SPOT_BIN_MAPPINGS(PREPARE_INPUT_VISIUM_HD.out.barcode_mappings, ch_barcodes_002um, PREPARE_INPUT_VISIUM_HD.out.sample_name)

        // pair each bin's tissue positions with its spot mappings, keyed on resolution
        ch_bins = PREPARE_INPUT_VISIUM_HD.out.tissue_positions_bins.join(SPOT_BIN_MAPPINGS.out.csv) // [resolution, tissue_positions, spot_mappings]

        // aggregate the 2um SEs into each requested bin resolution (e.g., 8um/16um)
        AGGREGATE_BINS_VISIUM_HD(ch_bins, ch_unique_002um)

        if (params.quantification_mode == 'clusteredEM') {
            // the 2um counts come straight from transcript discovery, every coarser bin from the aggregated SEs
            ch_counts_002um = BAMBU_TRANSCRIPT_DISCOVERY_VISIUM_HD.out.gene_counts_002um
                .combine(BAMBU_TRANSCRIPT_DISCOVERY_VISIUM_HD.out.col_data_002um)
                .map { gene_counts, col_data -> ['002um', gene_counts, col_data] }
            ch_counts_bins = AGGREGATE_BINS_VISIUM_HD.out.gene_counts.join(AGGREGATE_BINS_VISIUM_HD.out.col_data)

            // perform clustering at the requested resolution only
            ch_clustering = ch_counts_002um.concat(ch_counts_bins)
                .filter { resolution, _gene_counts, _col_data -> resolution == requested_bin } // [resolution, gene_counts, col_data]
            SEURAT_VISIUM_HD(ch_clustering)

            // since quantData is built at 2um, if clustering is performed at the bin-level,
            // we need to convert clusters from bin-level barcode -> cluster label
            // to 2um spot/barcode -> cluster label
            if (requested_bin == '002um') {
                ch_2um_spot_clusters = SEURAT_VISIUM_HD.out.clusters
            } else {
                ch_spot_mappings = SPOT_BIN_MAPPINGS.out.csv
                    .filter { resolution, _csv -> resolution == requested_bin }
                    .map { _resolution, csv -> csv }
                MAP_CLUSTERS_TO_2UM_SPOTS(SEURAT_VISIUM_HD.out.clusters, ch_spot_mappings)
                ch_2um_spot_clusters = MAP_CLUSTERS_TO_2UM_SPOTS.out.clusters
            }

            BAMBU_CLUSTER_LEVEL_QUANTIFICATION(ch_2um_spot_clusters, ch_quant_data, ch_extended_anno, ch_genome)

        } else if (params.quantification_mode == 'EM') {
            // run spot level quantification on all resolutions
            // at 2um resolution, tissue_positions and spot_mappings are not required
            ch_resolution_002um = channel.of(['002um', [], []])
            ch_resolutions      = ch_resolution_002um.concat(ch_bins) // [resolution, tissue_positions, spot_mappings]
            BAMBU_SPOT_LEVEL_QUANTIFICATION(ch_resolutions, ch_quant_data, ch_extended_anno, ch_genome)
        }

    } else {
        // restaring pipeline from quantData and manual cluster mapping
        ch_clusters   = ch_rows.map { row -> file(row.clusters_path, checkIfExists: true) }
        ch_quant_data = ch_rows.map { row -> file(row.quant_data_path, checkIfExists: true) }

        if (requested_bin == '002um') {
            // clusters made at 2um already name the spots bambu quantified
            ch_2um_spot_clusters = ch_clusters
        } else {
            // rebuild the spot mapping file
            EXTRACT_SPOT_BARCODES(ch_quant_data)
            ch_sample = EXTRACT_SPOT_BARCODES.out.sample_name.map { name -> name.text.trim() }
            CONVERT_BARCODE_MAPPINGS(requested_bin, channel.value(params.barcode_mappings))
            SPOT_BIN_MAPPINGS(CONVERT_BARCODE_MAPPINGS.out.csv, EXTRACT_SPOT_BARCODES.out.barcodes, ch_sample)

            ch_spot_mappings = SPOT_BIN_MAPPINGS.out.csv.map { _resolution, csv -> csv }
            MAP_CLUSTERS_TO_2UM_SPOTS(ch_clusters, ch_spot_mappings)
            ch_2um_spot_clusters = MAP_CLUSTERS_TO_2UM_SPOTS.out.clusters
        }

        BAMBU_CLUSTER_LEVEL_QUANTIFICATION(ch_2um_spot_clusters, ch_quant_data, ch_annotation, ch_genome)
    }
}

workflow {
    Validation.validateParams(params, workflow)

    // load reference files
    ch_genome     = channel.value(params.genome)
    ch_annotation = channel.value(params.annotation)

    if (params.genome.extension == 'gz') {
        DECOMPRESS_GENOME(ch_genome)
        ch_genome = DECOMPRESS_GENOME.out
    }

    if (params.annotation.extension == 'gz') {
        DECOMPRESS_ANNOTATION(ch_annotation)
        ch_annotation = DECOMPRESS_ANNOTATION.out
    }

    // parsing samplesheet csv file
    ch_input = channel.of(params.input)
    ch_rows  = ch_input.splitCsv(header:true, sep:',')

    if (params.visium_hd) {
        VISIUM_HD(ch_rows, ch_genome, ch_annotation, params.ndr, params.manual_clustering)
    } else {
        STANDARD(ch_rows, ch_genome, ch_annotation, params.ndr, params.manual_clustering)
    }

    channel.topic('versions').collectFile(name: 'software_versions.yml', storeDir: "${params.output_dir}")
}
