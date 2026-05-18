#! /usr/bin/env nextflow

include { DECOMPRESS as DECOMPRESS_GENOME }     from './modules/decompress.nf'
include { DECOMPRESS as DECOMPRESS_ANNOTATION } from './modules/decompress.nf'
include { PREPARE_INPUT_STANDARD } from './subworkflows/prepare_input_standard.nf'
include { PREPROCESS_FASTQ } from './modules/preprocess_fastq.nf'
include { ALIGNMENT } from './subworkflows/alignment.nf'
include { BAMBU_CONSTRUCT_READ_CLASS } from './modules/bambu/construct_read_class.nf'
include { BAMBU_PREPARE_ANNOTATION } from './modules/bambu/prepare_annotation.nf'
include { BAMBU_TRANSCRIPT_DISCOVERY } from './modules/bambu/transcript_discovery.nf'
include { CLUSTERING } from './subworkflows/clustering.nf'
include { BAMBU_EM } from './modules/bambu/EM_quant.nf'

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
    resolution: Float
}

workflow {
    Validation.validateParams(params, workflow)

    def ndr = params.ndr ?: 'NULL'

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

    // load config files
    ch_barcode_coordinate_config = file("${projectDir}/assets/10x_config/barcode_coordinate_config.csv", checkIfExists: true)
    ch_adapter_seq_config = file("${projectDir}/assets/10x_config/adapter_seq_config.csv", checkIfExists: true)
    ch_flank_seq_config = file("${projectDir}/assets/10x_config/flank_seq_config.csv", checkIfExists: true)

    // parsing samplesheet csv file
    ch_input = channel.of(params.input)

    ch_standard  = ch_input.splitCsv(header:true, sep:',')
    ch_n_samples = ch_standard.count()

    PREPARE_INPUT_STANDARD(ch_standard, ch_barcode_coordinate_config)

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
        BAMBU_TRANSCRIPT_DISCOVERY(ch_rds_files_collect, ch_genome, BAMBU_PREPARE_ANNOTATION.out.annotation, ndr)

        if (params.quantification_mode != 'no_quant') {
            if (params.quantification_mode == 'EM_clusters') {
                CLUSTERING(BAMBU_TRANSCRIPT_DISCOVERY.out.se_gene_counts, BAMBU_TRANSCRIPT_DISCOVERY.out.sample_names, ch_n_samples)
                ch_clusters = CLUSTERING.out.clusters.map { clusters -> [true, clusters] } // flag to indicate that clustering was performed
            } else {
                ch_clusters = channel.value([false, []]) // flag to indicate that clustering was not performed
            }
            BAMBU_EM(BAMBU_TRANSCRIPT_DISCOVERY.out.quant_data, BAMBU_TRANSCRIPT_DISCOVERY.out.extended_annotations, ch_clusters, ch_genome)
        }
    }

    channel.topic('versions').collectFile(name: 'software_versions.yml', storeDir: "${params.output_dir}")
}