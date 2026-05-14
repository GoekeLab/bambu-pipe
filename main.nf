#! /usr/bin/env nextflow

nextflow.enable.dsl=2

include { PREPARE_INPUT_STANDARD } from './subworkflows/prepare_input_standard.nf'
include { PREPROCESS_FASTQ } from './modules/preprocess_fastq.nf'
include { ALIGNMENT } from './subworkflows/alignment.nf'
include { BAMBU_CONSTRUCT_READ_CLASS } from './modules/bambu/construct_read_class.nf'
include { BAMBU_PREPARE_ANNOTATION } from './modules/bambu/prepare_annotation.nf'
include { BAMBU_TRANSCRIPT_DISCOVERY } from './modules/bambu/transcript_discovery.nf'
include { CLUSTERING } from './subworkflows/clustering.nf'
include { BAMBU_EM } from './modules/bambu/EM_quant.nf'

workflow {
    Validation.validateParams(params, workflow)

    def ndr = params.ndr ?: 'NULL'

    // load reference files
    ch_genome =  Channel.value(file(params.genome, checkIfExists: true))
	ch_annotation =  Channel.value(file(params.annotation, checkIfExists: true))

    // load config files
    ch_barcode_coordinate_config = file("${projectDir}/assets/10x_config/barcode_coordinate_config.csv", checkIfExists: true)
    ch_adapter_seq_config = file("${projectDir}/assets/10x_config/adapter_seq_config.csv", checkIfExists: true)
    ch_flank_seq_config = file("${projectDir}/assets/10x_config/flank_seq_config.csv", checkIfExists: true)

    // parsing samplesheet csv file
    ch_input = Channel.fromPath(params.input, checkIfExists: true)
    .ifEmpty { error "Cannot find samplesheet file: ${params.input}" }
    .map { file ->
        if (file.extension != "csv") {
            error "Invalid samplesheet. Must be a CSV file."
        }
        return file
    }

    // TODO: Add Visium-HD and non-standard routing
    ch_standard  = ch_input.splitCsv(header:true, sep:',')
    ch_n_samples = ch_standard.count()

    PREPARE_INPUT_STANDARD(ch_standard, ch_barcode_coordinate_config)
    ch_versions = PREPARE_INPUT_STANDARD.out.versions

    // input files are split by type (fastq, bam)
    ch_input_fastq = PREPARE_INPUT_STANDARD.out.fastq
    ch_input_bam = PREPARE_INPUT_STANDARD.out.bam

    // process fastq samples
    ch_preprocess_fastq_in = ch_input_fastq.map { sample, path, meta -> [sample, path, meta, meta.barcode] } // add whitelist path to fastq input tuple
    PREPROCESS_FASTQ(ch_preprocess_fastq_in, ch_flank_seq_config, ch_adapter_seq_config)
    ch_versions = ch_versions.mix(PREPROCESS_FASTQ.out.versions.first())
    ALIGNMENT(PREPROCESS_FASTQ.out.fastq, ch_genome, ch_annotation)
    ch_versions = ch_versions.mix(ALIGNMENT.out.versions)

    if (!params.bam_only) {
        // process bam samples
        ch_bam_files = ALIGNMENT.out.bam.mix(ch_input_bam)
        BAMBU_PREPARE_ANNOTATION(ch_annotation)
        ch_versions = ch_versions.mix(BAMBU_PREPARE_ANNOTATION.out.versions)
        BAMBU_CONSTRUCT_READ_CLASS(ch_bam_files, ch_genome, BAMBU_PREPARE_ANNOTATION.out.annotation)

        ch_rds_files_collect = BAMBU_CONSTRUCT_READ_CLASS.out.rds
            .map { sample, path, meta -> [sample, path, meta, meta.spatial_metadata] }
            .collect(flat:false)
            .map { collected_tup -> collected_tup.transpose() }
        BAMBU_TRANSCRIPT_DISCOVERY(ch_rds_files_collect, ch_genome, BAMBU_PREPARE_ANNOTATION.out.annotation, ndr)

        if (params.quantification_mode != 'no_quant') {
            if (params.quantification_mode == 'EM_clusters') {
                CLUSTERING(BAMBU_TRANSCRIPT_DISCOVERY.out.se_gene_counts, BAMBU_TRANSCRIPT_DISCOVERY.out.sample_names, ch_n_samples)
                ch_versions = ch_versions.mix(CLUSTERING.out.versions)
                ch_clusters = CLUSTERING.out.clusters.map { clusters -> [true, clusters] } // flag to indicate that clustering was performed
            } else {
                ch_clusters = Channel.value([false, []]) // flag to indicate that clustering was not performed
            }
            BAMBU_EM(BAMBU_TRANSCRIPT_DISCOVERY.out.quant_data, BAMBU_TRANSCRIPT_DISCOVERY.out.extended_annotations, ch_clusters, ch_genome)
        }
    }

    ch_versions.collectFile(name: 'software_versions.yml', storeDir: "${params.output_dir}")
}