#! /usr/bin/env nextflow

nextflow.enable.dsl=2

include { PREPARE_INPUT_STANDARD } from './subworkflows/prepare_input_standard.nf'
include { PREPROCESS_FASTQ } from './modules/preprocess_fastq.nf'
include { ALIGNMENT } from './subworkflows/alignment.nf'
include { BAMBU_CONSTRUCT_READ_CLASS } from './modules/bambu_construct_read_class.nf'
include { BAMBU_PREPARE_ANNOTATION } from './modules/bambu_prepare_annotation.nf'
include { BAMBU } from './modules/bambu.nf'
include { BAMBU_EM } from './modules/bambu_EM.nf'

workflow {
    def ndr = params.ndr ?: 'NULL'
    def run_read_class_construction = params.early_stop_stage != 'bam'
    def run_bambu_discovery = params.early_stop_stage == null
    def run_clustering = params.quantification_mode == 'EM_clusters'
    def run_bambu_em = run_bambu_discovery && params.quantification_mode != 'no_EM' 
    
    // checking required params
    if (!params.input) {
        error "params.input is not set — please provide a path to a CSV samplesheet"
    }

    if (!params.genome) {
        error "params.genome is not set — please provide a path to the reference genome FASTA file"
    }

    if (!params.annotation) {
        error "params.annotation is not set — please provide a path to the reference annotation GTF file"
    }

    // load reference files
    ch_genome =  Channel.value(file(params.genome, checkIfExists: true))
	ch_annotation =  Channel.value(file(params.annotation, checkIfExists: true))

    // load config files
    ch_barcode_coordinate_config = file("${projectDir}/10x_config/barcode_coordinate_config.csv", checkIfExists: true)
    ch_adapter_seq_config = file("${projectDir}/10x_config/adapter_seq_config.csv", checkIfExists: true)
    ch_flank_seq_config = file("${projectDir}/10x_config/flank_seq_config.csv", checkIfExists: true)

    // parsing samplesheet csv file
    ch_input = Channel.fromPath(params.input, checkIfExists: true)
    .ifEmpty { error "Cannot find samplesheet file: ${params.input}" }
    .map { file ->
        if (file.extension != "csv") {
            error "Invalid samplesheet. Must be a CSV file."
        }
        return file
    }

    PREPARE_INPUT_STANDARD(ch_input, ch_barcode_coordinate_config)

    // input files are split by type (fastq, bam, rds)
    ch_input_fastq = PREPARE_INPUT_STANDARD.out.fastq
    ch_input_bam = PREPARE_INPUT_STANDARD.out.bam
    ch_input_rds = PREPARE_INPUT_STANDARD.out.rds

    // process fastq samples
    ch_preprocess_fastq_in = ch_input_fastq.map { sample, path, meta -> [sample, path, meta, meta.barcode] } // add whitelist path to fastq input tuple
    PREPROCESS_FASTQ(ch_preprocess_fastq_in, ch_flank_seq_config, ch_adapter_seq_config)
    ALIGNMENT(PREPROCESS_FASTQ.out.fastq, ch_genome, ch_annotation, ch_input_fastq.count()) // fastq count is used to ensure paftools and minimap build index are skipped when there are no fastq samples

    // process bam samples
    if (run_read_class_construction) {
        ch_bam_files = ALIGNMENT.out.bam.concat(ch_input_bam) // concatenate aligned bam files with input bam files
        BAMBU_PREPARE_ANNOTATION(ch_annotation) // prepare annotation once for all samples
        BAMBU_CONSTRUCT_READ_CLASS(ch_bam_files, ch_genome, BAMBU_PREPARE_ANNOTATION.out)
    }

    // process rds samples
    if (run_bambu_discovery) {
        ch_rds_files = BAMBU_CONSTRUCT_READ_CLASS.out.rds.concat(ch_input_rds) // concatenate constructed read class rds files with input rds files
        // reshape and collect rds file channel 
        ch_rds_files_collect = ch_rds_files
        .map { sample, path, meta -> [sample, path, meta, meta.spatial_metadata] }
        .collect(flat:false) 
        .map { it.transpose() } 
        BAMBU(ch_rds_files_collect, ch_genome, BAMBU_PREPARE_ANNOTATION.out, ndr, run_clustering)
    }

    if (run_bambu_em) {
        ch_rds_files_em = ch_rds_files_collect.map { samples, paths, metas, spatial -> [samples, paths, metas] } // remove spatial metadata from input tuple for EM step
        BAMBU_EM(ch_rds_files_em, BAMBU.out.quant_data, BAMBU.out.extended_annotations, BAMBU.out.clusters, ch_genome)
    }
}