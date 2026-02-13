#! /usr/bin/env nextflow

nextflow.enable.dsl=2

include { PARSE_SAMPLESHEET } from './subworkflows/parse_samplesheet.nf'
include { PREPROCESS_FASTQ } from './subworkflows/preprocess_fastq.nf'
include { ALIGNMENT } from './subworkflows/alignment.nf'
include { BAMBU_CONSTRUCT_READ_CLASS } from './modules/bambu_construct_read_class.nf'
include { BAMBU_PREPARE_ANNOTATION } from './modules/bambu_prepare_annotation.nf'
include { BAMBU } from './modules/bambu.nf'
include { BAMBU_EM } from './modules/bambu_EM.nf'

workflow {
    def barcode_map_default = true
    def ndr = params.ndr ?: 'NULL'
    def run_em = params.quantification_mode != 'no_EM'
    def run_clustering = params.quantification_mode == 'EM_clusters'

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
    PARSE_SAMPLESHEET(ch_input, ch_barcode_coordinate_config)

    // input files are split by type (fastq, bam, rds)
    ch_fastq_rows = PARSE_SAMPLESHEET.out.fastq
    ch_bam_rows = PARSE_SAMPLESHEET.out.bam
    ch_rds_rows = PARSE_SAMPLESHEET.out.rds

    // process fastq samples
    PREPROCESS_FASTQ(ch_fastq_rows, ch_adapter_seq_config,  ch_flank_seq_config)
    ALIGNMENT(PREPROCESS_FASTQ.out, ch_genome, ch_annotation)

    // process bam samples
    ch_bam_files = ALIGNMENT.out.concat(ch_bam_rows) // concatenate aligned bam files with input bam files
    ch_bambu_annotation = BAMBU_PREPARE_ANNOTATION(ch_annotation) // prepare annotation once for all samples
    BAMBU_CONSTRUCT_READ_CLASS(ch_bam_files, ch_genome, ch_bambu_annotation)

    // process rds samples
    ch_rds_files = BAMBU_CONSTRUCT_READ_CLASS.out.concat(ch_rds_rows) // concatenate constructed read class rds files with input rds files
    ch_rds_files_collect = ch_rds_files.collect(flat:false).map { it.transpose() } // collect all rds files into a single tuple
    BAMBU(ch_rds_files_collect, ch_genome, ch_bambu_annotation, ndr, run_clustering)
	
    if(run_em){
    BAMBU_EM(ch_rds_files_collect, BAMBU.out.quant_data, BAMBU.out.extended_annotations, BAMBU.out.clusters, ch_genome)
	}
}