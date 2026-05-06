#! /usr/bin/env nextflow

nextflow.enable.dsl=2

include { PREPARE_INPUT_STANDARD } from './subworkflows/prepare_input_standard.nf'
include { PREPROCESS_FASTQ } from './modules/preprocess_fastq.nf'
include { ALIGNMENT } from './subworkflows/alignment.nf'
include { BAMBU_CONSTRUCT_READ_CLASS } from './modules/bambu/construct_read_class.nf'
include { BAMBU_PREPARE_ANNOTATION } from './modules/bambu/prepare_annotation.nf'
include { BAMBU_TRANSCRIPT_DISCOVERY } from './modules/bambu/transcript_discovery.nf'
include { SEURAT_CLUSTERING } from './modules/bambu/seurat_clustering.nf'
include { BAMBU_EM } from './modules/bambu/EM_quant.nf'

workflow {
    Validation.validateParams(params, workflow)

    def ndr = params.ndr ?: 'NULL'
    def run_read_class_construction = params.early_stop_stage != 'bam'
    def run_bambu_discovery = params.early_stop_stage == null
    def run_clustering = params.quantification_mode == 'EM_clusters'
    def run_bambu_em = run_bambu_discovery && params.quantification_mode != 'no_quant'

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

    // TODO: Visium HD routing — restore when Visium HD support is added
    // ch_input.splitCsv(header:true, sep:',')
    //     .toList()
    //     .map { rows ->
    //         def has_visium_hd = rows.any { it.containsKey('technology') && it.technology == 'visium-hd' }
    //         if (has_visium_hd && rows.size() > 1)
    //             error "Visium HD samples cannot be mixed with other samples"
    //         rows
    //     }
    //     .flatMap { it }
    //     .branch {
    //         visium_hd: it.containsKey('technology') && it.technology == 'visium-hd'
    //         standard: true
    //     }.set { ch_branched }
    // PREPARE_INPUT_STANDARD(ch_branched.standard, ch_barcode_coordinate_config)

    ch_input.splitCsv(header:true, sep:',')
        .set { ch_standard }

    PREPARE_INPUT_STANDARD(ch_standard, ch_barcode_coordinate_config)
    ch_versions = PREPARE_INPUT_STANDARD.out.versions

    // input files are split by type (fastq, bam, rds)
    ch_input_fastq = PREPARE_INPUT_STANDARD.out.fastq
    ch_input_bam = PREPARE_INPUT_STANDARD.out.bam
    ch_input_rds = PREPARE_INPUT_STANDARD.out.rds

    // process fastq samples
    ch_preprocess_fastq_in = ch_input_fastq.map { sample, path, meta -> [sample, path, meta, meta.barcode] } // add whitelist path to fastq input tuple
    PREPROCESS_FASTQ(ch_preprocess_fastq_in, ch_flank_seq_config, ch_adapter_seq_config)
    ch_versions = ch_versions.mix(PREPROCESS_FASTQ.out.versions.first())
    ALIGNMENT(PREPROCESS_FASTQ.out.fastq, ch_genome, ch_annotation)
    ch_versions = ch_versions.mix(ALIGNMENT.out.versions)

    // process bam samples
    if (run_read_class_construction) {
        ch_bam_files = ALIGNMENT.out.bam.concat(ch_input_bam) // concatenate aligned bam files with input bam files
        BAMBU_PREPARE_ANNOTATION(ch_annotation) // prepare annotation once for all samples
        ch_versions = ch_versions.mix(BAMBU_PREPARE_ANNOTATION.out.versions)
        BAMBU_CONSTRUCT_READ_CLASS(ch_bam_files, ch_genome, BAMBU_PREPARE_ANNOTATION.out.annotation)
    }

    // process rds samples
    if (run_bambu_discovery) {
        ch_rds_files = BAMBU_CONSTRUCT_READ_CLASS.out.rds.concat(ch_input_rds) // concatenate constructed read class rds files with input rds files
        // reshape and collect rds file channel
        ch_rds_files_collect = ch_rds_files
        .map { sample, path, meta -> [sample, path, meta, meta.spatial_metadata] }
        .collect(flat:false) 
        .map { it.transpose() } 
        BAMBU_TRANSCRIPT_DISCOVERY(ch_rds_files_collect, ch_genome, BAMBU_PREPARE_ANNOTATION.out.annotation, ndr)
    }

    if (run_bambu_em) {
        SEURAT_CLUSTERING(BAMBU_TRANSCRIPT_DISCOVERY.out.gene_counts, BAMBU_TRANSCRIPT_DISCOVERY.out.sample_names, run_clustering)
        ch_versions = ch_versions.mix(SEURAT_CLUSTERING.out.versions)
        BAMBU_EM(BAMBU_TRANSCRIPT_DISCOVERY.out.quant_data, BAMBU_TRANSCRIPT_DISCOVERY.out.extended_annotations, SEURAT_CLUSTERING.out.clusters, ch_genome)
    }

    ch_versions.collectFile(name: 'software_versions.yml', storeDir: "${params.output_dir}")
}