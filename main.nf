#! /usr/bin/env nextflow

nextflow.enable.dsl=2

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
    if (!params.reads) {
        error "params.reads is not set — please provide a path to a CSV samplesheet"
    }

    if (!params.genome) {
        error "params.genome is not set — please provide a path to the reference genome FASTA file"
    }

    if (!params.annotation) {
        error "params.annotation is not set — please provide a path to the reference annotation GTF file"
    }

    // ensure samplesheet exists and is a CSV file
    def reads_file = file(params.reads, checkIfExists: true)
    if (reads_file.getExtension() != "csv") {
        error "params.reads must be a CSV samplesheet"
    }

    // parsing genome and annotation files
    ch_genome =  Channel.value(file(params.genome, checkIfExists: true))
	ch_annotation =  Channel.value(file(params.annotation, checkIfExists: true))

    // parsing samplesheet csv file
    ch_reads = Channel.fromPath(params.reads)
    | splitCsv(header:true, sep:',')
    | map { row ->
        def meta = [
            chemistry: row.containsKey("chemistry") ? row.chemistry : params.chemistry,
            technology: row.containsKey("technology") ? row.technology : params.technology,
            barcode_map: row.containsKey("barcode_map") && row.barcode_map ? row.barcode_map : barcode_map_default // For barcode_map if the column or value is missing use barcode_map_default
        ]
        
        [row.sample, file(row.path), meta]
    }
   
    // filtering input files by type (fastq, bam, rds)
    ch_fastq_rows = ch_reads.filter { sample, path, meta -> path.name.endsWith('.fastq') || path.name.endsWith('.fastq.gz') }
    ch_bam_rows = ch_reads.filter { sample, path, meta -> path.name.endsWith('.bam') }
    ch_rds_rows = ch_reads.filter { sample, path, meta -> path.name.endsWith('.rds') }

    // process fastq samples
    PREPROCESS_FASTQ(ch_fastq_rows)
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