process BAMBU_TRANSCRIPT_DISCOVERY_VISIUM_HD {
    publishDir "$params.output_dir", mode: 'copy', pattern: 'extended_annotations.gtf'
    publishDir "$params.output_dir/unique_counts", mode: 'copy', pattern: 'unique_counts_002um'
    publishDir "$params.output_dir/gene_counts", mode: 'copy', pattern: 'gene_counts_002um'
    publishDir "$params.output_dir/intermediate_R", mode: 'copy', pattern: '*.rds', enabled: params.save_intermediates
    label "r"
    label "medium_cpu"
    label "high_mem"
    label "medium"

    input:
    tuple val(sample), path(rds_files), val(meta)
    path(genome)
    path(bambu_annotation)
    val(ndr)
    path(tissue_positions)

    output:
    path ('quant_data.rds'), emit: quant_data
    path ('extended_annotations.rds'), emit: extended_annotations
    path ('extended_annotations.gtf')
    path ('unique_counts_002um')
    path ('gene_counts_002um')
    path ('unique_counts_002um/se_unique_counts_002um.rds'), emit: se_unique_002um
    path ('unique_counts_002um/barcodes.tsv.gz'), emit: barcodes_002um
    path "versions.yml", topic: 'versions'

    script:
    """
    #!/usr/bin/env Rscript
    if ("$params.bambu_path" == "null") { library("bambu") } else { library("devtools"); load_all("$params.bambu_path") }
    library(DropletUtils)
    source(Sys.which("bambu_discovery.R"))
    source(Sys.which("save_counts.R"))

    annotation <- readRDS("$bambu_annotation")
    readClassFile <- setNames("$rds_files", "$sample")

    # the 2um tissue_positions is passed as sampleData so bambu joins the per-spot spatial metadata into colData
    result <- bambuDiscovery(reads = readClassFile, annotation = annotation, genome = "$genome",
        ncore = $task.cpus, ndr = $ndr, sampleData = "$tissue_positions")
    saveRDS(result\$extendedAnno, "extended_annotations.rds")
    writeToGTF(result\$extendedAnno, "extended_annotations.gtf")
    saveRDS(result\$quantData, "quant_data.rds")

    unique002um <- result\$seDiscovery
    gene002um   <- transcriptToGeneExpression(unique002um)

    saveCounts(unique002um, "unique_counts_002um", "Transcript Expression")
    saveCounts(gene002um,   "gene_counts_002um")

    writeLines(c('"${task.process}":', paste0('    R: ', R.Version()\$version.string), paste0('    bambu: ', as.character(packageVersion("bambu")))), "versions.yml")
    """
}
