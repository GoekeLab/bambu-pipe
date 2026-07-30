process BAMBU_EM {
    publishDir "$params.output_dir", mode: 'copy', pattern: 'transcript_counts_singlecell'
    label "r"
    label "low_cpu"
    label "high_mem"
    label "long"

    input:
    path(quant_data)
    path(extended_annotation)
    path(genome)

    output:
    path ('transcript_counts_singlecell')
    path "versions.yml", topic: 'versions'

    script:
    """
    #!/usr/bin/env Rscript
    if ("$params.bambu_path" == "null") { library("bambu") } else { library("devtools"); load_all("$params.bambu_path") }
    library(DropletUtils)
    source(Sys.which("save_counts.R"))

    extendedAnno <- readRDS("$extended_annotation")
    quantData    <- readRDS("$quant_data")

    se <- bambu.singlecell(
        reads = quantData,
        output = "EM",
        annotations = extendedAnno,
        genome = "$genome",
        ncore = $task.cpus,
        verbose = FALSE,
        opt.em = list(degradationBias = FALSE)
    )

    saveCounts(se, "transcript_counts_singlecell", "Transcript Expression")

    writeLines(c('"${task.process}":', paste0('    R: ', R.Version()\$version.string), paste0('    bambu: ', as.character(packageVersion("bambu")))), "versions.yml")
    """
}
