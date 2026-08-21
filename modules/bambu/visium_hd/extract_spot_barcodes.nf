process EXTRACT_SPOT_BARCODES {
    label "bambu"
    label "low_cpu"
    label "medium_mem"
    label "short"

    input:
    path(quant_data)

    output:
    path ('barcodes.tsv.gz'), emit: barcodes
    path ('sample_name.txt'), emit: sample_name
    path "versions.yml", topic: 'versions'

    script:
    """
    #!/usr/bin/env Rscript
    if ("$params.bambu_path" == "null") { library("bambu") } else { library("devtools"); load_all("$params.bambu_path") }

    # a Visium HD run holds a single sample, whose ids name the 2um spots bambu quantified
    quantData  <- readRDS("$quant_data")
    sampleData <- quantData[[1]]@sampleData

    write.table(sampleData\$id, gzfile("barcodes.tsv.gz"), quote = FALSE, row.names = FALSE, col.names = FALSE)
    writeLines(unique(sampleData\$sampleName), "sample_name.txt")

    writeLines(c('"${task.process}":', paste0('    R: ', R.Version()\$version.string), paste0('    bambu: ', as.character(packageVersion("bambu")))), "versions.yml")
    """
}
