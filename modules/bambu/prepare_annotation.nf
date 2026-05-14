process BAMBU_PREPARE_ANNOTATION{
    label "r"
    label "low_cpu"
    label "low_mem"
    label "short"

    input:
    path(annotation)

    output:
    path("bambu_annotation.rds"), emit: annotation
    path "versions.yml", topic: 'versions'

    script:
    """
    #!/usr/bin/env Rscript
    if ("$params.bambu_path" == "null") { library("bambu") } else { library("devtools"); load_all("$params.bambu_path") }

    annotation <- prepareAnnotations("$annotation")
    saveRDS(annotation, "bambu_annotation.rds")
    
    writeLines(c('"${task.process}":', paste0('    bambu: ', as.character(packageVersion("bambu")))), "versions.yml")
    """
}