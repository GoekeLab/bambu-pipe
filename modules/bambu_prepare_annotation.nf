process BAMBU_PREPARE_ANNOTATION{
    container "ghcr.io/ch99l/bambu-pipe-r:latest"
    label "low_cpu"
    label "low_mem"
    label "short"

    input:
    path(annotation)

    output:
    path("bambu_annotation.rds")

    script:
    """
    #!/usr/bin/env Rscript
    if ("$params.bambu_path" == "null") { library("bambu") } else { library("devtools"); load_all("$params.bambu_path") }

    annotation <- prepareAnnotations("$annotation")
    saveRDS(annotation, "bambu_annotation.rds")
    """
}