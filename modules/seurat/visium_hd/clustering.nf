process SEURAT_VISIUM_HD {
    publishDir "$params.output_dir", mode: 'copy', pattern: 'seurat_obj.rds'
    label "seurat"
    label "medium_cpu"
    label "high_mem"
    label "long"

    input:
    tuple val(resolution), path(gene_counts), path(col_data)

    output:
    path ("clusters.rds"), emit: clusters
    path ("seurat_obj.rds")
    path "versions.yml", topic: 'versions'

    script:
    """
    #!/usr/bin/env Rscript
    library(Seurat)
    source(Sys.which("create_seurat_object.R"))
    source(Sys.which("visium_hd_cluster.R"))

    banksy <- as.logical("$params.banksy")
    assay  <- "Spatial.$resolution"

    # Create Seurat object with bambu's colData as its metadata
    colData      <- readRDS("$col_data")
    binMetadata  <- colData[, colnames(colData) != "barcodes"] # drop barcodes column (list) to prevent downstream issues with Seurat
    object       <- createSeuratObject("$gene_counts", assay = assay, metadata = binMetadata)

    DefaultAssay(object) <- assay
    object <- NormalizeData(object, verbose = FALSE)

    # Cluster spatially with Banksy, or on gene expression alone
    if (banksy) {
        library(SeuratWrappers)
        library(Banksy)
        object <- clusterBanksy(object, assay, $params.banksy_lambda, $params.banksy_k_geom, $params.cluster_resolution)
    } else {
        object <- clusterExpression(object, assay, $params.cluster_resolution)
    }

    # Save the cluster assignment at the clustered resolution, plus the Seurat object
    clusters <- setNames(paste0("cluster_", object\$clusters), names(object\$clusters))
    saveRDS(clusters, "clusters.rds")
    saveRDS(object, "seurat_obj.rds")

    writeLines(c(
        '"${task.process}":',
        paste0('    R: ',      R.Version()\$version.string),
        paste0('    seurat: ', as.character(packageVersion("Seurat"))),
        paste0('    banksy: ', as.character(packageVersion("Banksy")))
    ), "versions.yml")
    """
}
