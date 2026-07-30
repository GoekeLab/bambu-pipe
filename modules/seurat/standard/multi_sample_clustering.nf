process SEURAT_MULTI_SAMPLE {
    publishDir "$params.output_dir", mode: 'copy', pattern: 'seurat_obj.rds'
    publishDir "$params.output_dir/intermediate_R", mode: 'copy', pattern: 'clusters.rds', enabled: params.save_intermediates
    label "seurat"
    label "medium_cpu"
    label "high_mem"
    label "medium"

    input:
    tuple path(gene_counts), path(col_data)

    output:
    path ('clusters.rds'), emit: clusters
    path ('seurat_obj.rds'), emit: seurat_obj
    path "versions.yml", topic: 'versions'

    script:
    """
    #!/usr/bin/env Rscript
    library(Seurat)
    source(Sys.which("create_seurat_object.R"))

    metadata <- readRDS("$col_data")
    dim      <- $params.seurat_dim_multi

    # Create Seurat object with bambu's colData as its metadata
    cellMix <- createSeuratObject("$gene_counts", metadata = metadata, minCells = 1, project = "cellMix")
    cellMix\$orig.ident <- cellMix\$sampleName

    # scRNA-seq multi-sample integration using Harmony adapted from https://satijalab.org/seurat/articles/seurat5_integration
    cellMix[["RNA"]] <- split(cellMix[["RNA"]], f = cellMix\$sampleName)
    cellMix <- NormalizeData(cellMix)
    cellMix <- FindVariableFeatures(cellMix)
    cellMix <- ScaleData(cellMix)
    cellMix <- RunPCA(cellMix)

    cellMix <- IntegrateLayers(
        object         = cellMix,
        method         = HarmonyIntegration,
        orig.reduction = "pca",
        new.reduction  = "harmony",
        group.by.vars  = c("technology", "chemistry"),
        verbose        = FALSE
    )

    dim <- min(dim, ncol(cellMix[["harmony"]]))
    cellMix <- FindNeighbors(cellMix, reduction = "harmony", dims = 1:dim)
    cellMix <- FindClusters(cellMix, resolution = $params.cluster_resolution, cluster.name = "harmony_clusters")
    saveRDS(cellMix, "seurat_obj.rds")

    clusters <- setNames(paste0("cluster_", cellMix\$harmony_clusters), names(cellMix\$harmony_clusters))
    saveRDS(clusters, "clusters.rds")

    writeLines(c(
        '"${task.process}":',
        paste0('    R: ',       R.Version()\$version.string),
        paste0('    seurat: ',  as.character(packageVersion("Seurat"))),
        paste0('    harmony: ', as.character(packageVersion("harmony")))
    ), "versions.yml")
    """
}
