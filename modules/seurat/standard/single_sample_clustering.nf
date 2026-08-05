process SEURAT_SINGLE_SAMPLE {
    publishDir "$params.output_dir", mode: 'copy', pattern: 'seurat_obj.rds'
    publishDir "$params.output_dir/intermediate_R", mode: 'copy', pattern: 'clusters.rds', enabled: params.save_intermediates
    label "seurat"
    label "medium_cpu"
    label "medium_mem"
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
    dim      <- $params.seurat_dim_single

    # Create Seurat object with bambu's colData as its metadata
    cellMix <- createSeuratObject("$gene_counts", metadata = metadata, minCells = 1, project = "cellMix")

    # Single sample scRNA-seq clustering adapted from https://satijalab.org/seurat/articles/pbmc3k_tutorial.html
    cellMix <- NormalizeData(cellMix, normalization.method = "LogNormalize", scale.factor = 10000)
    cellMix <- FindVariableFeatures(cellMix, selection.method = "vst", nfeatures = 2500)
    cellMix <- ScaleData(cellMix)
    npcs    <- ifelse(ncol(cellMix) > 50, 50, ncol(cellMix) - 1)
    cellMix <- RunPCA(cellMix, features = VariableFeatures(object = cellMix), npcs = npcs)
    dim     <- ifelse(dim >= dim(cellMix@reductions\$pca)[2], dim(cellMix@reductions\$pca)[2], dim)
    cellMix <- FindNeighbors(cellMix, dims = 1:dim)
    cellMix <- FindClusters(cellMix, resolution = $params.seurat_resolution, cluster.name = "clusters")
    
    saveRDS(cellMix, "seurat_obj.rds")

    clusters <- setNames(paste0("cluster_", cellMix\$clusters), names(cellMix\$clusters))
    saveRDS(clusters, "clusters.rds")
    
    writeLines(c(
        '"${task.process}":',
        paste0('    R: ',      R.Version()\$version.string),
        paste0('    seurat: ', as.character(packageVersion("Seurat")))
    ), "versions.yml")
    """
}
