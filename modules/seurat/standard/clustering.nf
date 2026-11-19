process SEURAT_CLUSTERING {
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
    library(harmony)
    source(Sys.which("create_seurat_object.R"))

    metadata <- readRDS("$col_data")
    dim      <- $params.seurat_dims

    # Create Seurat object with bambu's colData as its metadata
    cellMix <- createSeuratObject("$gene_counts", metadata = metadata, minCells = 1, project = "cellMix")
    cellMix\$orig.ident <- cellMix\$sampleName

    # Clustering adapted from https://satijalab.org/seurat/articles/pbmc3k_tutorial.html
    cellMix <- NormalizeData(cellMix)
    cellMix <- FindVariableFeatures(cellMix)
    cellMix <- ScaleData(cellMix)
    cellMix <- RunPCA(cellMix, npcs = min(50, ncol(cellMix) - 1))

    # Correct for technology and chemistry with Harmony (skipped when neither varies, e.g. a single sample)
    covariates <- c("technology", "chemistry")
    covariates <- covariates[vapply(covariates, function(v) length(unique(cellMix[[v, drop = TRUE]])) > 1, logical(1))]
    reduction  <- if (length(covariates) > 0) "harmony" else "pca"
    if (length(covariates) > 0) {
        set.seed($params.harmony_seed)
        cellMix <- RunHarmony(cellMix, group.by.vars = covariates, reduction.use = "pca", reduction.save = "harmony")
    }

    dim     <- min(dim, ncol(cellMix[[reduction]]))
    cellMix <- FindNeighbors(cellMix, reduction = reduction, dims = 1:dim)
    cellMix <- FindClusters(cellMix, resolution = $params.seurat_resolution, cluster.name = "clusters")
    saveRDS(cellMix, "seurat_obj.rds")

    clusters <- setNames(paste0("cluster_", cellMix\$clusters), names(cellMix\$clusters))
    saveRDS(clusters, "clusters.rds")

    writeLines(c(
        '"${task.process}":',
        paste0('    R: ',       R.Version()\$version.string),
        paste0('    seurat: ',  as.character(packageVersion("Seurat"))),
        paste0('    harmony: ', as.character(packageVersion("harmony")))
    ), "versions.yml")
    """
}
