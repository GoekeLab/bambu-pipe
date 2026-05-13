process SEURAT_SINGLE_SAMPLE {
    publishDir "$params.output_dir", mode: 'copy', pattern: 'seurat_obj.rds'
    publishDir "$params.output_dir/intermediate_R", mode: 'copy', pattern: 'clusters.rds', enabled: params.save_intermediates
    label "r"
    label "medium_cpu"
    label "medium_mem"
    label "medium"

    input:
    path(se)

    output:
    path ('clusters.rds'), emit: clusters
    path ('seurat_obj.rds'), emit: seurat_obj
    path "versions.yml", emit: versions

    script:
    """
    #!/usr/bin/env Rscript
    if ("$params.bambu_path" == "null") { library("bambu") } else { library("devtools"); load_all("$params.bambu_path") }
    library(Seurat)

    se     <- readRDS("$se")
    counts <- assays(se)\$counts
    dim    <- $params.seurat_dim_single

    # Single sample scRNA-seq clustering adapted from https://satijalab.org/seurat/articles/pbmc3k_tutorial.html
    cellMix <- CreateSeuratObject(counts = counts, project = "cellMix", min.cells = 1)
    cellMix <- NormalizeData(cellMix, normalization.method = "LogNormalize", scale.factor = 10000)
    cellMix <- FindVariableFeatures(cellMix, selection.method = "vst", nfeatures = 2500)
    cellMix <- ScaleData(cellMix)
    npcs    <- ifelse(ncol(counts) > 50, 50, ncol(counts) - 1)
    cellMix <- RunPCA(cellMix, features = VariableFeatures(object = cellMix), npcs = npcs)
    dim     <- ifelse(dim >= dim(cellMix@reductions\$pca)[2], dim(cellMix@reductions\$pca)[2], dim)
    cellMix <- FindNeighbors(cellMix, dims = 1:dim)
    cellMix <- FindClusters(cellMix, resolution = $params.resolution)
    saveRDS(cellMix, "seurat_obj.rds")

    x <- setNames(names(cellMix@active.ident), cellMix@active.ident)
    clusters <- list(splitAsList(unname(x), paste0("cluster", names(x))))

    saveRDS(cellMix, "cell_mix.rds")
    saveRDS(clusters, "clusters.rds")
    writeLines(c('"${task.process}":', paste0('    seurat: ', as.character(packageVersion("Seurat")))), "versions.yml")
    """
}
