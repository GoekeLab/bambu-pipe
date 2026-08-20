# Seurat clustering of a single Visium HD bin resolution, adapted from
# https://satijalab.org/seurat/articles/visiumhd_analysis_vignette
# Assumes Seurat is loaded by the caller, plus SeuratWrappers and Banksy
# when the Nextflow parameter--banksy is set to True.

# clusterBanksy() - spatially aware clustering using Banksy (refer to the vignette)
#   object            - Seurat object of one bin resolution, with the fullres pixel
#                       coordinates in meta.data
#   assay             - name of the assay to cluster on
#   lambda            - Banksy spatial weight
#   kGeom             - number of spatial neighbours Banksy builds its features from
#   clusterResolution - resolution passed to FindClusters()
#   npcs              - number of PCs, capped at ncol(object) - 1
# Returns the Seurat object with a 'clusters' column in meta.data
clusterBanksy <- function(object, assay, lambda, kGeom, clusterResolution, npcs = 30) {
    # dimx/dimy name the meta.data columns carrying each bin's coordinates
    object <- RunBanksy(object, lambda = lambda, assay = assay, slot = "data",
                        features = "variable", k_geom = kGeom,
                        dimx = "pxl_col_in_fullres", dimy = "pxl_row_in_fullres", verbose = FALSE)

    DefaultAssay(object) <- "BANKSY"
    npcs   <- min(npcs, ncol(object) - 1)
    object <- RunPCA(object, assay = "BANKSY", reduction.name = "pca.banksy",
                     features = rownames(object), npcs = npcs, verbose = FALSE)
    object <- FindNeighbors(object, reduction = "pca.banksy", dims = 1:npcs, verbose = FALSE)
    FindClusters(object, cluster.name = "clusters", resolution = clusterResolution, verbose = FALSE)
}

# clusterExpression() - expression-only clustering (refer to the vignette)
#   object            - Seurat object of one bin resolution
#   assay             - name of the assay to cluster on
#   clusterResolution - resolution passed to FindClusters()
#   dims              - number of PCs, capped at ncol(object) - 1
#   ncells            - number of cells SketchData() samples; sketch-based analysis is
#                       only performed when the object has more cells than this
# Returns the Seurat object with a 'clusters' column in meta.data
clusterExpression <- function(object, assay, clusterResolution, dims = 15, ncells = 50000) {
    object <- FindVariableFeatures(object, verbose = FALSE)
    object <- ScaleData(object, verbose = FALSE)

    sketched <- ncol(object) > ncells
    if (sketched) {
        object <- SketchData(object, ncells = ncells, method = "LeverageScore", sketched.assay = "sketch")
        DefaultAssay(object) <- "sketch"
        object <- FindVariableFeatures(object, verbose = FALSE)
        object <- ScaleData(object, verbose = FALSE)
    }

    reduction <- if (sketched) "pca.sketch" else "pca"
    npcs      <- min(dims, ncol(object) - 1)
    object    <- RunPCA(object, reduction.name = reduction, npcs = npcs, verbose = FALSE)
    object    <- FindNeighbors(object, reduction = reduction, dims = 1:npcs, verbose = FALSE)
    object    <- FindClusters(object, cluster.name = "clusters", resolution = clusterResolution, verbose = FALSE)

    if (sketched) {
        # project the sketched labels onto the full set of bins, then adopt them as the
        # cluster assignment so both paths emit the same 'clusters' column
        object <- ProjectData(object, assay = assay, full.reduction = "full.pca.sketch",
                              sketched.assay = "sketch", sketched.reduction = "pca.sketch",
                              dims = 1:npcs, refdata = list(clusters.projected = "clusters"))
        DefaultAssay(object) <- assay
        object$clusters      <- object$clusters.projected
    }

    object
}
