# Seurat clustering of a single Visium HD bin resolution, adapted from
# https://satijalab.org/seurat/articles/visiumhd_analysis_vignette
# Assumes Seurat is loaded by the caller, plus SeuratWrappers and Banksy
# when the spatially aware path is used.

# Seurat stores the map of bin-level barcode (e.g., 8um) -> cluster. To allow Bambu,
# to run quantification, we have to generate the 2um barocde -> cluster map
# as quantData is generated at the 2um resolution
mapSpotsToClusters <- function(binToClustersMap, spotMappings) {
    spotToClustersMap <- setNames(unname(binToClustersMap[spotMappings$sample_bin]), spotMappings$sample_barcode)
    spotToClustersMap[!is.na(spotToClustersMap)]
}

# Spatially aware clustering using Banksy (Refer to vignette for more information)
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

# Expression-only clustering (Refer to vignette for more information)
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
