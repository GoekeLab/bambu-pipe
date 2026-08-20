# Cluster assignments arrive either from Seurat or from the user, and bambu quantifies
# at 2um, so coarser bin labels have to be expanded down to the spots it quantified.
# Assumes data.table is available to the caller.

# readClusters() - load a cluster assignment as a named vector of id -> cluster label
#   path - a .rds holding a named vector, or a .csv/.tsv/.txt with id and cluster columns
readClusters <- function(path) {
    if (endsWith(path, ".rds")) return(readRDS(path))

    clusterDf <- data.table::fread(path, data.table = FALSE)
    setNames(as.character(clusterDf$cluster), clusterDf$id)
}

# mapSpotsToClusters() - build the 2 um barcode -> cluster map from a coarser one
#   binToClustersMap - named vector from readClusters(), bin barcode -> cluster
#   spotMappings     - data.frame mapping each 2 um sample_barcode to its sample_bin
# When clustering is performed at a coarser resolution, the original barcode to cluster
# map is at the bin-level (e.g., 8 um). However, quantData is generated at 2 um, so we
# need to generate a mapping from the 2 um barcodes to clusters.
mapSpotsToClusters <- function(binToClustersMap, spotMappings) {
    spotToClustersMap <- setNames(unname(binToClustersMap[spotMappings$sample_bin]), spotMappings$sample_barcode)
    spotToClustersMap[!is.na(spotToClustersMap)]
}
