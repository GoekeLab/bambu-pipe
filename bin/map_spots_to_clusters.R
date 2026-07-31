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

# Seurat stores the map of bin-level barcode (e.g., 8um) -> cluster. To allow Bambu,
# to run quantification, we have to generate the 2um barocde -> cluster map
# as quantData is generated at the 2um resolution
mapSpotsToClusters <- function(binToClustersMap, spotMappings) {
    spotToClustersMap <- setNames(unname(binToClustersMap[spotMappings$sample_bin]), spotMappings$sample_barcode)
    spotToClustersMap[!is.na(spotToClustersMap)]
}
