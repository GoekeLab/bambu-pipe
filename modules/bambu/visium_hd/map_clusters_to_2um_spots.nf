process MAP_CLUSTERS_TO_2UM_SPOTS {
    label "bambu"
    label "low_cpu"
    label "medium_mem"
    label "short"

    input:
    path(clusters)
    path(spot_mappings)

    output:
    path ('spot_clusters.rds'), emit: clusters
    path "versions.yml", topic: 'versions'

    script:
    """
    #!/usr/bin/env Rscript
    source(Sys.which("map_spots_to_clusters.R"))

    binToClustersMap <- readClusters("$clusters")
    spotMappings     <- read.csv("$spot_mappings", stringsAsFactors = FALSE)

    # bambu quantifies at 2um, so map the bin level labels down to the 2um spots
    spotToClustersMap <- mapSpotsToClusters(binToClustersMap, spotMappings)
    if (length(spotToClustersMap) == 0) {
        stop("no cluster id matches a bin in $spot_mappings -- check that the clusters were made at the ${params.clustering_bin}um resolution")
    }

    saveRDS(spotToClustersMap, "spot_clusters.rds")

    writeLines(c('"${task.process}":', paste0('    R: ', R.Version()\$version.string)), "versions.yml")
    """
}
