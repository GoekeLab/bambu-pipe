process SEURAT_CLUSTERING {
    container "ghcr.io/ch99l/bambu-pipe-r:latest"
    label "medium_cpu"
    label "high_mem"
    label "medium"

    input:
    path(quant_data)
    val(run_clustering)

    output:
    path ('clusters.rds'), emit: clusters

    script:
    """
    #!/usr/bin/env Rscript
    library("bambu")

    se = readRDS("$quant_data")

    # Seurat Clustering (if no clustering provided, automatically cluster)
    path <- Sys.getenv("PATH") |> strsplit(":")
    bin_path <- tail(path[[1]], n=1)
    clusters = NULL
    if(as.logical("$run_clustering")){
        clusters = list()
        cellMixs = list()
        source(file.path(bin_path,"/utilityFunctions.R"))
        for(quantData in se){
            quantData.gene = transcriptToGeneExpression(quantData)
            for(sample in unique(colData(quantData)\$sampleName)){
                i = which(colData(quantData)\$sampleName == sample)
                counts = assays(quantData.gene)\$counts[,i]
                cellMix = clusterCells(counts, resolution = $params.resolution)
                x = setNames(names(cellMix@active.ident), cellMix@active.ident)
                names(x) = paste0(sample,"_",names(x))
                clusters = c(clusters, splitAsList(unname(x), names(x)))
                cellMixs = c(cellMixs, cellMix)
            }
        }
        saveRDS(cellMixs, "cellMixs.rds")
    }

    saveRDS(clusters, "clusters.rds")
    """
}
