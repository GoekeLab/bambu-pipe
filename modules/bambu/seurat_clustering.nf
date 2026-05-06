process SEURAT_CLUSTERING {
    publishDir "$params.output_dir/intermediate_R", mode: 'copy', pattern: '*.rds', enabled: params.save_intermediates
    label "r"
    label "medium_cpu"
    label "high_mem"
    label "medium"

    input:
    path(gene_counts)
    path(sample_names)
    val(run_clustering)

    output:
    path ('clusters.rds'), emit: clusters
    path ('cell_mix.rds'), optional: true, emit: cell_mix
    path "versions.yml", emit: versions

    script:
    """
    #!/usr/bin/env Rscript
    if ("$params.bambu_path" == "null") { library("bambu") } else { library("devtools"); load_all("$params.bambu_path") }
    library(Seurat)

    # Seurat Clustering 
    clusterCells <- function(counts, resolution, dim = 15){
        cellMix <- CreateSeuratObject(counts = counts, project = "cellMix", min.cells = 1)
        cellMix <- NormalizeData(cellMix, normalization.method = "LogNormalize", scale.factor = 10000)
        cellMix <- FindVariableFeatures(cellMix, selection.method = "vst", nfeatures = 2500)
        cellMix <- ScaleData(cellMix)
        npcs <- ifelse(ncol(counts)>50, 50, ncol(counts)-1)
        cellMix <- RunPCA(cellMix, features = VariableFeatures(object = cellMix), npcs = npcs)
        dim <- ifelse(dim >= dim(cellMix@reductions\$pca)[2], dim(cellMix@reductions\$pca)[2],dim) # if data dimension is small, otherwise, cap dimension at 15
        cellMix <- FindNeighbors(cellMix, dims = 1:dim)
        cellMix <- FindClusters(cellMix, resolution = resolution)
        cellMix <- RunUMAP(cellMix, dims = 1:dim)

        return(cellMix)
    }

    clustersTemp <- NULL
    if ("$run_clustering" == "true") {
        # Joint clustering across all samples
        counts <- readRDS("$gene_counts")
        cellMix <- clusterCells(counts, $params.resolution)
        saveRDS(cellMix, "cell_mix.rds")

        # Extract barcodes and clusters from Seurat object
        allBarcodes   <- names(cellMix@active.ident)
        clusterLabels <- paste0("cluster", as.character(cellMix@active.ident))

        # Extract sample name for each barcode by stripping the trailing underscore
        sampleNames <- readRDS("$sample_names")
        sampleLabels <- sub("_[^_]+\$", "", allBarcodes)

        # Build an ordered list of CompressedCharacterLists, one per sample, in quantData order
        # Each CCL maps cluster labels to the barcodes belonging to that cluster
        clustersTemp <- setNames(lapply(sampleNames, function(s) {
            idx           <- which(sampleLabels == s)       # indices of barcodes belonging to sample s
            sampleBarcode <- allBarcodes[idx]               # barcodes for this sample in {sampleName}_{barcode} format
            clusterLabel  <- paste0(s, "_", clusterLabels[idx]) # cluster label per barcode: {sampleName}_cluster{N}
            splitAsList(sampleBarcode, clusterLabel)        # split barcodes into a CompressedCharacterList
        }), sampleNames)
    }

    saveRDS(clustersTemp, "clusters.rds")
    writeLines(c('"${task.process}":', paste0('    seurat: ', as.character(packageVersion("Seurat")))), "versions.yml")
    """
}
