process SEURAT_MULTI_SAMPLE {
    publishDir "$params.output_dir/intermediate_R", mode: 'copy', pattern: '*.rds', enabled: params.save_intermediates
    label "r"
    label "medium_cpu"
    label "high_mem"
    label "medium"

    input:
    path(se)
    path(sample_names)

    output:
    path ('clusters.rds'), emit: clusters
    path ('cell_mix.rds'), optional: true, emit: cell_mix
    path "versions.yml", emit: versions

    script:
    """
    #!/usr/bin/env Rscript
    if ("$params.bambu_path" == "null") { library("bambu") } else { library("devtools"); load_all("$params.bambu_path") }
    library(Seurat)

    # Extract gene count matrix and colData metadata
    se           <- readRDS("$se")
    counts       <- assays(se)\$counts
    dim          <- $params.seurat_dim_multi
    chemistry    <- as.character(colData(se)\$chemistry)
    technology   <- as.character(colData(se)\$technology)
    sampleLabels <- as.character(colData(se)\$sampleName)

    # Create Seurat object and append metadata
    cellMix <- CreateSeuratObject(counts = counts, project = "cellMix", min.cells = 1)
    cellMix\$sample     <- sampleLabels
    cellMix\$chemistry  <- chemistry
    cellMix\$technology <- technology

    # scRNA-seq multi-sample integration using Harmony adapted from https://satijalab.org/seurat/articles/seurat5_integration
    cellMix[["RNA"]] <- split(cellMix[["RNA"]], f = cellMix\$sample)
    cellMix <- NormalizeData(cellMix)
    cellMix <- FindVariableFeatures(cellMix)
    cellMix <- ScaleData(cellMix)
    cellMix <- RunPCA(cellMix)

    cellMix <- IntegrateLayers(
        object         = cellMix,
        method         = HarmonyIntegration,
        orig.reduction = "pca",
        new.reduction  = "harmony",
        group.by.vars  = c("technology", "chemistry"),
        verbose        = FALSE
    )

    cellMix <- FindNeighbors(cellMix, reduction = "harmony", dims = 1:dim)
    cellMix <- FindClusters(cellMix, resolution = $params.resolution, cluster.name = "harmony_clusters")
    saveRDS(cellMix, "cell_mix.rds")

    # Build ordered list of CompressedCharacterLists, one per sample, in quantData order
    allBarcodes   <- names(cellMix\$harmony_clusters)
    clusterLabels <- paste0("cluster", as.character(cellMix\$harmony_clusters))
    sampleNames   <- readRDS("$sample_names") # sampleNames contain the order of the samples in quantData 

    clusters <- setNames(lapply(sampleNames, function(s) {
        idx           <- which(cellMix\$sample == s)
        sampleBarcode <- allBarcodes[idx]
        clusterLabel  <- paste0(s, "_", clusterLabels[idx])
        splitAsList(sampleBarcode, clusterLabel)
    }), sampleNames)
    saveRDS(clusters, "clusters.rds")

    writeLines(c('"${task.process}":', paste0('    seurat: ', as.character(packageVersion("Seurat")))), "versions.yml")
    """
}
