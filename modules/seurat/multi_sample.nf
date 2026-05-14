process SEURAT_MULTI_SAMPLE {
    publishDir "$params.output_dir", mode: 'copy', pattern: 'seurat_obj.rds'
    publishDir "$params.output_dir/intermediate_R", mode: 'copy', pattern: 'clusters.rds', enabled: params.save_intermediates
    label "r"
    label "medium_cpu"
    label "high_mem"
    label "medium"

    input:
    path(se)
    path(sample_names)

    output:
    path ('clusters.rds'), emit: clusters
    path ('seurat_obj.rds'), emit: seurat_obj
    path "versions.yml", topic: 'versions'

    script:
    """
    #!/usr/bin/env Rscript
    library(SummarizedExperiment)
    library(IRanges)
    library(Seurat)

    # Extract gene count matrix and colData metadata
    se     <- readRDS("$se")
    counts <- assays(se)\$counts
    dim    <- $params.seurat_dim_multi

    # Create Seurat object and append metadata
    cellMix <- CreateSeuratObject(counts = counts, project = "cellMix", min.cells = 1)
    cellMix\$sample     <- setNames(colData(se)\$sampleName,  colnames(se))
    cellMix\$chemistry  <- setNames(colData(se)\$chemistry,   colnames(se))
    cellMix\$technology <- setNames(colData(se)\$technology,  colnames(se))
    cellMix\$orig.ident <- cellMix\$sample

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

    dim <- min(dim, ncol(cellMix[["harmony"]]))
    cellMix <- FindNeighbors(cellMix, reduction = "harmony", dims = 1:dim)
    cellMix <- FindClusters(cellMix, resolution = $params.resolution, cluster.name = "harmony_clusters")
    saveRDS(cellMix, "seurat_obj.rds")

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

    writeLines(c(
        '"${task.process}":',
        paste0('    seurat: ',              as.character(packageVersion("Seurat"))),
        paste0('    IRanges: ',             as.character(packageVersion("IRanges"))),
        paste0('    SummarizedExperiment: ', as.character(packageVersion("SummarizedExperiment")))
    ), "versions.yml")
    """
}
