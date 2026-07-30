# saveCounts() - write a SummarizedExperiment as a count directory: one
# MTX file per assay, plus a single shared barcodes.tsv.gz / features.tsv.gz,
# and the SummarisedExperiment object.
# Assumes DropletUtils and bambu have been loaded

saveCounts <- function(se, dir, geneType = "Gene Expression") {
    write10xCounts(dir, assays(se)$counts, version = "3", gene.type = geneType)
    file.rename(file.path(dir, "matrix.mtx.gz"), file.path(dir, "counts.mtx.gz"))
    extraAssays <- setdiff(assayNames(se), "counts")
    for (assayName in extraAssays) {
        mat     <- as(assays(se)[[assayName]], "CsparseMatrix")
        mtxPath <- file.path(dir, paste0(assayName, ".mtx"))
        Matrix::writeMM(mat, mtxPath)
        R.utils::gzip(mtxPath, overwrite = TRUE)
    }
    saveRDS(se, file.path(dir, paste0("se_", dir, ".rds")))
}
