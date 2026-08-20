# Assumes DropletUtils and bambu are loaded by the caller.

# saveCounts() - write a SummarizedExperiment as a 10x-style count directory, so
# downstream tools can read the matrices without loading the whole object.
#   se       - SummarizedExperiment object from bambu
#   dir      - name of the output directory
#   geneType - feature type written to the third column of features.tsv.gz,
#              "Transcript Expression" for transcript-level se, otherwise the
#              "Gene Expression" default.
# Writes into dir:
#   <assayName>.mtx.gz - one MatrixMarket file per assay, (e.g. counts, CPM)
#   barcodes.tsv.gz    - colnames(se), one per line, shared by every matrix
#   features.tsv.gz    - rownames(se) and geneType, shared by every matrix
#   se_<dir>.rds       - SummarizedExperiment object

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
