# save_counts() - write a SummarizedExperiment as a count directory: one
# MTX file per assay, plus a single shared barcodes.tsv.gz / features.tsv.gz,
# and the SummarisedExperiment object. 
# Assumes DropletUtils and bambu have been loaded

save_counts <- function(se, dir, gene.type = "Gene Expression") {
    write10xCounts(dir, assays(se)$counts, version = "3", gene.type = gene.type)
    file.rename(file.path(dir, "matrix.mtx.gz"), file.path(dir, "counts.mtx.gz"))
    extra_assays <- setdiff(assayNames(se), "counts")
    for (assay_name in extra_assays) {
        mat      <- as(assays(se)[[assay_name]], "CsparseMatrix")
        mtx_path <- file.path(dir, paste0(assay_name, ".mtx"))
        Matrix::writeMM(mat, mtx_path)
        R.utils::gzip(mtx_path, overwrite = TRUE)
    }
    saveRDS(se, file.path(dir, paste0("se_", dir, ".rds")))
}
