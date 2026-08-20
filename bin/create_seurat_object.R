# Assumes Seurat is loaded by the caller.

# createSeuratObject() - build a Seurat object from a count directory holding
# counts.mtx.gz, barcodes.tsv.gz and features.tsv.gz
#   countsDir - path to a count directory written by saveCounts()
#   metadata  - data.frame of cell metadata, one row per barcode with matching row names
#   assay, minCells, project - passed through to CreateSeuratObject()

createSeuratObject <- function(countsDir, assay = "RNA", metadata = NULL, minCells = 0, project = "SeuratProject") {
    counts <- ReadMtx(mtx      = file.path(countsDir, "counts.mtx.gz"),
                      cells    = file.path(countsDir, "barcodes.tsv.gz"),
                      features = file.path(countsDir, "features.tsv.gz"))

    CreateSeuratObject(counts, assay = assay, meta.data = metadata, min.cells = minCells, project = project)
}
