# Assumes SummarizedExperiment, Matrix and dplyr are loaded by the caller.
# A spot is a 2um barcode, i.e. one column of se; a bin is the coarser square
# (e.g. 8um) that several spots fall into.
# Build a bin-resolution SummarizedExperiment from a 2um-resolution SE by summing
# counts across the spots that share each bin. spotMappings gives the spot -> bin
# assignment for one resolution; tissuePositions supplies the new colData,
# alongside the id and sampleName carried over from the 2um SE.
# rowRanges/rowData and metadata are carried over.
aggregateResolution <- function(se, tissuePositions, spotMappings) {
    # Every column of se is one spot, so look up the bin each spot falls in. The
    # prefixed forms are used so the aggregated columns are named like bambu ids.
    matchIdx     <- match(colData(se)$id, spotMappings$sample_barcode)
    spotToBinMap <- factor(spotMappings$sample_bin[matchIdx])

    # The bins become the columns of the aggregated matrix, in level order
    binIds      <- levels(spotToBinMap) # bambu ids, e.g. <sample>_s_008um_<array_row>_<array_col>-1
    spotsPerBin <- t(fac2sparse(spotToBinMap)) # spots x bins, 1 where a spot falls in a bin

    # Sum each bin's spots: (features x spots) %*% (spots x bins) = features x bins
    spotCounts <- assays(se)$counts
    binCounts  <- as(spotCounts %*% spotsPerBin, "CsparseMatrix")

    ## build new colData for the aggregated SE object
    
    # bin id -> barcode (without sampleName prefix), e.g. sample_s_008um_00247_00090-1 -> s_008um_00247_00090-1
    binIdToBarcodeMap <- distinct(spotMappings, sample_bin, barcode = bin)

    # bin id -> the 2um barcodes the bin contains
    binToSpotsMap <- data.frame(sample_bin = as.character(spotToBinMap), barcode = colData(se)$barcode) %>%
        group_by(sample_bin) %>%
        summarise(barcodes = list(barcode), .groups = "drop")

    # creates dataframe containing one row per bin, each row contains metadata information (e.g., barcode, in_tissue, array and pixel coordinates)
    metadataPerBin <- data.frame(sample_bin = binIds) %>%
        left_join(binIdToBarcodeMap, by = "sample_bin") %>%
        left_join(binToSpotsMap,     by = "sample_bin") %>%
        left_join(tissuePositions,   by = "barcode") %>%
        select(-sample_bin)

    sampleName <- unique(colData(se)$sampleName)
    binColData <- DataFrame(id = binIds, sampleName = sampleName, metadataPerBin, row.names = binIds)

    binSe <- SummarizedExperiment(assays = SimpleList(counts = binCounts),
                                  rowRanges = rowRanges(se),
                                  colData = binColData)

    # Carry the bambu metadata matrices over, summed the same way
    incompatiblePerSpot <- metadata(se)$incompatibleCounts
    nonuniquePerSpot    <- metadata(se)$nonuniqueCounts

    metadata(binSe)$incompatibleCounts <- as(incompatiblePerSpot %*% spotsPerBin, "CsparseMatrix")
    metadata(binSe)$nonuniqueCounts    <- as(nonuniquePerSpot %*% spotsPerBin, "CsparseMatrix")
    metadata(binSe)$seType             <- metadata(se)$seType
    binSe
}
