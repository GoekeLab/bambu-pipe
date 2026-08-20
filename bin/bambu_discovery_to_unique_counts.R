# Shared by the standard and Visium HD transcript-discovery Nextflow modules.
# Assumes bambu is loaded by the caller.

# bambuDiscoveryToUniqueCounts() - run transcript discovery, then quantify without EM
#   reads      - named character vector of readClass .rds paths (names = sample names),
#                or a single path for a one-sample run
#   annotation - the bambu transcript annotation object (a GRangesList)
#   genome     - path to the reference genome FASTA
#   ncore      - integer, number of cores passed to bambu.singlecell
#   ndr        - numeric, the novel discovery rate (NDR) threshold for discovery
#   sampleData - optional path(s) to spatial / bin metadata 
# Returns a named list with extendedAnno, quantData and seDiscovery (the
# transcript-level unique-counts SE)
bambuDiscoveryToUniqueCounts <- function(reads, annotation, genome, ncore, ndr, sampleData = NULL) {
    # Transcript discovery
    extendedAnno <- bambu.singlecell(reads = reads, output = "extendedAnnotations",
        annotations = annotation, genome = genome, ncore = ncore,
        verbose = FALSE, NDR = ndr)

    # Read to transcript assignment
    quantData <- bambu.singlecell(reads = reads, output = "quantData",
        annotations = extendedAnno, genome = genome, ncore = ncore,
        verbose = FALSE, sampleData = sampleData)

    # Quantification without EM
    seDiscovery <- bambu.singlecell(reads = quantData, output = "uniqueCounts",
        annotations = extendedAnno)

    list(extendedAnno = extendedAnno, quantData = quantData, seDiscovery = seDiscovery)
}
