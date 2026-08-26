class Validation {

    static def validateParams(params, workflow) {
        // Container engine check
        if (!(workflow.containerEngine in ['docker', 'singularity']))
            throw new Exception("A container engine is required — please run with '-profile docker' or '-profile singularity'")

        // File existence checks
        if (!params.input.exists())
            throw new Exception("params.input '${params.input}' does not exist")

        if (params.input.extension != 'csv')
            throw new Exception("params.input '${params.input}' must be a CSV file")

        if (!params.genome.exists())
            throw new Exception("params.genome '${params.genome}' does not exist")

        if (!params.annotation.exists())
            throw new Exception("params.annotation '${params.annotation}' does not exist")

        // Enum checks
        if (!params.valid_quantification_modes.contains(params.quantification_mode))
            throw new Exception("Invalid params.quantification_mode '${params.quantification_mode}' — must be one of: ${params.valid_quantification_modes.join(', ')}")

        // Numeric range checks
        if (params.seurat_resolution <= 0)
            throw new Exception("Invalid params.seurat_resolution '${params.seurat_resolution}' — must be a positive number")

        if (params.ndr != null && (params.ndr < 0 || params.ndr > 1))
            throw new Exception("Invalid params.ndr '${params.ndr}' — must be a float between 0 and 1")

        // Visium HD checks
        if (params.visium_hd) {
            def needsBarcodeMappings

            if (params.manual_clustering) {
                // the manual path only expands bins down to 2um when the clusters were made at a coarser bin
                needsBarcodeMappings = params.clustering_bin != 2
            }
            else {
                // the bins samplesheet only drives the automatic path
                if (params.bins == null)
                    throw new Exception("params.bins is required when params.visium_hd is true")

                if (!params.bins.exists())
                    throw new Exception("params.bins '${params.bins}' does not exist")

                if (params.bins.extension != 'csv')
                    throw new Exception("params.bins '${params.bins}' must be a CSV file")

                validateVisiumHDBins(params.bins)

                // only the clustering mode reads the clustering bin, so only it has to resolve
                if (params.quantification_mode == 'clusteredEM')
                    validateClusteringBin(params.bins, params.clustering_bin)

                // every bin coarser than 2um is aggregated from the 2um spots, which needs the mappings
                needsBarcodeMappings = readBinResolutions(params.bins).any { res -> res != "2" }
            }

            if (needsBarcodeMappings) {
                if (params.barcode_mappings == null)
                    throw new Exception("params.barcode_mappings is required to aggregate 2 um spots into bins")

                if (!params.barcode_mappings.exists())
                    throw new Exception("params.barcode_mappings '${params.barcode_mappings}' does not exist")

                if (params.barcode_mappings.extension != 'parquet')
                    throw new Exception("params.barcode_mappings '${params.barcode_mappings}' must be a .parquet file")
            }
        }

        // Loupe alignment checks (standard Visium only)
        if (params.loupe_alignment != null) {
            if (params.visium_hd)
                throw new Exception("params.loupe_alignment is not used by the Visium HD workflow — tissue filtering uses the tissue_positions parquet files instead")

            if (!params.loupe_alignment.exists())
                throw new Exception("params.loupe_alignment '${params.loupe_alignment}' does not exist")

            if (params.loupe_alignment.extension != 'json')
                throw new Exception("params.loupe_alignment '${params.loupe_alignment}' must be a .json file")
        }
    }

    static def validateVisiumHDRows(rows) {
        if (rows.size() != 1)
            throw new Exception("Visium HD requires exactly 1 sample in the samplesheet, but found ${rows.size()}")

        def row = rows[0]
        ["sample", "path"].each { col ->
            if (!row.containsKey(col))
                throw new Exception("Samplesheet is missing a required '${col}' column")
            if (!row[col])
                throw new Exception("A row in the samplesheet has an empty '${col}' value")
        }

        if (!row.path.endsWith('.bam'))
            throw new Exception("Visium HD sample '${row.sample}' must point to a pre-aligned, barcode-tagged BAM file")
    }

    static def validateVisiumHDBins(bins) {
        if (!readBinResolutions(bins).contains("2"))
            throw new Exception("params.bins '${bins}' must include a row for the native 2 um resolution (resolution = 2)")
    }

    static def readBinResolutions(bins) {
        def lines = bins.text.readLines().findAll { line -> line.trim() }
        if (lines.size() < 2) return [] // no data rows to read

        def header = lines[0].split(',').collect { col -> col.trim() }
        def resIdx = header.indexOf("resolution")
        lines[1..-1].collect { line -> line.split(',')[resIdx].trim() }
    }

    static def validateClusteringBin(bins, clusteringBin) {
        def resolutions = readBinResolutions(bins)
        if (!resolutions.contains(clusteringBin.toString()))
            throw new Exception("params.clustering_bin '${clusteringBin}' is not listed in params.bins '${bins}' — available bins: ${resolutions.join(', ')}")
    }

    static def validateVisiumSampleCount(samples, loupeAlignment) {
        def has_visium = samples.any { sample, path, meta -> meta.chemistry.startsWith('visium') }
        if (has_visium && samples.size() > 1)
            throw new Exception("Visium chemistry requires exactly 1 sample, but found ${samples.size()}")

        if (loupeAlignment != null && !has_visium)
            throw new Exception("params.loupe_alignment was provided but no sample uses a visium chemistry")
    }

    static def validateRow(row, params, log) {
        ["sample", "path", "chemistry", "technology"].each { col ->
            if (!row.containsKey(col))
                throw new Exception("Samplesheet is missing a required '${col}' column")
            if (!row[col])
                throw new Exception("A row in the samplesheet has an empty '${col}' value")
        }

        if (!params.valid_chemistries.contains(row.chemistry)) {
            if (row.path.endsWith('.bam'))
                log.warn "Sample '${row.sample}' has custom chemistry '${row.chemistry}' — please check that this is intentional."
            else
                throw new Exception("Sample '${row.sample}' has invalid chemistry '${row.chemistry}' — must be one of: ${params.valid_chemistries.join(', ')}")
        }

        if (!params.valid_technologies.contains(row.technology))
            throw new Exception("Sample '${row.sample}' has invalid technology '${row.technology}' — must be one of: ${params.valid_technologies.join(', ')}")

        if (row.chemistry.startsWith('visium') && params.loupe_alignment == null)
            throw new Exception("Visium sample '${row.sample}' requires params.loupe_alignment — a Loupe manual alignment .json with the tissue selection")
    }

}
