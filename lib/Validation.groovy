class Validation {

    // FASTQ extensions the pipeline routes to preprocessing (any other sample must be a .bam)
    static final List<String> FASTQ_EXTENSIONS = ['.fastq', '.fq', '.fastq.gz', '.fq.gz']

    static def validateParams(params, workflow) {
        // Container engine check
        if (!(workflow.containerEngine in ['docker', 'singularity']))
            throw new Exception("A container engine is required — please run with '-profile docker' or '-profile singularity'")

        // a new run and a restart (--manual_clustering) take different inputs; a mix is rejected
        def newRunInputs  = [input: params.input, genome: params.genome, annotation: params.annotation]
        def restartInputs = [clusters: params.clusters, quant_data: params.quant_data, extended_annotations: params.extended_annotations]
        def (required, notUsed) = params.manual_clustering ? [restartInputs, newRunInputs] : [newRunInputs, restartInputs]

        def extra = notUsed.findAll { name, path -> path != null }.keySet().collect { name -> "--${name}" }
        if (extra) {
            if (params.manual_clustering)
                throw new Exception("A --manual_clustering restart takes --clusters, --quant_data and --extended_annotations instead — remove ${extra.join(', ')}")
            else
                throw new Exception("${extra.join(', ')} only apply when restarting with --manual_clustering")
        }

        required.each { name, path ->
            if (path == null)
                throw new Exception("params.${name} is required" + (params.manual_clustering ? " when params.manual_clustering is true" : ""))
            if (!path.exists())
                throw new Exception("params.${name} '${path}' does not exist")
        }

        // file type checks
        if (params.manual_clustering) {
            if (!['.csv', '.rds'].any { ext -> params.clusters.name.endsWith(ext) })
                throw new Exception("params.clusters '${params.clusters}' must be a .csv file with 'id' and 'cluster' columns")
            ['quant_data', 'extended_annotations'].each { name ->
                if (!params[name].name.endsWith('.rds'))
                    throw new Exception("params.${name} '${params[name]}' must be the ${name}.rds written to intermediate_R/ by the first run")
            }
        }
        else {
            if (params.input.extension != 'csv')
                throw new Exception("params.input '${params.input}' must be a CSV file")
            // bambu reads the annotation as GTF; a GFF3 file is read without an error into a wrong annotation
            if (!['.gtf', '.gtf.gz'].any { ext -> params.annotation.name.endsWith(ext) })
                throw new Exception("params.annotation '${params.annotation}' must be a .gtf or .gtf.gz file")
        }

        // Enum checks
        if (!params.valid_quantification_modes.contains(params.quantification_mode))
            throw new Exception("Invalid params.quantification_mode '${params.quantification_mode}' — must be one of: ${params.valid_quantification_modes.join(', ')}")

        // Numeric range checks
        if (params.seurat_resolution <= 0)
            throw new Exception("Invalid params.seurat_resolution '${params.seurat_resolution}' — must be a positive number")

        // clustering computes at most 50 principal components
        if (params.seurat_dims < 2 || params.seurat_dims > 50)
            throw new Exception("Invalid params.seurat_dims '${params.seurat_dims}' — must be an integer between 2 and 50")

        if (params.ndr != null && (params.ndr < 0 || params.ndr > 1))
            throw new Exception("Invalid params.ndr '${params.ndr}' — must be a float between 0 and 1")

        if (params.min_read_fraction_by_gene < 0 || params.min_read_fraction_by_gene > 1)
            throw new Exception("Invalid params.min_read_fraction_by_gene '${params.min_read_fraction_by_gene}' — must be a float between 0 and 1")

        // Visium HD checks
        if (params.visium_hd) {
            // banksy_lambda weighs the neighbourhood against the bin's own expression
            if (params.banksy_lambda < 0 || params.banksy_lambda > 1)
                throw new Exception("Invalid params.banksy_lambda '${params.banksy_lambda}' — must be a float between 0 and 1")

            if (params.banksy_k_geom < 1)
                throw new Exception("Invalid params.banksy_k_geom '${params.banksy_k_geom}' — must be a positive integer")

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

    // Validates a whole standard samplesheet: the checks on each row, then the checks across rows.
    // Returns the rows, so that samples are only processed once the samplesheet has passed.
    static def validateSamplesheet(rows, params, log) {
        if (rows.isEmpty())
            throw new Exception("Samplesheet has no samples")

        // checks on each row
        rows.each { row -> validateRow(row, params, log) }

        // checks across rows
        // duplicate names overwrite each other's outputs and collide in transcript discovery
        def duplicated = rows.collect { row -> row.sample }.countBy { name -> name }.findAll { name, n -> n > 1 }.keySet()
        if (duplicated)
            throw new Exception("Sample names must be unique — duplicated: ${duplicated.join(', ')}")

        def hasVisium = rows.any { row -> row.chemistry.startsWith('visium') }
        if (hasVisium && rows.size() > 1)
            throw new Exception("Visium chemistry requires exactly 1 sample, but found ${rows.size()}")

        if (params.loupe_alignment != null && !hasVisium)
            throw new Exception("params.loupe_alignment was provided but no sample uses a visium chemistry")

        rows
    }

    static def validateRow(row, params, log) {
        ["sample", "path", "chemistry", "technology"].each { col ->
            if (!row.containsKey(col))
                throw new Exception("Samplesheet is missing a required '${col}' column")
            if (!row[col])
                throw new Exception("A row in the samplesheet has an empty '${col}' value")
        }

        // samples are routed by extension; a file with any other extension would be silently dropped
        def supported = FASTQ_EXTENSIONS + ['.bam']
        if (!supported.any { ext -> row.path.endsWith(ext) })
            throw new Exception("Sample '${row.sample}' has an unsupported file '${row.path}' — must end in one of: ${supported.join(', ')}")

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
