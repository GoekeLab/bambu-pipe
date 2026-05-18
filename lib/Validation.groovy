class Validation {

    static def validateParams(params, workflow) {
        // Container engine check
        if (!(workflow.containerEngine in ['docker', 'singularity']))
            throw new Exception("A container engine is required — please run with '-profile docker' or '-profile singularity'")

        // Enum checks
        if (!params.valid_quantification_modes.contains(params.quantification_mode))
            throw new Exception("Invalid params.quantification_mode '${params.quantification_mode}' — must be one of: ${params.valid_quantification_modes.join(', ')}")

        // Numeric range checks
        if (params.resolution <= 0)
            throw new Exception("Invalid params.resolution '${params.resolution}' — must be a positive number")

        if (params.ndr != null && (params.ndr < 0 || params.ndr > 1))
            throw new Exception("Invalid params.ndr '${params.ndr}' — must be a float between 0 and 1")
    }

    static def validateVisiumSampleCount(samples) {
        def has_visium = samples.any { sample, path, meta -> meta.chemistry.startsWith('visium') }
        if (has_visium && samples.size() > 1)
            throw new Exception("Visium chemistry requires exactly 1 sample, but found ${samples.size()}")
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
    }

}
