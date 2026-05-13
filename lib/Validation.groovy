class Validation {

    static def validateParams(params, workflow) {
        // Container engine check
        if (!(workflow.containerEngine in ['docker', 'singularity']))
            error "A container engine is required — please run with '-profile docker' or '-profile singularity'"

        // Required inputs
        if (!params.input)
            error "params.input is not set — please provide a path to a CSV samplesheet"
        if (!params.genome)
            error "params.genome is not set — please provide a path to the reference genome FASTA file"
        if (!params.annotation)
            error "params.annotation is not set — please provide a path to the reference annotation GTF file"

        // Enum checks
        if (!params.valid_quantification_modes.contains(params.quantification_mode))
            error "Invalid params.quantification_mode '${params.quantification_mode}' — must be one of: ${params.valid_quantification_modes.join(', ')}"

        if (!(params.bam_only instanceof Boolean))
            error "Invalid params.bam_only '${params.bam_only}' — must be true or false"

        // Numeric range checks
        if (params.resolution <= 0)
            error "Invalid params.resolution '${params.resolution}' — must be a positive number"

        if (params.ndr != null && (params.ndr < 0 || params.ndr > 1))
            error "Invalid params.ndr '${params.ndr}' — must be a float between 0 and 1"
    }

}
