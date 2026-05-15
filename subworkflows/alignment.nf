include { MINIMAP_BUILD_INDEX } from '../modules/alignment/build_index.nf'
include { PAFTOOLS_GFF2BED   } from '../modules/alignment/gff2bed.nf'
include { MINIMAP_ALIGNMENT  } from '../modules/alignment/align.nf'

workflow ALIGNMENT {
    take:
    ch_unaligned_fastq
    ch_genome
    ch_annotation

    main:
    // ch_gate emits one item if fastq channel is non-empty, else emits nothing. Used to prevent MINIMAP_BUILD_INDEX and PAFTOOLS_GFF2BED
    // from running when there are no fastq samples to process
    ch_gate = ch_unaligned_fastq.first().map { _x -> true }

    // Build minimap2 index based on reference genome
    MINIMAP_BUILD_INDEX(ch_genome.combine(ch_gate).map { g, _gate -> g })

    // Convert gtf/gff annotation to bed format
    PAFTOOLS_GFF2BED(ch_annotation.combine(ch_gate).map { a, _gate -> a })

    // Minimap alignment
    MINIMAP_ALIGNMENT(ch_unaligned_fastq, MINIMAP_BUILD_INDEX.out.index.first(), PAFTOOLS_GFF2BED.out.first())

    emit:
    bam = MINIMAP_ALIGNMENT.out.bam
    versions = MINIMAP_BUILD_INDEX.out.versions.mix(MINIMAP_ALIGNMENT.out.versions.first())
}
