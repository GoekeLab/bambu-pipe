include { SEURAT_SINGLE_SAMPLE } from '../modules/seurat/single_sample.nf'
include { SEURAT_MULTI_SAMPLE  } from '../modules/seurat/multi_sample.nf'

workflow CLUSTERING {
    take:
    ch_se_gene_counts
    ch_sample_names
    ch_n_samples

    main:
    ch_se_branched = ch_se_gene_counts
        .combine(ch_n_samples)
        .branch { se, n ->
            single: n == 1
            multi:  n > 1
        }

    SEURAT_SINGLE_SAMPLE(ch_se_branched.single.map { se, n -> se })
    SEURAT_MULTI_SAMPLE(ch_se_branched.multi.map { se, n -> se }, ch_sample_names)

    emit:
    clusters = SEURAT_SINGLE_SAMPLE.out.clusters.mix(SEURAT_MULTI_SAMPLE.out.clusters)
    versions = SEURAT_SINGLE_SAMPLE.out.versions.mix(SEURAT_MULTI_SAMPLE.out.versions)
}
