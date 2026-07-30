include { SEURAT_SINGLE_SAMPLE } from '../modules/seurat/standard/single_sample_clustering.nf'
include { SEURAT_MULTI_SAMPLE  } from '../modules/seurat/standard/multi_sample_clustering.nf'

workflow CLUSTERING {
    take:
    ch_se_gene_counts
    ch_n_samples

    main:
    ch_se_branched = ch_se_gene_counts
        .combine(ch_n_samples)
        .branch { _se, n ->
            single: n == 1
            multi:  n > 1
        }

    SEURAT_SINGLE_SAMPLE(ch_se_branched.single.map { se, _n -> se })
    SEURAT_MULTI_SAMPLE(ch_se_branched.multi.map { se, _n -> se })

    emit:
    clusters = SEURAT_SINGLE_SAMPLE.out.clusters.mix(SEURAT_MULTI_SAMPLE.out.clusters)
}
