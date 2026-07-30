include { SEURAT_SINGLE_SAMPLE } from '../modules/seurat/standard/single_sample_clustering.nf'
include { SEURAT_MULTI_SAMPLE  } from '../modules/seurat/standard/multi_sample_clustering.nf'

workflow CLUSTERING {
    take:
    ch_gene_counts
    ch_col_data
    ch_n_samples

    main:
    ch_branched = ch_gene_counts
        .combine(ch_col_data)
        .combine(ch_n_samples)
        .branch { _gene_counts, _col_data, n ->
            single: n == 1
            multi:  n > 1
        }

    SEURAT_SINGLE_SAMPLE(ch_branched.single.map { gene_counts, col_data, _n -> [gene_counts, col_data] })
    SEURAT_MULTI_SAMPLE(ch_branched.multi.map { gene_counts, col_data, _n -> [gene_counts, col_data] })

    emit:
    clusters = SEURAT_SINGLE_SAMPLE.out.clusters.mix(SEURAT_MULTI_SAMPLE.out.clusters)
}
