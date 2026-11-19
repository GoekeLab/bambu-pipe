# Changelog

This file contains all notable changes to Bambu-Pipe. 

---

## [v0.10.4] - 2026-11-19

### Added
- `--flexiplex_full_rank_search` for 10x Single Cell chemistries, which makes flexiplex-filter search the whole barcode rank curve for the inflection point (`-u 0`); the default search stops at the 95th percentile of barcode counts and can remove real cells in some samples
- `single_cpu` process label (1 CPU) for single-threaded processes
- `max_forks_rc` developer parameter (default 4; 100 with the `hpc` profile, where each task is its own job), the number of `BAMBU_CONSTRUCT_READ_CLASS` tasks that run at once; each task loads its whole BAM into memory, so running every sample at once multiplied the peak memory
- `--min_read_fraction_by_gene` (default 0.05), the minimum fraction of its gene's reads a novel transcript needs in at least one sample to be kept during transcript discovery
- `--seurat_dims` (default 30), the number of principal components used for clustering in every workflow
- `harmony_seed` developer parameter (default 42), the random seed for Harmony, so multi-sample clusters are reproducible

### Changed
- With the developer parameter `bambu_path` set, its directory is now mounted into the containers (singularity and docker profiles), so the local Bambu can be loaded
- Seurat container `bambu-pipe-seurat:1.0.1`: R packages are now installed at fixed versions with `remotes::install_version`, with the same versions as 1.0.0 (Seurat 5.4.0, SeuratObject 5.3.0, harmony 1.2.4, R.utils 2.13.0, Banksy 1.6.0 from Bioconductor 3.22, SeuratWrappers at `9586c29`)
- The filtered BAM produced for Visium and Visium HD samples (`<sample>_filtered.bam` and `.bai`) is now always published to `bam/`, instead of `intermediate_bam/` only with `--save_intermediates`
- Transcript discovery (`BAMBU_TRANSCRIPT_DISCOVERY`, `BAMBU_TRANSCRIPT_DISCOVERY_VISIUM_HD`) now runs on a fixed 2 cores (`discovery_ncore`), instead of the `medium_cpu` allocation (16–64)
- CPU labels no longer grow with retries: `low_cpu` is 4, `medium_cpu` 16 and `high_cpu` 32 CPUs on every attempt; retries still double memory and extend time
- `BAMBU_CONSTRUCT_READ_CLASS`, `BAMBU_PREPARE_ANNOTATION`, `PAFTOOLS_GFF2BED`, `EXTRACT_SPOT_BARCODES`, `AGGREGATE_BINS_VISIUM_HD`, `MAP_CLUSTERS_TO_2UM_SPOTS`, `VISIUM_BUILD_TISSUE_POSITIONS`, `CONVERT_BARCODE_MAPPINGS`, `CONVERT_TISSUE_POSITIONS` and `SPOT_BIN_MAPPINGS` now use `single_cpu`, instead of `low_cpu` (4–16), as they run on one core
- `BAMBU_CONSTRUCT_READ_CLASS` now runs bambu at its default `ncore` (1), instead of `ncore = task.cpus`; each task processes one BAM, so the extra cores were unused
- Failed steps are retried at most twice (was 3 times), and R errors (exit code 1) in Bambu and Seurat steps only once, since an ordinary R error exits with the same code as a BiocParallel worker running out of memory; the largest memory request drops from 1024 GB to 512 GB
- Single- and multi-sample Seurat clustering are now one process, `SEURAT_CLUSTERING` (replacing `SEURAT_SINGLE_SAMPLE`, `SEURAT_MULTI_SAMPLE` and the `CLUSTERING` subworkflow). All runs now use Seurat's default 2,000 variable genes (single-sample runs used 2,500), 30 principal components by default (single-sample runs used 15) and the `high_mem` label, and `seurat_obj.rds` names the cluster column `clusters` (multi-sample runs used `harmony_clusters`)
- **Breaking:** a `--manual_clustering` restart now takes `--clusters`, `--quant_data` and `--extended_annotations`, instead of a samplesheet (`clusters_path`, `quant_data_path`) and an `extended_annotations.rds` given to `--annotation`. A restart takes no `--input`, `--genome` or `--annotation`; passing them, or the restart flags without `--manual_clustering`, stops the pipeline with an error. `--quantification_mode` is ignored on a restart
- `SEURAT_VISIUM_HD` is renamed `SEURAT_CLUSTERING_VISIUM_HD`; Visium HD clustering with `--banksy false` now uses 30 principal components by default (was 15)
- The quantification processes (`BAMBU_CLUSTER_LEVEL_QUANTIFICATION`, `BAMBU_SINGLE_CELL_QUANTIFICATION`, `BAMBU_SPOT_LEVEL_QUANTIFICATION`) no longer take the reference genome, which Bambu does not use when quantifying

### Removed
- `--deduplicate_umis`, which has had no effect since v0.9: Bambu always deduplicates UMIs in single-cell mode
- `--chemistry` and `--technology`, which have had no effect since v0.9; set `chemistry` and `technology` in the samplesheet
- `seurat_dim_single` and `seurat_dim_multi` developer parameters, replaced by `--seurat_dims`
- The unused fusion detection code (`bin/JAFFAL.groovy`, `bin/fusion_detection*.R`, `jaffal/`) and the `fusion_mode`, `jaffal_ref_dir` and `jaffal_code_dir` developer parameters

### Fixed
- The pipeline failed to compile on Nextflow 26.10 ("Process `…` cannot be used as a variable"): `main.nf` no longer enables static typing (`nextflow.enable.types`), whose rules Nextflow 26.10 now enforces; the typed `params` block, and the conversion of command-line values to each parameter's type, are unchanged
- A sample whose file had an unsupported extension (e.g. `.FASTQ.GZ`, `.fq.bz2`, `.cram`) was silently dropped from the run; the pipeline now stops with an error naming the sample and the supported extensions (`.fastq`, `.fq`, `.fastq.gz`, `.fq.gz`, `.bam`)
- `--annotation` now accepts GTF only (`.gtf`/`.gtf.gz`). A GFF3 file was read by Bambu without an error into a wrong annotation (each exon a separate transcript), and the README listed `.gff`/`.gff.gz` as supported
- A `--manual_clustering` restart samplesheet with several rows ran several quantifications whose outputs overwrote each other, and a missing column failed with an unclear error; restarts now take one file per flag
- Duplicate sample names are now rejected before any sample is processed; previously the samples were aligned and processed twice before transcript discovery failed on a file name collision
- The samplesheet is now validated as a whole before any sample is processed (each row, duplicate names, the one-sample rule for Visium, `--loupe_alignment` without a Visium sample, and an empty samplesheet); previously the Visium checks ran beside the first alignment jobs
- `--banksy_lambda` (0–1) and `--banksy_k_geom` (at least 1) are now validated at launch
- Multi-sample Harmony integration now corrects for sequencing technology and capture chemistry, as documented: Seurat's `IntegrateLayers(method = HarmonyIntegration)` ignored `group.by.vars` and corrected by sample. Harmony now runs through `RunHarmony` on whichever of technology and chemistry differ between samples (it is skipped when neither does), and is seeded, so reruns give the same clusters
- The example Seurat objects (`examples/seurat_obj_single_sample.rds`, `examples/seurat_obj_multi_sample.rds`) are regenerated, and the README's visualisation examples now use the `clusters` and `sampleName` columns (the multi-sample example used a `sample` column that runs no longer produce)
- Rerunning into an existing `--output_dir` kept the earlier run's `pipeline_info/` reports ("Failed to render execution report"); the timeline, report, trace and DAG are now overwritten
- The README described the cluster-level count columns as `clusterId` for single-sample runs; they are always `sampleName_clusterId` (e.g. `sample1_cluster_0`), and barcodes are always `sampleName_barcode`

## [v0.10.3] - 2026-10-01

### Changed
- Increased the minimum adapter overlap in the `PREPROCESS_FASTQ` adapter re-search from the cutadapt default (3 bp) to 10 bp, reducing the number of valid reads discarded due to short chance matches at read ends
- Shortened the polyT right flank for 3' and Visium chemistries (`flexiplex_config.csv`) from 30 T to 9 T, matching the flexiplex presets; ONT reads rarely call the full 30 T homopolymer, so few reads matched the flank during barcode discovery
- Visium samples (`visium-v*`) are now demultiplexed directly against the spot whitelist, as recommended by flexiplex; barcode discovery and flexiplex-filter are skipped, replacing the `-u 0` knee detection workaround
- Updated the supported chemistry names in the README to the official 10x Genomics assay names, with a link to the 10x Genomics long-read compatibility article
- Updated the flexiplex flank (`-f`) and barcode (`-e`) edit distances to `-f 8` and `-e 2`, matching the shortened 9 T polyT flank and standardising with the flexiplex documentation; these are now set per chemistry in `flexiplex_config.csv` (`flank_max_edit_distance`, `barcode_max_edit_distance`), replacing the `flexiplex_f_5prime`, `flexiplex_f_3prime` and `flexiplex_e` developer parameters
- Moved the 10x config files out of `assets/10x_config/` into `assets/`, and renamed `adapter_seq_config.csv` to `cutadapt_config.csv` and `flank_seq_config.csv` to `flexiplex_config.csv`

### Removed
- `visium-v4` and `visium-v5` (Visium CytAssist Spatial Gene Expression) chemistries, as these probe-based assays do not produce full-length transcripts

### Fixed
- Corrected the 10x3v3 and 10x3v4 end adapter (`rev_primer_f`/`rev_primer_r` in `cutadapt_config.csv`) to `CCCATGTACTCTGCGTTGATACCACTGCTT`, matching 10x3v2 and Visium; the previous sequence does not occur in these reads, so the adapter was rarely trimmed
- Resolved an edge case where `PREPROCESS_FASTQ` failed if the input FASTQ filename had the same length as the barcode (16) or UMI (10/12) pattern; the unquoted `?` wildcards in the flexiplex pattern were expanded by the shell to match the filename, which was not intended

## [v0.10.2] - 2026-09-23

### Fixed
- `--output_dir` is now typed as `String` instead of `Path`; typed `Path` params must already exist, so a new output directory failed validation when the trace, timeline, report and DAG reports (which create it first) were disabled

## [v0.10.1] - 2026-08-26

### Added
- `--loupe_alignment` for Visium Spatial Gene Expression samples, taking the manual alignment `.json` exported after fiducial alignment and tissue detection in Loupe Browser. Required for `visium-v*` samples
  - Out-of-tissue barcodes are filtered from the BAM before transcript discovery and quantification
  - Introduce the `VISIUM_BUILD_TISSUE_POSITIONS` module, which builds the tissue positions file for `visium-v*` samples; the spatial metadata in this file is attached to the `colData` of the `SummarizedExperiment` objects
  - Loupe alignment example (`examples/loupe_alignment_visium_example.json`), used by the `test_visium` smoke test
- `CB`/`UB` tags in aligned BAM files (minimap2 `-y`), carrying the barcode and UMI from the FASTQ header comments

### Changed
- The spatial metadata attached to Visium `SummarizedExperiment` objects now follows the Space Ranger tissue positions format (`barcode`, `in_tissue`, `array_row`, `array_col`, `pxl_row_in_fullres`, `pxl_col_in_fullres`), replacing the `x_coordinate`/`y_coordinate` columns
- `FILTER_BARCODED_BAM` moved to `modules/prepare_input/shared/` and is used by both the standard Visium and Visium HD workflows; it now fails when no reads remain after filtering
- User-supplied BAM files must carry the barcode and UMI in the `CB`/`UB` tags; barcodes encoded in the read name are no longer supported

### Fixed
- flexiplex-filter's knee detection could discard most in-tissue barcodes for Visium samples; the inflection search now covers the whole barcode rank curve (`-u 0`) for `visium-v*` chemistries

## [v0.10.0] - 2026-08-17

### Added
- Visium HD workflow (`--visium_hd`), run as a single sample from a Spaceranger-aligned, barcode-tagged BAM
  - Transcript discovery and read-to-transcript assignment at the native 2 µm resolution, with counts aggregated to every bin listed in `--bins`
  - `--barcode_mappings` for the Spaceranger `barcode_mappings.parquet`, used to assign 2 µm spots to bins
  - Out-of-tissue reads filtered from the BAM using the 2 µm `tissue_positions.parquet`
  - Spatially aware clustering with Banksy (`--banksy`, `--banksy_lambda`, `--banksy_k_geom`), or gene expression alone
  - `--clustering_bin` to select the resolution to cluster at; cluster labels are expanded back to 2 µm spots for quantification
  - Spot-level quantification at every resolution under `--quantification_mode EM`
  - `test_visium_hd` smoke test profile with synthetic example data
- `--manual_clustering` to restart the pipeline from cluster assignments generated outside the pipeline, for both standard and Visium HD runs
- `test_sc_quant_data` and `test_visium_hd_quant_data` smoke test profiles covering the manual clustering restart
- Self-hosted `bambu` and `seurat` container images published to `ghcr.io/goekelab`, built by the `build_container.yml` GitHub Actions workflow
- Shared R helpers in `bin/` for transcript discovery, Seurat object creation, count saving, and cluster mapping

### Changed
- Renamed `--resolution` to `--seurat_resolution`
- Renamed the `--quantification_mode` option `EM_clusters` to `clusteredEM`, matching the `bambu.singlecell` API
- `quant_data.rds` and `extended_annotations.rds` are now always published to `intermediate_R/`, so a manual clustering run can restart from them
- Seurat objects are built from the published count directories and Bambu's `colData` instead of the `SummarizedExperiment`
- `clusters.rds` is now a named vector of `id -> cluster` label, replacing the per-sample list of `CompressedCharacterList`
- Cluster-level quantification moved into a single module shared by the standard and Visium HD workflows
- Restructured modules into `standard/`, `visium_hd/`, and `shared/` directories
- Smoke tests now run on pull request and manual dispatch only, with in-progress runs cancelled on a new push

## [v0.9.1] - 2026-05-20

### Added
- Harmony batch correction for multi-sample Seurat clustering
- Processing of CB/UB tagged custom BAM files
- GitHub Actions workflow to run smoke test on push and pull request to `main` and `devel` branches

### Changed
- Upgraded pipeline to support Nextflow version `26.04.0` and above

## [v0.9-beta] - 2026-05-11

### Added
- Quality score filtering with Chopper
- Primer removal with Cutadapt
- Reverse complement FASTQ utility script (`bin/reverse_complement_fastq.py`) to enable stranded alignment in minimap2
- Automatic extraction of 10x barcodes and spatial coordinates from the Spaceranger container
- Support for multiple sample analysis using Nextflow parallelisation
- Modularised codebase into discrete modules and subworkflows (`modules/bambu/`, `modules/alignment/`, `modules/prepare_input_standard/`)
- External 10x config asset files for barcode coordinates, adapter sequences, and flank sequences
- `params` block centralising all pipeline parameters (previously defined in main.nf)
- `process` block with dynamic retry strategy
- Resource labels for CPU, memory, and time
- HPC execution profile (`conf/`) to support parallelisation on high performance computing systems
- Minimal end-to-end smoke test (`conf/test.config`)
- Manifest block with author and version metadata
- Emit software versions in a .yml file
- Input validation via `lib/Validation.groovy`
- `quantification_mode` parameter to control quantification strategy (`no_quant`, `EM`, `EM_clusters`)
- Seurat clustering as a dedicated process (`SEURAT_CLUSTERING`) for cluster-based EM quantification
- Joint clustering across all samples on a combined gene counts matrix (previously per-sample)
- Cluster output restructured to an ordered list of `CompressedCharacterList`, one per sample in `quantData` order (previously a flat single CCL mixing all samples)
- `SEURAT_CLUSTERING` now takes gene counts matrix and sample names as inputs instead of the full `quantData` object
- `clusterCells` helper inlined into the process (previously sourced from `bin/utilityFunctions.R`)
- `early_stop_stage` parameter to terminate the pipeline after BAM or RDS generation

### Changed
- Migration to Wave community containers (previously root-level `Dockerfile`)
- Removed deprecated parameters
- Removed hardcoded values and redundant code
- Simplified input logic using a single samplesheet
- Enhanced input validation check

---

## [v0.1-beta] - 2025-05-19

### Added
- Initial pipeline release
