# Changelog

This file contains all notable changes to Bambu-Pipe. 

---

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
