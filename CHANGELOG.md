# Changelog

This file contains all notable changes to Bambu-Pipe. 

---

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
