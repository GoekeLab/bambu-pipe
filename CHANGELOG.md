# Changelog

This file contains all notable changes to Bambu-Pipe. 

---

## [1.0.0] - 2026-MM-DD

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
- Cluster output restructured to an ordered list of `CompressedCharacterList`s, one per sample in `quantData` order (previously a flat single CCL mixing all samples)
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

## [BETA] - 2023-05-03

### Added
- Initial pipeline release
