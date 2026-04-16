# Test Cases for bambu-singlecell-spatial Pipeline

## Input Entry Point Tests

| # | Name | Input type | Description | Expected |
|---|------|-----------|-------------|----------|
| 1 | fastq-basic | FASTQ | Single sample, ONT, 10x3v3 | Full pipeline runs end-to-end, produces `*_se.rds` |
| 2 | bam-basic | BAM | Pre-aligned, demultiplexed BAM | Skips alignment, starts from read class construction |
| 3 | rds-basic | RDS | Pre-computed read class file | Skips alignment + read class, starts from BAMBU |
| 4 | mixed-inputs | FASTQ + BAM + RDS | All three types in one samplesheet | Each branch runs correctly, all samples aggregated into combined output |
| 5 | bam-with-barcode-map | BAM | Non-demultiplexed BAM + `barcode_map` column | Demultiplexing applied correctly before read class |

---

## Chemistry / Technology Tests

| # | Name | Chemistry | Technology | Description | Expected |
|---|------|----------|-----------|-------------|----------|
| 6 | 3prime-ont | 10x3v3 | ONT | Standard 3' chemistry | Reads reverse complemented, flexiplex tolerance=13 |
| 7 | 5prime-ont | 10x5v2 | ONT | 5' chemistry | No reverse complement, TSO trimming, tolerance=8 |
| 8 | pacbio | 10x3v3 | PacBio | PacBio technology | minimap2 uses `splice:hq` preset |
| 9 | visium-standard | visium-v1 | ONT | Spatial Visium | Spatial metadata passed to BAMBU, coordinates in colData |
| 10 | visium-mixed | visium-v1 + 10x3v3 | ONT | Spatial + non-spatial samples combined | Non-visium sampleData set to NA, visium retains coordinates |
| 11 | chemistry-from-param | — | — | Chemistry/technology via `--chemistry`/`--technology` instead of samplesheet column | Correctly overrides per-row values |

---

## `early_stop_stage` Tests

| # | Name | `early_stop_stage` | Expected |
|---|------|-------------------|----------|
| 12 | stop-bam | `"bam"` | Pipeline stops after alignment; `*_demultiplexed.bam` published, no RDS, no GTF |
| 13 | stop-rds | `"rds"` | Pipeline stops after read class construction; `*_readClassFile.rds` published, no GTF |
| 14 | full-run | `null` (default) | Full pipeline completes; GTF + SE object published |

---

## `quantification_mode` Tests

| # | Name | `quantification_mode` | Expected |
|---|------|----------------------|----------|
| 15 | no-quant | `"no_quant"` | BAMBU runs (discovery + no-EM quant), SEURAT_CLUSTERING and BAMBU_EM skipped |
| 16 | em-only | `"EM"` | BAMBU + BAMBU_EM run, SEURAT_CLUSTERING skipped; `degBias=FALSE` in BAMBU_EM |
| 17 | em-clusters | `"EM_clusters"` (default) | Full pipeline; SEURAT_CLUSTERING runs, clusters passed to BAMBU_EM; `degBias=TRUE` |

---

## Multi-Sample Tests

| # | Name | Samples | Expected |
|---|------|---------|----------|
| 18 | single-sample | 1 | Output named `{sample}_se.rds` (not `combined_`) |
| 19 | multi-sample | 2+ | Output named `combined_se.rds`, all samples aggregated |
| 20 | multi-sample-mixed-chemistry | 2 samples, different chemistries | Each sample preprocessed with its own chemistry config |

---

## `bambu_path` Tests

| # | Name | `bambu_path` | Expected |
|---|------|-------------|----------|
| 21 | default-bambu | `null` | `library("bambu")` used in all 4 modules |
| 22 | custom-bambu | valid path to local bambu source | `devtools::load_all()` called in all 4 modules |
| 23 | invalid-bambu-path | non-existent path | `load_all()` fails with clear error |

---

## Parameter Validation Tests

| # | Name | Parameter | Input | Expected |
|---|------|-----------|-------|----------|
| 24 | missing-input | `--input` | omitted | Pipeline exits with error: "input is required" |
| 25 | missing-genome | `--genome` | omitted | Pipeline exits with error |
| 26 | missing-annotation | `--annotation` | omitted | Pipeline exits with error |
| 27 | invalid-quant-mode | `--quantification_mode` | `"foo"` | Validation error listing valid options |
| 28 | invalid-early-stop | `--early_stop_stage` | `"foo"` | Validation error |
| 29 | invalid-ndr | `--ndr` | `1.5` | Validation error: NDR must be between 0 and 1 |
| 30 | invalid-resolution | `--resolution` | `0` or negative | Validation error: resolution must be > 0 |
| 31 | missing-samplesheet-columns | samplesheet | `sample` column missing | Error indicating required column |

---

## NDR / Deduplication Tests

| # | Name | Description | Expected |
|---|------|-------------|----------|
| 32 | ndr-set | `--ndr 0.1` | NDR passed to BAMBU discovery call |
| 33 | ndr-null | `--ndr null` (default) | `NDR = NULL` in BAMBU R call |
| 34 | dedup-off | `--deduplicate_umis false` | `dedupUMI = FALSE` in BAMBU_CONSTRUCT_READ_CLASS |

---

## Regression / Edge Cases

| # | Name | Description | Expected |
|---|------|-------------|----------|
| 35 | qscore-filtering-off | `--qscore_filtering false` | `chopper` step skipped in PREPROCESS_FASTQ |
| 36 | process-by-chromosome-off | `--process_by_chromosome false` | `processByChromosome = FALSE` in read class construction |
| 37 | save-intermediates | `--save_intermediates true` | Intermediate files written to disk |
| 38 | retry-on-oom | Trigger OOM (exit 137) | Process retries with doubled memory |
| 39 | all-visium | All samples are Visium | `sampleData` non-null for all, spatial coordinates in output |
| 40 | all-non-visium | No Visium samples | `sampleData = NULL` passed to BAMBU |

---

## Notes

- The existing chr9 example data in `examples/` can serve as the test fixture for most cases.
- Highest priority cases to implement first: #1 (end-to-end FASTQ), #12–14 (early stop stages), #15–17 (quantification modes), and #24–31 (validation) — these cover the main branching logic with minimal data requirements.
