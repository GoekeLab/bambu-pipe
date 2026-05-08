# **Context-Aware Transcript Quantification from Long-Read Single-Cell and Spatial Transcriptomics Data**
### **Content**
- [Overview](#overview)
- [Installation](#installation)
- [General Usage](#general-usage)
- [Parameters](#parameters)
- [Output](#output)
- [Spatial Analysis](#spatial-analysis)
- [Advanced Usage](#advanced-usage)
- [Additional Information](#additional-information)
- [Release History](#release-history)
- [Citation](#citation)
- [Contributors](#contributors)


### **Overview**
This pipeline performs context-aware transcript discovery and quantification from long-read single-cell and spatial transcriptomics data. The workflow is divided into three stages:

**Preprocessing**
1. (Optional) Quality score filtering with [Chopper](https://github.com/wdecoster/chopper)
2. Barcode/UMI identification and demultiplexing with [Flexiplex](https://davidsongroup.github.io/flexiplex/)
3. Primer removal with [Cutadapt](https://cutadapt.readthedocs.io/en/stable/)

**Alignment**

4. Genome alignment with [Minimap2](https://lh3.github.io/minimap2/minimap2.html)

**Transcript Discovery and Quantification**

5. Read class construction and transcript discovery with [Bambu](https://github.com/GoekeLab/bambu/tree/BambuDev) (performed jointly across all samples)
6. (Optional) Transcript quantification with Bambu using one of two modes:
   - Cluster-level EM: Gene expression-based cell clustering with [Seurat](https://github.com/satijalab/seurat) across all samples, followed by per-sample cluster-level transcript quantification
   - Single-cell EM: Per-cell transcript quantification 

### **Installation** 
Install the following dependencies before running the pipeline:
- [Nextflow](https://www.nextflow.io/docs/latest/install.html) 
- [Docker](https://docs.docker.com/engine/install/ubuntu/) (or [Singularity](https://docs.sylabs.io/guides/3.0/user-guide/installation.html) if you do not have user permissions for Docker). 

The latest version for each dependency is recommended. 

### **General Usage** 
To run the pipeline, you must provide a samplesheet, reference genome, and reference annotation file as input. The pipeline performs transcript discovery and quantification on either a single sample or multiple samples based on the number of samples specified in the samplesheet. Refer to the [Parameters](#parameters) and Samplesheet (CSV) sections below for more details. 

**Running the pipeline**

Use the command below to run the pipeline on the test data provided in `examples/`
``` 
nextflow run main.nf \
  --input examples/samplesheet_test_fastq.csv \
  --genome examples/GRCh38.primary_assembly.genome.chr9_1_1000000.fa \
  --annotation examples/gencode.v49.primary_assembly.annotation.chr9_1_1000000.gtf \
  -profile singularity,hpc
``` 

**Samplesheet (CSV)**

The pipeline requires a `.csv` formatted samplesheet to define the input data. This file is mandatory, regardless of the number of samples being processed. Each row in the samplesheet represents a single sample and its corresponding file path and metadata. 

*Required Columns*

The samplesheet must include the following columns:
- `sample`: sample name (no spaces or non-alphanumeric characters)
- `path`: path to the input file (FASTQ, BAM, or RDS)
- `chemistry`: 10x library chemistry (see Supported 10x Library Chemistries below)
- `technology`: sequencing technology (`ONT` or `PacBio`)

Note: The first row of the samplesheet must be a header containing the exact column names: `sample`, `path`, `chemistry`, and `technology`. 

*Supported Input Formats*

The pipeline is designed to be flexible. Depending on your starting point in the workflow, the `path` column can point to the following file types:
- **FASTQ**: Raw reads (compressed `.gz` or uncompressed)
- **BAM**: Demultiplexed and aligned reads
- **RDS**: Pre-processed bambu read class objects

For more details on starting the pipeline from specific stages, please refer to the [Advanced Usage](#advanced-usage) section. 

*Example Samplesheet (Single Sample)*
```csv
sample,path,chemistry,technology
sample1,path/to/sample1_fastq.gz,10x3v2,ONT
```

*Example Samplesheet (Multiple Samples)*
```csv
sample,path,chemistry,technology
sample1,path/to/sample1_fastq.gz,10x3v2,ONT
sample2,path/to/sample2_fastq.gz,10x3v3,PacBio
sample3,path/to/sample3_fastq.gz,10x3v4,ONT
```

Note: Example samplesheets are provided in `examples/`.
If all samples share the same library chemistry and/or sequencing technology, you may omit the `chemistry` and `technology` columns and use the `--chemistry` and `--technology` flags instead.


*Supported 10x Library Chemistries*

The following single cell and spatial library chemistries are supported. Please specify the sample chemistry in the samplesheet as shown:
- `10x3v2` (Single Cell 3' v2)
- `10x3v3` (Single Cell 3' v3 & Next GEM Single Cell 3' v3.1)
- `10x3v4` (GEM-X Single Cell 3' v4)
- `10x5v2` (Single Cell 5' v2)
- `10x5v3` (GEM-X Single Cell 5' v3)
- `visium-v1` (Visium Spatial Gene Expression Slide 6.5 mm; serial prefix V1)
- `visium-v2` (Visium Spatial Gene Expression Slide 6.5 mm; serial prefix V2)
- `visium-v3` (Visium Spatial Gene Expression Slide 6.5 mm; serial prefix V3)
- `visium-v4` (Visium CytAssist Spatial Gene Expression Slide 6.5 mm; serial prefix V4)
- `visium-v5` (Visium CytAssist Spatial Gene Expression Slide 11mm; serial prefix V5)

**Pipeline Configuration**

*Nextflow Profiles*

To configure the executor and container, pass profile types via the `-profile` argument.

- Container profiles:
  - `singularity`: use Singularity images (recommended on HPC systems)
  - `docker`: use Docker images

- Executor profiles:
  - `hpc`: execute on an HPC system (default executor: `slurm`; edit `process.executor` in `nextflow.config` to switch to `pbs`, `sge`, etc.)
  - `local`: execute on a local machine with reduced resource limits — not recommended for full-size datasets

### **Parameters**

**Mandatory**
- `--input` [string]: Path to the samplesheet .csv file 
- `--genome` [string]: Path to the reference genome .fa or .fasta file 
- `--annotation` [string]: Path to the reference annotation .gtf or .gff file 

**Optional**
- `--output_dir` [string, default: 'output']: Path to the output directory
- `--chemistry` [string, default: null]: Specify if all samples in the samplesheet share the same library chemistry 
- `--technology` [string, default: null]: Specify if all samples in the samplesheet share the same sequencing technology
- `--early_stop_stage` [string, default: null]: Stop the pipeline early and output intermediate files (see Advanced Usage section). Options:
  - "bam": Stops pipeline after minimap2 alignment
  - "rds": Stops pipeline after Bambu read class construction
- `--qscore_filtering` [boolean, default: true]: Enable or disable quality score filtering of reads
- `--ndr` [float, default: null]: NDR threshold for Bambu transcript discovery. If not set, Bambu will recommend a suitable value
- `--deduplicate_umis` [boolean, default: true]: If true, Bambu will perform UMI deduplication 
- `--quantification_mode` [string, default: "EM_clusters"]: Quantification mode for transcript counts. Available options are:
  - "no_quant": Transcript quantification is not performed
  - "EM": Performs transcript quantification for each cell/spatial coordinate
  - "EM_clusters": Performs gene expression-based cell clustering using [Seurat](https://satijalab.org/seurat/), followed by transcript quantification at the cluster level
- `--resolution` [float, default: 0.8]: Seurat clustering resolution

> **Warning:** We currently recommend processing one sample at a time if `--quantification_mode` is set to `EM_clusters`. Mixing multiple samples is not advised, as batch effect correction across samples during clustering has not yet been implemented and will be available in a future release.

### **Output**
All outputs from the pipeline are written to the directory specified by the `--output_dir` parameter. The pipeline produces per-sample alignment files, per-sample read class files used by Bambu, and the combined transcript discovery and quantification results. The examples below show the output directory structure for both single and multi-sample runs:

*Output Structure (Single Sample)*
```
output/
├── bam/                                
│   ├── sample1_demultiplexed.bam
│   └── sample1_demultiplexed.bam.bai
│
├── read_class/                                
│   └── sample1_read_class.rds
│
├── extended_annotations.gtf
├── se_unique_counts.rds
├── se_gene_counts.rds
│
│   # single-cell EM:
├── se_transcript_counts_singlecell.rds
│
│   # clustered EM:
├── se_transcript_counts_clusters.rds
├── se_gene_counts_clusters.rds
│
└── software_versions.yml
```

*Output Structure (Multiple Samples)*
```
output/
├── bam/                                
│   ├── sample1_demultiplexed.bam
│   ├── sample1_demultiplexed.bam.bai
│   ├── sample2_demultiplexed.bam
│   └── sample2_demultiplexed.bam.bai
│
├── read_class/                                
│   ├── sample1_read_class.rds
│   └── sample2_read_class.rds
│
├── extended_annotations.gtf
├── se_unique_counts.rds
├── se_gene_counts.rds
│
│   # single-cell EM:
├── se_transcript_counts_singlecell.rds
│
│   # clustered EM:
├── se_transcript_counts_clusters.rds
├── se_gene_counts_clusters.rds
│
└── software_versions.yml
``` 

**Description of the Output Files**
| File | Description 
|---|---
| <sample_name>_demultiplexed.bam | BAM file containing demultiplexed, trimmed and aligned reads
| <sample_name>_demultiplexed.bam.bai | BAM index for the corresponding BAM file
| <sample_name>_read_class.rds |  An intermediate metadata file used by Bambu that contains the constructed read classes. This file can be used as input in subsequent runs to bypass the initial preprocessing and alignment steps. 
| extended_annotations.gtf | A `.gtf` file containing the novel transcripts discovered by Bambu as well as the reference annotations provided by the user.
| se_unique_counts.rds | A [RangedSummarizedExperiment](https://www.rdocumentation.org/packages/SummarizedExperiment/versions/1.2.3/topics/RangedSummarizedExperiment-class) object containing transcript-level unique counts at single-cell resolution, produced prior to EM quantification.
| se_gene_counts.rds | A RangedSummarizedExperiment object containing gene-level counts at single-cell resolution
| se_transcript_counts_singlecell.rds | A RangedSummarizedExperiment object containing per-cell transcript counts after EM quantification. Only produced when EM is run without cluster assignments.
| se_transcript_counts_clusters.rds | A RangedSummarizedExperiment object containing cluster-level transcript counts after EM quantification. Columns follow the `sampleName_clusterId` naming convention. Only produced when EM is run with clustering.
| se_gene_counts_clusters.rds | A RangedSummarizedExperiment object containing cluster-level gene counts. Only produced when EM is run with clustering.
| software_versions.yml | A YAML file listing the versions of all software tools used during the pipeline run.

**Count Matrices**

The [RangedSummarizedExperiment](https://www.rdocumentation.org/packages/SummarizedExperiment/versions/1.2.3/topics/RangedSummarizedExperiment-class) object contains four distinct types of count matrices, which can be accessed in R using the `assays()` function. Depending on your analysis requirements you can choose from the following:
- `counts`: expression estimates
- `CPM`: sequencing depth normalised estimates
- `fullLengthCounts`: estimates of read counts mapped as full length reads for each transcript
- `uniqueCounts`: counts of reads that are uniquely mapped to each transcript 

Note: In `se_unique_counts.rds`, unique counts are stored under the `counts` assay, not `uniqueCounts`.


### **Spatial Analysis**
The pipeline applies the same processing steps to both single-cell and spatial samples. However, for spatial data, the generated `SummarizedExperiment` object is appended with spatial mapping information, which is stored in `colData`.  

**Example - Spatial Mapping Information (`visium-v*`)**:

For `visium-v*` samples, `colData` contains the spatial barcode and the corresponding X and Y spatial coordinates. 

| Barcode            | X coordinate | Y coordinate| 
|:---|:---|:---|
| AAACAACGAATAGTTC | 17 | 1 |
| AAACAAGTATCTCCCA  | 103 | 51 |
| AAACAATCTACTAGCA | 44 | 4 |


### **Visium HD Spatial Analysis (Under Development)**

This feature is still under development and will be released in a future update.


### **Fusion Transcript Analysis (Under Development)**
This feature is still under development and will be released in a future update.


### **Advanced Usage**

**Minimal End-to-End Smoke Test**

Example data and pre-configured profiles are provided in `examples/` to run the pipeline end-to-end automatically without preparing your own data. The commands below must be run from the project's root directory. Combine the profile `test_base` with one of the profiles below and a container profile (`singularity` or `docker`).

| Profile | Description |
|---|---|
| `test_fastq` | Single-sample ONT run from raw reads |
| `test_bam` | Single-sample ONT run from demultiplexed BAM |
| `test_rds` | Single-sample ONT run from pre-computed read class object |
| `test_multi` | Multi-sample run with ONT and PacBio samples |

```bash
# Test from FASTQ input
nextflow run . -profile test_base,test_fastq,singularity

# Test from BAM input
nextflow run . -profile test_base,test_bam,singularity

# Test from a pre-computed RDS
nextflow run . -profile test_base,test_rds,singularity

# Test with multiple samples (ONT + PacBio)
nextflow run . -profile test_base,test_multi,singularity
```

The output files from the smoke tests are written to `.smoke_test/<profile>/output/`.

**Stopping the Pipeline Early**

The `--early_stop_stage` parameter allows you to stop the pipeline at an intermediate stage and save the outputs for later use. This is useful when you want to inspect intermediate files or when you plan to re-run downstream steps separately.

- `--early_stop_stage bam`: Stops after genome alignment. BAM files are saved to `output/bam/`.
- `--early_stop_stage rds`: Stops after Bambu read class construction. Read class `.rds` files are saved to `output/read_class/`.

```bash
# Stop after read class construction
nextflow run main.nf \
  --input samplesheet.csv \
  --genome reference.fa \
  --annotation reference.gtf \
  --early_stop_stage rds \
  -profile singularity,hpc
```

**Restarting from a Specific Stage**

Because the pipeline accepts FASTQ, BAM, or RDS files as input, you can restart from any intermediate stage by providing the corresponding files in your samplesheet. This avoids re-running expensive preprocessing and alignment steps when they have already been completed.

*Example: Incremental sample addition*

A common use case is to process an initial set of samples through to read class `.rds` files, then re-run the pipeline once additional samples are available. 

*Step 1* — Run the first batch of samples from FASTQ to `.rds`:
```bash
nextflow run main.nf \
  --input samplesheet_batch1.csv \
  --genome reference.fa \
  --annotation reference.gtf \
  --early_stop_stage rds \
  -profile singularity,hpc
```

This produces `output/read_class/sample1_read_class.rds`, `output/read_class/sample2_read_class.rds`, etc.

*Step 2* — When new samples are ready, run all samples together. Point the `path` column at the existing `.rds` files for the original samples and at the new FASTQ/BAM files for the new samples:

```csv
sample,path,chemistry,technology
sample1,output/read_class/sample1_read_class.rds,10x3v3,ONT
sample2,output/read_class/sample2_read_class.rds,10x3v3,ONT
sample3,path/to/sample3.fastq.gz,10x3v3,ONT
```

```bash
nextflow run main.nf \
  --input samplesheet_all.csv \
  --genome reference.fa \
  --annotation reference.gtf \
  -profile singularity,hpc
```

The pipeline will skip preprocessing and alignment for `sample1` and `sample2`, process `sample3` from FASTQ through to `.rds`, and then perform transcript discovery and quantification jointly across all three samples.

**Manual Clustering (Under Development)**

Currently, cell clustering is performed automatically as part of the pipeline. In a future release, a tutorial will be provided that allows users to stop the pipeline after transcript discovery, perform their own custom clustering, and then resume the pipeline to run Bambu transcript quantification using their cluster assignments.

### **Additional Information**
UMI correction is done at the barcode level. The longest read for each unique barcode-UMI combination is kept for analysis.

### **Release History** 

- v0.1-beta: 2025-May-19
- v0.9-beta: 2026-May-11


### **Citation**
If you use this pipeline, please cite our paper:

Sim, A., Ling, M. H., Chen, Y., Lu, H., See, Y. X., Perrin, A., Leng Agnes, O. B., Cao, E. Y., Chia, B., Liu, J., Wüstefeld, T., Shin, J. W., & Göke, J. (2025). Isoform-level discovery, quantification and fusion analysis from single-cell and spatial long-read RNA-seq data with Bambu-Clump. https://doi.org/10.1101/2024.12.30.630828

The following are citations for the other tools used in this pipeline:

#### Chopper
De Coster Wouter, & Rademakers, R. (2023). NanoPack2: Population scale evaluation of long-read sequencing data. Bioinformatics, 39(5). https://doi.org/10.1093/bioinformatics/btad311

#### Cutadapt
Martin, M. (2011). Cutadapt removes adapter sequences from high-throughput sequencing reads. EMBnet.journal, 17(1), 10. https://doi.org/10.14806/ej.17.1.200

#### Flexiplex
Cheng, O., Ling, M. H., Wang, C., Wu, S., Ritchie, M. E., Göke, J., Amin, N., & Davidson, N. M. (2024). Flexiplex: a versatile demultiplexer and search tool for omics data. Bioinformatics, 40(3). https://doi.org/10.1093/bioinformatics/btae102

#### Minimap2
Li, H. (2021). New strategies to improve minimap2 alignment accuracy. Bioinformatics, 37(23), 4572–4574. https://doi.org/10.1093/bioinformatics/btab705

#### Samtools
Danecek, P., Bonfield, J. K., Liddle, J., Marshall, J., Ohan, V., Pollard, M. O., Whitwham, A., Keane, T., McCarthy, S. A., Davies, R. M., & Li, H. (2021). Twelve years of SAMtools and BCFtools. GigaScience, 10(2). https://doi.org/10.1093/gigascience/giab008

#### Seurat
Hao, Y., Stuart, T. A., Kowalski, M. H., Choudhary, S., Hoffman, P., Hartman, A., Srivastava, A., Molla, G., Shaista Madad, Fernandez-Granda, C., & Rahul Satija. (2023). Dictionary learning for integrative, multimodal and scalable single-cell analysis. Nature Biotechnology. https://doi.org/10.1038/s41587-023-01767-y

### **Contributors**
This package is developed and maintained by [Andre Sim](https://github.com/andredsim), [Chin Hao Lee](https://github.com/ch99l), [Min Hao Ling](https://github.com/lingminhao), and [Jonathan Goeke](https://github.com/jonathangoeke) at the Genome Institute of Singapore. If you wish to contribute, please leave an issue. Thank you.
