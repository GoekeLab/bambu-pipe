# **BETA**
This is a pre-release version of the pipeline intended for testers only. Please use with caution.

# **Context-Aware Transcript Quantification from Long Read Single-Cell and Spatial Transcriptomics data**
This pipeline performs context-aware transcript discovery and quantification from long read single-cell and spatial transcriptomics data. The workflow consists of: 
1. (Optional) Quality score filtering with [chopper](https://github.com/wdecoster/chopper)
2. Barcode/UMI identification and demultiplexing with [flexiplex](https://davidsongroup.github.io/flexiplex/)
3. Primer removal with [cutadapt](https://cutadapt.readthedocs.io/en/stable/)
4. Genome alignment with [minimap2](https://lh3.github.io/minimap2/minimap2.html)
5. Transcript discovery and quantification with [Bambu](https://github.com/GoekeLab/bambu/tree/BambuDev).

The final output includes novel transcripts found in the sample and transcript level count matrices for each barcode/spatial coordinate. 
#

### **Content** 
- [Installation](#installation)
- [General Usage](#General-Usage)
- [Release History](#Release-History)
- [Citation](#Citation)
- [Contributors](#Contributors)
#

### **Installation** 
Install the following dependencies before running the pipeline:
- [Nextflow](https://www.nextflow.io/docs/latest/install.html) 
- [Docker](https://docs.docker.com/engine/install/ubuntu/) (or [Singularity](https://docs.sylabs.io/guides/3.0/user-guide/installation.html_) if you do not have user permissions for Docker). 

The latest version for each dependency is recommended. 
#

### **General Usage** 
To run the pipeline, you must provide a samplesheet, reference genome, and reference annotation file as input. The pipeline performs transcript discovery and quantification on either a single sample or on multiple samples based on the number of samples specified in the samplesheet. Refer to the [Parameters](#parameters) and [Samplesheet (CSV)](#samplesheet-csv) sections below for more details. 

**Running the pipeline**

Use the command below to run the pipeline on the test data provided in `examples/`
``` 
nextflow run $PWD/bambu-singlecell-spatial \
  --input $PWD/examples/samplesheet.csv \
  --genome $PWD/examples/Homo_sapiens.GRCh38.dna_sm.primary_assembly_chr9_1_1000000.fa \
  --annotation $PWD/examples/Homo_sapiens.GRCh38.91_chr9_1_1000000.gtf \
  -profile singularity,hpc
``` 

**Nextflow Profiles**

To configure the executor and container used by the pipeline, pass the following profile types through the `-profile` argument in Nextflow.  

Container profiles:
- `singularity`: use Sigularity images (recommended on HPC systems)
- `docker`: use Docker images

Executor profiles:
- `local`: execute pipeline on a local machine (suitable for small datasets)
- `hpc`: execute pipeline on an HPC system

Note: By default, the executor for the `hpc` profile is set to 'slurm'. To change this, modify the nextflow.config file. 

#

### **Samplesheet (CSV)**

The pipeline requires a `.csv` formatted samplesheet to define the input data. This file is mandatory, regardless of the number of samples being processed. Each row in the samplesheet must represent a single sample and its corresponding file paths and metadata. 

**Required Columns**

The following columns must be present in the CSV:

- `sample`: sample name
- `path`: path to the input file (FASTQ, BAM, or RDS)
- `chemistry`: 10x library chemistry (see [Supported 10x Library Chemistries](#supported-10x-library-chemistries))
- `technology`: sequencing technology (see [Supported Sequencing Technologies](#supported-sequencing-technologies))

**Supported Input Formats**

The pipeline is designed to be flexible. Depending on your starting point in the workflow, the `path` column can point to the following file types:
- **FASTQ**: Raw reads (compressed `.gz` or uncompressed)
- **BAM**: Demultiplexed and aligned reads
- **RDS**: Pre-processed bambu read class objects

For more details on starting the pipeline from specific stages, please refer to the [Advanced Usage]() section. 

**Example Samplesheet (Single Sample)**
| sample | path | chemistry | technology |
|---|---|---|---|
| sample1 | path/to/sample1_fastq.gz | 10x3v2 | ONT

**Example Samplesheet (Multiple Samples)**
| sample | path | chemistry | technology |
|---|---|---|---|
| sample1 | path/to/sample1_fastq.gz | 10x3v2 | ONT
| sample2 | path/to/sample2_fastq.gz | 10x3v3 | PacBio
| sample3 | path/to/sample3_fastq.gz | 10x3v4 | ONT

Note: A samplesheet template is provided at `examples/samplesheet.csv`.
If all samples share the same library chemistry and/or sequencing technology, you may omit the `chemistry` and `technology` columns and use the `--chemistry` and `--technology` flags instead.
#

### **Supported 10x Library Chemistries**
The following library chemistries are supported. Please specify the sample chemistry in the samplesheet as shown:
- `10x3v2` (Single Cell 3' v2)
- `10x3v3` (Single Cell 3' v3 & Next GEM Single Cell 3' v3.1)
- `10x3v4` (GEM-X Single Cell 3' v4)
- `10x5v2` (Single Cell 5' v2)
- `10x5v3` (GEM-X Single Cell 5' v3)
- `visium-v1` (Visium Spatial Gene Expression Slide 6.5 mm; serial prefix V1)
- `visium-v2` (Visium Spatial Gene Expression Slide 6.5 mm; serial prefix V2)
- `visium-v3` (Visium Spatial Gene Expression Slide 6.5mm; serial prefix V3)
- `visium-v4` (Visium CytAssist Spatial Gene Expression Slide 6.5mm; serial prefix V4)
- `visium-v5` (Visium CytAssist Spatial Gene Expression Slide 11mm; serial prefix V5)
#

### **Supported Sequencing Technologies**
The following sequencing technologies are supported. Please specify the technology for each sample in the samplesheet as shown:
- `ONT` (Oxford Nanopore Technologies)
- `PacBio` (Pacific Biosciences)
#

### **Parameters**

**Mandatory**
- `--input` [string]: Path to the samplesheet .csv file 
- `--genome` [string]: Path to the reference genome .fa or .fasta file 
- `--annotation` [string]: Path to the reference annotation .gtf or .gff file 

**Optional**
- `--outdir` [string, default: 'output']: Path to the output directory
- `--chemistry` [string, default: null]: Specify if all samples in the samplesheet share the same library chemistry 
- `--technology` [string, default: null]: Specify if all samples in the samplesheet share the same sequencing technology
- `early_stop_stage` [string, default: null]: Stop the pipeline early and output intermediate files (see Advanced Usage section). Options:
  - "bam": Stops pipeline after minimap2 alignment
  - "rds": Stops pipeline after Bambu read class construction
- `--qscore_filtering` [boolean, default: true]: Enable or disable quality score filtering of reads
- `--ndr` [float, default: null]: NDR threshold for Bambu transcript discovery. If not set, Bambu will recommend a suitable value
- `--process_by_chromosome` [boolean, default: true]: If true, run Bambu steps separately for each chromosome to reduce memory usage
- `--deduplicate_umis` [boolean, default: true]: If true, Bambu will perform UMI deduplication 
- `--quantification_mode` [string, default: "EM_clusters"]: Quantification mode for transcript counts. Available options are:
  - "no_EM": No EM quantification
  - "EM": Perform EM quantification across all cells/spatial coordinates
  - "EM_clusters": Performs pseudo-bulk clustering using [Seurat](https://satijalab.org/seurat/) followed by EM quantification at the cluster level
- `--resolution` [float, default: 0.8]: Seurat clustering resolution
#

### **Output** ###
All outputs from the pipeline are written to the directory specified by the `--outdir` parameter. The pipeline produces per-sample alignment files, per-sample read class files used by Bambu, and the combined transcript discovery and quantification results. The examples below show the output directory structure for both single and multi-sample runs:

**Output Structure (Single Sample)**
```
output/
├── bam/                                
│   ├── sample1_demultiplexed.bam
│   └── sample1_demultiplexed.bam.bai
│
├── read_class/                                
│   └── sample1_readClassFile.rds
│
├── sample1_extended_annotations.gtf
└── sample1_se.rds
```

**Output Structure (Multiple Samples)**
```
output/
├── bam/                                
│   ├── sample1_demultiplexed.bam
│   ├── sample1_demultiplexed.bam.bai
│   ├── sample2_demultiplexed.bam
│   └── sample2_demultiplexed.bam.bai
│
├── read_class/                                
│   ├── sample1_readClassFile.rds
│   └── sample2_readClassFile.rds
│
├── combined_extended_annotations.gtf
└── combined_se.rds
```

Note: For single sample runs, the `extended_annotations.gtf` and `se.rds` are prefixed with the `sample_name`. For multi-sample runs, the `combined` prefix is used instead. 

**Description of the Output Files**
| File | Description 
|---|---
| <sample_name>_demultiplexed.bam | BAM file containing demultiplexed, trimmed and aligned reads
| <sample_name>_demultiplexed.bam.bai | BAM index for the corresponding BAM file
| <sample_name>_readClassFile.rds |  An intermediate metadata file used by Bambu the constructed read classes. This file can be used as input to the pipeline. 
| *_extended_annotations.gtf | A `.gtf` file containing the novel transcripts discovered by Bambu as well as the reference annotations provided by the user.
| *_se.rds | A RangedSummarizedExperiment object containing count matrices (`.mtx`) from transcript quantification by Bambu. Depending on the `quantification_mode`, the matrices are provided at either pseudobulk or single-cell level. The rows of the matrices represent transcript names, while the columns follow the `sampleName_cellBarcode` or `sampleName_clusterId` naming convention.

**Count Matrices**

The `SummarizedExperiment` object contains four distrinct types of count matrices, which can be accessed in R using the `assays()` function. Depending on your analysis requirements you can choose from the following:
- `counts`: expression estimates
- `CPM`: seqencing depth normalised estimates
- `fullLengthCounts`: estimates of read counts mapped as full length reads for each transcript
- `uniqueCounts`: counts of reads that are uniquely mapped to each transcript 

#

### **Spatial Analysis** ##
The pipeline applies the same processing steps to both single-cell and spatial samples. However, for spatial data, the generated `SummarizedExperiment` object is appended with spatial mapping information, which is stored in `colData`.  

**Example - Spatial Mapping Information (Non-Visium HD)**:

For non-Visium HD samples, the SummarizedExperiment's colData contains the spatial barcode and the corresponding X and Y spatial coordinates. 

| Barcode            | X coordinate | Y coordinate| 
|:---|:---|:---|
| AAACAACGAATAGTTC | 17 | 1 |
| AAACAAGTATCTCCCA  | 103 | 51 |
| AAACAATCTACTAGCA | 44 | 4 |
#

### **Visium HD Spatial Analysis (Under Development)**

This feature is still under development and will be released in a future update.
#

### **Fusion Transcript Analysis (Under Development)** ##
This feature is still under development and will be released in a future update.
# 

### **Additional Information**
UMI correction is done at the barcode level. The longest read for each unique barcode-UMI combination is kept for analysis.

### **Release History** 

Beta Release: 2023-May-03

### **Citation**
Bambu
Chen, Y., Sim, A., Wan, Y.K. et al. Context-aware transcript quantification from long-read RNA-seq data with Bambu. Nat Methods (2023). https://doi.org/10.1038/s41592-023-01908-w

Chopper

Cutadapt

Flexiplex
Oliver Cheng, Min Hao Ling, Changqing Wang, Shuyi Wu, Matthew E Ritchie, Jonathan Göke, Noorul Amin, Nadia M Davidson, Flexiplex: a versatile demultiplexer and search tool for omics data, Bioinformatics, Volume 40, Issue 3, March 2024, btae102, https://doi.org/10.1093/bioinformatics/btae102

Minimap2
Li, H. (2021). New strategies to improve minimap2 alignment accuracy. Bioinformatics, 37:4572-4574.

Samtools

Seurat

### **Contributors**
This package is developed and maintained by [Andre Sim](https://github.com/andredsim), [Chin Hao Lee](https://github.com/ch99l), [Min Hao Ling](https://github.com/lingminhao), and [Jonathan Goeke](https://github.com/jonathangoeke) at the Genome Institute of Singapore. If you wish to contribute, please leave an issue. Thank you.
