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
To run this pipeline, you will need the following dependencies:
- [Nextflow](https://www.nextflow.io/docs/latest/install.html) 
- [Docker](https://docs.docker.com/engine/install/ubuntu/) (or [Singularity](https://docs.sylabs.io/guides/3.0/user-guide/installation.html_) if you do not have user permissions for Docker). 

The latest version for each dependency is recommended. 
#

### **General Usage** 
To run the pipeline, you must provide a CSV samplesheet, a reference genome, and a reference annotation file (see [Parameters](#parameters) section). The pipeline performs transcript discovery and quantification on either a single sample or multiple samples, as specified in the samplesheet. See the [Samplesheet (CSV)](#samplesheet-csv) section below for more details. 

**Running the pipeline**

The example below shows how to run the pipeline on a test dataset provided in `examples/`
``` 
nextflow run $PWD/bambu-singlecell-spatial \
  --input $PWD/examples/samplesheet.csv \
  --genome $PWD/examples/Homo_sapiens.GRCh38.dna_sm.primary_assembly_chr9_1_1000000.fa \
  --annotation $PWD/examples/Homo_sapiens.GRCh38.91_chr9_1_1000000.gtf \
  -profile singularity,hpc
``` 

**Profiles**

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

The samplesheet is always required, regardless of the number of samples you wish to perform transcript discovery and quantification on. Each row represents one sample and its associated input file and metadata. The sample can either be a compressed or uncompressed FASTQ, BAM, or RDS file depending on your starting point in the workflow (see Advanced Usage section for more information). 

The following columns are required for each sample:

- `sample`: sample name
- `path`: path to the input file (FASTQ, BAM, or RDS)
- `chemistry`: 10x library chemistry (see [Supported 10x Library Chemistries](#supported-10x-library-chemistries))
- `technology`: sequencing technology (see [Supported Sequencing Technologies](#supported-sequencing-technologies))

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

A samplesheet template is also provided at `examples/samplesheet.csv`

Note: If all samples share the same chemistry and/or technology, you can omit the chemistry and technology columns in the samplesheet and use the `--chemistry` and `--technology` parameters instead.
#

### **Supported 10x Library Chemistries**
The following library chemistries are supported. Please specify the sample chemistry in the samplesheet as shown:
- `10x3v2` (Single Cell 3' v2)
- `10x3v3` (Single Cell 3' v3 & Next GEM Single Cell 3' v3.1)
- `10x3v4` (GEM-X Single Cell 3' v4)
- `10x5v2` (Single Cell 5' v2)
- `10x5v3` (GEM-X Single Cell 5' v3)
- `visium-v1` (Visium Spatial Gene Expression Slide 6.5 mm; serial prefix **V1**)
- `visium-v2` (Visium Spatial Gene Expression Slide 6.5 mm; serial prefix **V2**)
- `visium-v3` (Visium Spatial Gene Expression Slide 6.5mm; serial prefix **V3**)
- `visium-v4` (Visium CytAssist Spatial Gene Expression Slide 6.5mm; serial prefix **V4**)
- `visium-v5` (Visium CytAssist Spatial Gene Expression Slide 11mm; serial prefix **V5**)
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
After the run, all outputs are written to the directory specified by the `--outdir` parameter. The pipeline produces per-sample alignment files, per-sample read class files used by Bambu, and the combined transcript discovery and quantification results. An example output structure for a single sample and multi-sample run are shown below:

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

**Description of the Output Files**
| File | Description 
|---|---
| <sample_name>_demultiplexed.bam |
| <sample_name>_demultiplexed.bam.bai | BAM index for the corresponding BAM file
| <sample_name>_readClassFile.rds |  
| *_extended_annotations.gtf | 
| *_se.rds 

#

### **Spatial Analysis** ##
The pipeline applies the same processing steps to both single-cell and spatial samples. The primary difference is that, for spatial samples, the resulting `SummarisedExperiment` object includes additional spatial mapping information stored in the `colData` field. 

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
