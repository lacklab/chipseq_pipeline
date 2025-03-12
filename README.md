# ChIP-Seq Analysis Pipeline

This repository contains a Snakemake pipeline for ChIP-Seq data processing and analysis, covering quality control, alignment, peak calling, and downstream visualization.

---

## 📖 Overview

This pipeline automates the essential steps involved in ChIP-Seq data analysis, including data preparation, alignment, peak calling, and quality control reporting.


### **1. Data Preparation and Quality Control**
- Adapter trimming using `trim_galore` to remove sequencing adapters and low-quality bases.
- Quality control of raw and trimmed reads using `FastQC`.

### **2. Alignment and BAM Processing**
- Read alignment with `BWA`.
- BAM file sorting, indexing, and duplicate marking/removal using `samtools`.
- Quality metrics generation with `samtools idxstats`, `flagstat`, and `stats`.

### **3. Peak Calling and Signal Processing**
- Peak calling with `MACS3` for identifying narrow peak regions.
- Signal normalization and BigWig generation using `bedtools` and `bedGraphToBigWig`.
- Integration with external datasets via ChIP-Atlas (BED and BigWig files).

### **4. Annotation and Quality Control**
- Peak annotation using HOMER’s `annotatePeaks.pl`.
- Fraction of Reads in Peaks (FRiP) calculation for quality metrics.


### **5. Reporting and Visualization**
Centralized quality control report using MultiQC, including:
- **Fingerprint Plots:** Assess library complexity and duplication.
- **PCA (Principal Component Analysis):** Visualize global sample similarity.
- **Correlation Matrix:** Heatmap representation of sample correlations.
- **FastQC Reports:** Evaluates sequence quality, GC content, and potential contaminants for both raw and trimmed reads.
- **Trimgalore Reports:** Provides adapter trimming efficiency metrics and post-trimming quality scores.
- **Samtools Metrics**:
    - **flagstat:** Provides summary alignment statistics for raw, filtered, and deduplicated BAM files.
    - **idxstats:** Reports per-chromosome read distribution.
    - **stats:** Generates detailed alignment metrics across multiple processing stages.
- **FRIP Score Calculation:** Fraction of reads overlapping identified peaks.
- **Number of Peaks:** Total number of peaks identified during peak calling.
- **Median Fragment Length:** Estimates fragment size distribution in the BAM files.

---

## 📦 Requirements

### Tools:
- Python (3.x)
- Snakemake
- BWA
- Samtools
- TrimGalore
- MACS3
- HOMER
- Bedtools
- MultiQC
- ChIP-Atlas dependencies (aria2, wget)
- SRA-tools
- paralle-fastq-dump
- deeptools
- bedGraphToBigWig

### Python Packages:
- `pandas`
- `pysam`

**You don't need to install these individually. As long as you have Python, Snakemake and conda in your path, new environments will be created.**

---

## 📂 Folder Structure

```plaintext
├── assets/
│   ├── annotatepeaks.asset # MultiQC custom asset for ploting
├── config/
│   ├── config.yaml         # Configuration file for the pipeline
│   ├── references.yaml     # A yaml file for the paths for the references used
│   └── samples.tsv         # A table that the samples to be run
├── link/                   # Will be created to softlink the raw data
├── logs/                   # Will be created for every rules to printout the logs
├── profile/
│   ├── config.yaml         # Configuration file for the cluster setup
├── qc/                     # QC reporting
├── workflow/                 
│   ├── envs/               # Conda environment files
│   ├── scripts/            # Scripts used in rules
│   ├── rules/              # Snakemake rules for each step
│   │   ├── annot.smk
│   │   ├── bw.smk
│   │   ├── chipatlas.smk
│   │   ├── common.smk
│   │   ├── deeptools.smk
│   │   ├── map.smk
│   │   ├── peak.smk
│   │   ├── qc.smk
│   │   ├── sra.smk
│   │   └── trim.smk
│   └── Snakefile           # Main entry point for the pipeline
├── raw-data/               # Raw sequencing data files (not included)
├── sra-data/               # Raw sequencing data files from SRA (not included)
├── results_{ref}/          # Output results directory
├── trimmed/                # Trimmed reads directory
└── run_snakemake.sh        # Script to run the pipeline
```

---


##  📜 Input File
The pipeline requires a **sample metadata file** in TSV format. Below is an example:



| Name           | Unit | Control     | Library | Fastq1                                                                 | Fastq2                                                                 | GSM       | Genome |
|-------------- |------|------------|---------|------------------------------------------------------------------------|------------------------------------------------------------------------|----------|--------|
| LNCaP_AR      | r1   | LNCaP_input | Paired  | /groups/lackgrp/raw_data/LNCP/ChIP-seq/AR_4h_R1_1.fq.gz               | /groups/lackgrp/raw_data/LNCP/ChIP-seq/AR_4h_R1_2.fq.gz               | -        | hg38   |
| LNCaP_H3K27ac | r1   | LNCaP_input | Paired  | /groups/lackgrp/raw_data/LNCP/ChIP-seq/H3K27ac_IP_4h_R1_1.fq.gz       | /groups/lackgrp/raw_data/LNCP/ChIP-seq/H3K27ac_IP_4h_R1_2.fq.gz       | -        | hg38   |
| LNCaP_input   | r1   | -          | Paired  | /groups/lackgrp/raw_data/LNCP/ChIP-seq/INPUT_IP_4h_R1_1.fq.gz         | /groups/lackgrp/raw_data/LNCP/ChIP-seq/INPUT_IP_4h_R1_2.fq.gz         | -        | hg38   |
| PDX35_HOXB13  | r1   | -          | Single  | SRR11856307                                                            | -                                                                      | GSM4569938 | hg38  |



- **Name**: Group name.
- **Unit**: Unit or replicate identifier.
- **Control**: \<Name\> of the sample.
- **Library**: Single or Paired ended experiment.
- **Fastq1**: Path to raw FASTQ file (Forward Pair). **PS**: To trigger SRA fetch put the SRR number, if multiple SRR for the experiments, they will be merged under Name column.
- **Fastq2**: Path to raw FASTQ file (Reverse Pair).
- **GSM**: GSM identifier of publised experiment to allow ChIP-Atlas fetch.
- **Genome**: The reference built that you want to run the sample in. See (`config/references.yaml`)

---

## 🎛️ Configurations
The configuration file (`config/config.yaml`) specifies pipeline parameters. Below is an example:

```yaml
SAMPLES: config/samples.tsv

OUTPUT:
    RUN:
        QC: True
        PEAKS: True
        BWS: True
        CHIPATLASBED: False
        CHIPATLASBIGWIG: False
    BW_NORMALIZATIONS:
        - rawcount
        - FPM
    MACS_THRESHOLD:
        - 0.01
        - 0.05
        - 0.1
    BAMPROCESS_PARAMS: -q 30
    CHIPATLASBED_THRESHOLD: '05'

CUT_ADAPTERS: True
UPDATE_CHIPATLAS: False
```

The reference path are specified in (`config/references.yaml`). Below is an example:

```yaml
hg19:
    FA: /groups/lackgrp/genomeAnnotations/hg19/hg19.fa
    BWA_IDX: /groups/lackgrp/genomeAnnotations/hg19/hg19.bwa.idx
    CHROM_SIZES: /groups/lackgrp/genomeAnnotations/hg19/hg19.chrom.sizes
    BLACKLIST: /groups/lackgrp/genomeAnnotations/hg19/hg19-blacklist.v2.bed
hg38:
    FA: /groups/lackgrp/genomeAnnotations/hg38/hg38.fa
    BWA_IDX: /groups/lackgrp/genomeAnnotations/hg38/hg38.bwa.idx
    CHROM_SIZES: /groups/lackgrp/genomeAnnotations/hg38/hg38.chrom.sizes
    BLACKLIST: /groups/lackgrp/genomeAnnotations/hg38/hg38-blacklist.v2.bed
```

The related conda environments are located in `workflow/envs/*.yaml`. No configuration is needed

---



## ▶️ Running the Pipeline

### 1. Configure the Pipeline
Edit the `config/config.yaml` and `config/references.yaml` files to specify:
- Output directory structure.
- Flags to enable or disable specific steps (e.g., QC, peak calling, BigWig normalization).
- Reference genome details.

Ensure that the `config/samples.tsv` file is properly formatted with the sample information.

---

### 2. Slurm Profile
The pipeline uses a Slurm cluster via the `profile/` directory. The `config.yaml` for the Slurm profile should include the following:

```yaml
executor: slurm
default-resources:
  slurm_partition: "express,normal,big-mem,long"
  mem_mb_per_cpu: 8096
  mem_mb: None
  disk_mb: 1024000
  runtime: 600
restart-times: 1
max-jobs-per-second: 10
max-status-checks-per-second: 1
local-cores: 1
latency-wait: 60
jobs: 12
keep-going: True
rerun-incomplete: True
printshellcmds: True
scheduler: greedy
use-conda: True
```
- **Cluster Resources**: Adjust memory (mem_mb), disk (disk_mb), and partition names according to your Slurm setup.

---
### 3. Submission Script

Use the following run_pipeline.sh script to submit the pipeline to the Slurm cluster. The script activates the required conda environment and runs Snakemake with the specified profile.

```bash
#!/bin/bash

snakemake --profile profile/

```
---
### 4. Submit the Pipeline

Run the following command to execute the pipeline:

```bash
sbatch run_pipeline.sh
```
This will:
- Automatically submit jobs to the Slurm cluster.
- Use the configuration specified in the profile/config.yaml file.
- Execute all defined rules in the pipeline.
