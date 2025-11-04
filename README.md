# Viral Variant Calling Pipeline from Metagenome (VirVarCall)

**VirVarCall** is a comprehensive, robust, and flexible bioinformatics pipeline designed for viral variant calling and consensus sequence generation from paired-end sequencing data. The pipeline processes raw FASTQ files through quality control, host filtering, alignment, variant calling, and generates multiple output formats for downstream analysis.

## 🚀 Features

- **End-to-end Processing**: Quality trimming, host read removal, alignment, variant calling, and consensus generation
- **Flexible Host Filtering**: Supports both bbduk and bowtie2 with automatic fallback
- **Comprehensive Output**: Consensus sequences, SNP reports, coverage statistics, and frequency profiles
- **Robust Error Handling**: Comprehensive logging and error recovery mechanisms
- **Resource Optimization**: Configurable memory and thread usage with automatic resource scaling
- **Quality Control**: Multiple QC checkpoints and validation steps

## 📋 Requirements

### **Bioinformatics Tools**

```bash
# Core tools
trimmomatic
bwa
samtools
bcftools
freebayes

# Host filtering (one of these)
bbduk.sh  # from BBMap
bowtie2

# Standard Unix tools
awk
sed
gzip
```

### **Reference Genomes**
- **Viral target reference** (FASTA format)
- **Optional host reference** for filtering (human, mouse, etc.)

## 🛠️ Installation

```bash
# Clone repository
git clone https://github.com/yourusername/viral-variant-pipeline
cd viral-variant-pipeline

# Make script executable
chmod +x VirVarCall.sh

# Install dependencies (example with conda)
conda create -n virvarcall trimmomatic bwa samtools bcftools freebayes bbmap bowtie2
conda activate virvarcall
```

## 🎯 Quick Start

### Basic Usage
```bash
./VirVarCall.sh
```

### Advanced Usage
```bash
# Run with custom resources
./VirVarCall.sh --threads 32 --memory 80g

# Skip trimming and keep intermediates
./VirVarCall.sh --skip-trim --keep-intermediates

# Disable host filtering
./VirVarCall.sh --no-host-filter

# Use custom host reference
./VirVarCall.sh --host-ref /path/to/host.fasta
```

## 📁 Input Requirements

### **File Naming Convention**
```bash
# Forward reads
{sample}_R1.fastq.gz

# Reverse reads  
{sample}_R2.fastq.gz
```

### **Required Files**
```bash
# Viral reference genome (configurable)
NC_003977.fasta

# Host reference genome (optional, for filtering)
/path/to/host_reference.fasta
```

## 📊 Output Files

For each sample `{sample}`, the pipeline generates:

| **File** | **Description** |
|----------|----------------|
| **`{sample}_consensus.fasta`** | Consensus sequence with low-coverage masking |
| **`{sample}_SNPs.tsv`** | SNP report with variant typing |
| **`{sample}_coverage.tsv`** | Coverage statistics across genome |
| **`{sample}_frequency_report.tsv`** | Base-by-base frequency analysis |
| **`{sample}_filtered_variants.tsv`** | Filtered variant calls |
| **`{sample}_host_filtering_report.txt`** | Host filtering statistics |
| **Various log files** | Processing logs for each step |

## ⚙️ Configuration

Edit the configuration section in `VirVarCall.sh`:

```bash
#!/bin/bash
set -eo pipefail

#############################################
# Configuration Section                    #
#############################################

# Reference genome file (must be in FASTA format)
REF_GENOME="NC_003977.fasta"

# Number of CPU threads to use for parallel processing
THREADS=80

# Minimum depth threshold for coverage calculations
MIN_DEPTH=1

# FreeBayes parameters
MIN_ALT_COUNT=1
MIN_ALT_FRACTION=0.1
MIN_MAPPING_QUAL=20
MIN_BASE_QUAL=20

# Output and temporary directories
OUTDIR="Results"
TRIMDIR="Trimmed"
HOSTFILTERDIR="Host_Filtered"

# Memory settings for host filtering
BBDUK_MEMORY="100g"
```

## 🔧 Pipeline Steps

```bash
# 1. Quality Trimming - Trimmomatic for adapter removal
trimmomatic PE -threads $THREADS input_R1.fastq.gz input_R2.fastq.gz ...

# 2. Host Filtering - bbduk or bowtie2 to remove host reads
bbduk.sh -Xmx${BBDUK_MEMORY} in1=... in2=... ref=$HOST_REFERENCE

# 3. Alignment - BWA MEM alignment to viral reference  
bwa mem -t $THREADS -Y "$REF_GENOME" R1.fastq R2.fastq

# 4. Coverage Analysis - samtools coverage statistics
samtools coverage sample.sorted.bam

# 5. Variant Calling - FreeBayes for variant detection
freebayes -f "$REF_GENOME" --ploidy 1 sample.sorted.bam

# 6. Variant Filtering - bcftools for quality filtering
bcftools view -i "QUAL>20 & INFO/DP>2" variants.vcf.gz

# 7. Consensus Generation - bcftools consensus
bcftools consensus -f "$REF_GENOME" -m lowcov.bed filtered.vcf.gz

# 8. Report Generation - Comprehensive output files
```

## 🧬 Consensus Generation

The pipeline provides flexible consensus generation options using `bcftools consensus`:

| **Goal** | **Keep Reference** | **Mask Low-Coverage** | **Remove Reference** |
|----------|-------------------|----------------------|---------------------|
| **Command** | `-H 1` | `-H 1 -m lowcov.bed` | `-H A --missing '-'` |
| **Output** | full consensus | consensus with Ns | variants-only |
| **Ref. bases kept?** | ✅ Yes | ✅ except masked | ❌ No |

### **Consensus Modes Explained:**

```bash
# 1. Full Consensus (-H 1)
# Purpose: Generate complete consensus sequence
# Behavior: Uses reference bases where no variants are called
bcftools consensus -f reference.fasta -H 1 filtered.vcf.gz

# 2. Masked Consensus (-H 1 -m lowcov.bed)  
# Purpose: Generate consensus while masking low-coverage regions
# Behavior: Replaces low-coverage positions with 'N'
bcftools consensus -f reference.fasta -H 1 -m lowcov.bed filtered.vcf.gz

# 3. Variants-Only (-H A --missing '-')
# Purpose: Extract only variant positions
# Behavior: Outputs only positions with called variants  
bcftools consensus -f reference.fasta -H A --missing '-' filtered.vcf.gz
```

### **Low-Coverage Masking:**
```bash
# Threshold: Configurable via MIN_DEPTH parameter (default: 1)
MIN_DEPTH=1

# Implementation: Positions with depth < MIN_DEPTH are masked with 'N'
awk -v min=$MIN_DEPTH '$3 < min {print $1 "\t" $2-1 "\t" $2}' sample_depth.txt > sample_lowcov.bed

# File: Generated as {sample}_lowcov.bed during pipeline execution
```

## 🐛 Troubleshooting

### **Common Issues**

```bash
# 1. Memory errors - Use --memory to adjust allocation
./VirVarCall.sh --memory 120g

# 2. Host filtering failures - Check available tools
which bbduk.sh
which bowtie2

# 3. No variants called - Check input BAM file
samtools flagstat sample.sorted.bam

# 4. Empty outputs - Check error logs
cat pipeline_errors.log
```

### **Check Logs**
```bash
# Main error log
cat pipeline_errors.log

# Sample-specific logs
ls Results/*.log

# Trimmomatic logs  
cat Results/*_trimmomatic.log

# Alignment logs
cat Results/*_alignment.log
```

## 📝 Command Line Options

```bash
./VirVarCall.sh [options]

Options:
  --skip-trim          Skip trimming step (use existing trimmed files)
  --keep-intermediates Keep intermediate files
  --no-host-filter     Disable host filtering
  --host-ref FILE      Specify custom host reference file
  --threads N          Number of threads to use (default: 80)
  --memory MEM         Memory for host filtering (default: 100g)
  -h, --help           Show help message

Examples:
  ./VirVarCall.sh --skip-trim --keep-intermediates
  ./VirVarCall.sh --no-host-filter
  ./VirVarCall.sh --host-ref /path/to/host_genome.fasta --memory 120g
  ./VirVarCall.sh --threads 40 --memory 80g
```

## 🎓 Example Workflow

```bash
# 1. Prepare your data
# Place paired-end FASTQ files in working directory
ls -1 *.fastq.gz
# Sample1_R1.fastq.gz
# Sample1_R2.fastq.gz
# Sample2_R1.fastq.gz
# Sample2_R2.fastq.gz

# 2. Prepare reference files
cp your_virus.fasta NC_003977.fasta

# 3. Run pipeline with optimal settings
./VirVarCall.sh --threads 32 --memory 80g --keep-intermediates

# 4. Check results
ls -la Results/

# 5. Examine host filtering report
cat Results/Sample1_host_filtering_report.txt

# 6. Check consensus sequence
head -20 Results/Sample1_consensus.fasta

# 7. View SNP report
head Results/Sample1_SNPs.tsv
```

## 🤝 Contributing

Contributions are welcome! Please feel free to submit pull requests or open issues for bugs and feature requests.

## 📄 License

This project is licensed under the MIT License - see the LICENSE file for details.

## 📞 Support

For issues and questions:

```bash
# 1. Check the main error log
cat pipeline_errors.log

# 2. Review sample-specific logs
ls Results/*.log

# 3. Verify dependencies are installed
which trimmomatic bwa samtools bcftools freebayes

# 4. Check input file formats
file *.fastq.gz
file NC_003977.fasta
```

---

**Note**: This pipeline is designed for viral genomics research and requires appropriate computational resources for optimal performance.
