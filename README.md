# Viral Variant Calling Pipeline from Metagenome-virvarcall
VirVarCall is a comprehensive, robust, and flexible bioinformatics pipeline designed for viral variant calling and consensus sequence generation from paired-end sequencing data. The pipeline processes raw FASTQ files through quality control, host filtering, alignment, variant calling, and generates multiple output formats for downstream analysis.
 Features
*End-to-end Processing: Quality trimming, host read removal, alignment, variant calling, and consensus generation

Flexible Host Filtering: Supports both bbduk and bowtie2 with automatic fallback

Comprehensive Output: Consensus sequences, SNP reports, coverage statistics, and frequency profiles

Robust Error Handling: Comprehensive logging and error recovery mechanisms

Resource Optimization: Configurable memory and thread usage with automatic resource scaling

Quality Control: Multiple QC checkpoints and validation steps

📋 Requirements
Bioinformatics Tools
Core: trimmomatic, bwa, samtools, bcftools, freebayes

Host Filtering: bbduk.sh (from BBMap) or bowtie2

Standard: awk, sed, gzip

Reference Genomes
Viral target reference (FASTA format)

Optional host reference for filtering (human, mouse, etc.)

🛠️ Installation
bash
# Clone repository
git clone https://github.com/yourusername/viral-variant-pipeline
cd viral-variant-pipeline

# Make script executable
chmod +x pipeline.sh

# Install dependencies (example with conda)
conda create -n virvarcall trimmomatic bwa samtools bcftools freebayes bbmap bowtie2
conda activate virvarcall
🎯 Quick Start
Basic Usage
bash
./pipeline.sh
Advanced Usage
bash
# Run with custom resources
./pipeline.sh --threads 32 --memory 80g

# Skip trimming and keep intermediates
./pipeline.sh --skip-trim --keep-intermediates

# Disable host filtering
./pipeline.sh --no-host-filter

# Use custom host reference
./pipeline.sh --host-ref /path/to/host.fasta
📁 Input Requirements
File Naming Convention
Forward reads: {sample}_R1.fastq.gz

Reverse reads: {sample}_R2.fastq.gz

Required Files
Viral reference genome: NC_003977.fasta (configurable)

Host reference genome (optional, for filtering)

📊 Output Files
For each sample {sample}, the pipeline generates:

#File	Description
{sample}_consensus.fasta	Consensus sequence with low-coverage masking
{sample}_SNPs.tsv	SNP report with variant typing
{sample}_coverage.tsv	Coverage statistics across genome
{sample}_frequency_report.tsv	Base-by-base frequency analysis
{sample}_filtered_variants.tsv	Filtered variant calls
{sample}_host_filtering_report.txt	Host filtering statistics
Various log files	Processing logs for each step
⚙️ Configuration
Edit the configuration section in pipeline.sh:

bash
# Reference genome
REF_GENOME="NC_003977.fasta"

# Computational resources
THREADS=80
BBDUK_MEMORY="100g"

# Quality thresholds
MIN_DEPTH=1
MIN_ALT_COUNT=1
MIN_ALT_FRACTION=0.1

# Directories
OUTDIR="Results"
TRIMDIR="Trimmed"
HOSTFILTERDIR="Host_Filtered"
🔧 Pipeline Steps
Quality Trimming - Trimmomatic for adapter removal and quality filtering

Host Filtering - bbduk or bowtie2 to remove host reads

Alignment - BWA MEM alignment to viral reference

Coverage Analysis - samtools coverage statistics

Variant Calling - FreeBayes for variant detection

Variant Filtering - bcftools for quality filtering

Consensus Generation - bcftools consensus with low-coverage masking

Report Generation - Comprehensive output files

🐛 Troubleshooting
Common Issues
Memory errors: Use --memory to adjust allocation

Host filtering failures: Pipeline automatically uses fallback methods

No variants called: Pipeline creates consensus from reference

Empty outputs: Check pipeline_errors.log and input quality

Check Logs
bash
# Main error log
cat pipeline_errors.log

# Sample-specific logs
ls Results/*.log
📝 Command Line Options
bash
./pipeline.sh [options]

Options:
  --skip-trim          Skip trimming step (use existing trimmed files)
  --keep-intermediates Keep intermediate files
  --no-host-filter     Disable host filtering
  --host-ref FILE      Specify custom host reference file
  --threads N          Number of threads to use (default: 80)
  --memory MEM         Memory for host filtering (default: 100g)
  -h, --help           Show help message
🎓 Example Workflow
bash
# 1. Prepare your data
# Place paired-end FASTQ files in working directory
# Sample1_R1.fastq.gz, Sample1_R2.fastq.gz, etc.

# 2. Prepare reference files
cp your_virus.fasta NC_003977.fasta

# 3. Run pipeline with optimal settings
./pipeline.sh --threads 32 --memory 80g --keep-intermediates

# 4. Check results
ls -la Results/
cat Results/Sample1_host_filtering_report.txt
