## 🧪 Test Run Example

### **Quick Test Setup**

For a quick test run with minimal data, you can use the following setup:

```bash
# Sample files needed for test run:
JUST1_R1.fastq.gz    # Forward reads
JUST1_R2.fastq.gz    # Reverse reads  
NC_003977.fasta      # Viral reference genome
```

### **Test Run Command**

```bash
# Quick test with host filtering disabled and minimal resources
./pipeline.sh --no-host-filter --threads 8 --memory 16g --keep-intermediates
```

### **Expected Output for Test Run**

After the test completes, you should see these output files in the `Results/` directory:

```bash
Results/
├── JUST1_consensus.fasta              # Consensus sequence
├── JUST1_SNPs.tsv                     # SNP calls
├── JUST1_coverage.tsv                 # Coverage statistics
├── JUST1_frequency_report.tsv         # Base frequency report
├── JUST1_filtered_variants.tsv        # Filtered variants
├── JUST1_trimmomatic.log              # Trimming log
├── JUST1_alignment.log                # Alignment log
├── JUST1_freebayes.log                # Variant calling log
└── JUST1_filtering.log                # Filtering log
```

### **Quick Verification Commands**

```bash
# Check if pipeline completed successfully
echo "Pipeline exit code: $?"

# Verify main output files were created
ls -la Results/JUST1_*

# Check consensus sequence
head -n 2 Results/JUST1_consensus.fasta

# View coverage summary
cat Results/JUST1_coverage.tsv

# Check SNP calls
cat Results/JUST1_SNPs.tsv
```

### **Test Run Configuration**

The test run will use these modified settings:
- **Threads**: 8 (reduced from default 80)
- **Memory**: 16g (reduced from default 100g)  
- **Host filtering**: Disabled with `--no-host-filter`
- **Intermediate files**: Kept for debugging with `--keep-intermediates`

### **Expected Runtime**
- **Small dataset (≤1GB)**: 10-30 minutes
- **Medium dataset (1-5GB)**: 30-90 minutes
- **Large dataset (>5GB)**: 2+ hours

### **Troubleshooting Test Runs**

If the test run fails, check these common issues:

```bash
# 1. Check file permissions
chmod +x pipeline.sh

# 2. Verify input files exist and are not empty
ls -lh JUST1_R1.fastq.gz JUST1_R2.fastq.gz NC_003977.fasta

# 3. Check if dependencies are installed
which trimmomatic bwa samtools bcftools freebayes

# 4. View detailed error log
cat pipeline_errors.log

# 5. Check individual step logs
cat Results/JUST1_*.log
```

This test configuration is ideal for:
- ✅ Testing pipeline functionality
- ✅ Debugging installation issues  
- ✅ Processing small datasets quickly
- ✅ Learning the pipeline workflow
- ✅ Validating input file formats
