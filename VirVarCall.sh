#!/bin/bash
set -eo pipefail

#############################################
# Configuration Section with Explanations   #
#############################################

# Reference genome file (must be in FASTA format)
REF_GENOME="NC_003977.fasta"

# Number of CPU threads to use for parallel processing
THREADS=80

# Minimum depth threshold for coverage calculations and consensus generation
MIN_DEPTH=1

# FreeBayes parameters
MIN_ALT_COUNT=1
MIN_ALT_FRACTION=0.1
MIN_MAPPING_QUAL=20
MIN_BASE_QUAL=20

# Trimmomatic adapters file path
ADAPTERS="/home/user/miniconda3/share/trimmomatic/adapters/TruSeq3-PE.fa"

# Output and temporary directories
OUTDIR="Results"
TRIMDIR="Trimmed"  # Directory for storing trimmed files
HOSTFILTERDIR="Host_Filtered"  # Directory for host-filtered files
KEEP_INTERMEDIATES=false
SKIP_TRIMMING=true  # Set to true to skip trimming if files exist

# Host filtering parameters
HOST_FILTER=true  # Set to false to disable host filtering
HOST_REFERENCE="/home/user/SR/db/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"

# Variant filtering criteria
FILTER_CRITERIA="QUAL>20 & INFO/DP>2 & INFO/AO>2 & INFO/SAF+INFO/SAR>2"

# Memory settings for host filtering
BBDUK_MEMORY="100g"  # Increased memory for bbduk
BOWTIE2_MEMORY="50g"

#############################################
# Pipeline Execution                        #
#############################################

# Initialize error log
ERROR_LOG="pipeline_errors.log"
echo "Pipeline error log - $(date)" > $ERROR_LOG

#####################################################
# Basic usage with increased memory
# ./pipeline.sh --memory 120g

# Skip trimming and keep intermediates
# ./pipeline.sh --skip-trim --keep-intermediates --memory 100g

# Disable host filtering
# ./pipeline.sh --no-host-filter

# Use custom host reference
# ./pipeline.sh --host-ref /path/to/host.fasta --memory 120g --threads 60

# Show help
# ./pipeline.sh --help
##############################################################################
# Function to show usage
usage() {
  cat << EOF
Usage: $0 [options]
Options:
  --skip-trim          Skip the trimming step (use existing files in $TRIMDIR/)
  --keep-intermediates Keep intermediate files
  --no-host-filter     Disable host filtering
  --host-ref FILE      Specify custom host reference file
  --threads N          Number of threads to use (default: $THREADS)
  --memory MEM         Memory for host filtering (default: $BBDUK_MEMORY)
  -h, --help           Show this help message

Examples:
  $0 --skip-trim --keep-intermediates
  $0 --no-host-filter
  $0 --host-ref /path/to/host_genome.fasta --memory 120g
  $0 --threads 40 --memory 80g

Output files per sample:
  - Consensus sequence: Results/{sample}_consensus.fasta
  - SNPs: Results/{sample}_SNPs.tsv
  - Coverage: Results/{sample}_coverage.tsv
  - Frequency report: Results/{sample}_frequency_report.tsv
  - Filtered variants: Results/{sample}_filtered_variants.tsv
  - Host filtering report: Results/{sample}_host_filtering_report.txt
EOF
  exit 0
}

# Parse command line arguments
while [[ $# -gt 0 ]]; do
  case $1 in
    --skip-trim)
      SKIP_TRIMMING=true
      shift
      ;;
    --keep-intermediates)
      KEEP_INTERMEDIATES=true
      shift
      ;;
    --no-host-filter)
      HOST_FILTER=false
      shift
      ;;
    --host-ref)
      HOST_REFERENCE="$2"
      shift 2
      ;;
    --threads)
      THREADS="$2"
      shift 2
      ;;
    --memory)
      BBDUK_MEMORY="$2"
      shift 2
      ;;
    -h|--help)
      usage
      ;;
    *)
      echo "Unknown option: $1"
      usage
      ;;
  esac
done

# Create directories
mkdir -p "$OUTDIR"
mkdir -p "$TRIMDIR"
mkdir -p "$HOSTFILTERDIR"

# Check for required tools
check_tool() {
  if ! command -v "$1" &> /dev/null; then
    echo "Error: $1 not found in PATH" | tee -a $ERROR_LOG
    return 1
  fi
}

echo "Checking required tools..."
for tool in trimmomatic bwa samtools bcftools freebayes; do
  if ! check_tool "$tool"; then
    exit 1
  fi
done

# Check for host filtering tools
if [ "$HOST_FILTER" = true ]; then
  if command -v bbduk.sh &> /dev/null; then
    HOST_FILTER_METHOD="bbduk"
    echo "Using bbduk for host filtering"
  elif command -v bowtie2 &> /dev/null; then
    HOST_FILTER_METHOD="bowtie2"
    echo "Using bowtie2 for host filtering"
  else
    echo "WARNING: Neither bbduk nor bowtie2 found. Disabling host filtering." | tee -a $ERROR_LOG
    HOST_FILTER=false
  fi
fi

# Index reference if needed
if [ ! -f "${REF_GENOME}.bwt" ]; then
  echo "Indexing reference genome..."
  bwa index "$REF_GENOME" 2>>$ERROR_LOG || {
    echo "ERROR: bwa index failed" | tee -a $ERROR_LOG
    exit 1
  }
fi

# Index host reference if host filtering is enabled
if [ "$HOST_FILTER" = true ] && [ ! -f "${HOST_REFERENCE}.bwt" ]; then
  echo "Indexing host reference genome..."
  bwa index "$HOST_REFERENCE" 2>>$ERROR_LOG || {
    echo "ERROR: Host reference genome indexing failed" | tee -a $ERROR_LOG
    exit 1
  }
fi

# Function to count reads
count_reads() {
  local file=$1
  if [[ $file == *.gz ]]; then
    zcat "$file" 2>/dev/null | awk 'END {print NR/4}'
  else
    awk 'END {print NR/4}' "$file" 2>/dev/null
  fi
}

# Function for host filtering with better memory management
run_host_filtering() {
  local sample_name=$1
  local skip_trim=$2
  
  echo "Step 2/8: Host filtering for $sample_name..."
  
  # Determine input files based on whether trimming was skipped
  if [ "$skip_trim" = true ]; then
    # Use original input files
    R1="${sample_name}_R1.fastq.gz"
    R2="${sample_name}_R2.fastq.gz"
  else
    # Use trimmed files
    R1="${TRIMDIR}/${sample_name}_R1_trimmed.fastq.gz"
    R2="${TRIMDIR}/${sample_name}_R2_trimmed.fastq.gz"
  fi
  
  # Check if input files exist
  if [ ! -f "$R1" ] || [ ! -f "$R2" ]; then
    echo "Error: Input files not found for ${sample_name}" | tee -a $ERROR_LOG
    echo "Looking for: $R1 and $R2" | tee -a $ERROR_LOG
    return 1
  fi
  
  # Count input reads
  echo "Counting input reads..."
  input_reads_R1=$(count_reads "$R1")
  input_reads_R2=$(count_reads "$R2")
  total_input_reads=$((input_reads_R1 + input_reads_R2))
  echo "Total input reads: $total_input_reads"
  
  # Use selected method for host filtering
  case $HOST_FILTER_METHOD in
    "bbduk")
      run_bbduk_filtering "$sample_name" "$R1" "$R2" "$total_input_reads"
      ;;
    "bowtie2")
      run_bowtie2_filtering "$sample_name" "$R1" "$R2" "$total_input_reads"
      ;;
    *)
      echo "ERROR: No host filtering method available" | tee -a $ERROR_LOG
      return 1
      ;;
  esac
}

# Function for bbduk host filtering with memory optimization
run_bbduk_filtering() {
  local sample_name=$1
  local R1=$2
  local R2=$3
  local total_input_reads=$4
  
  echo "Using bbduk for host filtering with ${BBDUK_MEMORY} memory..."
  
  # Try different k-mer sizes if default fails
  for kmer in 23 27 31; do
    echo "Trying k-mer size: $kmer"
    
    bbduk.sh -Xmx${BBDUK_MEMORY} in1="$R1" in2="$R2" \
      out1="${HOSTFILTERDIR}/${sample_name}_R1_filtered.fastq" \
      out2="${HOSTFILTERDIR}/${sample_name}_R2_filtered.fastq" \
      ref="$HOST_REFERENCE" \
      k=$kmer \
      threads=$THREADS \
      stats="${OUTDIR}/${sample_name}_bbduk_stats.txt" \
      2>>"${OUTDIR}/${sample_name}_bbduk.log" && {
      echo "bbduk successful with k=$kmer"
      break
    } || {
      echo "bbduk failed with k=$kmer, trying next size..." | tee -a $ERROR_LOG
      cat "${OUTDIR}/${sample_name}_bbduk.log" >> $ERROR_LOG
      
      # If this was the last kmer, try with reduced memory and threads
      if [ $kmer -eq 31 ]; then
        echo "All kmer sizes failed, trying with reduced resources..."
        bbduk.sh -Xmx80g in1="$R1" in2="$R2" \
          out1="${HOSTFILTERDIR}/${sample_name}_R1_filtered.fastq" \
          out2="${HOSTFILTERDIR}/${sample_name}_R2_filtered.fastq" \
          ref="$HOST_REFERENCE" \
          k=23 \
          threads=40 \
          stats="${OUTDIR}/${sample_name}_bbduk_stats.txt" \
          2>>"${OUTDIR}/${sample_name}_bbduk.log" || {
          echo "ERROR: bbduk host filtering failed for $sample_name with all parameters" | tee -a $ERROR_LOG
          return 1
        }
      fi
    }
  done
  
  # Count filtered reads
  filtered_reads_R1=$(count_reads "${HOSTFILTERDIR}/${sample_name}_R1_filtered.fastq")
  filtered_reads_R2=$(count_reads "${HOSTFILTERDIR}/${sample_name}_R2_filtered.fastq")
  total_filtered_reads=$((filtered_reads_R1 + filtered_reads_R2))
  mapped_reads=$((total_input_reads - total_filtered_reads))
  
  write_host_filtering_report "$sample_name" "$total_input_reads" "$mapped_reads" "$total_filtered_reads"
}

# Function for bowtie2 host filtering
run_bowtie2_filtering() {
  local sample_name=$1
  local R1=$2
  local R2=$3
  local total_input_reads=$4
  
  echo "Using bowtie2 for host filtering..."
  
  # Build bowtie2 index if needed
  if [ ! -f "${HOST_REFERENCE}.1.bt2" ]; then
    echo "Building bowtie2 index for host reference..."
    bowtie2-build --threads $THREADS "$HOST_REFERENCE" "$HOST_REFERENCE" 2>>$ERROR_LOG || {
      echo "ERROR: bowtie2-build failed" | tee -a $ERROR_LOG
      return 1
    }
  fi
  
  # Align and extract unmapped reads
  bowtie2 -x "$HOST_REFERENCE" -1 "$R1" -2 "$R2" \
    --threads $THREADS \
    --very-sensitive \
    --un-conc-gz "${HOSTFILTERDIR}/${sample_name}_filtered.fastq.gz" \
    -S "${HOSTFILTERDIR}/${sample_name}_host_aligned.sam" 2>>"${OUTDIR}/${sample_name}_bowtie2.log" || {
    echo "ERROR: bowtie2 host filtering failed for $sample_name" | tee -a $ERROR_LOG
    cat "${OUTDIR}/${sample_name}_bowtie2.log" >> $ERROR_LOG
    return 1
  }
  
  # Rename files to consistent naming
  mv "${HOSTFILTERDIR}/${sample_name}_filtered.1.fastq.gz" "${HOSTFILTERDIR}/${sample_name}_R1_filtered.fastq.gz" 2>/dev/null || true
  mv "${HOSTFILTERDIR}/${sample_name}_filtered.2.fastq.gz" "${HOSTFILTERDIR}/${sample_name}_R2_filtered.fastq.gz" 2>/dev/null || true
  
  # Decompress for consistency
  if [ -f "${HOSTFILTERDIR}/${sample_name}_R1_filtered.fastq.gz" ]; then
    gunzip -c "${HOSTFILTERDIR}/${sample_name}_R1_filtered.fastq.gz" > "${HOSTFILTERDIR}/${sample_name}_R1_filtered.fastq"
    gunzip -c "${HOSTFILTERDIR}/${sample_name}_R2_filtered.fastq.gz" > "${HOSTFILTERDIR}/${sample_name}_R2_filtered.fastq"
  fi
  
  # Count reads
  filtered_reads_R1=$(count_reads "${HOSTFILTERDIR}/${sample_name}_R1_filtered.fastq")
  filtered_reads_R2=$(count_reads "${HOSTFILTERDIR}/${sample_name}_R2_filtered.fastq")
  total_filtered_reads=$((filtered_reads_R1 + filtered_reads_R2))
  mapped_reads=$(samtools view -c -@ $THREADS "${HOSTFILTERDIR}/${sample_name}_host_aligned.sam" 2>>$ERROR_LOG || echo "0")
  
  write_host_filtering_report "$sample_name" "$total_input_reads" "$mapped_reads" "$total_filtered_reads"
}

# Function to write host filtering report
write_host_filtering_report() {
  local sample_name=$1
  local total_input_reads=$2
  local mapped_reads=$3
  local total_filtered_reads=$4
  
  # Calculate percentages
  mapped_percent=$(awk -v m="$mapped_reads" -v t="$total_input_reads" 'BEGIN {if(t>0) printf "%.2f", (m/t)*100; else print "0.00"}')
  unmapped_percent=$(awk -v u="$total_filtered_reads" -v t="$total_input_reads" 'BEGIN {if(t>0) printf "%.2f", (u/t)*100; else print "0.00"}')
  
  # Check if we have any unmapped reads left
  if [ "$total_filtered_reads" -eq 0 ]; then
    echo "WARNING: No unmapped reads found after host filtering!" | tee -a $ERROR_LOG
    echo "Using original input files instead of host-filtered files" | tee -a $ERROR_LOG
    
    # Use original files instead
    if [[ "$R1" == *.gz ]]; then
      gunzip -c "$R1" > "${HOSTFILTERDIR}/${sample_name}_R1_filtered.fastq"
      gunzip -c "$R2" > "${HOSTFILTERDIR}/${sample_name}_R2_filtered.fastq"
    else
      cp "$R1" "${HOSTFILTERDIR}/${sample_name}_R1_filtered.fastq"
      cp "$R2" "${HOSTFILTERDIR}/${sample_name}_R2_filtered.fastq"
    fi
    
    # Update counts to reflect that we're using all reads
    total_filtered_reads=$total_input_reads
    unmapped_percent=100.00
    mapped_reads=0
    mapped_percent=0.00
  fi
  
  # Write status report
  {
    echo "Host Filtering Status for ${sample_name}"
    echo "====================================="
    echo "Input Reads: $total_input_reads"
    echo "Host Reads: $mapped_reads ($mapped_percent%)"
    echo "Unmapped Reads: $total_filtered_reads ($unmapped_percent%)"
    echo "Filtering Method: $HOST_FILTER_METHOD"
    if [ "$total_filtered_reads" -eq 0 ]; then
      echo "WARNING: No reads passed host filtering, using original files"
    fi
    echo "Filtered Files:"
    echo "- R1: ${HOSTFILTERDIR}/${sample_name}_R1_filtered.fastq"
    echo "- R2: ${HOSTFILTERDIR}/${sample_name}_R2_filtered.fastq"
    echo "====================================="
    echo "Generated on: $(date)"
  } > "${OUTDIR}/${sample_name}_host_filtering_report.txt"
  
  # Cleanup intermediate files
  if [ "$KEEP_INTERMEDIATES" = false ]; then
    rm -f "${HOSTFILTERDIR}/${sample_name}_host_aligned.sam" \
          "${HOSTFILTERDIR}/${sample_name}_R1_filtered.fastq.gz" \
          "${HOSTFILTERDIR}/${sample_name}_R2_filtered.fastq.gz" 2>/dev/null || true
  fi
  
  echo "Host filtering completed for $sample_name"
  echo "Status report: ${OUTDIR}/${sample_name}_host_filtering_report.txt"
}

# Process each sample
for R1_FILE in *_R1.fastq.gz; do
  # Skip if no files found
  [ -e "$R1_FILE" ] || { echo "No FASTQ files found!"; exit 1; }
  
  SAMPLE_ID=${R1_FILE%_R1.fastq.gz}
  R2_FILE=${SAMPLE_ID}_R2.fastq.gz
  
  echo "========================================"
  echo "Processing sample: $SAMPLE_ID"
  echo "========================================"

  # Check if R2 file exists
  if [ ! -f "$R2_FILE" ]; then
    echo "ERROR: R2 file not found: $R2_FILE" | tee -a $ERROR_LOG
    continue
  fi

  # Initialize cleanup trap
  cleanup() {
    if [ "$KEEP_INTERMEDIATES" = false ]; then
      echo "Cleaning up intermediate files for $SAMPLE_ID..."
      rm -f \
        "${SAMPLE_ID}.sam" "${SAMPLE_ID}.bam" \
        "${SAMPLE_ID}.namesorted.bam" "${SAMPLE_ID}.fixmate.bam" \
        "${SAMPLE_ID}.sorted.bam" "${SAMPLE_ID}.sorted.bam.bai" \
        "${SAMPLE_ID}.dedup.bam" "${SAMPLE_ID}.dedup.bam.bai" \
        "${SAMPLE_ID}.mpileup.txt" "${SAMPLE_ID}_depth.txt" \
        "${SAMPLE_ID}_lowcov.bed" "${SAMPLE_ID}_temp_cons.fa" \
        "${SAMPLE_ID}.raw_variants.vcf.gz" "${SAMPLE_ID}.raw_variants.vcf.gz.csi" \
        "${SAMPLE_ID}.filtered.vcf.gz" "${SAMPLE_ID}.filtered.vcf.gz.csi" \
        2>/dev/null || true
    else
      echo "Keeping intermediate files for $SAMPLE_ID as requested"
    fi
  }
  trap cleanup EXIT

  # 1. Trim reads (with skip option)
  echo "Step 1/8: Trimming reads..."
  if [ "$SKIP_TRIMMING" = true ] && \
     [ -f "${TRIMDIR}/${SAMPLE_ID}_R1_trimmed.fastq.gz" ] && \
     [ -f "${TRIMDIR}/${SAMPLE_ID}_R2_trimmed.fastq.gz" ]; then
    echo "Using existing trimmed files from $TRIMDIR"
  else
    echo "Running Trimmomatic..."
    trimmomatic PE -threads $THREADS \
      "$R1_FILE" "$R2_FILE" \
      "${TRIMDIR}/${SAMPLE_ID}_R1_trimmed.fastq.gz" "${TRIMDIR}/${SAMPLE_ID}_R1_unpaired.fastq.gz" \
      "${TRIMDIR}/${SAMPLE_ID}_R2_trimmed.fastq.gz" "${TRIMDIR}/${SAMPLE_ID}_R2_unpaired.fastq.gz" \
      ILLUMINACLIP:${ADAPTERS}:2:30:10:3:TRUE \
      LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:50 2> "${OUTDIR}/${SAMPLE_ID}_trimmomatic.log" || {
      echo "ERROR: Trimmomatic failed for $SAMPLE_ID" | tee -a $ERROR_LOG
      continue
    }
    # Remove unpaired files
    rm -f "${TRIMDIR}/${SAMPLE_ID}_R1_unpaired.fastq.gz" "${TRIMDIR}/${SAMPLE_ID}_R2_unpaired.fastq.gz" 2>/dev/null || true
  fi

  # 2. Host filtering (if enabled)
  if [ "$HOST_FILTER" = true ]; then
    if ! run_host_filtering "$SAMPLE_ID" "$SKIP_TRIMMING"; then
      echo "ERROR: Host filtering failed for $SAMPLE_ID" | tee -a $ERROR_LOG
      # Continue without host filtering
      echo "Continuing without host filtering..."
      R1_INPUT="${TRIMDIR}/${SAMPLE_ID}_R1_trimmed.fastq.gz"
      R2_INPUT="${TRIMDIR}/${SAMPLE_ID}_R2_trimmed.fastq.gz"
    else
      # Use host-filtered files for subsequent steps
      R1_INPUT="${HOSTFILTERDIR}/${SAMPLE_ID}_R1_filtered.fastq"
      R2_INPUT="${HOSTFILTERDIR}/${SAMPLE_ID}_R2_filtered.fastq"
    fi
  else
    # Use trimmed files directly
    R1_INPUT="${TRIMDIR}/${SAMPLE_ID}_R1_trimmed.fastq.gz"
    R2_INPUT="${TRIMDIR}/${SAMPLE_ID}_R2_trimmed.fastq.gz"
  fi

  # Check if input files exist and have content
  if [ ! -s "$R1_INPUT" ] || [ ! -s "$R2_INPUT" ]; then
    echo "ERROR: Input files are empty or missing for $SAMPLE_ID" | tee -a $ERROR_LOG
    continue
  fi

  # 3. Alignment to target reference
  echo "Step 3/8: Aligning reads to target reference..."
  if [[ "$R1_INPUT" == *.gz ]]; then
    # For gzipped files
    bwa mem -t $THREADS -Y "$REF_GENOME" \
      <(zcat "$R1_INPUT") <(zcat "$R2_INPUT") \
      2>>"${OUTDIR}/${SAMPLE_ID}_alignment.log" \
      | samtools view -@ $THREADS -Sb - \
      2>>"${OUTDIR}/${SAMPLE_ID}_alignment.log" \
      | samtools sort -@ $THREADS -o "${SAMPLE_ID}.sorted.bam" \
      2>>"${OUTDIR}/${SAMPLE_ID}_alignment.log" || {
      echo "ERROR: Alignment failed for $SAMPLE_ID" | tee -a $ERROR_LOG
      cat "${OUTDIR}/${SAMPLE_ID}_alignment.log" >> $ERROR_LOG
      continue
    }
  else
    # For uncompressed files
    bwa mem -t $THREADS -Y "$REF_GENOME" \
      "$R1_INPUT" "$R2_INPUT" \
      2>>"${OUTDIR}/${SAMPLE_ID}_alignment.log" \
      | samtools view -@ $THREADS -Sb - \
      2>>"${OUTDIR}/${SAMPLE_ID}_alignment.log" \
      | samtools sort -@ $THREADS -o "${SAMPLE_ID}.sorted.bam" \
      2>>"${OUTDIR}/${SAMPLE_ID}_alignment.log" || {
      echo "ERROR: Alignment failed for $SAMPLE_ID" | tee -a $ERROR_LOG
      cat "${OUTDIR}/${SAMPLE_ID}_alignment.log" >> $ERROR_LOG
      continue
    }
  fi
  
  # Check if alignment produced any reads
  if [ ! -s "${SAMPLE_ID}.sorted.bam" ]; then
    echo "ERROR: No reads aligned to target reference for $SAMPLE_ID" | tee -a $ERROR_LOG
    continue
  fi
  
  samtools index -@ $THREADS "${SAMPLE_ID}.sorted.bam" 2>>$ERROR_LOG || {
    echo "ERROR: samtools index failed for $SAMPLE_ID" | tee -a $ERROR_LOG
    continue
  }

  # 4. Coverage statistics
  echo "Step 4/8: Calculating coverage..."
  samtools coverage "${SAMPLE_ID}.sorted.bam" > "${OUTDIR}/${SAMPLE_ID}_coverage.tsv" 2>>$ERROR_LOG || {
    echo "ERROR: samtools coverage failed for $SAMPLE_ID" | tee -a $ERROR_LOG
    continue
  }

  # 5. Generate base frequency report
  echo "Step 5/8: Generating base frequency report..."
  samtools mpileup -aa -A -Q 0 -d 0 -f "$REF_GENOME" "${SAMPLE_ID}.sorted.bam" \
    > "${SAMPLE_ID}.mpileup.txt" 2>>$ERROR_LOG || {
    echo "ERROR: samtools mpileup failed for $SAMPLE_ID" | tee -a $ERROR_LOG
    continue
  }
  
  awk -F'\t' -v OFS='\t' '
    BEGIN {print "Position", "Ref", "A", "C", "G", "T", "Insertions", "Deletions", "TotalDepth"}
    {
      pos = $2; ref = $3; cov = $4; bases = toupper($5)
      
      # Count standard bases
      a = gsub(/A/, "", bases)
      c = gsub(/C/, "", bases)
      g = gsub(/G/, "", bases)
      t = gsub(/T/, "", bases)
      
      # Initialize INDEL counts
      ins = "0"
      del = "0"
      
      # Simple INDEL detection
      if (match(bases, /\+[0-9]+[ACGTN]+/)) {
        ins_len = substr(bases, RSTART+1, RLENGTH-1)
        ins_seq = substr(bases, RSTART+1+length(ins_len), ins_len)
        ins = ins_len ":" ins_seq
      }
      if (match(bases, /-[0-9]+[ACGTN]+/)) {
        del_len = substr(bases, RSTART+1, RLENGTH-1)
        del_seq = substr(bases, RSTART+1+length(del_len), del_len)
        del = del_len ":" del_seq
      }
      
      print pos, ref, a, c, g, t, ins, del, cov
    }' "${SAMPLE_ID}.mpileup.txt" > "${OUTDIR}/${SAMPLE_ID}_frequency_report.tsv"

  # 6. Variant calling with FreeBayes
  echo "Step 6/8: Calling variants with FreeBayes..."
  freebayes \
    -f "$REF_GENOME" \
    --ploidy 1 \
    --min-alternate-count $MIN_ALT_COUNT \
    --min-alternate-fraction $MIN_ALT_FRACTION \
    --min-mapping-quality $MIN_MAPPING_QUAL \
    --min-base-quality $MIN_BASE_QUAL \
    --report-monomorphic \
    "${SAMPLE_ID}.sorted.bam" \
    2>"${OUTDIR}/${SAMPLE_ID}_freebayes.log" \
    | bcftools sort -Oz -o "${SAMPLE_ID}.raw_variants.vcf.gz" 2>>$ERROR_LOG || {
    echo "ERROR: FreeBayes failed for $SAMPLE_ID" | tee -a $ERROR_LOG
    cat "${OUTDIR}/${SAMPLE_ID}_freebayes.log" >> $ERROR_LOG
    continue
  }

  # Check if variant calling produced any results
  if [ ! -s "${SAMPLE_ID}.raw_variants.vcf.gz" ]; then
    echo "WARNING: No variants called for $SAMPLE_ID, creating empty files" | tee -a $ERROR_LOG
    # Create empty output files
    touch "${OUTDIR}/${SAMPLE_ID}_filtered_variants.tsv"
    touch "${OUTDIR}/${SAMPLE_ID}_SNPs.tsv"
    # Create consensus from reference
    sed "s/>/>${SAMPLE_ID}_consensus /" "$REF_GENOME" > "${OUTDIR}/${SAMPLE_ID}_consensus.fasta"
    continue
  fi

  bcftools index "${SAMPLE_ID}.raw_variants.vcf.gz" 2>>$ERROR_LOG || {
    echo "ERROR: bcftools index failed for raw variants" | tee -a $ERROR_LOG
    continue
  }

  # 7. Filter variants
  echo "Step 7/8: Filtering variants..."
  echo "Using filtering criteria: $FILTER_CRITERIA"
  
  bcftools view -i "$FILTER_CRITERIA" \
    "${SAMPLE_ID}.raw_variants.vcf.gz" \
    2>>"${OUTDIR}/${SAMPLE_ID}_filtering.log" \
    | bcftools norm -f "$REF_GENOME" -m - \
    2>>"${OUTDIR}/${SAMPLE_ID}_filtering.log" \
    | bcftools sort \
    2>>"${OUTDIR}/${SAMPLE_ID}_filtering.log" \
    | bcftools view -Oz -o "${SAMPLE_ID}.filtered.vcf.gz" \
    2>>"${OUTDIR}/${SAMPLE_ID}_filtering.log" || {
    echo "ERROR: Variant filtering failed for $SAMPLE_ID" | tee -a $ERROR_LOG
    cat "${OUTDIR}/${SAMPLE_ID}_filtering.log" >> $ERROR_LOG
    continue
  }
  
  bcftools index "${SAMPLE_ID}.filtered.vcf.gz" 2>>$ERROR_LOG || {
    echo "ERROR: bcftools index failed for filtered variants" | tee -a $ERROR_LOG
    continue
  }

  # Generate filtered variants report
  bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%DP\t%AF\n' "${SAMPLE_ID}.filtered.vcf.gz" \
    > "${OUTDIR}/${SAMPLE_ID}_filtered_variants.tsv" 2>>$ERROR_LOG || {
    echo "ERROR: bcftools query failed for $SAMPLE_ID" | tee -a $ERROR_LOG
    continue
  }

  # 8. Generate consensus with low-coverage masking
  echo "Step 8/8: Generating consensus sequence..."
  
  # Generate depth file
  samtools depth -aa "${SAMPLE_ID}.sorted.bam" > "${SAMPLE_ID}_depth.txt" 2>>$ERROR_LOG || {
    echo "ERROR: samtools depth failed for $SAMPLE_ID" | tee -a $ERROR_LOG
    continue
  }
  
  # Create low-coverage BED file
  awk -v min=$MIN_DEPTH '$3 < min {print $1 "\t" $2-1 "\t" $2}' "${SAMPLE_ID}_depth.txt" > "${SAMPLE_ID}_lowcov.bed" 2>>$ERROR_LOG
  
  bcftools consensus -f "$REF_GENOME" \
    -m "${SAMPLE_ID}_lowcov.bed" \
    -H 1 \
    "${SAMPLE_ID}.filtered.vcf.gz" \
    2>>"${OUTDIR}/${SAMPLE_ID}_consensus.log" \
    | sed "s/>/>${SAMPLE_ID}_consensus /" \
    > "${OUTDIR}/${SAMPLE_ID}_consensus.fasta" || {
    echo "Error generating consensus sequence" | tee -a $ERROR_LOG
    cat "${OUTDIR}/${SAMPLE_ID}_consensus.log" >> $ERROR_LOG
    continue
  }

  # Generate detailed SNP report
  if [ -s "${OUTDIR}/${SAMPLE_ID}_filtered_variants.tsv" ]; then
    awk -F'\t' -v OFS='\t' '
      BEGIN {print "Position", "Ref", "Alt", "Depth", "Frequency", "VariantType"}
      NR > 1 {
        type = (length($3) == length($4)) ? "SNP" : (length($3) < length($4)) ? "INS" : "DEL"
        print $2, $3, $4, $5, $6, type
      }' "${OUTDIR}/${SAMPLE_ID}_filtered_variants.tsv" > "${OUTDIR}/${SAMPLE_ID}_SNPs.tsv"
  else
    echo "No variants to report" > "${OUTDIR}/${SAMPLE_ID}_SNPs.tsv"
  fi

  # Validate output
  echo "Verifying output files..."
  for outfile in \
    "${OUTDIR}/${SAMPLE_ID}_consensus.fasta" \
    "${OUTDIR}/${SAMPLE_ID}_SNPs.tsv" \
    "${OUTDIR}/${SAMPLE_ID}_coverage.tsv" \
    "${OUTDIR}/${SAMPLE_ID}_frequency_report.tsv" \
    "${OUTDIR}/${SAMPLE_ID}_filtered_variants.tsv"; do
    if [ ! -s "$outfile" ]; then
      echo "WARNING: Output file $outfile is empty!" | tee -a $ERROR_LOG
    fi
  done

  echo "Successfully processed $SAMPLE_ID"
  echo "Output files:"
  ls -lh "${OUTDIR}/${SAMPLE_ID}"_* 2>/dev/null | awk '{print $9}'
  echo ""

done

echo "Pipeline completed!"
echo "Error summary:"
if [ -s $ERROR_LOG ]; then
  cat $ERROR_LOG
else
  echo "No errors reported"
fi

exit 0
