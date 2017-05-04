#!/usr/bin/env bash
########
# Script to execute pipeline on specified set of data
#

set -e
set -u

## initial setup
genome="genome/genome"

fname_pref="$data_dir/tmp"
tmp_reads_file="${fname_pref}.fastq"
bam="${fname_pref}.bam"
sam="${fname_pref}.sam"

if [[ $# -ne 1 ]]; then
    echo "Usage: $0 <path to data-dir>"
    exit
fi

cd "$1"

> "$log_file" # clear log file

# check if needed files exist
if [[ ! -f "./input/$genome_file" ]]; then
    echo "No genome found (./input/$genome_file)"
    exit 1
fi
if [[ ! -f "./input/$reads_file" ]]; then
    echo "No reads found (./input/$reads_file)"
    exit 1
fi

if [[ ! -d "$data_dir" ]]; then
    mkdir -p "$data_dir"
fi

## begin of pipeline

# extract reads of specified length
echo "Filter fastq ($read_min_len <= |seq| <= $read_max_len)"
echo -ne "${BACKGROUND}"
cutadapt \
    -m $read_min_len \
    -M $read_max_len \
    "./input/$reads_file" \
    -o "$tmp_reads_file" \
|& tee -a "$log_file"
echo -ne "${RESET}"

# quality assessment
if [[ ! -d "$fastqc_dir" ]]; then
    echo "FastQC quality analysis"
    mkdir -p "$fastqc_dir"

    echo -ne "${BACKGROUND}"
    FASTQC_START_TIME=$SECONDS
    fastqc \
        --outdir="$fastqc_dir" \
        "$tmp_reads_file" \
    |& tee -a "$log_file"
    FASTQC_ELAPSED_TIME=$(($SECONDS - $FASTQC_START_TIME))
    echo -ne "${RESET}"
else
    echo "Use existing FastQC result"
    FASTQC_ELAPSED_TIME=0
fi

# generate genome if needed
if [[ ! -d "$(dirname "$genome")" ]]; then
    echo "Generate genome"
    mkdir genome

    echo -ne "${BACKGROUND}"
    BOWTIEBUILD_START_TIME=$SECONDS
    bowtie2-build "./input/$genome_file" "$genome" > /dev/null
    BOWTIEBUILD_ELAPSED_TIME=$(($SECONDS - $BOWTIEBUILD_START_TIME))
    echo -ne "${RESET}"
    echo "[bowtie2-build stdout truncated]" >> "$log_file"
else
    echo "Use existing genome"
    BOWTIEBUILD_ELAPSED_TIME=0
fi

# map reads
if [[ ! -f "$sam" ]]; then
    echo "Generate mapping"

    echo -ne "${BACKGROUND}"
    BOWTIE_START_TIME=$SECONDS
    bowtie2 \
        --very-sensitive \
        -x "$genome" \
        -U "$tmp_reads_file" \
        -S "$sam" \
    |& tee -a "$log_file"
    BOWTIE_ELAPSED_TIME=$(($SECONDS - $BOWTIE_START_TIME))
    echo -ne "${RESET}"

    echo "Convert sam to bam"
    samtools view -b "$sam" > "$bam"
else
    echo "Use existing mapping"
    BOWTIE_ELAPSED_TIME=0
fi

## post-process bam
BAM_START_TIME=$SECONDS
# sort bam into different directory
echo "Sort bam"
samtools sort -o "$data_dir/sorted.bam" "$bam"

# store only aligned reads
echo "Generating additional data"
samtools view -F 4 -b -o "aligned_reads.bam" "$data_dir/sorted.bam"
samtools view -f 4 -b -o "not_aligned_reads.bam" "$data_dir/sorted.bam"
BAM_ELAPSED_TIME=$(($SECONDS - $BAM_START_TIME))

# (re)create results using scripts
rm -rf "$result_dir"
mkdir -p "$result_dir"

SCRIPT_START_TIME=$SECONDS
for script in "./scripts/"*; do
    ext="${script##*.}"
    if [ "$ext" == "py" ]; then
        echo -e "${SUCCESS}[Executing python script \"$script\"]${RESET}"

        echo -ne "${BACKGROUND}"
        python "$script" "$data_dir/sorted" |& tee -a "$log_file"
        echo -ne "${RESET}"
    elif [ "$ext" == "sh" ]; then
        echo -e "${SUCCESS}[Executing shell script \"$script\"]${RESET}"

        echo -ne "${BACKGROUND}"
        "$script" "$data_dir/sorted" |& tee -a "$log_file"
        echo -ne "${RESET}"
    else
        echo -e "${WARNING}[Unknown extension \"$ext\" for \"$script\"]${RESET}"
    fi
done
SCRIPT_ELAPSED_TIME=$(($SECONDS - $SCRIPT_START_TIME))

# print timing information
print_seconds() {
    printf '%02dh:%02dm:%02ds\n' $(($1/3600)) $(($1%3600/60)) $(($1%60))
}

echo "Runtime statistics" |& tee -a "$log_file"
echo " > fastQC: $(print_seconds $FASTQC_ELAPSED_TIME)" |& tee -a "$log_file"
echo " > bowtie-build: $(print_seconds $BOWTIEBUILD_ELAPSED_TIME)" |& tee -a "$log_file"
echo " > bowtie: $(print_seconds $BOWTIE_ELAPSED_TIME)" |& tee -a "$log_file"
echo " > BAM conversions: $(print_seconds $BAM_ELAPSED_TIME)" |& tee -a "$log_file"
echo " > scripts: $(print_seconds $SCRIPT_ELAPSED_TIME)" |& tee -a "$log_file"
