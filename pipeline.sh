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
    fastqc \
        --outdir="$fastqc_dir" \
        "$tmp_reads_file" \
    |& tee -a "$log_file"
    echo -ne "${RESET}"
else
    echo "Use existing FastQC result"
fi

# generate genome if needed
if [[ ! -d "$(dirname "$genome")" ]]; then
    echo "Generate genome"
    mkdir genome

    echo -ne "${BACKGROUND}"
    bowtie2-build "./input/$genome_file" "$genome" > /dev/null
    echo -ne "${RESET}"
    echo "[bowtie2-build stdout truncated]" >> "$log_file"
else
    echo "Use existing genome"
fi

# map reads
if [[ ! -f "$sam" ]]; then
    echo "Generate mapping"

    echo -ne "${BACKGROUND}"
    bowtie2 \
        --very-sensitive \
        -x "$genome" \
        -U "$tmp_reads_file" \
        -S "$sam" \
    |& tee -a "$log_file"
    echo -ne "${RESET}"

    echo "Convert sam to bam"
    samtools view -b "$sam" > "$bam"
else
    echo "Use existing mapping"
fi

## post-process bam
# sort bam into different directory
echo "Sort bam"
samtools sort -o "$data_dir/sorted.bam" "$bam"

# store only aligned reads
echo "Generating additional data"
samtools view -F 4 -b -o "aligned_reads.bam" "$data_dir/sorted.bam"
samtools view -f 4 -b -o "not_aligned_reads.bam" "$data_dir/sorted.bam"

# (re)create results using scripts
rm -rf "$result_dir"
mkdir -p "$result_dir"

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
