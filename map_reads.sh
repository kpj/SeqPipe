########
# Script to execute pipeline on specified set of data
#

set -e
set -u
. config.sh

## initial setup
genome="genome/genome"

fname_pref="out"
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
if [[ ! -f "$genome_file" ]]; then
    echo "No genome found ($genome_file)"
    exit 1
fi
if [[ ! -f "$reads_file" ]]; then
    echo "No reads found ($reads_file)"
    exit 1
fi

## begin of pipeline

# extract reads of specified length
echo "Filter fastq ($read_min_len < |seq| < $read_max_len)"
cutadapt \
    -m $read_min_len \
    -M $read_max_len \
    "$reads_file" \
    -o "$tmp_reads_file" \
| tee -a "$log_file"

# quality assessment
if [[ ! -d "$fastqc_dir" ]]; then
    fastqc \
        --outdir="$fastqc_dir" \
        "$tmp_reads_file" \
    | tee -a "$log_file"
fi

# generate genome if needed
if [[ ! -d "$(dirname "$genome")" ]]; then
    echo "Generate genome"
    mkdir genome

    bowtie2-build "$genome_file" "$genome" | tee -a "$log_file"
else
    echo "Use existing genome"
fi

# map reads
if [[ ! -f "$sam" ]]; then
    echo "Generate mapping"
    bowtie2 \
        --very-sensitive \
        -x "$genome" \
        -U "$tmp_reads_file" \
        -S "$sam" \
    | tee -a "$log_file"

    echo "Convert sam to bam"
    samtools view -b "$sam" > "$bam"
else
    echo "Use existing mapping"
fi

# post-process bam
if [[ ! -d "$data_dir" ]]; then
    mkdir "$data_dir"

    # sort bam into different directory
    echo "Sort bam"
    samtools sort -o "$data_dir/sorted.bam" "$bam"
fi

# (re)create results using scripts
rm -rf "$result_dir"
mkdir -p "$result_dir"

for script in "./scripts/"*; do
    ext="${script##*.}"
    if [ "$ext" == "py" ]; then
        echo "[Executing python script \"$script\"]"
        python "$script" "$data_dir/sorted" | tee -a "$log_file"
    elif [ "$ext" == "sh" ]; then
        echo "[Executing shell script \"$script\"]"
        "$script" "$data_dir/sorted" | tee -a "$log_file"
    else
        echo "[Unknown extension \"$ext\"]"
    fi
done
