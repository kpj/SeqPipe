########
# Script to execute pipeline on specified set of data
#

set -e
set -u

## initial setup

genome_fasta="genome.fa"
inp_reads="reads.fastq"

log_file="log.txt"
pyscript="coverage.py"
genome="genome/genome"
data_dir="data/"
image_dir="images/"
fastqc_dir="fastqc/"

fname_pref="out"
read_file="${fname_pref}.fastq"
bam="${fname_pref}.bam"
sam="${fname_pref}.sam"

read_len_thres="25"

if [[ $# -ne 1 ]] ; then
    echo "Usage: $0 <path to data-dir>"
    exit
fi

cd "$1"

> "$log_file"

# check if needed files exist
if [[ ! -f "$genome_fasta" ]] ; then
    echo "No genome found ($genome_fasta)"
    exit 1
fi
if [[ ! -f "$inp_reads" ]] ; then
    echo "No reads found ($inp_reads)"
    exit 1
fi

if [[ ! -f "$pyscript" ]] ; then
    echo "No python script found ($pyscript)"
    exit 1
fi

## begin of pipeline

# extract reads of certain length
echo "Filter fastq"
cutadapt -m 0 -M $read_len_thres "$inp_reads" -o "$read_file" | tee -a "$log_file"

# quality assessment
if [[ ! -d "$fastqc_dir" ]] ; then
    fastqc --outdir="$fastqc_dir" "$read_file" | tee -a "$log_file"
fi

# generate genome if needed
if [[ ! -d "$(dirname "$genome")" ]] ; then
    echo "Generate genome"
    mkdir genome

    bowtie2-build "$genome_fasta" "$genome" | tee -a "$log_file"
else
    echo "Use existing genome"
fi

# map reads
if [[ ! -f "$sam" ]]; then
    echo "Generate mapping"
    bowtie2 --very-sensitive -x "$genome" -U "$read_file" -S "$sam" | tee -a "$log_file"

    echo "Convert sam to bam"
    samtools view -b "$sam" > "$bam"
else
    echo "Use existing mapping"
fi

# post-process bam
if [[ ! -d "$data_dir" ]] ; then
    mkdir "$data_dir"

    # sort bam into different directory
    echo "Sort bam"
    samtools sort -o "$data_dir/sorted.bam" "$bam"
fi

# create coverage plots
if [[ ! -d "$image_dir" ]] ; then
    mkdir "$image_dir"

    echo "Creating final plots"
    python "$pyscript" "$data_dir/sorted" | tee -a "$log_file"
fi
