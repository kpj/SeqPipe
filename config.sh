#########
# General variable definitions

## pipeline options
ts=$(date +"%Y-%m-%d_%H-%M-%S")
default_output_dir="${ts}_mapping_result"

read_min_len=0
read_max_len=100

## per analysis options
reads_file="reads.fastq"
genome_file="genome.fa"
log_file="log.txt"

data_dir="data/"
result_dir="results/"
fastqc_dir="fastqc/"

## internal helpers
SUCCESS="\033[32m"
WARNING="\033[33m"
ERROR="\033[31m"
BACKGROUND="\033[90m"
INFO="\033[36m"
RESET="\033[0m"
