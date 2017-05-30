#!/usr/bin/env bash
#########
# General variable definitions

## pipeline options
ts=$(date +"%Y-%m-%d_%H-%M-%S")
export default_output_dir="mapping_result_${ts}"

export read_min_len=0
export read_max_len=99999999999
export bowtie_params=""

export disable_scripts=false
export show_prefix=false

## per analysis options
export meta_info_file="info.json"

export reads_file="reads.fastq"
export genome_file="genome.fa"
export log_file="log.txt"

export data_dir="data/"
export result_dir="results/"
export fastqc_dir="fastqc/"

## internal helpers
export SUCCESS="\033[32m"
export WARNING="\033[33m"
export ERROR="\033[31m"
export BACKGROUND="\033[90m"
export INFO="\033[36m"
export RESET="\033[0m"
