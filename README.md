# RNAseq Analysis Pipeline

## Usage

Execute `run.sh` as follows:
```
$ run.sh <path to read directory> <path to reference genome file>
```
This will create a `mapping_results` directory which will contain all results.

Note: make sure to delete old results before creating new ones.

## Overview

* needed
    * reference genome: genome/reference.fa
    * reads: data/reads.fastq
* process reference genome
    * bowtie-build genome/reference.fa genome/genome
* assess read quality
    * fastqc --outdir=fastqc_data data/reads.fastq
    * http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/
* extract needed features
    * cutadapt -m $min_len -M $max_len "$in_fastq" -o "$out_fastq"
* map reads
    * bowtie -a -v 1 --sam genome/genome data/reads.fastq out.sam
    * samtools view -b out.sam > out.bam
    * [samtools sort -o sorted.bam out.bam]
