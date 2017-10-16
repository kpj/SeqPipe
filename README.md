# SeqPipe (snakemake edition)

A sequencing pipeline based on Snakemake.

## Installation

Install [Snakemake](https://snakemake.readthedocs.io/en/stable/) in order to run the pipeline:

```bash
$ pip install snakemake
```

Furthermore, the following software is used:

* FastQC
* BWA

## Usage

Set up `config.yaml` to your liking and execute the pipeline as follows:

```bash
$ snakemake -p
```

Create an overview of the pipeline:

```bash
$ snakemake --dag | dot -Tpng > overview.png
```
