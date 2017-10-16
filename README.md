# SeqPipe (snakemake edition)

A sequencing pipeline based on Snakemake.

## Installation

```
$ pip install snakemake
```

## Usage

Execute as follows:

```bash
$ snakemake -p
```

Create an overview of the pipeline:

```bash
$ snakemake --dag | dot -Tpng > overview.png
```
