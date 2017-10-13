"""
Sequencing workflow
"""

import os


###
# setup
configfile: 'config.yaml'

input_dir = config['directories']['input']
output_dir = config['directories']['output']

###
# file gathering
sample_list = []

# quality control
for entry in os.scandir(input_dir):
    tmp = entry.name.split('.')[0]
    sample_list.append(tmp)

###
# rule definitions
ruleorder: copy > gunzip

rule all:
    input:
        expand(
            os.path.join(
                output_dir, 'quality_control', '{sample}'), sample=sample_list)

rule gunzip:
    input:
        fname = os.path.join(input_dir, '{sample}.fastq.gz')
    output:
        os.path.join(output_dir, 'input', '{sample}.fastq')
    shell:
        """
        gunzip -c {input.fname} > {output}
        """

rule copy:
    input:
        fname = os.path.join(input_dir, '{sample}.fastq')
    output:
        os.path.join(output_dir, 'input', '{sample}.fastq')
    shell:
        """
        cp {input.fname} {output}
        """

rule quality_control:
    input:
        read_file = os.path.join(output_dir, 'input', '{sample}.fastq')
    output:
        os.path.join(output_dir, 'quality_control', '{sample}')
    shell:
        """
        mkdir -p {output}
        fastqc --outdir={output} {input.read_file}
        """
