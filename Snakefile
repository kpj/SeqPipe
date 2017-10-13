"""
Sequencing workflow
"""

import os


###
# setup
configfile: 'config.yaml'

p = lambda p: os.path.relpath(p, start=workflow.basedir)
input_dir = p(config['directories']['reads'])
ref_dir = p(config['directories']['references'])
output_dir = p(config['directories']['results'])

###
# file gathering
sample_list = []
for entry in os.scandir(input_dir):
    tmp = entry.name.split('.')[0]
    sample_list.append(tmp)

reference_list = []
for entry in os.scandir(ref_dir):
    tmp = entry.name.split('.')[0]
    reference_list.append(tmp)

qc_files = expand(
    os.path.join(
        output_dir, 'quality_control', '{sample}'), sample=sample_list)
map_files = expand(
    os.path.join(
        output_dir, 'read_mapping', '{reference}', '{sample}.sam'),
    sample=sample_list, reference=reference_list)
all_result_files = qc_files + map_files

###
# rule definitions
ruleorder: copy > gunzip

rule all:
    input:
        all_result_files

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
    benchmark:
        os.path.join(
            output_dir, 'benchmark', 'quality_control',
            '{sample}.txt')
    shell:
        """
        mkdir -p {output}
        fastqc --outdir={output} {input.read_file}
        """

rule read_mapping:
    input:
        read_file = os.path.join(output_dir, 'input', '{sample}.fastq'),
        reference_file = srcdir(os.path.join(ref_dir, '{reference}.fa'))
    output:
        os.path.join(output_dir, 'read_mapping', '{reference}', '{sample}.sam')
    params:
        cwd = os.path.join(
            output_dir, 'read_mapping', '{reference}', '{sample}_tmp')
    benchmark:
        os.path.join(
            output_dir, 'benchmark', 'read_mapping',
            '{sample}_{reference}.txt')
    shell:
        """
        mkdir -p {params.cwd}
        mkdir -p $(dirname {workflow.basedir}/{output})

        cd {params.cwd}
        cp {input.reference_file} reference.fa

        bwa index reference.fa
        bwa mem \
            reference.fa {workflow.basedir}/{input.read_file} \
            > {workflow.basedir}/{output}
        """
