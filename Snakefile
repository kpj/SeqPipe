"""
Sequencing workflow
"""

import os
import re
import sys

from typing import Any, Dict


###
# setup
configfile: os.path.join(workflow.basedir, 'config.yaml')

p = lambda p: os.path.relpath(p, start=workflow.basedir)
input_dir = p(config['directories']['reads'])
ref_dir = p(config['directories']['references'])
output_dir = p(config['directories']['results'])

os.chdir(workflow.basedir)

do_adapter_trimming = os.path.exists(config['conditions']['adapter_file'])
premapsuffix = 'noadapters' if do_adapter_trimming else 'raw'

paired_read_file_identifier = 'R2.fastq'

###
# helper functions
def clean_sample_name(sname: str) -> str:
    """ Remove filename suffixes
    """
    suffixes = ['.fastq', '.fastq.gz']
    for suf in suffixes:
        if sname.endswith(suf):
            return sname[:-len(suf)]
    return sname

def generate_read_mapping_input(wildcards: Any) -> Dict[str, str]:
    """ Generate input files needed to handle both single- and paired-end reads
    """
    # order R1 and R2 correctly
    if wildcards.sample.endswith('R1'):
        sample1 = wildcards.sample
        sample2 = wildcards.sample.replace('R1', 'R2')
    elif wildcards.sample.endswith('R2'):
        sample1 = wildcards.sample.replace('R2', 'R1')
        sample2 = wildcards.sample
    else:
        sample1 = wildcards.sample
        sample2 = ''

    # assemble filenames
    read_file1 = os.path.join(
        output_dir, 'input', f'{sample1}.{premapsuffix}.fastq')
    read_file2 = os.path.join(
        output_dir, 'input', f'{sample2}.{premapsuffix}.fastq') \
        if len(sample2) > 0 else read_file1

    # remove superfluous files
    if sample1 not in sample_list:
        read_file1 = read_file2
    elif sample2 not in sample_list:
        read_file2 = read_file1

    return {
        'read_file1': read_file1,
        'read_file2': read_file2
    }

###
# file gathering
sample_list = set()
secondary_sample_map = {}

print('Parsing input', file=sys.stderr)
for entry in os.scandir(input_dir):
    # remember sample
    cur_sample = clean_sample_name(entry.name)
    sample_list.add(cur_sample)

    # check for paired-end reads
    if paired_read_file_identifier in entry.name:
        pair_fname = entry.path.replace(
            paired_read_file_identifier, 'R1.fastq')

        if os.path.exists(pair_fname):
            print(
                f' > Assuming that "{entry.path}" '
                f'is paired with "{pair_fname}".', file=sys.stderr)

            clean_paired_fname = clean_sample_name(os.path.basename(pair_fname))
            assert clean_paired_fname not in secondary_sample_map
            secondary_sample_map[clean_paired_fname] = \
                clean_sample_name(entry.name)

secondary_samples = set(secondary_sample_map.values())
primary_samples = sample_list - secondary_samples

reference_list = []
for entry in os.scandir(ref_dir):
    tmp = entry.name.split('.')[0]
    reference_list.append(tmp)

# print pipeline execution summary
print('Mapping-pipeline overview', file=sys.stderr)
for i, ref in enumerate(sorted(reference_list)):
    ref_sep = '│' if i < len(reference_list)-1 else ' '
    ref_knick = '├' if i < len(reference_list)-1 else '└'

    print(f' {ref_knick}── {ref}', file=sys.stderr)
    for j, sample in enumerate(sorted(primary_samples)):
        if sample in secondary_sample_map:
            paired_sample = secondary_sample_map[sample]
            msg = f'"{sample}" -- "{paired_sample}"'
        else:
            msg = f'"{sample}"'

        sam_knick = '├' if j < len(primary_samples)-1 else '└'
        print(f' {ref_sep}   {sam_knick}── {msg}', file=sys.stderr)

# figure out which result files to generate
qc_files = expand(
    os.path.join(
        output_dir, 'quality_control', '{sample}'), sample=sample_list) \
    if config['conditions']['assess_read_quality'] else []
rdist_overview = [os.path.join(
    output_dir, 'results', 'read_distribution_overview.tsv')]
all_result_files = qc_files + rdist_overview

# only consider primary samples for read_mapping (paired-end reads)
sample_wildcard_constraint = f'({")|(".join(map(re.escape, primary_samples))})'

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
        os.path.join(output_dir, 'input', '{sample}.raw.fastq')
    shell:
        """
        gunzip -c {input.fname} > {output}
        """

rule copy:
    input:
        fname = os.path.join(input_dir, '{sample}.fastq')
    output:
        os.path.join(output_dir, 'input', '{sample}.raw.fastq')
    shell:
        """
        cp {input.fname} {output}
        """

rule quality_control:
    input:
        read_file = os.path.join(
            output_dir, 'input', '{sample}.'+premapsuffix+'.fastq')
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

rule adapter_trimming:
    input:
        fname = os.path.join(output_dir, 'input', '{sample}.raw.fastq')
    output:
        os.path.join(output_dir, 'input', '{sample}.noadapters.fastq')
    shell:
        """
        cutadapt \
            -a $(cat {config[conditions][adapter_file]}) \
            -o {output} \
            {input.fname}
        """

rule convert_reference_suffix:
    input:
        fname = '{prefix}.fa'
    output:
        temp('{prefix}.fasta')
    shell:
        """
        cp {input.fname} {output}
        """

rule reference_indexing:
    input:
        reference_file = srcdir(os.path.join(ref_dir, '{reference}.fasta'))
    output:
        os.path.join(output_dir, 'read_mapping', '{reference}', 'index')
    shell:
        """
        mkdir -p {output}
        cd {output}

        cp {input.reference_file} reference.fa
        bwa index reference.fa
        """

rule read_mapping:
    input:
        unpack(generate_read_mapping_input),
        reference_dir = os.path.join(
            output_dir, 'read_mapping', '{reference}', 'index')
    output:
        os.path.join(
            output_dir, 'read_mapping', '{reference}', '{sample}', 'reads.sam')
    wildcard_constraints:
        sample = sample_wildcard_constraint
    threads:
        config['runtime']['threads']
    benchmark:
        os.path.join(
            output_dir, 'benchmark', 'read_mapping',
            '{sample}_{reference}.txt')
    shell:
        """
        mkdir -p $(dirname {workflow.basedir}/{output})

        ref_file="{workflow.basedir}/{input.reference_dir}/reference.fa"
        read_file1="{workflow.basedir}/{input.read_file1}"
        read_file2="{workflow.basedir}/{input.read_file2}"

        # handle single-end reads
        if [ "$read_file1" == "$read_file2" ]; then
            read_file2=""
        fi

        bwa mem \
            -t {threads} \
            "$ref_file" \
            "$read_file1" $read_file2 \
            > {workflow.basedir}/{output}
        """

rule samtools_filter:
    input:
        '{sample}.sam'
    output:
        '{sample}.filtered.bam'
    shell:
        """
        # filters: read quality, secondary reads, supplementary reads
        samtools view \
            -b \
            -q config[parameters][minimal_mapping_quality] \
            -F 256 \
            -F 2048 \
            {input} > {output}
        """

rule samtools_sort:
    input:
        '{sample}.filtered.bam'
    output:
        '{sample}.filtered.sorted.bam'
    threads: config['runtime']['threads']
    wrapper:
        '0.2.0/bio/samtools/sort'

rule samtools_index:
    input:
        '{sample}.filtered.sorted.bam'
    output:
        '{sample}.filtered.sorted.bam.bai'
    wrapper:
        '0.2.0/bio/samtools/index'

rule read_distribution_overview:
    input:
        bam_files = expand(
            os.path.join(
                output_dir, 'read_mapping',
                '{reference}', '{sample}', 'reads.filtered.sorted.bam'),
            sample=primary_samples, reference=reference_list),
        index_files = expand(
            os.path.join(
                output_dir, 'read_mapping',
                '{reference}', '{sample}', 'reads.filtered.sorted.bam.bai'),
            sample=primary_samples, reference=reference_list),
    output:
        report = os.path.join(output_dir, 'results', 'read_distribution_overview.tsv')
    script:
        'scripts/rdist_overview.py'
