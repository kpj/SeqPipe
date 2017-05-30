"""
Handle execution of pipeline
"""

import io
import os
import time
import shlex
import datetime
from contextlib import contextmanager

from typing import Union, Dict, Generator

import sh


@contextmanager
def cwd(path: str) -> Generator:
    """ Temporarily change working directory
    """
    oldpwd = os.getcwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(oldpwd)

def seconds2string(seconds: int) -> str:
    """ Convert seconds to nicely readable string
    """
    return str(datetime.timedelta(seconds=seconds))

def execute(out_stream: Union[io.StringIO, str], cmd: sh.Command) -> int:
    """ Execute command and return various pieces of meta-information
    """
    start_time = time.time()
    cmd(_out=out_stream, _err_to_out=True)
    duration = round(time.time() - start_time)

    return duration

def pipeline(
    read_path: str, genome_path: str,
    output_dir: str, param_obj: Dict,
    out_stream: io.StringIO
) -> Dict:
    """ Execute whole sequencing pipeline
    """
    # trim reads
    if param_obj['min_read_len'] == -1 and param_obj['max_read_len'] == -1:
        print('Skipping read trimming', file=out_stream)
        trimmed_read_path = read_path
        cutadapt_duration = 0
    else:
        trimmed_read_path = f'{os.path.splitext(read_path)[0]}_trimmed.fastq'
        cutadapt = sh.cutadapt.bake(
            read_path, o=trimmed_read_path,
            m=param_obj['min_read_len'], M=param_obj['max_read_len'])
        cutadapt_duration = execute(out_stream, cutadapt)

    # quality control
    fastqc_dir = os.path.join(output_dir, 'fastqc')
    os.makedirs(fastqc_dir)

    fastqc = sh.fastqc.bake(trimmed_read_path, outdir=fastqc_dir)
    fastqc_duration = execute(out_stream, fastqc)

    # build reference genome
    genome_base = os.path.basename(genome_path)
    genome_prefix = os.path.join(output_dir, 'genome', genome_base)
    os.makedirs(os.path.dirname(genome_prefix))

    bowtiebuild = sh.bowtie_build.bake(genome_path, genome_prefix)
    bowtiebuild_duration = execute(os.devnull, bowtiebuild)

    # map reads
    mapfile_base = os.path.join(output_dir, 'data', 'mapping')
    os.makedirs(os.path.dirname(mapfile_base))

    sam_path = os.path.join(output_dir, f'{mapfile_base}.sam')

    bowtie = sh.bowtie.bake('-S', genome_prefix, trimmed_read_path, sam_path)
    if param_obj['bowtie_args'] is not None:
        bowtie = bowtie.bake(*shlex.split(param_obj['bowtie_args']))
    bowtie_duration = execute(out_stream, bowtie)

    # bam post-processing
    bam_duration = 0
    bam_path = f'{mapfile_base}.bam'
    bam_sorted_path = f'{mapfile_base}_sorted.bam'
    bam_aligned_path = os.path.join(output_dir, 'aligned_reads.bam')
    bam_notaligned_path = os.path.join(output_dir, 'not_aligned_reads.bam')

    st_view = sh.samtools.bake('view', '-b', sam_path, o=bam_path)
    bam_duration += execute(out_stream, st_view)

    st_sort = sh.samtools.bake('sort', bam_path, o=bam_sorted_path)
    bam_duration += execute(out_stream, st_sort)

    st_alig = sh.samtools.bake(
        'view', '-b', bam_sorted_path, F=4, o=bam_aligned_path)
    bam_duration += execute(out_stream, st_alig)

    st_noalig = sh.samtools.bake(
        'view', '-b', bam_sorted_path, f=4, o=bam_notaligned_path)
    bam_duration += execute(out_stream, st_noalig)

    # execute scripts
    script_duration = 0

    if param_obj['exec_scripts']:
        script_dir = os.path.join(output_dir, 'scripts')
        for _, __, filenames in os.walk(script_dir):
            for fn in filenames:
                try:
                    cmd = sh.Command(fn, search_paths=[script_dir]).bake(
                        os.path.splitext(bam_sorted_path)[0])

                    print(f'Executing "{fn}"', file=out_stream)
                    with cwd(output_dir):
                        script_duration += execute(out_stream, cmd)
                except sh.CommandNotFound:
                    print(f'Skipping "{fn}"', file=out_stream)
    else:
        print('Skipping script execution', file=out_stream)

    # final timings
    print('Runtime statistics', file=out_stream)
    print(
        f' > cutadapt: {seconds2string(cutadapt_duration)}', file=out_stream)
    print(
        f' > fastQC: {seconds2string(fastqc_duration)}', file=out_stream)
    print(
        f' > bowtie-build: {seconds2string(bowtiebuild_duration)}',
        file=out_stream)
    print(
        f' > bowtie: {seconds2string(bowtie_duration)}', file=out_stream)
    print(
        f' > BAM processing: {seconds2string(bam_duration)}',
        file=out_stream)
    print(
        f' > Script execution: {seconds2string(script_duration)}',
        file=out_stream)

    return {
        'read_name': os.path.basename(read_path),
        'results_path': os.path.join(output_dir, 'results')
    }