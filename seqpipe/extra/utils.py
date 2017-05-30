"""
Utilities for various cases
"""

import os
import sys
import json
import subprocess

from typing import Sized, Dict

import pandas as pd

import click
import pysam
from tqdm import tqdm
from joblib import Parallel, delayed, cpu_count


def count_fastq_sequences(fname: str) -> int:
    """ Count entries in given FASTQ file
    """
    #records = SeqIO.parse(fname, 'fastq')
    #return len(list(records))

    # hacky (but faster) way to count sequences in fastq
    out = subprocess.Popen(
        ['wc', '-l', fname],
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT
    ).communicate()[0]
    return int(int(out.partition(b' ')[0]) / 4)

def count_bam_reads(fname: str) -> int:
    """ Count number of reads in given BAM file
    """
    cproc = subprocess.run(
        ['samtools', 'view', '-F', '0x904', '-c', fname],
        check=True, stdout=subprocess.PIPE)
    res = cproc.stdout.decode('utf-8').rstrip()
    return int(res)

def count_bam_reads_per_sub(bam: str) -> Dict[str, int]:
    """ Count number of reads in given BAM file for each sub-reference
    """
    pysam.index(bam)
    bamf = pysam.AlignmentFile(bam, 'rb')
    assert bamf.has_index()

    counts = {}
    for ref in bamf.references:
        counts[ref] = bamf.count(reference=ref)
    return counts

def parse_single_mapping_result(fname_base: str, split: bool) -> pd.DataFrame:
    """ Compute read-count statistics
    """
    dir_pref = os.path.dirname(fname_base) or '.'
    fname_result = dir_pref + '/statistics_' + os.path.basename(fname_base)
    fname_cache = os.path.join(
        fname_result, f'stats_{"split" if split else "nosplit"}.csv')

    if os.path.exists(fname_cache):
        print('Cached', fname_cache)
        return pd.read_csv(fname_cache, index_col=0)

    if not os.path.exists(fname_result):
        os.makedirs(fname_result)

    # compute statistics
    df = pd.DataFrame()

    print('> Getting counts for each input read-file')
    for read_entry in os.scandir(os.path.join(fname_base, 'runs')):
        with open(os.path.join(read_entry.path, 'meta.json')) as fd:
            cur_meta = json.load(fd)

        if not split:
            mapped_count = count_bam_reads(
                os.path.join(read_entry.path, 'aligned_reads.bam'))
            total_count = count_fastq_sequences(
                os.path.join(
                    read_entry.path, 'data', cur_meta['trimmed_read_path']))

            df = df.append({
                'read_name': cur_meta['read_base'],
                'reference': cur_meta['genome_base'],
                'sub_reference': None,
                'mapped_count': mapped_count,
                'total_count': total_count
            }, ignore_index=True)
        else:
            total_count = count_fastq_sequences(
                os.path.join(
                    read_entry.path, 'data', cur_meta['trimmed_read_path']))

            bam = os.path.join(read_entry.path, 'aligned_reads.bam')
            counts = count_bam_reads_per_sub(bam)

            for sref, count in counts.items():
                df = df.append({
                    'read_name': cur_meta['read_base'],
                    'reference': cur_meta['genome_base'],
                    'sub_reference': sref,
                    'mapped_count': count,
                    'total_count': total_count
                }, ignore_index=True)

    df.to_csv(fname_cache)
    return df

def compute_statistics(fnames: Sized, split: bool = False) -> pd.DataFrame:
    """ Aggregate statistics over all given mappings
    """
    core_num = int(cpu_count() * 4/5)
    result = Parallel(n_jobs=core_num)(
        delayed(parse_single_mapping_result)(fn, split) for fn in tqdm(fnames))
    return pd.concat(result).reset_index(drop=True)

@click.command()
@click.option('--split/--no-split', default=False)
@click.argument('files', nargs=-1, type=click.Path(exists=True))
def main(split: bool, files: Sized) -> None:
    if len(files) == 0:
        print('No files provided')
        return

    df = compute_statistics(files, split=split)
    import ipdb; ipdb.set_trace()

if __name__ == '__main__':
    main()
