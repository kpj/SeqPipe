"""
Utilities for various cases
"""

import os
import sys
import json
import subprocess

import pandas as pd

from tqdm import tqdm
from joblib import Parallel, delayed, cpu_count


def count_fastq_sequences(fname):
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

def count_bam_reads(fname):
    """ Count number of reads in given BAM file
    """
    cproc = subprocess.run(
        ['samtools', 'view', '-F', '0x904', '-c', fname],
        check=True, stdout=subprocess.PIPE)
    res = cproc.stdout.decode('utf-8').rstrip()
    return int(res)

def parse_single_mapping_result(fname_base):
    """ Compute read-count statistics
    """
    fname_result = fname_base.rstrip('/') + '_statistics'
    fname_cache = os.path.join(fname_result, 'stats.csv')

    if os.path.exists(fname_cache):
        print('Cached', fname_cache)
        return pd.read_csv(fname_cache, index_col=0)

    if not os.path.exists(fname_result):
        os.makedirs(fname_result)

    # compute statistics
    df = pd.DataFrame()

    with open(os.path.join(fname_base, 'info.json')) as fd:
        meta_info = json.load(fd)

    print('> Getting counts for each input read-file')
    for read_entry in os.scandir(os.path.join(fname_base, 'runs')):
        mapped_count = count_bam_reads(
            os.path.join(read_entry.path, 'aligned_reads.bam'))
        total_count = count_fastq_sequences(
            os.path.join(read_entry.path, 'data', 'tmp.fastq'))

        df = df.append({
            'read_name': read_entry.name,
            'reference': meta_info['reference'],
            'mapped_count': mapped_count,
            'total_count': total_count
        }, ignore_index=True)

    df.to_csv(fname_cache)
    return df

def compute_statistics(fnames):
    """ Aggregate statistics over all given mappings
    """
    core_num = int(cpu_count() * 4/5)
    result = Parallel(n_jobs=core_num)(
        delayed(parse_single_mapping_result)(fn) for fn in tqdm(fnames))
    return pd.concat(result)

def main():
    df = compute_statistics(sys.argv[1:])
    import ipdb; ipdb.set_trace()

if __name__ == '__main__':
    main()
