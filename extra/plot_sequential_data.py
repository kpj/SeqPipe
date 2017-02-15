"""
Visualize fraction of reads mapped to intermediate results
"""

import os
import sys
import shutil
import subprocess

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from Bio import SeqIO
from tqdm import tqdm


def analyze_entry(path, fname='aligned_reads.fastq'):
    """ Count entries in given FASTQ file
    """
    #records = SeqIO.parse(os.path.join(path, fname), 'fastq')
    #return len(list(records))

    # hacky (but faster) way to count sequences in fastq
    if os.path.getsize(os.path.join(path, fname)) == 0:
        return 0

    cproc = subprocess.run(
        ['grep', '-c', '^@', os.path.join(path, fname)],
        check=True, stdout=subprocess.PIPE)
    res = cproc.stdout.decode('utf-8').rstrip()
    return int(res)

def compute_statistics(fname_base):
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
    total_counts = {}
    read_file_list = []
    df = pd.DataFrame()

    print('> Getting absolute counts for each input read-file')
    read_dir = '{}/initial/input/'.format(fname_base)
    for read_file in tqdm(os.listdir(read_dir)):
        read_file_list.append(read_file)
        total_counts[read_file] = analyze_entry(
            os.path.join(
                '{}/initial/mapping/runs/'.format(fname_base),
                read_file, 'data'
            ), fname='tmp.fastq')

    for stage in tqdm(os.listdir(fname_base)):
        out_dir = os.path.join(fname_base, stage, 'output')
        if os.path.exists(out_dir):
            with open(os.path.join(fname_base, stage, 'info.txt')) as fd:
                reference = fd.read()

            for fname_read in tqdm(read_file_list):
                cur_dir = os.path.join(out_dir, fname_read + '_extra')
                count = analyze_entry(cur_dir)

                df = df.append({
                    'read_name': fname_read,
                    'reference': reference.strip(),
                    'count': count
                }, ignore_index=True)

    # count no-match entries
    nm_data = []
    for read_name, group in df.groupby('read_name'):
        cur_count = group['count'].sum()
        tot_count = total_counts[read_name]

        nm_data.append({
            'read_name': read_name,
            'reference': 'no match',
            'count': tot_count-cur_count
        })
    df = df.append(nm_data, ignore_index=True)

    df.to_csv(fname_cache)
    return df

def plot_piecharts(df, fname_base):
    """ Visualize read-distribution over sources and references
    """
    fname_result = fname_base.rstrip('/') + '_statistics'

    print('\n> Creating plots')
    for read_name, group in df.groupby('read_name'):
        names, values = [], []
        for i, row in group.iterrows():
            names.append(row['reference'])
            values.append(row['count'])

        plt.figure()
        plt.pie(
            values, labels=names, explode=(.03,)*len(values),
            autopct='%1.1f%%', shadow=True)
        plt.axis('equal')
        plt.title(read_name)

        #plt.tight_layout()
        plt.savefig(os.path.join(fname_result, read_name.replace(' ', '_')) + '.pdf')
        plt.close()

def main(fname):
    df = compute_statistics(fname)
    plot_piecharts(df, fname)


if __name__ == '__main__':
    if len(sys.argv) != 2:
        print('Usage: {} <seq. data dir>'.format(sys.argv[0]))
        exit(-1)

    sns.set_style('white')
    plt.style.use('seaborn-poster')

    main(sys.argv[1])
