"""
Visualize fraction of reads mapped to intermediate results
"""

import os
import sys
import shutil
import subprocess

from Bio import SeqIO
from tqdm import tqdm
import matplotlib.pyplot as plt


def analyze_entry(path, fname='aligned_reads.fastq'):
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

def main(fname):
    data = {}
    totals = {}

    # check given read-files
    print('> Getting absolute counts for each input read-file')
    read_dir = '{}/initial/input/'.format(fname)
    for read_file in tqdm(os.listdir(read_dir)):
        data[read_file] = {}
        totals[read_file] = analyze_entry(
            os.path.join(
                '{}/initial/mapping/runs/'.format(fname),
                read_file, 'data'
            ), fname='tmp.fastq')

    # acquire data
    print('> Counting mapped reads per stage')
    for stage in tqdm(os.listdir(fname)):
        out_dir = os.path.join(fname, stage, 'output')
        if os.path.exists(out_dir):
            with open(os.path.join(fname, stage, 'info.txt')) as fd:
                stage_name = fd.read()

            for read in tqdm(data.keys()):
                cur_dir = os.path.join(out_dir, read + '_extra')
                res = analyze_entry(cur_dir)
                data[read][stage_name] = res

    # visualize data
    out_dir = fname.rstrip('/') + '_images'
    if os.path.exists(out_dir):
        shutil.rmtree(out_dir)
    os.makedirs(out_dir)

    print('\n> Creating plots')
    for read_name, dat in tqdm(data.items()):
        names, values = [], []
        for n, v in dat.items():
            names.append(n)
            values.append(v)

        tot = sum(values)
        names.append('no match')
        values.append(totals[read_name]-tot)

        plt.figure()
        plt.pie(
            values, labels=names, explode=(.03,)*len(values),
            autopct='%1.1f%%', shadow=True)
        plt.axis('equal')
        plt.title(read_name)

        plt.tight_layout()
        plt.savefig(os.path.join(out_dir, read_name.replace(' ', '_')) + '.pdf')
        plt.close()

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print('Usage: {} <seq. data dir>'.format(sys.argv[0]))
        exit(-1)

    main(sys.argv[1])
