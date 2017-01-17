"""
Visualize fraction of reads mapped to intermediate results
"""

import os
import sys
import subprocess

import matplotlib.pyplot as plt
from Bio import SeqIO


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
    read_dir = '{}/initial/input/'.format(fname)
    for read_file in os.listdir(read_dir):
        data[read_file] = {}
        totals[read_file] = analyze_entry(read_dir, fname=read_file)

    # acquire data
    for stage in os.listdir(fname):
        out_dir = os.path.join(fname, stage, 'output')
        if os.path.exists(out_dir):
            with open(os.path.join(fname, stage, 'info.txt')) as fd:
                stage_name = fd.read()

            for read in data.keys():
                cur_dir = os.path.join(out_dir, read + '_extra')
                res = analyze_entry(cur_dir)
                data[read][stage_name] = res

    # visualize data
    fig, axarr = plt.subplots(1, len(data), figsize=(20,8))

    for (read_name, dat), ax in zip(data.items(), axarr):
        names, values = [], []
        for n, v in dat.items():
            names.append(n)
            values.append(v)

        tot = sum(values)
        names.append('no match')
        values.append(totals[read_name]-tot)

        ax.pie(
            values, labels=names, explode=(.03,)*len(values),
            autopct='%1.1f%%', shadow=True)
        ax.axis('equal')
        ax.set_title(read_name)

    plt.savefig('output.pdf')

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print('Usage: {} <seq. data dir>'.format(sys.argv[0]))
        exit(-1)

    main(sys.argv[1])
