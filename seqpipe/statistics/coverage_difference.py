"""
Plot coverage difference plots
"""

import os
import json
import itertools
import collections

from typing import List

import numpy as np
import pandas as pd

import seaborn as sns
import matplotlib.pyplot as plt

import pysam
from tqdm import tqdm


Conf = collections.namedtuple(
    'config', ['read1', 'read2', 'ref', 'sub_ref'])

def gather_files(dir_: str) -> List:
    """ Find all alignment files
    """
    result = []
    for entry in tqdm(os.scandir(dir_), total=len(os.listdir(dir_))):
        aligned_bam_file = os.path.join(entry.path, 'aligned_reads.bam')

        pysam.index(aligned_bam_file)
        samfile = pysam.AlignmentFile(aligned_bam_file, 'rb')
        assert samfile.has_index()

        with open(os.path.join(entry.path, 'meta.json')) as fd:
            meta = json.load(fd)

        result.append({
            'read': meta['read_base'],
            'reference': meta['genome_base'],
            'sam': samfile
        })
    return result

def compute_coverage(data: List) -> pd.DataFrame:
    """ Compute per-base coverage
    """
    result = []
    for entry in tqdm(data):
        sam = entry['sam']
        tmp = {}
        for ref, slen in zip(sam.references, sam.lengths):
            cov_tmp = sam.count_coverage(ref, 0, slen)
            cov = np.asarray(cov_tmp).sum(axis=0).astype(int)

            assert len(cov_tmp[0]) == cov.shape[0] == slen
            tmp[ref] = cov

        result.append({
            'reference': entry['reference'],
            'read': entry['read'],
            'coverage': tmp
        })
    return pd.DataFrame(result)

def plot_entry(
    cov1: np.ndarray, cov2: np.ndarray,
    conf: Conf,
    output_dir: str
) -> None:
    """ Plot a particular configuration
    """
    max_cov = max(cov1.max(), cov2.max()) or 1

    plt.figure()

    # individual plots
    plt.subplot(311)
    plt.plot(cov1)
    plt.title(conf.read1)
    plt.ylim((0, max_cov))

    plt.subplot(312)
    plt.plot(cov1-cov2)
    plt.title(f'Difference: {conf.read1}-{conf.read2}')
    plt.ylim((-max_cov, max_cov))

    plt.subplot(313)
    plt.plot(cov2)
    plt.title(conf.read2)
    plt.ylim((0, max_cov))

    # surrounding graphics
    plt.suptitle(f'reference: {conf.ref}, sub-reference: {conf.sub_ref}')

    with sns.axes_style('white'):
        plt.gcf().add_subplot(111, frameon=False)
        plt.tick_params(
            labelcolor='none',
            top='off', bottom='off', left='off', right='off')
        plt.xlabel('base position')
        plt.ylabel('base coverage')

    # save result
    plt.tight_layout()
    plt.subplots_adjust(top=0.85)

    fname = os.path.join(
        output_dir,
        f'covdiff_{conf.ref}_{conf.sub_ref}_{conf.read1}_{conf.read2}.pdf')
    plt.savefig(fname)

def plot_coverage_differences(df: pd.DataFrame, output_dir: str) -> None:
    """ Create difference plots
    """
    for ref, group in tqdm(df.groupby('reference')):
        all_reads = group['read'].unique()
        for read1, read2 in itertools.combinations(all_reads, 2):
            tmp1 = group[group['read']==read1]
            tmp2 = group[group['read']==read2]
            assert tmp1.shape[0] == 1 == tmp2.shape[0]

            cov1_d = tmp1.iloc[0]['coverage']
            cov2_d = tmp2.iloc[0]['coverage']
            assert cov1_d.keys() == cov2_d.keys()

            for sub in cov1_d.keys():
                conf = Conf(read1=read1, read2=read2, ref=ref, sub_ref=sub)
                plot_entry(cov1_d[sub], cov2_d[sub], conf, output_dir)

def main(files: List, output_dir: str) -> None:
    if len(files) == 0:
        print('No files provided')
        return

    data = []
    for fname in tqdm(files):
        fn = os.path.join(fname, 'runs')
        tmp = gather_files(fn)
        data.append((os.path.dirname(fname), tmp))

    df_list = []
    for idx, tmp in tqdm(data):
        df_tmp = compute_coverage(tmp)
        df_list.append(df_tmp)
    df = pd.concat(df_list)

    plot_coverage_differences(df, output_dir)
