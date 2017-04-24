"""
Visualize expression differences
"""

import os
import sys
import itertools
import collections

import pysam
import numpy as np

import seaborn as sns
import matplotlib.pyplot as plt

from tqdm import tqdm


def gather_files(stage_dir):
    """ Find all alignment files
    """
    result = []
    for stage_entry in os.scandir(stage_dir):
        read_file_dir = os.path.join(stage_entry.path, 'mapping', 'runs')

        if not os.path.exists(read_file_dir):
            continue
        print('>', stage_entry.name)

        for read_entry in os.scandir(read_file_dir):
            print(' >>', read_entry.name)
            aligned_bam_file = os.path.join(read_entry.path, 'aligned_reads.bam')

            pysam.index(aligned_bam_file)
            samfile = pysam.AlignmentFile(aligned_bam_file, 'rb')
            assert samfile.has_index()

            with open(os.path.join(stage_entry.path, 'info.txt')) as fd:
                stage = fd.read().strip()

            result.append({
                'stage': stage,
                'read': read_entry.name,
                'sam': samfile
            })
    return result

def compute_coverage(alig_data):
    """ Compute per-base coverage
    """
    result = []
    for entry in tqdm(alig_data):
        sam = entry['sam']
        tmp = {}
        for ref, slen in zip(sam.references, sam.lengths):
            cov_tmp = sam.count_coverage(ref, 0, slen)
            cov = np.asarray(cov_tmp).sum(axis=0).astype(int)

            assert len(cov_tmp[0]) == cov.shape[0] == slen
            tmp[ref] = cov

        result.append({
            'stage': entry['stage'],
            'read': entry['read'],
            'coverage': tmp
        })
    return result

def plot_result(cov_data):
    """ Visualize coverage results
    """
    def plot(stage, ref, read1, read2, cov1, cov2):
        max_cov = max(cov1.max(), cov2.max())

        plt.figure()

        # individual plots
        plt.subplot(311)
        plt.plot(cov1)
        plt.title(read1)
        plt.ylim((0,max_cov))

        plt.subplot(312)
        plt.plot(cov1-cov2)
        plt.title(f'Difference: {read1}-{read2}')
        plt.ylim((-max_cov,max_cov))

        plt.subplot(313)
        plt.plot(cov2)
        plt.title(read2)
        plt.ylim((0,max_cov))

        # surrounding graphics
        plt.suptitle(f'Stage: {stage}, reference: {ref}')

        with sns.axes_style('white'):
            plt.gcf().add_subplot(111, frameon=False)
            plt.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
            plt.xlabel('base position')
            plt.ylabel('base coverage')

        # save result
        plt.tight_layout()
        plt.subplots_adjust(top=0.85)
        plt.savefig(f'images/diffexpr_{stage}_{ref}_{read1}_{read2}.pdf')

    # aggregate data
    cov_tmp = collections.defaultdict(lambda: collections.defaultdict(list))
    for entry in cov_data:
        cov_tmp[entry['read']][entry['stage']].append(entry['coverage'])

    cov_tmp = dict(cov_tmp)
    all_reads = list(cov_tmp.keys())

    # find read-pairs and plot
    for read1,read2 in tqdm(itertools.combinations(all_reads, 2)):
        stages1 = cov_tmp[read1]
        stages2 = cov_tmp[read2]
        assert stages1.keys() == stages2.keys()

        for stage in tqdm(stages1):
            refs1 = stages1[stage]
            refs2 = stages2[stage]
            assert len(refs1) == len(refs2) == 1

            refs1 = refs1[0]
            refs2 = refs2[0]
            assert refs1.keys() == refs2.keys()

            for ref in refs1:
                assert len(refs1[ref]) == len(refs2[ref])
                plot(stage, ref, read1, read2, refs1[ref], refs2[ref])

def main(dname):
    alig_data = gather_files(dname)
    cov_data = compute_coverage(alig_data)

    plot_result(cov_data)


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print(f'Usage: {sys.argv[0]} <seq. data dir>')
        exit(-1)

    main(sys.argv[1])
