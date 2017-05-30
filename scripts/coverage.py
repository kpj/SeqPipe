#!/usr/bin/env python3

"""
Investigate read coverage

Generate coverage, etc with
  genomeCoverageBed -ibam foo.bam -d > foo.cov

  bedtools genomecov -ibam foo.bam -bga > foo.bed
  annotatePeaks.pl foo.bed none -gff3 foo.gff > foo.ann
"""

import os
import sys
import tempfile
import itertools
import collections

import numpy as np
import pandas as pd

import seaborn as sns
import matplotlib.pylab as plt
import matplotlib.ticker as mtick
from matplotlib.font_manager import FontProperties

from tqdm import tqdm


def compute_base_enrichment(sequences):
    """
    Compute relative base enrichment at each position.

    Arguments
        sequences
            List of sequences to be analyzed
    """
    dic = collections.defaultdict(lambda: collections.defaultdict(int))
    for seq in sequences:
        for pos, base in enumerate(seq):
            dic[pos][base] += 1
    return dict({k: dict(v) for k,v in dic.items()})

def extract_name(attribute):
    """
    Extract name from attribute field of GFF.

    Arguments
        attribute
            Attribute string from GFF
    """
    dattr = dict((k.strip(), v.strip()) for k,v in (item.split('=') for item in attribute.split(';')))
    return dattr['Name']

def annotate_with_gff(references, gff_fname, pattern):
    """
    Annotate `pattern` matches in GFF.

    Arguments
        references
            Reference sequences to find annotations for
        gff_fname
            Filename of GFF data
        pattern
            String to search for in GFF attributes
    """
    # read DataFrame
    df_gff = pd.read_table(
        gff_fname, comment='#', header=None,
        names=['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes'])
    df_gff.start = df_gff.start.astype(int)
    df_gff.end = df_gff.end.astype(int)
    print('Done')

    # process data
    data = {}
    print(f'Pattern: "{pattern}"')
    for ref in references:
        df = df_gff[(df_gff.seqid == ref) & df_gff.attributes.str.contains(pattern)]

        data[ref] = []
        for index, row in df.iterrows():
            data[ref].append({
                'start': row.start,
                'end': row.end,
                'name': extract_name(row.attributes)
            })

        print(' >', ref, len(data[ref]))
    return data

def annotate_peak_positions(references, df_cov, depth_thres_frac=.5):
    """
    Annotate peaks with exact bp positions.

    Arguments
        references
            Reference sequences to find annotations for
        df_cov
            DataFrame with coverage information
        depth_thres
            Threshold on how many reads constitute a peak
    """
    data = {}
    for ref in references:
        df = df_cov[df_cov.reference == ref]
        max_depth = df['depth'].max()
        depth_thres = max_depth * depth_thres_frac
        print(f'[{ref}] Threshold: {depth_thres} ({depth_thres_frac} * {max_depth})')
        sub = df[df.depth > depth_thres]

        data[ref] = []
        start = None
        cur = None
        strand = None
        for index, row in sub.iterrows():
            if start is None:
                start = row.position
                cur = start
                strand = row.strand
            else:
                if cur+1 == row.position:
                    cur += 1
                else:
                    data[ref].append({
                        'start': start,
                        'end': cur,
                        'name': f'{start}.\n.{cur}',
                        'strand': strand
                    })

                    start = row.position
                    cur = start
                    strand = row.strand

        if not start is None and start != cur:
            data[ref].append({
                'start': start,
                'end': cur,
                'name': f'{start}.\n.{cur}',
                'strand': strand
            })

        print(' >', ref, len(data[ref]))
    return data

def annotate_read_pattern(
    references, df_cov, bam_fname,
    cache_fname='data/coverage_annos.txt'
):
    """
    Annotate each peak with corresponding read pattern.

    Arguments
        references
            Reference sequences to find annotations for
        df_cov
            DataFrame with coverage information
        bam_fname
            Filename of bam with needed alignment information
        cache_fname
            Where to save annotations to
    """
    def get_reads(ref, start, end):
        sequences = []
        with tempfile.NamedTemporaryFile(suffix='sam') as tf:
            os.system(f'samtools view -o {tf.name} {bam_fname} {ref}:{start}-{end}')

            for line in tf:
                parts = line.split()
                pos = int(parts[3])
                seq = parts[9].decode('utf-8')

                if pos == start:
                    sequences.append(seq)
        return sequences

    def extract_enriched_bases(base_freqs):
        enriched_bases = ''
        score_string = ''

        for pos, bases in base_freqs.items():
            best_base = max(bases, key=lambda b: bases[b])

            total_count = sum(bases.values())
            score = bases[best_base] / total_count

            enriched_bases += best_base

            cur = str(int(abs(score * 10) % 10))
            score_string += '\u2714' if cur == '0' else cur

        return enriched_bases, score_string

    peak_data = annotate_peak_positions(references, df_cov)
    open(cache_fname, 'w').close() # clear old cache

    data = {}
    for ref, peak_list in tqdm(peak_data.items()):
        with open(cache_fname, 'a') as fd:
            fd.write(ref + '\n')

        data[ref] = []
        for peak in sorted(peak_list, key=lambda x: x['strand']):
            start, end = peak['start'], peak['end']

            reads = get_reads(ref, start, end)
            base_freqs = compute_base_enrichment(reads)

            ebases, score = extract_enriched_bases(base_freqs)
            pos_range = f'{start}{" "*(len(ebases)-len(str(start))-len(str(end)))}{end}'

            annotation = f'{pos_range}\n{ebases}\n{score}'
            data[ref].append({
                'start': start,
                'end': end,
                'name': annotation
            })

            with open(cache_fname, 'a') as fd:
                fd.write(f'{peak["strand"]} strand\n' + annotation + '\n\n')
    return data

def file_checks(base_fname, two_strands=True):
    """
    Assure that all needed files are in place.

    Arguments
        *similar to main*
    """
    coverage_fname = f'{base_fname}.cov'
    bam_fname = f'{base_fname}.bam'

    assert os.path.isfile(bam_fname), bam_fname

    # check that required files are available
    bam_index_fname = f'{bam_fname}.bai'
    if not os.path.isfile(bam_index_fname):
        print('No bam index found, generating one...')
        cmd = f'samtools index {bam_fname}'
        print(f' > "{cmd}"...', end=' ', flush=True)
        os.system(cmd)
        print('Done')
    else:
        print('Using cached bam index')

    if not os.path.isfile(coverage_fname):
        if two_strands:
            bam_p = f'{base_fname}_pstrand.bam'
            bam_m = f'{base_fname}_mstrand.bam'

            cov_p = f'{base_fname}_pstrand.cov'
            cov_m = f'{base_fname}_mstrand.cov'

            cmds = [
                f'samtools view -F 16 {bam_fname} -b > {bam_p}',
                f'samtools view -f 16 {bam_fname} -b > {bam_m}',

                f'genomeCoverageBed -ibam {bam_p} -d > {cov_p}',
                f'genomeCoverageBed -ibam {bam_m} -d > {cov_m}',

                f'sed -e "s/$/\t+/" -i {cov_p}',
                f'sed -e "s/$/\t-/" -i {cov_m}',

                f'cat {cov_p} {cov_m} > {coverage_fname}'
            ]
        else:
            cmds = [
                f'genomeCoverageBed -ibam {bam_fname} -d > {coverage_fname}',
                f'sed -e "s/$/\t+-/" -i {coverage_fname}'
            ]

        print('No coverage file found, generating one...')
        for cmd in cmds:
            print(f' > "{cmd}"...', end=' ', flush=True)
            ret_code = os.system(cmd)
            print('Done' if ret_code == 0 else 'Fail')
    else:
        print('Using cached coverage files')

    return bam_fname, coverage_fname

def main(base_fname):
    """
    Main interface.

    Arguments
        base_fname
            Base filename of data to be investigated
    """
    # assemble filenames
    bam_fname, coverage_fname = file_checks(base_fname)

    # get data
    print('Reading data...', end=' ', flush=True)
    df_cov = pd.read_table(
        coverage_fname,
        header=None, names=['reference', 'position', 'depth', 'strand'])
    df_cov.position = df_cov.position.astype(int)
    df_cov.depth = df_cov.depth.astype(int)

    # find interesting regions from annotations
    references = df_cov.reference.unique()

    # find annotations
    #annos = collections.defaultdict(list)
    #annos = annotate_peak_positions(references, df_cov)
    #annos = annotate_with_gff(references, 'data/gff3/chromosome_pure.gff', 'snoRNA')
    annos = annotate_read_pattern(references, df_cov, bam_fname)

    for ref in tqdm(references):
        pos = df_cov[df_cov.reference==ref].position
        depth = df_cov[(df_cov.reference==ref) & (df_cov.strand == '+-')].depth

        pos_p = df_cov[(df_cov.reference==ref) & (df_cov.strand == '+')].position
        depth_p = df_cov[(df_cov.reference==ref) & (df_cov.strand == '+')].depth
        pos_m = df_cov[(df_cov.reference==ref) & (df_cov.strand == '-')].position
        depth_m = df_cov[(df_cov.reference==ref) & (df_cov.strand == '-')].depth

        # plot result
        fig = plt.figure(figsize=(33,11))

        font = FontProperties()
        font.set_family('monospace')

        #plt.gca().xaxis.set_major_formatter(mtick.FormatStrFormatter('%.e'))
        #plt.gca().yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))
        #plt.yscale('log')

        frac = 4/5
        if len(depth) == 0:
            plt.plot(pos_p, depth_p, label='+ strand')
            plt.plot(pos_m, -depth_m, label='- strand')

            min_d = -frac * max(depth_m)
            max_d = frac * max(depth_p)

            read_num = depth_m.sum() + depth_p.sum()
        else:
            plt.plot(pos, depth, label='+- strand')

            min_d = (1-frac) * max(depth)
            max_d = frac * max(depth)
            read_num = depth.sum()

        annotation_pos_cycler = itertools.cycle(
            np.linspace(int(min_d), int(max_d), 10))
        for e in annos[ref]:
            plt.axvspan(
                xmin=e['start'], xmax=e['end'],
                alpha=0.2, color='red')

            plt.annotate(e['name'],
                xy=(e['start'], next(annotation_pos_cycler)), xycoords='data',
                xytext=(-50, 20), textcoords='offset points',
                arrowprops=dict(arrowstyle='->'),
                fontsize=3, fontproperties=font)

        plt.title(f'{ref} (total base hits: {read_num})')
        plt.xlabel('bp position')
        plt.ylabel('read count')

        plt.xlim((0, max(pos)))
        plt.legend(loc='best')

        plt.tight_layout()
        fig.savefig(f'results/coverage_{ref}.pdf')

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print(f'Usage: {sys.argv[0]} <bam without prefix>')
        exit(-1)

    main(sys.argv[1])
