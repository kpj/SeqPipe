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
    print('Pattern: "{}"'.format(pattern))
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

def annotate_peak_positions(references, df_cov, depth_thres=1e4):
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
    print('Threshold:', depth_thres)
    for ref in references:
        df = df_cov[(df_cov.reference == ref) & (df_cov.depth > depth_thres)]

        data[ref] = []
        start = None
        cur = None
        strand = None
        for index, row in df.iterrows():
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
                        'name': '{}.\n.{}'.format(start, cur),
                        'strand': strand
                    })

                    start = row.position
                    cur = start
                    strand = row.strand

        if not start is None and start != cur:
            data[ref].append({
                'start': start,
                'end': cur,
                'name': '{}.\n.{}'.format(start, cur),
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
            os.system('samtools view -o {} {} {}:{}-{}'.format(
                tf.name, bam_fname, ref, start, end))

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
            score_string += str(int(abs(score * 10) % 10))

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
            pos_range = '{}{}{}'.format(
                start, ' '*(len(ebases)-len(str(start))-len(str(end))), end)

            annotation = '{}\n{}\n{}'.format(pos_range, ebases, score)
            data[ref].append({
                'start': start,
                'end': end,
                'name': annotation
            })

            with open(cache_fname, 'a') as fd:
                fd.write('{} strand\n'.format(peak['strand']) + annotation + '\n\n')
    return data

def file_checks(base_fname, two_strands=True):
    """
    Assure that all needed files are in place.

    Arguments
        *similar to main*
    """
    coverage_fname = '{}.cov'.format(base_fname)
    bam_fname = '{}.bam'.format(base_fname)

    assert os.path.isfile(bam_fname), bam_fname

    # check that required files are available
    bam_index_fname = '{}.bai'.format(bam_fname)
    if not os.path.isfile(bam_index_fname):
        print('No bam index found, generating one...')
        cmd = 'samtools index {inp}'.format(inp=bam_fname)
        print(' > "{}"...'.format(cmd), end=' ', flush=True)
        os.system(cmd)
        print('Done')
    else:
        print('Using cached bam index')

    if not os.path.isfile(coverage_fname):
        if two_strands:
            bam_p = '{}_pstrand.bam'.format(base_fname)
            bam_m = '{}_mstrand.bam'.format(base_fname)

            cov_p = '{}_pstrand.cov'.format(base_fname)
            cov_m = '{}_mstrand.cov'.format(base_fname)

            cmds = [
                'samtools view -F 16 {inp} -b > {outp}'.format(inp=bam_fname, outp=bam_p),
                'samtools view -f 16 {inp} -b > {outp}'.format(inp=bam_fname, outp=bam_m),

                'genomeCoverageBed -ibam {inp} -d > {outp}'.format(inp=bam_p, outp=cov_p),
                'genomeCoverageBed -ibam {inp} -d > {outp}'.format(inp=bam_m, outp=cov_m),

                'sed -e "s/$/\t+/" -i {fname}'.format(fname=cov_p),
                'sed -e "s/$/\t-/" -i {fname}'.format(fname=cov_m),

                'cat {inp1} {inp2} > {outp}'.format(inp1=cov_p, inp2=cov_m, outp=coverage_fname)
            ]
        else:
            cmds = [
                'genomeCoverageBed -ibam {inp} -d > {outp}'.format(inp=bam_fname, outp=coverage_fname),
                'sed -e "s/$/\t+-/" -i {fname}'.format(fname=coverage_fname)
            ]

        print('No coverage file found, generating one...')
        for cmd in cmds:
            print(' > "{}"...'.format(cmd), end=' ', flush=True)
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
    df_cov = pd.read_table(coverage_fname, header=None, names=['reference', 'position', 'depth', 'strand'])
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

        plt.title('{} (total base hits: {})'.format(ref, read_num))
        plt.xlabel('bp position')
        plt.ylabel('read count')

        plt.xlim((0, max(pos)))
        plt.legend(loc='best')

        plt.tight_layout()
        fig.savefig('results/coverage_{}.pdf'.format(ref))

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print('Usage: {} <bam without prefix>'.format(sys.argv[0]))
        exit(-1)

    main(sys.argv[1])
