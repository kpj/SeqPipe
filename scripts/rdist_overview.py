"""
Generate read distribution overview
"""

import os
import json
import collections

from typing import List, Dict

import pysam

import matplotlib.pyplot as plt
from matplotlib_venn import venn2


def main(bam_files: List[str], report_fname: str) -> None:
    # get reads per bam and group accordingly
    read_ids = collections.defaultdict(dict)  # type: Dict[str, Dict[str, List]]
    for fname in bam_files:
        alf = pysam.AlignmentFile(fname)
        read_names = [seg.query_name for seg in alf.fetch()]

        dname = os.path.dirname(fname)
        read_name = dname.split('/')[-1]
        ref_name = dname.split('/')[-2]

        read_ids[read_name][ref_name] = read_names
    read_ids = dict(read_ids)

    # plot result
    for read_name, ref_map in read_ids.items():
        if len(ref_map) != 2:
            print(f'[{read_name}] Skipping Venn diagrams as not exactly two references are present')
            continue

        label_list, data_list = [], []
        for ref_name, reads in ref_map.items():
            label_list.append(ref_name)
            data_list.append(set(reads))

        plt.figure()
        venn2(subsets=data_list, set_labels=label_list)
        plt.title(f'Readfile: "{read_name}"')
        plt.savefig(f'{os.path.dirname(report_fname)}/{read_name}_venn.pdf')

    # save raw data
    with open(report_fname, 'w') as fd:
        for read_name, ref_map in read_ids.items():
            for ref_name, reads in ref_map.items():
                fd.write(f'{read_name}\t{ref_name}\t{len(reads)}\n')

if __name__ == '__main__':
    main(snakemake.input.bam_files, snakemake.output.report)
