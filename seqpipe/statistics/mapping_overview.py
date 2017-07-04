"""
Generate various plots for general mapping overviews
"""

import os
import sys

from typing import Dict, List

import pandas as pd

import matplotlib as mpl
mpl.use('Agg')

import seaborn as sns
import matplotlib.pyplot as plt

from .utils import compute_statistics


def compute_relative_counts(
    row: pd.Series,
    total_align_num: Dict[str, int]
) -> pd.Series:
    """ Compute fraction of total number of alignments
    """
    row['rel_count'] = row['mapped_count'] / total_align_num[row['read_name']]
    return row

def read_distribution(
    df: pd.DataFrame,
    absolute: bool, split: bool, output_dir: str
) -> None:
    """ Plot read fractions mapped to references
    """
    total_align_num = df.groupby('read_name').sum()['mapped_count'].to_dict()
    df = df.apply(compute_relative_counts, axis=1, args=(total_align_num,))
    df.to_csv(f'{output_dir}/read_distribution.csv')

    count_ref = 'mapped_count' if absolute else 'rel_count'
    ref_key = 'sub_reference' if split else 'reference'

    sub = df.pivot('read_name', ref_key)[count_ref]
    sub.plot(kind='bar', stacked=True)

    plt.tight_layout()
    plt.savefig(f'{output_dir}/read_distribution.pdf')

def main(
    files: List, output_dir: str,
    absolute: bool, split: bool = False
) -> None:
    if len(files) == 0:
        print('No files provided')
        return

    df = compute_statistics(files, split=split)

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    read_distribution(df, absolute, split, output_dir)
