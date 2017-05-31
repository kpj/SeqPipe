"""
Generate various plots for general mapping overviews
"""

import os
import sys

from typing import Dict, Sized

import pandas as pd

import seaborn as sns
import matplotlib.pyplot as plt

from .utils import compute_statistics


def ensure_sanity() -> None:
    """ Check that environment looks usable
    """
    if not os.path.isdir('images'):
        os.mkdir('images')

def compute_relative_counts(
    row: pd.Series,
    total_align_num: Dict[str, int]
) -> pd.Series:
    """ Compute fraction of total number of alignments
    """
    row['rel_count'] = row['mapped_count'] / total_align_num[row['read_name']]
    return row

def read_distribution(df: pd.DataFrame, absolute: bool) -> None:
    """ Plot read fractions mapped to references
    """
    total_align_num = df.groupby('read_name').sum()['mapped_count'].to_dict()
    df = df.apply(compute_relative_counts, axis=1, args=(total_align_num,))
    df.to_csv('images/read_distribution.csv')

    count_ref = 'mapped_count' if absolute else 'rel_count'
    sub = df.pivot('read_name', 'reference')[count_ref]
    sub.plot(kind='bar', stacked=True)

    plt.tight_layout()
    plt.savefig('images/read_distribution.pdf')

def main(files: Sized, absolute: bool = False) -> None:
    if len(files) == 0:
        print('No files provided')
        return

    ensure_sanity()
    df = compute_statistics(files)

    read_distribution(df, absolute)

if __name__ == '__main__':
    main(sys.argv[1:])
