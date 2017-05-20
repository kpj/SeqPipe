"""
Generate various plots for general mapping overviews
"""

import os
import sys

import seaborn as sns
import matplotlib.pyplot as plt

from utils import compute_statistics


def ensure_sanity():
    """ Check that environment looks usable
    """
    if not os.path.isdir('images'):
        os.mkdir('images')

def read_distribution(df):
    """ Plot read fractions mapped to references
    """
    total_align_num = df.groupby('read_name').sum()['mapped_count'].to_dict()
    def func(row):
        row['rel_count'] = row['mapped_count'] / total_align_num[row['read_name']]
        return row
    df = df.apply(func, axis=1)
    df.to_csv('images/read_distribution.csv')

    sub = df.pivot('read_name', 'reference')['rel_count']
    sub.plot(kind='bar', stacked=True)

    plt.tight_layout()
    plt.savefig('images/read_distribution.pdf')

def main():
    ensure_sanity()
    df = compute_statistics(sys.argv[1:])

    read_distribution(df)

if __name__ == '__main__':
    main()
