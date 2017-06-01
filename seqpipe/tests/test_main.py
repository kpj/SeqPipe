import os
import json
import tempfile

import pysam
import numpy as np
import pandas as pd

import pytest
import pandas.util.testing as pdt
from click.testing import CliRunner

import seqpipe


DATA_ROOT = 'seqpipe/tests/data/'

def test_file_gathering() -> None:
    fnames = seqpipe.gather_files([
        f'{DATA_ROOT}/references/00-ref.fa',
        f'{DATA_ROOT}/reads',
        f'{DATA_ROOT}/references/20-ref.fa'])

    assert fnames == [
        f'{DATA_ROOT}/references/00-ref.fa',
        f'{DATA_ROOT}/reads/bar.fastq.gz',
        f'{DATA_ROOT}/reads/baz.fastq.gz',
        f'{DATA_ROOT}/reads/foo.fastq.gz',
        f'{DATA_ROOT}/references/20-ref.fa'
    ]

@pytest.mark.skipif(
    'TRAVIS' in os.environ and os.environ['TRAVIS'] == 'true',
    reason='Skip on Travis CI.')
def test_whole_pipeline() -> None:
    runner = CliRunner()

    with tempfile.TemporaryDirectory() as tmpdirname:
        output_dir = f'{tmpdirname}/output/'

        ## test mapping
        result = runner.invoke(seqpipe.cli, [
            'map',
            '--read', f'{DATA_ROOT}/reads/foo.fastq.gz',
            '--genome', f'{DATA_ROOT}/references/20-ref.fa',
            '-o', output_dir, '-b', '-a -v0'
        ])
        assert result.exit_code == 0

        # test metadata
        with open(f'{output_dir}/info.json') as fd:
            meta = json.load(fd)
        assert meta['read_path_list'] == [
            os.path.abspath(f'{DATA_ROOT}/reads/foo.fastq.gz')]
        assert meta['genome_path_list'] == [
            os.path.abspath(f'{DATA_ROOT}/references/20-ref.fa')]

        # test mapping result
        res_path = os.path.join(output_dir, 'runs', '20-ref.fa___foo.fastq.gz')

        bam_path = os.path.join(res_path, 'data', 'mapping_sorted.bam')
        pysam.index(bam_path)
        bamf = pysam.AlignmentFile(bam_path, 'rb')
        assert bamf.has_index()

        assert bamf.references == (
            'checkAAA', 'checkCCC',
            'special', 'anotherSpecial')

        assert bamf.count(reference='checkAAA', start=0, end=1) == 3
        assert bamf.count(reference='checkAAA', start=0, end=8) == 18

        assert bamf.count(reference='checkCCC', start=0, end=1) == 2
        assert bamf.count(reference='checkCCC', start=0, end=8) == 12

        assert bamf.count(reference='special', start=0, end=7) == 1
        assert bamf.count(reference='anotherSpecial', start=0, end=15) == 0


        ## test statistics
        stats_dir = os.path.join(tmpdirname, 'statistics/')

        result = runner.invoke(seqpipe.cli, [
            'stats', 'plot_rdist',
            '-o', stats_dir,
            output_dir
        ])
        assert result.exit_code == 0

        df = pd.read_csv(
            os.path.join(stats_dir, 'read_distribution.csv'),
            index_col=0)
        pdt.assert_frame_equal(df.sort_index(axis=1), pd.DataFrame({
            'mapped_count': [31.],
            'read_name': ['foo.fastq'],
            'reference': ['20-ref.fa'],
            'sub_reference': [np.nan],
            'total_count': [6.],
            'rel_count': [1.]
        }).sort_index(axis=1))
