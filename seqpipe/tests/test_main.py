import os
import json
import tempfile

import pysam
from click.testing import CliRunner

import seqpipe


DATA_ROOT = 'seqpipe/tests/data/'

def test_file_gathering() -> None:
    fnames = seqpipe.main.gather_files([
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

def test_whole_pipeline() -> None:
    with tempfile.TemporaryDirectory() as tmpdirname:
        output_dir = f'{tmpdirname}/output/'

        runner = CliRunner()
        result = runner.invoke(seqpipe.main.run, [
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
