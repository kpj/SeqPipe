import seqpipe


DATA_ROOT = 'seqpipe/tests/data/'

def test_file_gathering():
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
