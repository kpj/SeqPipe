"""
Hub for organizing all sequencing related things
"""

import re
import io
import os
import gzip
import json
import shutil
import itertools
import collections

from typing import List, Tuple, Dict

import sh
import click
from joblib import Parallel, delayed

from .pipeline import pipeline, USED_TOOLS


def gather_files(path_list: List[str]) -> List[str]:
    """ Gather all files (eg from directories) in given list
    """
    out = []
    for path in path_list:
        if os.path.isfile(path):
            out.append(path)
        elif os.path.isdir(path):
            for dirpath, _, filenames in os.walk(path):
                out.extend([os.path.join(dirpath, fn)
                    for fn in sorted(filenames)])
        else:
            raise RuntimeError(f'Did not understand {path}')

    # return unique paths while preserving order
    return list(collections.OrderedDict.fromkeys(out))

def all_tools_available() -> bool:
    """ Check that all needed tools are available on startup
    """
    missing = []
    for tool in USED_TOOLS:
        try:
            sh.Command(tool)
        except sh.CommandNotFound:
            missing.append(tool)

    if len(missing) > 0:
        print(f'Cannot start pipeline, missing tools: {missing}')

    return len(missing) == 0

class PipelineStream(io.StringIO):
    """ Print and save output at the same time
    """
    def __init__(self, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)
        self.ansi_escape = re.compile(r'\x1b[^m]*m')

    def write(self, s: str) -> int:
        print(s, end='', flush=True)
        return super().write(self.ansi_escape.sub('', s))

class SequencingRun:
    """ Store information regarding a single sequencing run
    """
    def __init__(
        self,
        read_path: str, genome_path: str,
        output_dir: str
    ) -> None:
        self.read_path = read_path
        self.genome_path = genome_path
        self.output_dir = output_dir

    def __call__(self, param_obj: Dict) -> Dict:
        # prepare data/directory structures
        rp, gp, od = self._prepare_environment(
            self.read_path, self.genome_path, self.output_dir)

        # execute pipeline
        stream = PipelineStream()
        res = pipeline(rp, gp, od, param_obj, stream)
        log_value = stream.getvalue()

        # save log
        log_path = os.path.join(od, 'log.txt')
        with open(log_path, 'w') as fd:
            fd.write(log_value)

        return res

    @staticmethod
    def _prepare_environment(
        read_path: str, genome_path: str,
        output_dir: str
    ) -> Tuple[str, str, str]:
        """ Setup needed data structures
        """
        # set prefixes
        cur_dir = os.path.dirname(os.path.abspath(__file__))

        read_base = os.path.basename(read_path)
        genome_base = os.path.basename(genome_path)
        pipeline_dir = os.path.join(
            output_dir, 'runs',
            f'{genome_base}___{read_base}')

        read_remote = os.path.join(
            pipeline_dir, 'input', os.path.splitext(read_base)[0])
        genome_remote = os.path.join(pipeline_dir, 'input', genome_base)

        os.makedirs(os.path.join(pipeline_dir, 'input'))
        os.makedirs(os.path.join(pipeline_dir, 'results'))

        # copy needed files/directories
        shutil.copyfile(genome_path, genome_remote)

        read_base = os.path.basename(read_path)
        if read_base.endswith('.gz'):
            with gzip.open(read_path, 'rb') as fd_in:
                with open(read_remote, 'wb') as fd_out:
                    for line in fd_in:
                        fd_out.write(line)
        else:
            shutil.copyfile(read_path, read_remote)

        shutil.copytree(
            os.path.join(cur_dir, 'scripts'),
            os.path.join(pipeline_dir, 'scripts'))

        return read_remote, genome_remote, pipeline_dir

def run_pipeline(
    read_path_list: List[str],
    genome_path_list: List[str],
    output_dir: str,
    exec_scripts: bool,
    min_read_len: int, max_read_len: int,
    bowtie_args: str,
    threads: int
) -> None:
    """ Create read-genome matrix and compute all read alignments.
        Subsequently, apply various scripts and aggregate results.
    """
    if not all_tools_available():
        exit(-1)

    if len(read_path_list) == 0:
        print('No paths given...')
        exit(-1)
    if len(genome_path_list) == 0:
        print('No genomes given...')
        exit(-1)
    if os.path.exists(output_dir):
        print('Output directory already exists')
        exit(-1)

    # clean filenames
    read_path_list = gather_files(read_path_list)
    genome_path_list = gather_files(genome_path_list)

    # prepare sequencing runs
    run_list = []
    for read, genome in itertools.product(read_path_list, genome_path_list):
        run_list.append(SequencingRun(read, genome, output_dir))

    # save meta-information
    param_obj = click.get_current_context().params

    os.makedirs(output_dir)
    param_path = os.path.join(output_dir, 'info.json')
    with open(param_path, 'w') as fd:
        json.dump(param_obj, fd)

    # commence pipelines
    results = Parallel(n_jobs=threads)(
        delayed(run)(param_obj) for run in run_list)

    # aggregate results
    result_dir = os.path.join(output_dir, 'results/')
    os.makedirs(result_dir)

    print('Gathering results')
    found_result = False
    for res in results:
        idx = f'{res["genome_base"]}___{res["read_base"]}'
        cur_dir = os.path.join(result_dir, idx)
        os.makedirs(cur_dir)

        for entry in os.scandir(res['results_path']):
            print(f' - {entry.path}')
            shutil.copy(entry.path, cur_dir)
            found_result = True
    if not found_result:
        print('  -- no results found --')
