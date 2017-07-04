"""
CLI entry point
"""

import datetime

import click
from joblib import cpu_count


##################
# Helper functions

def get_timestamp() -> str:
    """ Return current time as formatted string
    """
    return datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")


###########
# Setup CLI

@click.group()
def cli():
    pass

@cli.command(name='map', help='Read alignment pipeline.')
@click.option(
    '--read', '-r', 'read_path_list',
    multiple=True, type=click.Path(exists=True, resolve_path=True),
    help='Path to read file/directory.', required=True)
@click.option(
    '--genome', '-g', 'genome_path_list',
    multiple=True, type=click.Path(exists=True, resolve_path=True),
    help='Path to genome file/directory.', required=True)
@click.option(
    '--output', '-o', 'output_dir',
    type=click.Path(file_okay=False, resolve_path=True),
    help='Directory to save results to.',
    default=f'mapping_result_{get_timestamp()}')
@click.option(
    '--scripts/--no-scripts', 'exec_scripts', default=True,
    help='Whether to execute scripts or not.')
@click.option(
    '--min-read-len', '-m', default=-1,
    help='Minimal read length.')
@click.option(
    '--max-read-len', '-M', default=-1,
    help='Maximal read length.')
@click.option(
    '--bowtie-args', '-b',
    help='Extra arguments for bowtie.')
@click.option(
    '--threads', '-t', default=int(cpu_count() * 4/5),
    help='How many threads to run in.')
def cmd1(*args, **kwargs) -> None:
    from .mapping.env_setup import run_pipeline
    run_pipeline(*args, **kwargs)

@cli.group(help='Post-processing tools.')
def stats() -> None:
    pass

@stats.command(name='plot_rdist', help='Plot read alignment distributions.')
@click.option(
    '--absolute/--relative', default=False,
    help='Plot absolute or relative mapped read counts.')
@click.option(
    '--split/--no-split', default=False,
    help=('Split along entries of each reference, '
          'instead of references themselves.'))
@click.option(
    '--output', '-o', 'output_dir',
    type=click.Path(file_okay=False, resolve_path=True),
    help='Directory to save results to.',
    default=f'results_{get_timestamp()}')
@click.argument('files', nargs=-1, type=click.Path(exists=True))
def cmd2(*args, **kwargs) -> None:
    from .statistics import plot_mapping_overview
    plot_mapping_overview(*args, **kwargs)

@stats.command(name='generate', help='Create CSV containing basic statistics.')
@click.option(
    '--split/--no-split', default=False,
    help=('Split along entries of each reference, '
          'instead of references themselves.'))
@click.argument('files', nargs=-1, type=click.Path(exists=True))
def cmd3(*args, **kwargs) -> None:
    from .statistics import generate_csv
    generate_csv(*args, **kwargs)

@stats.command(name='plot_covdiff', help='Plot coverage difference plots.')
@click.option(
    '--reference', 'references',
    multiple=True, help='References to consider.')
@click.option(
    '--sub-reference', 'sub_references',
    multiple=True, help='Sub-references to consider.')
@click.option(
    '--output', '-o', 'output_dir',
    type=click.Path(file_okay=False, resolve_path=True),
    help='Directory to save results to.',
    default=f'results_{get_timestamp()}')
@click.argument('files', nargs=-1, type=click.Path(exists=True))
def cmd4(*args, **kwargs) -> None:
    from .statistics import plot_coverage_difference
    plot_coverage_difference(*args, **kwargs)


if __name__ == '__main__':
    cli()
