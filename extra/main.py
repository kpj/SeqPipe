"""
Main dashboard for the application
"""

from typing import Sized

import click

from mapping_overview import main as plot_distribution


@click.group()
def cli() -> None:
    """ Toolkit to analyse sequencing data
    """
    pass

@cli.command()
@click.option('--absolute/--relative', default=False, help='Plot absolute or relative mapped read counts')
@click.argument('files', nargs=-1, type=click.Path(exists=True))
def rdist(absolute: bool, files: Sized) -> None:
    """ Plot read alignment distribution
    """
    if len(files) == 0:
        print('No files provided')
        return

    plot_distribution(files, absolute)


if __name__ == '__main__':
    cli()
