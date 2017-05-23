"""
Main dashboard for the application
"""

from typing import Sized

import click

from mapping_overview import main as plot_distribution


@click.group()
def cli() -> None:
    pass

@cli.command()
@click.option('--absolute/--no-absolute', default=False)
@click.argument('files', nargs=-1, type=click.Path(exists=True))
def read_distribution(absolute: bool, files: Sized) -> None:
    if len(files) == 0:
        print('No files provided')
        return

    plot_distribution(files, absolute)


if __name__ == '__main__':
    cli()
