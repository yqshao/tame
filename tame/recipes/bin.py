"""Main entry point of tame"""
import tame
import click
from tame.recipes.jacf    import jacf_cmd
from tame.recipes.mdc     import mdc_cmd
from tame.recipes.onsager import onsager_cmd
from tame.recipes.persist import persist_cmd
from tame.recipes.pmsd    import pmsd_cmd
from tame.recipes.rdf     import rdf_cmd

@click.group()
def main():
    """TAME CLI - Command line interface of TAME"""
    pass

main.add_command(jacf_cmd)
main.add_command(mdc_cmd)
main.add_command(onsager_cmd)
main.add_command(persist_cmd)
main.add_command(pmsd_cmd)
main.add_command(rdf_cmd)

@click.command()
def version():
    click.echo(f'TAME version: {tame.__version__}')

if __name__ == "__main__":
    main()

