"""Main entry point of tame"""

import argparse

def main():
    from importlib import import_module
    parser = argparse.ArgumentParser(prog='tame')
    subpss = parser.add_subparsers(dest='command', title='subcommands', required=True)
    subps = subpss.add_parser('cond', help='Conductivity Analysess')

    subp = subps.add_subparsers(title='methods', dest='command', required=True)
    p = subp.add_parser('jacf', help='Green-Kubo method '\
                    'with current density autocorrelation function')
    import_module('tame.recipes.cond_jacf').set_parser(p)
    p = subp.add_parser('pmsd', help='Green-Kubo method '\
                            'with polarization mean square displacement')
    import_module('tame.recipes.cond_pmsd').set_parser(p)


    subps = subpss.add_parser('diff', help='Diffusion Coefficients')
    subp = subps.add_subparsers(title='methods', dest='command', required=True)
    p = subp.add_parser('dcf', help='with displacement correlations')
    import_module('tame.recipes.diff_msd').set_parser(p)

    args = parser.parse_args()
    args.func(args);

if __name__ == "__main__":
    main()

