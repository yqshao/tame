"""Main entry point of mdppp"""

import argparse

def main():
    from importlib import import_module
    parser = argparse.ArgumentParser(prog='mdppp')
    subps = parser.add_subparsers(dest='command', title='subcommands')
    subp = subps.add_parser('cond', help='Conductivity Analysess')
    
    subps = subp.add_subparsers(title='methods', dest='command')
    subp = subps.add_parser('jacf', help='Green-Kubo method '\
                            'with current density autocorrelation function')
    import_module('mdppp.recipes.cond_jacf').set_parser(subp)
    subp = subps.add_parser('pmsd', help='Green-Kubo method '\
                            'with polarization mean square displacement')
    import_module('mdppp.recipes.cond_pmsd').set_parser(subp)

    args = parser.parse_args()
    args.func(args);

if __name__ == "__main__":
    main()

