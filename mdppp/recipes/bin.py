import argparse

def main():
    from mdppp.recipes.cond_jacf import set_parser
    parser = argparse.ArgumentParser(prog='mdppp')
    subps = parser.add_subparsers(dest='command', title='subcommands', required=True)
    subp = subps.add_parser('cond', help='Conductivity Analysess')
    subps = subp.add_subparsers(title='method', dest='command', required=True)
    subp = subps.add_parser('jacf', help='Green-Kubo method with J-J autocorrelation function')
    set_parser(subp)
    args = parser.parse_args()
    args.func(args);

if __name__ == "__main__":
    main()

