import argparse
import sys

from chemonym import __version__, smiles_to_iupac


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="chemonym",
        description=(
            "Convert SMILES strings to IUPAC names. "
            "Prints just the name if one SMILES is given, "
            "or '<smiles>\\t<name>' per line if multiple are given."
        ),
    )
    parser.add_argument("smiles", nargs="+", help="one or more SMILES strings")
    parser.add_argument(
        "--version",
        action="version",
        version=__version__,
    )
    args = parser.parse_args(argv)

    single = len(args.smiles) == 1
    had_error = False
    for smiles in args.smiles:
        try:
            name = smiles_to_iupac(smiles)
        except (ValueError, NotImplementedError) as e:
            had_error = True
            print(f"{smiles}: {e}", file=sys.stderr)
            continue
        print(name if single else f"{smiles}\t{name}")

    return 1 if had_error else 0


if __name__ == "__main__":
    sys.exit(main())
