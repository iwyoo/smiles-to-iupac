import argparse
import sys

from smiles_to_iupac import __version__, smiles_to_alternative_names, smiles_to_iupac


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="smiles-to-iupac",
        description=(
            "Convert SMILES strings to IUPAC names. "
            "Prints just the name if one SMILES is given, "
            "or '<smiles>\\t<name>' per line if multiple are given."
        ),
    )
    parser.add_argument("smiles", nargs="+", help="one or more SMILES strings")
    parser.add_argument(
        "--alternatives",
        action="store_true",
        help="print the general-nomenclature alternatives to the preferred name as '<kind>\\t<name>' lines",
    )
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
            if args.alternatives:
                lines = [f"{alt.kind}\t{alt.name}" for alt in smiles_to_alternative_names(smiles)]
            else:
                lines = [smiles_to_iupac(smiles)]
        except (ValueError, NotImplementedError) as e:
            had_error = True
            print(f"{smiles}: {e}", file=sys.stderr)
            continue
        for line in lines:
            print(line if single else f"{smiles}\t{line}")

    return 1 if had_error else 0


if __name__ == "__main__":
    sys.exit(main())
