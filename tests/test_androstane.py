from smiles_to_iupac import smiles_to_iupac


def test_androstane():
    # Androstane's PIN per the 1989 IUPAC steroid nomenclature (Rule 3S-2.3,
    # https://iupac.qmul.ac.uk/steroid/3S02a.html): "the hydrocarbon with
    # methyl groups at C-10 and C-13 but without a side chain at C-17".
    # Structure (C19H32) cross-checked against PubChem CID 6857536's
    # ConnectivitySMILES, which canonicalizes identically; also
    # independently confirmed that removing this SMILES's two methyl
    # carbons reproduces `_gonane.py`'s own canonical gonane skeleton
    # exactly (see module docstring), proving the methyls sit at the real
    # C10/C13 angular positions, not some other pair of ring-fusion atoms.
    assert smiles_to_iupac("CC12CCCC1C3CCC4CCCCC4(C3CC2)C") == "androstane"


def test_stereo_specified_androstane_falls_through_to_von_baeyer():
    # Androstane's own definition (Rule 3S-2.3) doesn't specify ring-fusion
    # or angular-methyl stereochemistry, so this module deliberately leaves
    # stereo-marked input unmatched, falling through to the general von
    # Baeyer engine instead -- this exact SMILES is real 5alpha-androstane
    # (PubChem CID 94144's own IsomericSMILES).
    assert (
        smiles_to_iupac("C[C@@]12CCC[C@H]1[C@@H]3CC[C@H]4CCCC[C@@]4([C@H]3CC2)C")
        == "2,15-dimethyltetracyclo[8.7.0.0^2,7.0^11,15]heptadecane"
    )


def test_gonane_is_not_androstane():
    # unsubstituted gonane itself must not be misrecognized as androstane
    # (the two modules' canonical keys are structurally distinct).
    assert smiles_to_iupac("C1CCCC2CCC3C(C12)CCC4C3CCC4") == "gonane"


def test_extra_methyl_on_androstane_falls_through_to_von_baeyer():
    # a third methyl beyond the two angular ones changes the canonical
    # SMILES, so this module correctly leaves it unmatched, falling
    # through to the general von Baeyer engine as a plain substituted
    # hydrocarbon instead.
    assert (
        smiles_to_iupac("CC12CCCC1C3CC(C)C4CCCCC4(C3CC2)C")
        == "2,8,15-trimethyltetracyclo[8.7.0.0^2,7.0^11,15]heptadecane"
    )
