from smiles_to_iupac import smiles_to_iupac


def test_estrane():
    # Estrane's PIN per the 1989 IUPAC steroid nomenclature (Rule 3S-2.2,
    # https://iupac.qmul.ac.uk/steroid/3S02a.html): "the hydrocarbon with a
    # methyl group at C-13 but without a methyl group at C-10 and without a
    # side chain at C-17". Structure (C18H30) cross-checked against PubChem
    # CID 5460658's IsomericSMILES with stereo markers stripped, which
    # canonicalizes identically; also independently confirmed that removing
    # `_androstane.py`'s C10 methyl (the one at the six-six-ring-fusion
    # atom, not C13's five-six-ring-fusion atom) reproduces this exact
    # canonical SMILES (see module docstring).
    assert smiles_to_iupac("CC12CCCC1C1CCC3CCCCC3C1CC2") == "estrane"


def test_stereo_specified_estrane_falls_through_to_von_baeyer():
    # Estrane's own definition (Rule 3S-2.2) doesn't specify ring-fusion or
    # angular-methyl stereochemistry, so this module deliberately leaves
    # stereo-marked input unmatched, falling through to the general von
    # Baeyer engine instead -- this exact SMILES is PubChem CID 5460658's
    # own IsomericSMILES.
    assert (
        smiles_to_iupac("C[C@@]12CCC[C@H]1[C@@H]3CCC4CCCC[C@@H]4[C@H]3CC2")
        == "15-methyltetracyclo[8.7.0.0^2,7.0^11,15]heptadecane"
    )


def test_gonane_and_androstane_are_not_estrane():
    # unsubstituted gonane and androstane must not be misrecognized as
    # estrane (all three modules' canonical keys are structurally distinct).
    assert smiles_to_iupac("C1CCCC2CCC3C(C12)CCC4C3CCC4") == "gonane"
    assert smiles_to_iupac("CC12CCCC1C3CCC4CCCCC4(C3CC2)C") == "androstane"


def test_extra_methyl_on_estrane_falls_through_to_von_baeyer():
    # a second methyl beyond the one C13 angular methyl changes the
    # canonical SMILES, so this module correctly leaves it unmatched,
    # falling through to the general von Baeyer engine as a plain
    # substituted hydrocarbon instead.
    assert (
        smiles_to_iupac("CC12CCCC1C1CC(C)C3CCCCC3C1CC2")
        == "8,15-dimethyltetracyclo[8.7.0.0^2,7.0^11,15]heptadecane"
    )
