from smiles_to_iupac import smiles_to_iupac


def test_pregnane():
    # Pregnane's PIN per the 1989 IUPAC steroid nomenclature (Rule 3S-2.4,
    # https://iupac.qmul.ac.uk/steroid/3S02a.html, Table 1): methyl groups
    # at both C-10 and C-13 (like androstane) plus a plain ethyl side chain
    # at C-17. Structure (C21H36) cross-checked against PubChem CID
    # 439513's ConnectivitySMILES, which canonicalizes identically; also
    # independently confirmed that removing that side chain's two carbons
    # reproduces `_androstane.py`'s own canonical SMILES exactly (see
    # module docstring).
    assert smiles_to_iupac("CCC1CCC2C1(CCC3C2CCC4C3(CCCC4)C)C") == "pregnane"


def test_stereo_specified_pregnane_falls_through_to_von_baeyer():
    # Pregnane's own definition (Rule 3S-2.4, first row of Table 1) doesn't
    # specify ring-fusion stereochemistry (side-chain C20 isn't itself a
    # stereocenter here, unlike the longer side chains further down Table
    # 1), so this module deliberately leaves stereo-marked input unmatched,
    # falling through to the general von Baeyer engine instead -- this
    # exact SMILES is PubChem CID 439513's own IsomericSMILES.
    assert (
        smiles_to_iupac("CC[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC[C@H]4[C@@]3(CCCC4)C)C")
        == "14-ethyl-2,15-dimethyltetracyclo[8.7.0.0^2,7.0^11,15]heptadecane"
    )


def test_gonane_androstane_and_estrane_are_not_pregnane():
    # unsubstituted gonane/androstane/estrane must not be misrecognized as
    # pregnane (all four modules' canonical keys are structurally distinct).
    assert smiles_to_iupac("C1CCCC2CCC3C(C12)CCC4C3CCC4") == "gonane"
    assert smiles_to_iupac("CC12CCCC1C3CCC4CCCCC4(C3CC2)C") == "androstane"
    assert smiles_to_iupac("CC12CCCC1C1CCC3CCCCC3C1CC2") == "estrane"


def test_extra_methyl_on_pregnane_falls_through_to_von_baeyer():
    # a third methyl beyond the two angular methyls (and the C17 ethyl
    # side chain) changes the canonical SMILES, so this module correctly
    # leaves it unmatched, falling through to the general von Baeyer
    # engine as a plain substituted hydrocarbon instead.
    assert (
        smiles_to_iupac("CCC1CC(C)C2C1(CCC3C2CCC4C3(CCCC4)C)C")
        == "14-ethyl-2,12,15-trimethyltetracyclo[8.7.0.0^2,7.0^11,15]heptadecane"
    )
