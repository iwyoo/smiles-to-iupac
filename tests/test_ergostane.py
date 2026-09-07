from smiles_to_iupac import smiles_to_iupac


def test_ergostane():
    # Ergostane's PIN per the 1989 IUPAC steroid nomenclature (Rule
    # 3S-2.4, https://iupac.qmul.ac.uk/steroid/3S02a.html, Table 1): methyl
    # groups at both C-10 and C-13 (like the rest of this family) plus a
    # nine-carbon C-17 side chain (C20-C28) -- cholestane's own eight-carbon
    # side chain with one extra methyl branch at C24. Structure (C28H50)
    # cross-checked against PubChem CID 6857535's ConnectivitySMILES, which
    # canonicalizes identically; also independently confirmed that
    # removing that side chain's nine carbons reproduces `_androstane.py`'s
    # own canonical SMILES exactly (see module docstring).
    assert smiles_to_iupac("CC(C)C(C)CCC(C)C1CCC2C1(CCC3C2CCC4C3(CCCC4)C)C") == "ergostane"


def test_stereo_specified_ergostane_falls_through_to_von_baeyer():
    # C20 and C24 are both genuine stereocenters, so this module
    # deliberately leaves stereo-marked input unmatched, falling through to
    # the general von Baeyer engine instead -- this exact SMILES is
    # PubChem CID 6857535's own IsomericSMILES.
    assert (
        smiles_to_iupac("C[C@H](CC[C@H](C)C(C)C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CCC4[C@@]3(CCCC4)C)C")
        == "14-(5,6-dimethylheptan-2-yl)-2,15-dimethyltetracyclo[8.7.0.0^2,7.0^11,15]heptadecane"
    )


def test_earlier_steroid_parents_are_not_ergostane():
    # unsubstituted gonane/androstane/estrane/pregnane/cholane/cholestane
    # must not be misrecognized as ergostane (all seven modules' canonical
    # keys are structurally distinct).
    assert smiles_to_iupac("C1CCCC2CCC3C(C12)CCC4C3CCC4") == "gonane"
    assert smiles_to_iupac("CC12CCCC1C3CCC4CCCCC4(C3CC2)C") == "androstane"
    assert smiles_to_iupac("CC12CCCC1C1CCC3CCCCC3C1CC2") == "estrane"
    assert smiles_to_iupac("CCC1CCC2C1(CCC3C2CCC4C3(CCCC4)C)C") == "pregnane"
    assert smiles_to_iupac("CCCC(C)C1CCC2C1(CCC3C2CCC4C3(CCCC4)C)C") == "cholane"
    assert smiles_to_iupac("CC(C)CCCC(C)C1CCC2C1(CCC3C2CCC4C3(CCCC4)C)C") == "cholestane"


def test_extra_methyl_on_ergostane_falls_through_to_von_baeyer():
    # a further methyl beyond the two angular methyls (and the C17 side
    # chain's own methyls) changes the canonical SMILES, so this module
    # correctly leaves it unmatched, falling through to the general von
    # Baeyer engine as a plain substituted hydrocarbon instead.
    assert (
        smiles_to_iupac("CC(C)C(C)CC(C)C1CCC2C1(CCC3C2CCC4C3(CCCC4)C)C")
        == "14-(4,5-dimethylhexan-2-yl)-2,15-dimethyltetracyclo[8.7.0.0^2,7.0^11,15]heptadecane"
    )
