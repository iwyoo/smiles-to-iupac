from smiles_to_iupac import smiles_to_iupac


def test_cholestane():
    # Cholestane's PIN per the 1989 IUPAC steroid nomenclature (Rule
    # 3S-2.4, https://iupac.qmul.ac.uk/steroid/3S02a.html, Table 1): methyl
    # groups at both C-10 and C-13 (like the rest of this family) plus an
    # eight-carbon C-17 side chain (C20-C27) with a methyl branch at C20
    # and a symmetric gem-dimethyl terminus at C25. Structure (C27H48)
    # cross-checked against PubChem CID 6857534's ConnectivitySMILES, which
    # canonicalizes identically; also independently confirmed that
    # removing that side chain's eight carbons reproduces `_androstane.py`'s
    # own canonical SMILES exactly (see module docstring).
    assert smiles_to_iupac("CC(C)CCCC(C)C1CCC2C1(CCC3C2CCC4C3(CCCC4)C)C") == "cholestane"


def test_stereo_specified_cholestane_falls_through_to_von_baeyer():
    # Like cholane's C20, cholestane's C20 (bearing a methyl branch) is a
    # genuine stereocenter, so this module deliberately leaves stereo-marked
    # input unmatched, falling through to the general von Baeyer engine
    # instead -- this exact SMILES is PubChem CID 6857534's own
    # IsomericSMILES.
    assert (
        smiles_to_iupac("C[C@H](CCCC(C)C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CCC4[C@@]3(CCCC4)C)C")
        == "2,15-dimethyl-14-(6-methylheptan-2-yl)tetracyclo[8.7.0.0^2,7.0^11,15]heptadecane"
    )


def test_earlier_steroid_parents_are_not_cholestane():
    # unsubstituted gonane/androstane/estrane/pregnane/cholane must not be
    # misrecognized as cholestane (all six modules' canonical keys are
    # structurally distinct).
    assert smiles_to_iupac("C1CCCC2CCC3C(C12)CCC4C3CCC4") == "gonane"
    assert smiles_to_iupac("CC12CCCC1C3CCC4CCCCC4(C3CC2)C") == "androstane"
    assert smiles_to_iupac("CC12CCCC1C1CCC3CCCCC3C1CC2") == "estrane"
    assert smiles_to_iupac("CCC1CCC2C1(CCC3C2CCC4C3(CCCC4)C)C") == "pregnane"
    assert smiles_to_iupac("CCCC(C)C1CCC2C1(CCC3C2CCC4C3(CCCC4)C)C") == "cholane"


def test_extra_methyl_on_cholestane_falls_through_to_von_baeyer():
    # a third methyl beyond the two angular methyls (and the C17 side
    # chain's own methyls) changes the canonical SMILES, so this module
    # correctly leaves it unmatched, falling through to the general von
    # Baeyer engine as a plain substituted hydrocarbon instead.
    assert (
        smiles_to_iupac("CC(C)CCCC(C)C1CC(C)C2C1(CCC3C2CCC4C3(CCCC4)C)C")
        == "2,12,15-trimethyl-14-(6-methylheptan-2-yl)tetracyclo[8.7.0.0^2,7.0^11,15]heptadecane"
    )
