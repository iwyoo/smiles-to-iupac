from smiles_to_iupac import smiles_to_iupac


def test_all_seven_parent_hydrides():
    # Each canonical SMILES cross-checked against its PubChem CID's
    # ConnectivitySMILES (see module docstring for the CIDs and the
    # side-chain/methyl removal chain confirming each skeleton reduces to
    # the previous one).
    assert smiles_to_iupac("C1CCCC2CCC3C(C12)CCC4C3CCC4") == "gonane"
    assert smiles_to_iupac("C1CCC2C(C1)CCC3C2CCC4C3CCC4") == "gonane"
    assert smiles_to_iupac("CC12CCCC1C3CCC4CCCCC4(C3CC2)C") == "androstane"
    assert smiles_to_iupac("CC12CCCC1C1CCC3CCCCC3C1CC2") == "estrane"
    assert smiles_to_iupac("CCC1CCC2C1(CCC3C2CCC4C3(CCCC4)C)C") == "pregnane"
    assert smiles_to_iupac("CCCC(C)C1CCC2C1(CCC3C2CCC4C3(CCCC4)C)C") == "cholane"
    assert smiles_to_iupac("CC(C)CCCC(C)C1CCC2C1(CCC3C2CCC4C3(CCCC4)C)C") == "cholestane"
    assert smiles_to_iupac("CC(C)C(C)CCC(C)C1CCC2C1(CCC3C2CCC4C3(CCCC4)C)C") == "ergostane"


def test_stereo_specified_parent_hydrides_fall_through_to_von_baeyer():
    # None of Rule 2.1/3S-2.2/3S-2.3 specify ring-fusion stereochemistry,
    # and cholane/cholestane/ergostane's side-chain stereocenters (C20, and
    # C24 for ergostane, per Table 1) are likewise out of scope here, so a
    # stereo-specified input for any of the seven always falls through to
    # the general von Baeyer engine instead -- each SMILES below is its
    # PubChem CID's own IsomericSMILES (CID 6857523 gonane, 94144
    # androstane, 5460658 estrane, 439513 pregnane, 6857459 cholane, 6857534
    # cholestane, 6857535 ergostane).
    assert (
        smiles_to_iupac("C1CCC2[C@H](C1)CCC3C2CCC4C3CCC4")
        == "tetracyclo[8.7.0.0^2,7.0^11,15]heptadecane"
    )
    assert (
        smiles_to_iupac("C[C@@]12CCC[C@H]1[C@@H]3CC[C@H]4CCCC[C@@]4([C@H]3CC2)C")
        == "2,15-dimethyltetracyclo[8.7.0.0^2,7.0^11,15]heptadecane"
    )
    assert (
        smiles_to_iupac("C[C@@]12CCC[C@H]1[C@@H]3CCC4CCCC[C@@H]4[C@H]3CC2")
        == "15-methyltetracyclo[8.7.0.0^2,7.0^11,15]heptadecane"
    )
    assert (
        smiles_to_iupac("CC[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC[C@H]4[C@@]3(CCCC4)C)C")
        == "14-ethyl-2,15-dimethyltetracyclo[8.7.0.0^2,7.0^11,15]heptadecane"
    )
    assert (
        smiles_to_iupac("CCC[C@@H](C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC[C@H]4[C@@]3(CCCC4)C)C")
        == "2,15-dimethyl-14-(pentan-2-yl)tetracyclo[8.7.0.0^2,7.0^11,15]heptadecane"
    )
    assert (
        smiles_to_iupac("C[C@H](CCCC(C)C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CCC4[C@@]3(CCCC4)C)C")
        == "2,15-dimethyl-14-(6-methylheptan-2-yl)tetracyclo[8.7.0.0^2,7.0^11,15]heptadecane"
    )
    assert (
        smiles_to_iupac("C[C@H](CC[C@H](C)C(C)C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CCC4[C@@]3(CCCC4)C)C")
        == "14-(5,6-dimethylheptan-2-yl)-2,15-dimethyltetracyclo[8.7.0.0^2,7.0^11,15]heptadecane"
    )


def test_extra_methyl_beyond_each_parent_hydride_falls_through_to_von_baeyer():
    # A methyl (or, for gonane, both angular methyls) beyond what each
    # skeleton's own definition allows changes the canonical SMILES, so it
    # correctly falls through to the general von Baeyer engine as a plain
    # substituted hydrocarbon instead of being misrecognized as the next
    # skeleton up.
    assert (
        smiles_to_iupac("CC12CCCC1(C)CCC1C2CCC2C1CCCC2")
        == "11,15-dimethyltetracyclo[8.7.0.0^2,7.0^11,15]heptadecane"
    )
    assert (
        smiles_to_iupac("CC12CCCC1C3CC(C)C4CCCCC4(C3CC2)C")
        == "2,8,15-trimethyltetracyclo[8.7.0.0^2,7.0^11,15]heptadecane"
    )
    assert (
        smiles_to_iupac("CC12CCCC1C1CC(C)C3CCCCC3C1CC2")
        == "8,15-dimethyltetracyclo[8.7.0.0^2,7.0^11,15]heptadecane"
    )
    assert (
        smiles_to_iupac("CCC1CC(C)C2C1(CCC3C2CCC4C3(CCCC4)C)C")
        == "14-ethyl-2,12,15-trimethyltetracyclo[8.7.0.0^2,7.0^11,15]heptadecane"
    )
    assert (
        smiles_to_iupac("CCCC(C)C1CC(C)C2C1(CCC3C2CCC4C3(CCCC4)C)C")
        == "2,12,15-trimethyl-14-(pentan-2-yl)tetracyclo[8.7.0.0^2,7.0^11,15]heptadecane"
    )
    assert (
        smiles_to_iupac("CC(C)CCCC(C)C1CC(C)C2C1(CCC3C2CCC4C3(CCCC4)C)C")
        == "2,12,15-trimethyl-14-(6-methylheptan-2-yl)tetracyclo[8.7.0.0^2,7.0^11,15]heptadecane"
    )
    assert (
        smiles_to_iupac("CC(C)C(C)CC(C)C1CCC2C1(CCC3C2CCC4C3(CCCC4)C)C")
        == "14-(4,5-dimethylhexan-2-yl)-2,15-dimethyltetracyclo[8.7.0.0^2,7.0^11,15]heptadecane"
    )
