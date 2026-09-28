from smiles_to_iupac import smiles_to_iupac


def test_morphine():
    # PubChem CID 5288826's IsomericSMILES. Real PIN (Blue Book P-101.8,
    # `tmp/bluebook/P10.txt` lines 1512-1516): '(5betaH)-17-methyl-7,8-
    # didehydrofuro[2',3',4',5':4,12,13,5]morphinan-3,6alpha-diol' --
    # stereodescriptors omitted (see module docstring).
    smiles = "CN1CC[C@]23[C@@H]4[C@H]1CC5=C2C(=C(C=C5)O)O[C@H]3[C@H](C=C4)O"
    assert smiles_to_iupac(smiles) == "17-methyl-7,8-didehydrofuro[2′,3′,4′,5′:4,12,13,5]morphinan-3,6-diol"


def test_codeine():
    # PubChem CID 5284371's IsomericSMILES -- morphine's 3-O-methyl ether.
    smiles = "CN1CC[C@]23[C@@H]4[C@H]1CC5=C2C(=C(C=C5)OC)O[C@H]3[C@H](C=C4)O"
    assert smiles_to_iupac(smiles) == "3-methoxy-17-methyl-7,8-didehydrofuro[2′,3′,4′,5′:4,12,13,5]morphinan-6-ol"
