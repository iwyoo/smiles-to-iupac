from smiles_to_iupac import smiles_to_iupac


def test_morphine():
    # PubChem CID 5288826's IsomericSMILES. Real name (Blue Book P-101.8,
    # `tmp/bluebook/P1.txt` lines 556-557's bridge-prefix alternative):
    # '4,5alpha-epoxy-17-methyl-7,8-didehydromorphinan-3,6alpha-diol' --
    # stereodescriptors are out of scope here (see module docstring), so
    # the expected result omits them, matching this codebase's own
    # established policy for the sibling epoxy-bridge module.
    smiles = "CN1CC[C@]23[C@@H]4[C@H]1CC5=C2C(=C(C=C5)O)O[C@H]3[C@H](C=C4)O"
    assert smiles_to_iupac(smiles) == "4,5-epoxy-17-methyl-7,8-didehydromorphinan-3,6-diol"


def test_codeine():
    # PubChem CID 5284371's IsomericSMILES -- morphine's 3-O-methyl ether.
    smiles = "CN1CC[C@]23[C@@H]4[C@H]1CC5=C2C(=C(C=C5)OC)O[C@H]3[C@H](C=C4)O"
    assert smiles_to_iupac(smiles) == "4,5-epoxy-3-methoxy-17-methyl-7,8-didehydromorphinan-6-ol"
