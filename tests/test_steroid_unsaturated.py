from smiles_to_iupac import smiles_to_iupac


def test_androst_5_ene():
    # PubChem CID 21670825.
    assert smiles_to_iupac("C[C@@]12CCC[C@H]1[C@@H]3CC=C4CCCC[C@@]4([C@H]3CC2)C") == "androst-5-ene"


def test_androst_4_ene():
    # A double bond moved to the 4,5 locant pair via the same androstane
    # skeleton (constructed by substructure-locating locants 4/5 and
    # setting that ring bond to double).
    assert smiles_to_iupac("CC12CCCC1C1CCC3=CCCCC3(C)C1CC2") == "androst-4-ene"


def test_plain_androstane_unaffected():
    assert smiles_to_iupac("CC12CCCC1C3CCC4CCCCC4(C3CC2)C") == "androstane"


def test_ring_fusion_double_bond_still_falls_through():
    # A double bond between two ring-fusion locants (5,10 on the estrane
    # skeleton) needs the compound-locant convention (e.g.
    # 'estra-5(10)-ene'), which is separate, unimplemented machinery --
    # this must not be silently misnamed as a plain sequential '-ene'.
    result = smiles_to_iupac("CC12CCCC1C1CCC3=C(CCCC3)C1CC2")
    assert "estra" not in result
    assert "ene" in result
