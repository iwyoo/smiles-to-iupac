from smiles_to_iupac import smiles_to_iupac


def test_benzenide_name():
    # Blue Book P-72.2.2.1's own worked example: the phenyl anion names as
    # 'benzenide (PIN)', with no locant citation (every ring position is
    # equivalent by symmetry).
    assert smiles_to_iupac("[c-]1ccccc1") == "benzenide"


def test_benzene_neutral_parent_unaffected():
    assert smiles_to_iupac("c1ccccc1") == "benzene"
