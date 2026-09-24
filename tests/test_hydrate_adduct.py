from smiles_to_iupac import smiles_to_iupac


def test_dihydrate_adduct():
    # Blue Book P-14.8: the proportion (1/n) is always cited.
    assert smiles_to_iupac("OC(=O)C(=O)O.O.O") == "ethanedioic acid—water (1/2)"


def test_monohydrate_adduct():
    # (1/1) is cited even for a single water, never omitted.
    assert smiles_to_iupac("CC(=O)O.O") == "ethanoic acid—water (1/1)"
