from smiles_to_iupac import smiles_to_iupac


def test_carbonic_acid():
    # PubChem CID 767.
    assert smiles_to_iupac("OC(=O)O") == "carbonic acid"


def test_dimethyl_carbonate():
    # PubChem CID 12021.
    assert smiles_to_iupac("COC(=O)OC") == "dimethyl carbonate"


def test_ethyl_methyl_carbonate():
    # PubChem CID 522046 -- confirms alphanumeric word ordering.
    assert smiles_to_iupac("CCOC(=O)OC") == "ethyl methyl carbonate"


def test_methyl_hydrogen_carbonate():
    # Partial ester, PubChem CID 78579.
    assert smiles_to_iupac("COC(=O)O") == "methyl hydrogen carbonate"


def test_cyclic_carbonate_is_not_this_modules_territory():
    # Ethylene carbonate -- _ketone.py's own ring-ketone naming, not this
    # module's acyclic ester citation style.
    assert smiles_to_iupac("O=C1OCCO1") == "1,3-dioxolan-2-one"
