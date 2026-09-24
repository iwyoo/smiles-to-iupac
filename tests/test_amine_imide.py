from smiles_to_iupac import smiles_to_iupac


def test_trimethylhydrazinium_ide_amine_imide():
    # Blue Book P-74.2.1.3 worked example: 1,2,2-trimethylhydrazin-2-ium-1-ide.
    assert smiles_to_iupac("C[N-][NH+](C)C") == "1,2,2-trimethylhydrazin-2-ium-1-ide"


def test_unsubstituted_amine_imide():
    assert smiles_to_iupac("[NH-][NH3+]") == "hydrazin-2-ium-1-ide"


def test_single_substituent_per_nitrogen_amine_imide():
    # A substituent on each nitrogen, confirming locant citation splits
    # correctly across the fixed anion=1/cation=2 assignment.
    assert smiles_to_iupac("C[N-][NH2+]C") == "1,2-dimethylhydrazin-2-ium-1-ide"
