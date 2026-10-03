import pytest

from smiles_to_iupac import smiles_to_iupac


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("COCCOC", "1,2-dimethoxyethane"),
        ("COCCOCCOC", "1-methoxy-2-(2-methoxyethoxy)ethane"),
        ("COCCOCCOCCOC", "2,5,8,11-tetraoxadodecane"),
        ("COCCOCCOCCOCCOC", "2,5,8,11,14-pentaoxapentadecane"),
        ("CSCCSCCSCCSC", "2,5,8,11-tetrathiadodecane"),
        ("CNCCNCCNCCNC", "2,5,8,11-tetraazadodecane"),
        ("CSCCSC", "1,2-bis(methylsulfanyl)ethane"),
        ("NCCNCCN", "N1-(2-aminoethyl)ethane-1,2-diamine"),
        ("NCCN(C)CCN", "N1-(2-aminoethyl)-N1-methylethane-1,2-diamine"),
        ("NCCNCCNCCN", "N1-{2-[(2-aminoethyl)amino]ethyl}ethane-1,2-diamine"),
        ("C1COCCOCCOCCO1", "1,4,7,10-tetraoxacyclododecane"),
        ("C1CNCCNCCCNCCNC1", "1,4,8,11-tetraazacyclotetradecane"),
        ("C1CNCCNCCN1", "1,4,7-triazonane"),
        ("C1CNCCOCCN1", "1,4,7-oxadiazonane"),
    ],
)
def test_hetero_chain_and_macrocycle_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("[Cl][Pt]12<-[NH2]CC[NH]->1CC[NH2]->2", "chlorido[N1-(2-aminoethyl)ethane-1,2-diamine-κ3N,N',N'']platinum"),
        ("C[O]1CC[O](C)->[Pt]<-1([Cl])[Cl]", "dichlorido(1,2-dimethoxyethane-κ2O,O')platinum"),
        (
            "C1C[O]2->[K]3456<-[O]1CC[O]->3CC[O]->4CC[O]->5CC[O]->6CC2",
            "(1,4,7,10,13,16-hexaoxacyclooctadecane-κ6O,O',O'',O''',O'''',O''''')potassium",
        ),
        (
            "C1C[NH]2->[Cu]34<-[NH](C1)CC[NH]->3CCC[NH]->4CC2",
            "(1,4,8,11-tetraazacyclotetradecane-κ4N,N',N'',N''')copper",
        ),
    ],
)
def test_polydentate_ligand_complexes(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
