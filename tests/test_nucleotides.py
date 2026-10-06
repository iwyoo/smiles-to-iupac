import pytest
from smiles_to_iupac import NonPreferredNameWarning, smiles_to_iupac


# P-106.1 retained names, the primed locant of the ester position, P-106.3.1 derivatives and P-106.3.5 P-thio
@pytest.mark.parametrize(
    "smiles,expected",
    [
        ('Nc1ncnc2c1ncn2[C@@H]1O[C@H](COP(=O)(O)O)[C@@H](O)[C@H]1O', '5′-adenylic acid'),
        ('Nc1nc(=O)c2ncn([C@@H]3O[C@H](COP(=O)(O)O)[C@@H](O)[C@H]3O)c2[nH]1', '5′-guanylic acid'),
        ('O=c1nc[nH]c2c1ncn2[C@@H]1O[C@H](COP(=O)(O)O)[C@@H](O)[C@H]1O', '5′-inosinic acid'),
        ('O=c1[nH]c(=O)c2ncn([C@@H]3O[C@H](CO)[C@@H](OP(=O)(O)O)[C@H]3O)c2[nH]1', '3′-xanthylic acid'),
        ('Nc1ccn([C@@H]2O[C@H](COP(=O)(O)O)[C@@H](O)[C@H]2O)c(=O)n1', '5′-cytidylic acid'),
        ('Cc1cn([C@H]2C[C@H](O)[C@@H](COP(=O)(O)O)O2)c(=O)[nH]c1=O', '5′-thymidylic acid'),
        ('O=c1ccn([C@@H]2O[C@H](COP(=O)(O)O)[C@@H](O)[C@H]2O)c(=O)[nH]1', '5′-uridylic acid'),
        ('Nc1ncnc2c1ncn2[C@@H]1O[C@H](CO)[C@@H](O)[C@H]1OP(=O)(O)O', '2′-adenylic acid'),
        ('Nc1ncnc2c1ncn2[C@@H]1O[C@H](CO)[C@@H](OP(=O)(O)O)[C@H]1O', '3′-adenylic acid'),
        ('Cc1cn([C@H]2C[C@H](OP(=O)(O)O)[C@@H](CO)O2)c(=O)[nH]c1=O', '3′-thymidylic acid'),
        ('Nc1ncnc2c1ncn2[C@H]1C[C@H](O)[C@@H](COP(=O)(O)O)O1', '2′-deoxy-5′-adenylic acid'),
        ('CCCNC(=O)Nc1ncnc2c1ncn2[C@@H]1O[C@H](COP(=O)(O)O)[C@@H](O)[C@H]1O', 'N6-(propylcarbamoyl)-5′-adenylic acid'),
        ('Nc1ncnc2c1ncn2[C@@H]1O[C@H](COP(=O)(O)O)[C@H]2O[C@H](/C=C/c3ccccc3)O[C@H]21', '2′,3′-O-[(1S,2E)-3-phenylprop-2-ene-1,1-diyl]-5′-adenylic acid'),
        ('Nc1ncnc2c1ncn2[C@@H]1O[C@H](COP(=O)(O)O)[C@@H](O)[C@H]1OC', '2′-O-methyl-5′-adenylic acid'),
        ('Nc1ncnc2c1ncn2[C@@H]1O[C@H](COP(O)(O)=S)[C@@H](O)[C@H]1O', 'P-thio-5′-adenylic acid'),
        ('Nc1ncnc2c1ncn2[C@@H]1O[C@H](COP(=O)(O)S)[C@@H](O)[C@H]1O', 'P-thio-5′-adenylic acid'),
    ],
)
def test_nucleotide_retained_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


# P-106.2 di-, tri- and polyphosphates, P-106.3.1 derivatives of them
@pytest.mark.parametrize(
    "smiles,expected",
    [
        ('Nc1ncnc2c1ncn2[C@@H]1O[C@H](COP(=O)(O)OP(=O)(O)O)[C@@H](O)[C@H]1O', 'adenosine 5′-(trihydrogen diphosphate)'),
        ('Nc1ncnc2c1ncn2[C@@H]1O[C@H](COP(=O)(O)OP(=O)(O)OP(=O)(O)O)[C@@H](O)[C@H]1O', 'adenosine 5′-(tetrahydrogen triphosphate)'),
        ('Nc1nc(=O)c2ncn([C@@H]3O[C@H](COP(=O)(O)OP(=O)(O)O)[C@@H](O)[C@H]3O)c2[nH]1', 'guanosine 5′-(trihydrogen diphosphate)'),
        ('Nc1nc(=O)c2ncn([C@@H]3O[C@H](COP(=O)(O)OP(=O)(O)OP(=O)(O)O)[C@@H](O)[C@H]3O)c2[nH]1', 'guanosine 5′-(tetrahydrogen triphosphate)'),
        ('O=c1ccn([C@@H]2O[C@H](COP(=O)(O)OP(=O)(O)O)[C@@H](O)[C@H]2O)c(=O)[nH]1', 'uridine 5′-(trihydrogen diphosphate)'),
        ('O=c1ccn([C@@H]2O[C@H](COP(=O)(O)OP(=O)(O)OP(=O)(O)O)[C@@H](O)[C@H]2O)c(=O)[nH]1', 'uridine 5′-(tetrahydrogen triphosphate)'),
        ('Nc1ccn([C@@H]2O[C@H](COP(=O)(O)OP(=O)(O)O)[C@@H](O)[C@H]2O)c(=O)n1', 'cytidine 5′-(trihydrogen diphosphate)'),
        ('Nc1ccn([C@@H]2O[C@H](COP(=O)(O)OP(=O)(O)OP(=O)(O)O)[C@@H](O)[C@H]2O)c(=O)n1', 'cytidine 5′-(tetrahydrogen triphosphate)'),
        ('Cc1cn([C@H]2C[C@H](O)[C@@H](COP(=O)(O)OP(=O)(O)O)O2)c(=O)[nH]c1=O', 'thymidine 5′-(trihydrogen diphosphate)'),
        ('Cc1cn([C@H]2C[C@H](O)[C@@H](COP(=O)(O)OP(=O)(O)OP(=O)(O)O)O2)c(=O)[nH]c1=O', 'thymidine 5′-(tetrahydrogen triphosphate)'),
        ('Nc1ncnc2c1ncn2[C@H]1C[C@H](O)[C@@H](COP(=O)(O)OP(=O)(O)O)O1', '2′-deoxyadenosine 5′-(trihydrogen diphosphate)'),
        ('Nc1ncnc2c1ncn2[C@H]1C[C@H](O)[C@@H](COP(=O)(O)OP(=O)(O)OP(=O)(O)O)O1', '2′-deoxyadenosine 5′-(tetrahydrogen triphosphate)'),
        ('Nc1ncnc2c1ncn2[C@@H]1O[C@H](COP(=O)(O)OP(=O)(O)OP(=O)(O)OP(=O)(O)O)[C@@H](O)[C@H]1O', 'adenosine 5′-(pentahydrogen tetraphosphate)'),
        ('O=c1[nH]c(=O)c2ncn([C@@H]3O[C@H](CO)[C@@H](OP(=O)(O)OP(=O)(O)O)[C@H]3O)c2[nH]1', 'xanthosine 3′-(trihydrogen diphosphate)'),
        ('Nc1ncnc2c1ncn2[C@@H]1O[C@H](COP(=O)(O)O)[C@@H](OP(=O)(O)O)[C@H]1O', 'adenosine 3′,5′-bis(dihydrogen phosphate)'),
        ('CC(=O)OC[C@H]1O[C@@H](n2cnc3c2nc(N)[nH]c3=O)C[C@@H]1OP(=O)(O)OP(=O)(O)O', '5′-O-acetyl-2′-deoxyguanosine 3′-(trihydrogen diphosphate)'),
    ],
)
def test_nucleoside_polyphosphates(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


# P-106.3.2 functional replacement in the diphosphate/triphosphate chain, P-67.2.2 numbering from the esterified phosphorus
@pytest.mark.parametrize(
    "smiles,expected",
    [
        ('Nc1ncnc2c1ncn2[C@@H]1O[C@H](COP(O)(=S)OP(=O)(O)O)[C@@H](O)[C@H]1O', 'adenosine 5′-(trihydrogen 1-thiodiphosphate)'),
        ('Nc1ncnc2c1ncn2[C@@H]1O[C@H](COP(=O)(O)OP(O)(O)=S)[C@@H](O)[C@H]1O', 'adenosine 5′-(trihydrogen 3-thiodiphosphate)'),
        ('Nc1ncnc2c1ncn2[C@@H]1O[C@H](COP(=O)(O)SP(=O)(O)O)[C@@H](O)[C@H]1O', 'adenosine 5′-(trihydrogen 2-thiodiphosphate)'),
        ('Nc1ncnc2c1ncn2[C@@H]1O[C@H](COP(=O)(O)OP(=O)(O)OP(O)(O)=S)[C@@H](O)[C@H]1O', 'adenosine 5′-(tetrahydrogen 5-thiotriphosphate)'),
        ('Nc1ncnc2c1ncn2[C@@H]1O[C@H](COP(=O)(O)OP(=O)(O)NP(=O)(O)O)[C@@H](O)[C@H]1O', 'adenosine 5′-(tetrahydrogen 4-imidotriphosphate)'),
        ('Nc1nc(=O)c2ncn([C@@H]3O[C@H](COP(=O)(O)OP(=O)(O)NP(=O)(O)O)[C@@H](O)[C@H]3O)c2[nH]1', 'guanosine 5′-(tetrahydrogen 4-imidotriphosphate)'),
    ],
)
def test_nucleoside_polyphosphate_analogues(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


# P-106.3.3 the nucleotide acyl group on a group senior to the phosphoric acid residue
@pytest.mark.parametrize(
    "smiles,expected",
    [
        ('C1=NC2=C(N1[C@H]3[C@@H]([C@@H]([C@H](O3)COP(=O)(O)Oc4cccc(C(=O)O)c4)O)O)NC(=NC2=O)N', '3-(5′-guanylyloxy)benzoic acid'),
        ('C1=NC(=C2C(=N1)N(C=N2)[C@H]3[C@@H]([C@@H]([C@H](O3)COP(=O)(O)OS(=O)(=O)O)O)O)N', '5′-adenylyl hydrogen sulfate'),
        ('C1=NC(=C2C(=N1)N(C=N2)[C@H]3[C@@H]([C@@H]([C@H](O3)COP(=O)(O)OS(=O)(=O)O)OP(=O)(O)O)O)N', '3′-O-phosphono-5′-adenylyl hydrogen sulfate'),
        ('C1=NC(=O)C2=C(N1)N(C=N2)[C@H]3[C@@H]([C@@H]([C@H](O3)COP(=O)(O)OCC(=O)O)O)O', '(5′-inosinylyloxy)acetic acid'),
        ('C1=NC2=C(N1[C@H]3C[C@@H]([C@H](O3)COP(=S)(O)Oc4cccc(C(=O)O)c4)O)NC(=NC2=O)N', '3-(2′-deoxy-P-thio-5′-guanylyloxy)benzoic acid'),
        ('Nc1ncnc2c1ncn2[C@@H]1O[C@H](COP(=O)(O)OCP(=O)(O)O)[C@@H](O)[C@H]1O', '[(5′-adenylyloxy)methyl]phosphonic acid'),
        ('Nc1ncnc2c1ncn2[C@@H]1O[C@H](COP(=O)(O)OCC(=S)O)[C@@H](O)[C@H]1O', '2-(5′-adenylyloxy)ethanethioic O-acid'),
    ],
)
def test_nucleotidyl_groups(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


# P-106.3.4 oligonucleotides read from the end that gives the lower link locants, P-106.3.5 P-thio links
@pytest.mark.parametrize(
    "smiles,expected",
    [
        ('Nc1nc2c(ncn2[C@H]2C[C@H](OP(=O)(O)OC[C@H]3O[C@@H](n4ccc(=O)[nH]c4=O)C[C@@H]3OP(=O)(O)OC[C@H]3O[C@@H](n4cnc5c(=O)[nH]c(N)nc54)C[C@@H]3O)[C@@H](CO)O2)c(=O)[nH]1', '2′-deoxyguanylyl-(3′→5′)-2′-deoxyuridylyl-(3′→5′)-2′-deoxyguanosine'),
        ('Nc1nc2c(ncn2[C@H]2C[C@H](OP(O)(=S)OC[C@H]3O[C@@H](n4ccc(=O)[nH]c4=O)C[C@@H]3OP(O)(=S)OC[C@H]3O[C@@H](n4cnc5c(=O)[nH]c(N)nc54)C[C@@H]3O)[C@@H](CO)O2)c(=O)[nH]1', '2′-deoxy-P-thioguanylyl-(3′→5′)-2′-deoxy-P-thiouridylyl-(3′→5′)-2′-deoxyguanosine'),
        ('Cc1cn([C@H]2C[C@H](OP(=O)(O)OC[C@H]3O[C@@H](n4cnc5c(N)ncnc54)C[C@@H]3O)[C@@H](COP(=O)(O)O[C@H]3[C@@H](O)[C@H](n4ccc(=O)[nH]c4=O)O[C@@H]3COP(=O)(O)O[C@H]3[C@@H](O)[C@H](n4cnc5c(N)ncnc54)O[C@@H]3CO)O2)c(=O)[nH]c1=O', 'adenylyl-(3′→5′)-uridylyl-(3′→5′)-thymidylyl-(3′→5′)-2′-deoxyadenosine'),
        ('Nc1ncnc2c1ncn2[C@@H]1O[C@H](COP(=O)(O)O[C@H]2[C@@H](O)[C@H](n3cnc4c(N)ncnc43)O[C@@H]2CO)[C@@H](O)[C@H]1O', 'adenylyl-(3′→5′)-adenosine'),
        ('Nc1ncnc2c1ncn2[C@@H]1O[C@H](COP(=O)(O)O[C@@H]2[C@H](O)[C@@H](CO)O[C@H]2n2cnc3c(N)ncnc32)[C@@H](O)[C@H]1O', 'adenylyl-(2′→5′)-adenosine'),
        ('Nc1nc2c(ncn2[C@H]2C[C@H](O)[C@@H](COP(=O)(O)O[C@H]3C[C@H](n4ccc(=O)[nH]c4=O)O[C@@H]3CO)O2)c(=O)[nH]1', '2′-deoxyuridylyl-(3′→5′)-2′-deoxyguanosine'),
        ('Cc1cn([C@H]2C[C@H](OP(=O)(O)O[C@H]3C[C@H](n4cnc5c(N)ncnc54)O[C@@H]3CO)[C@@H](CO)O2)c(=O)[nH]c1=O', '2′-deoxyadenylyl-(3′→3′)-thymidine'),
        ('Nc1ncnc2c1ncn2[C@@H]1O[C@H](COP(=O)(O)O[C@H]2[C@@H](O)[C@H](n3cnc4c(N)ncnc43)O[C@@H]2COP(=O)(O)O)[C@@H](O)[C@H]1O', 'adenylyl-(5′→3′)-5′-adenylic acid'),
        ('Nc1ncnc2c1ncn2[C@@H]1O[C@H](COP(=O)(O)O[C@H]2[C@@H](OC)[C@H](n3cnc4c(N)ncnc43)O[C@@H]2COP(=O)(O)O[C@H]2[C@@H](O)[C@H](n3cnc4c(N)ncnc43)O[C@@H]2CO)[C@@H](O)[C@H]1O', 'adenylyl-(3′→5′)-2′-O-methyladenylyl-(3′→5′)-adenosine'),
    ],
)
def test_oligonucleotides(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


# P-106.3.2: three alternative names are listed for the carba analogue and none is marked preferred;
# P-106 gives no name for cyclic nucleoside phosphates
@pytest.mark.parametrize(
    "smiles,expected",
    [
        ('Nc1nc(=O)c2ncn([C@@H]3O[C@H](COP(=O)(O)CP(=O)(O)O)[C@@H](O)[C@H]3O)c2[nH]1', 'guanosine 5′-(trihydrogen 2-carbadiphosphate)'),
        ('Nc1ncnc2c1ncn2[C@@H]1O[C@@H]2COP(=O)(O)O[C@H]2[C@H]1O', '3′,5′-dideoxyadenosine-3′,5′-diyl hydrogen phosphate'),
        ('Nc1ncnc2c1ncn2[C@@H]1O[C@H](CO)[C@H]2OP(=O)(O)O[C@H]21', '2′,3′-dideoxyadenosine-2′,3′-diyl hydrogen phosphate'),
        ('Nc1nc(=O)c2ncn([C@@H]3O[C@@H]4COP(=O)(O)O[C@H]4[C@H]3O)c2[nH]1', '3′,5′-dideoxyguanosine-3′,5′-diyl hydrogen phosphate'),
    ],
)
def test_nucleotide_names_without_a_preferred_form_warn(smiles, expected):
    with pytest.warns(NonPreferredNameWarning):
        assert smiles_to_iupac(smiles) == expected
