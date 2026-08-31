import pytest

from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._common import UnsupportedStructure


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # oxa series (P-22.2.1 Table 22.1), formulas cross-checked:
        # oxirane C2H4O, oxetane C3H6O, oxolane C4H8O, oxane C5H10O,
        # oxepane C6H12O.
        ("C1CO1", "oxirane"),
        ("C1CCO1", "oxetane"),
        ("C1CCCO1", "oxolane"),
        ("C1CCCCO1", "oxane"),
        ("C1CCCCCO1", "oxepane"),
        # thia series, formulas cross-checked: thiirane C2H4S, thietane
        # C3H6S, thiolane C4H8S, thiane C5H10S, thiepane C6H12S.
        ("C1CS1", "thiirane"),
        ("C1CCS1", "thietane"),
        ("C1CCCS1", "thiolane"),
        ("C1CCCCS1", "thiane"),
        ("C1CCCCCS1", "thiepane"),
        # aza series: aziridine/azetidine coincide with the literal
        # Hantzsch-Widman stem construction; pyrrolidine/piperidine are
        # retained names that are themselves the PIN (not
        # azolidine/azinane) per P-22.2.1 Table 2.3; azepane is the stem
        # form again. Formulas cross-checked: aziridine C2H5N, azetidine
        # C3H7N, pyrrolidine C4H9N, piperidine C5H11N, azepane C6H13N.
        ("C1CN1", "aziridine"),
        ("C1CCN1", "azetidine"),
        ("C1CCCN1", "pyrrolidine"),
        ("C1CCCCN1", "piperidine"),
        ("C1CCCCCN1", "azepane"),
        # Mancude (aromatic) single-heteroatom monocycles, P-22.2.1 Table
        # 2.2 -- all confirmed as PubChem's IUPACName for the exact
        # SMILES: furan (CID 8029), thiophene (CID 8030), selenophene
        # (CID 136130), tellurophene (CID 136131), pyridine (CID 1049).
        ("c1ccoc1", "furan"),
        ("c1ccsc1", "thiophene"),
        ("c1cc[se]c1", "selenophene"),
        ("c1cc[te]c1", "tellurophene"),
        ("c1ccncc1", "pyridine"),
        # Pyrrole's PIN carries the indicated-hydrogen prefix (PubChem
        # CID 8027: "1H-pyrrole", not bare "pyrrole").
        ("c1cc[nH]c1", "1H-pyrrole"),
        # Mancude two-heteroatom monocycles, P-22.2.1 Table 2.2 -- all
        # confirmed as PubChem's IUPACName for the exact SMILES: imidazole
        # (CID 795), pyrazole (CID 1048), 1,3-oxazole (CID 9255), 1,2-oxazole
        # / isoxazole (CID 9254), 1,3-thiazole (CID 9256), 1,2-thiazole /
        # isothiazole (CID 67515), pyridazine (CID 9259), pyrimidine
        # (CID 9260), pyrazine (CID 9261).
        ("c1cnc[nH]1", "1H-imidazole"),
        ("c1cc[nH]n1", "1H-pyrazole"),
        ("c1cocn1", "1,3-oxazole"),
        ("c1ccon1", "1,2-oxazole"),
        ("c1cscn1", "1,3-thiazole"),
        ("c1ccsn1", "1,2-thiazole"),
        ("c1ccnnc1", "pyridazine"),
        ("c1ccncn1", "pyrimidine"),
        ("c1cnccn1", "pyrazine"),
        # Se/Te analogues of the N+S mancude pair, Table 2.2 -- 1,3- and
        # 1,2-selenazole confirmed as PubChem's IUPACName for the exact
        # SMILES (CID 11686913, CID 13224788); 1,2-tellurazole likewise
        # (CID 102212476). 1,3-tellurazole has no PubChem record for the
        # exact SMILES (CID 0), confirmed from Table 2.2's text alone
        # instead, matching the same symmetric O->S->Se->Te pattern as the
        # other three.
        ("c1cnc[se]1", "1,3-selenazole"),
        ("c1ccn[se]1", "1,2-selenazole"),
        ("c1cnc[te]1", "1,3-tellurazole"),
        ("c1ccn[te]1", "1,2-tellurazole"),
    ],
)
def test_smiles_to_iupac_hetero_monocyclic(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_substituted_hetero_monocyclic_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC1CCCCO1")


def test_ring_size_outside_scope_raises():
    # 8-membered ring: out of scope (this module covers 3-7).
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1CCCCCCO1")


def test_two_heteroatoms_raises():
    # 1,4-dioxane: two ring heteroatoms, out of scope here.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1COCCO1")


def test_unsupported_heteroatom_element_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1CCCCP1")


def test_substituted_mancude_monocyclic_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("Cc1ccncc1")


def test_substituted_two_heteroatom_mancude_ring_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("Cc1cnc[nH]1")


def test_three_heteroatom_mancude_ring_raises():
    # 1,2,4-triazole: three ring heteroatoms, out of scope here.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("c1nc[nH]n1")
