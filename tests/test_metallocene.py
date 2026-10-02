import pytest

from smiles_to_iupac import smiles_to_iupac


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("C1=C[CH]C=C1.C1=C[CH]C=C1.[Fe]", "ferrocene"),  # neutral biradical, PubChem CID 10219726
        ("C1=C[CH]C=C1.C1=C[CH]C=C1.[Ru]", "ruthenocene"),  # neutral biradical, CID 11020720
        ("C1=C[CH]C=C1.C1=C[CH]C=C1.[Os]", "osmocene"),  # neutral biradical form, CID 102601604
        ("[CH-]1C=CC=C1.[CH-]1C=CC=C1.[Ni+2]", "nickelocene"),  # anionic-dication, CID 62390
        ("[CH-]1C=CC=C1.[CH-]1C=CC=C1.[Cr+2]", "chromocene"),  # anionic-dication, CID 79154
        ("[CH-]1C=CC=C1.[CH-]1C=CC=C1.[Co+2]", "cobaltocene"),  # anionic-dication, CID 92884
        ("[CH-]1C=CC=C1.[CH-]1C=CC=C1.[V+2]", "vanadocene"),  # anionic-dication, CID 71311455
        ("[CH-]1C=CC=C1.[CH-]1C=CC=C1.[Fe+2]", "ferrocene"),  # opposite charge convention also matches
        ("C1=C[CH]C=C1.C1=C[CH]C=C1.[Ni]", "nickelocene"),  # opposite charge convention also matches
    ],
)
def test_metallocene_retained_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_metallocene_wrong_metal_is_not_matched():
    # no retained "-ocene" name applies; the general P-69.2 name is used
    assert smiles_to_iupac("C1=C[CH]C=C1.C1=C[CH]C=C1.[Zn]") == "bis(\u03b75-cyclopenta-2,4-dien-1-yl)zinc"


def test_metallocene_three_rings_is_not_matched():
    # no retained "-ocene" name applies; the general P-69.2 name is used
    assert smiles_to_iupac("C1=C[CH]C=C1.C1=C[CH]C=C1.C1=C[CH]C=C1.[Fe]") == "tris(\u03b75-cyclopenta-2,4-dien-1-yl)iron"


def test_metallocene_single_ring_is_not_matched():
    # no retained "-ocene" name applies; the general P-69.2 name is used
    assert smiles_to_iupac("C1=C[CH]C=C1.[Fe]") == "(\u03b75-cyclopenta-2,4-dien-1-yl)iron"
