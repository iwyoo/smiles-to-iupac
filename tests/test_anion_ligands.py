import pytest

from smiles_to_iupac import smiles_to_iupac

K = "κ"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("CC(=O)[O][Cu][O]C(C)=O", "diacetatocopper"),
        ("O=C([O][Zn][O]C(=O)c1ccccc1)c1ccccc1", "dibenzoatozinc"),
        ("O=C[O][Zn][O]C=O", "diformatozinc"),
        ("O=[N+]([O-])[O][Cu][O][N+](=O)[O-]", "dinitratocopper"),
        ("N#C[S][Hg][S]C#N", f"bis(thiocyanato-{K}S)mercury"),
        ("S=C=[N][Hg][N]=C=S", f"bis(thiocyanato-{K}N)mercury"),
        ("CC[S][Hg][S]CC", "diethanethiolatomercury"),
        ("C[N](C)[Ti]([N](C)C)([N](C)C)[N](C)C", "tetrakis(dimethylazanido)titanium"),
        ("C[P](C)[Zr][P](C)C", "bis(dimethylphosphanido)zirconium"),
        ("CC1=CC(C)=[O]->[Cu]<-[O-]1", f"(pentane-2,4-dionato-{K}2O,O')cuprate(1-)"),
        ("O=C1[O-]->[Fe]<-[O-]C1=O", f"(oxalato-{K}2O,O')ferrate(2-)"),
        ("O=C1[O][Cu][O]1", f"(carbonato-{K}2O,O')copper"),
    ],
)
def test_anionic_ligand_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
