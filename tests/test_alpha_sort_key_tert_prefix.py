from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._substituents import alpha_sort_key


def test_tert_butyl_sorts_under_b_not_t():
    assert sorted(["tert-butyl", "methyl"], key=alpha_sort_key) == ["tert-butyl", "methyl"]


def test_tert_butyl_cited_before_methyl_on_ring():
    # PubChem reference: 3-tert-butyl-3-methylcyclohexan-1-one.
    assert smiles_to_iupac("CC(C)(C)C1(C)CCCC(=O)C1") == "3-tert-butyl-3-methylcyclohexan-1-one"


def test_tert_butyl_wins_lowest_locant_tie_break():
    # Symmetric ring positions (1,4 either way): P-14.5.2 gives the lower
    # locant to the alphabetically-first substituent, 'tert-butyl' ("b").
    assert smiles_to_iupac("CC(C)(C)C1CCC(C)CC1") == "1-tert-butyl-4-methylcyclohexane"
