from smiles_to_iupac import smiles_to_iupac


def test_didehydrooxepane_matches_blue_book_worked_example():
    # tmp/bluebook/P1.txt line 422: oxepane (PIN) -> 2,3-didehydrooxepane,
    # "the prefix 'didehydro' indicates loss of 2 hydrogen atoms".
    assert smiles_to_iupac("C1CCCC=CO1") == "2,3-didehydrooxepane"


def test_didehydro_locants_pick_lower_of_the_two_ring_directions():
    # Same shape, double bond moved elsewhere on the ring -- the two walk
    # directions from O both land on the same locant pair here ({4,5}),
    # confirming the direction search picks the correct, lower one rather
    # than an arbitrary SMILES-atom-order-dependent pair.
    assert smiles_to_iupac("O1CCC=CCC1") == "4,5-didehydrooxepane"


def test_oxepane_itself_is_unaffected():
    assert smiles_to_iupac("C1CCCCCO1") == "oxepane"


def test_didehydropiperidine_nitrogen_ring():
    # Same mechanism on the saturated N-heterocycle table (piperidine);
    # the ring N keeps its own N-H, the double bond sits between two
    # carbons elsewhere on the ring.
    assert smiles_to_iupac("C1CC=CCN1") == "3,4-didehydropiperidine"
