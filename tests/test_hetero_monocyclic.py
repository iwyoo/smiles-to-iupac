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


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # One substituent on a mancude parent, fixed heteroatom locants
        # unchanged from the unsubstituted table above -- all confirmed as
        # PubChem's IUPACName for the exact SMILES. Symmetric parents
        # (single heteroatom, and the three 6-membered two-nitrogen rings)
        # pick the lower of two valid numbering directions automatically;
        # the heteroatom-asymmetric ones (thiazole/oxazole family) have
        # only one valid direction to begin with.
        ("Cc1ccoc1", "3-methylfuran"),  # CID 13587
        ("Cc1ccsc1", "3-methylthiophene"),  # CID 12024
        ("Cc1cc[se]c1", "3-methylselenophene"),  # CID 13022371
        ("Cc1cc[te]c1", "3-methyltellurophene"),  # CID 13022372
        ("Cc1ccccn1", "2-methylpyridine"),  # CID 7975
        ("Clc1ccncc1", "4-chloropyridine"),  # CID 12288
        ("Cc1ccnnc1", "4-methylpyridazine"),  # CID 136882
        ("Cc1ncccn1", "2-methylpyrimidine"),  # CID 78748
        ("Cc1cnccn1", "2-methylpyrazine"),  # CID 7976
        ("Cc1ccon1", "3-methyl-1,2-oxazole"),  # CID 96098
        ("Cc1cscn1", "4-methyl-1,3-thiazole"),  # CID 12748
        ("Cc1ccsn1", "3-methyl-1,2-thiazole"),  # CID 12747
        # Pyrrole/imidazole/pyrazole's own N-H position (locant 1) is a
        # real, unambiguous substitutable position -- substituting it
        # directly drops the parent's '1H-' indicated-hydrogen prefix
        # entirely, since the locant '1-' alone already pins it.
        ("Cn1cccc1", "1-methylpyrrole"),  # CID 7304
        ("Cn1ccnc1", "1-methylimidazole"),  # CID 1390
        ("Cn1cccn1", "1-methylpyrazole"),  # CID 70255
        # Two or more substituents, all confirmed as PubChem's IUPACName for the exact SMILES.
        # Symmetric parents pick the lowest locant *set* automatically;
        # mixed substituent kinds are cited alphabetically.
        ("Cc1ccc(C)o1", "2,5-dimethylfuran"),  # CID 12266
        ("Cc1ccc(Cl)s1", "2-chloro-5-methylthiophene"),  # CID 140208
        ("Clc1ccncc1Cl", "3,4-dichloropyridine"),  # CID 2736081
        ("Cc1nc(C)co1", "2,4-dimethyl-1,3-oxazole"),  # CID 138961
        ("Clc1nccnc1Cl", "2,3-dichloropyrazine"),  # CID 78575
        ("Cc1c(C)coc1", "3,4-dimethylfuran"),  # CID 34338
        # imidazole/pyrazole multi-substituent is only safe when one of
        # the substituents sits at the tautomer-fixing N-H locant (1) --
        # see `test_imidazole_multi_substituent_without_n1_raises` below
        # for the case that's still rejected.
        ("Cn1cc(Cl)nc1", "4-chloro-1-methylimidazole"),  # CID 12514200
        ("Cn1cc(Cl)cn1", "4-chloro-1-methylpyrazole"),  # CID 13844024
        # A pyrrole carbon substituent alongside the untouched N-H is
        # unambiguous (no tautomer axis involved, unlike imidazole/
        # pyrazole), so the '1H-' prefix is retained.
        ("Cc1cc(Cl)c[nH]1", "4-chloro-2-methyl-1H-pyrrole"),  # CID 57109453
    ],
)
def test_smiles_to_iupac_hetero_monocyclic_substituent(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_substituted_hetero_monocyclic_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("CC1CCCCO1")


def test_hetero_monocyclic_substituent_with_non_alkyl_branch_raises():
    # A branch containing a heteroatom (here, an amine nitrogen and a
    # carboxylic acid) must not be silently walked as if it were a plain
    # carbon chain and named as a fictitious alkyl substituent -- found
    # via real-data testing: this exact SMILES was misnamed
    # '3-chloro-2-(3,4-dimethylpentyl)pyridine' (PubChem PIN is
    # '2-amino-4-(3-chloro-2-pyridinyl)butanoic acid'), silently dropping
    # both the amino and carboxylic acid groups.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("NC(CCc1ncccc1Cl)C(=O)O")


def test_imidazole_multi_substituent_without_n1_raises():
    # Same tautomer ambiguity as the single-substituent case below, but
    # with two ring-carbon substituents and neither at the N-H-derived
    # locant 1: PubChem's own name for this exact SMILES ("2-chloro-
    # 5-methyl-1H-imidazole", CID 313195) doesn't match this input's own
    # literal, structurally-fixed N-H position (which would give locants
    # 2 and 4), confirming PubChem silently renormalizes to its own
    # canonical tautomer before naming.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("Cc1cnc(Cl)[nH]1")


def test_ring_size_outside_scope_raises():
    # 8-membered ring: out of scope (this module covers 3-7).
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1CCCCCCO1")


def test_two_heteroatoms_raises():
    # 1,4-related Se+Te pair: two ring heteroatoms, out of scope here --
    # unlike every other pair below (including the other Se/Te-containing
    # ones), this one specific combination has no name registered in
    # PubChem, so it stays unsupported.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("[Te]1CC[Se]CC1")


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # 1,4-related two-heteroatom saturated 6-membered rings with their
        # own retained/systematic name (P-22.2.1) -- PubChem-verified:
        # morpholine (CID 8083), piperazine (CID 4837), thiomorpholine
        # (CID 67164), 1,4-dioxane (CID 31275), 1,4-oxathiane (CID 27596),
        # 1,4-dithiane (CID 10452).
        ("C1COCCN1", "morpholine"),
        ("C1CNCCN1", "piperazine"),
        ("C1CSCCN1", "thiomorpholine"),
        ("C1COCCO1", "1,4-dioxane"),
        ("C1COCCS1", "1,4-oxathiane"),
        ("C1CSCCS1", "1,4-dithiane"),
    ],
)
def test_two_heteroatom_saturated_ring_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # 1,4-related two-heteroatom saturated 7-membered rings -- same
        # N/O/S axis as the 6-membered family above, one ring atom larger.
        # PubChem PUG REST auto-generated names match exactly for all six.
        ("C1CNCCNC1", "1,4-diazepane"),
        ("C1CNCCOC1", "1,4-oxazepane"),
        ("C1CNCCSC1", "1,4-thiazepane"),
        ("C1COCCOC1", "1,4-dioxepane"),
        ("C1COCCSC1", "1,4-oxathiepane"),
        ("C1CSCCSC1", "1,4-dithiepane"),
    ],
)
def test_seven_membered_two_heteroatom_saturated_ring_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # 1,3-related two-heteroatom saturated 7-membered rings -- same
        # N/O/S axis as the 1,4-seven-membered family above, heteroatoms
        # one carbon apart instead of two. PubChem PUG REST auto-generated
        # names match exactly for all six.
        ("N1CNCCCC1", "1,3-diazepane"),
        ("O1CNCCCC1", "1,3-oxazepane"),
        ("S1CNCCCC1", "1,3-thiazepane"),
        ("O1COCCCC1", "1,3-dioxepane"),
        ("C1CCCOCS1", "1,3-oxathiepane"),
        ("S1CSCCCC1", "1,3-dithiepane"),
    ],
)
def test_seven_membered_1_3_two_heteroatom_saturated_ring_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # 1,2-related two-heteroatom saturated 7-membered rings --
        # heteroatoms directly adjacent. PubChem's raw IUPACName drops the
        # locant for these six ('diazepane', 'oxazepane', ...), but that
        # bare stem clashes with the 1,3-axis above (P-22.2.2.1.7 only
        # omits locants when there's no ambiguity), so this module keeps
        # '1,2-' explicit -- see the module docstring.
        ("C1CCCCNN1", "1,2-diazepane"),
        ("C1CCCCON1", "1,2-oxazepane"),
        ("C1CCCCSN1", "1,2-thiazepane"),
        ("C1CCCCOO1", "1,2-dioxepane"),
        ("C1CCCCOS1", "1,2-oxathiepane"),
        ("C1CCCCSS1", "1,2-dithiepane"),
    ],
)
def test_seven_membered_1_2_two_heteroatom_saturated_ring_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # 1,4-related two-heteroatom saturated 6-membered rings whose pair
        # includes Se and/or Te -- PubChem-verified: selenomorpholine,
        # telluromorpholine, 1,4-oxaselenane, 1,4-oxatellurane,
        # 1,4-thiaselenane, 1,4-thiatellurane, 1,4-diselenane,
        # 1,4-ditellurane (Se+Te itself has no registered name -- see
        # `test_two_heteroatoms_raises` above).
        ("C1CNCC[Se]1", "selenomorpholine"),
        ("C1CNCC[Te]1", "telluromorpholine"),
        ("C1COCC[Se]1", "1,4-oxaselenane"),
        ("C1COCC[Te]1", "1,4-oxatellurane"),
        ("C1CSCC[Se]1", "1,4-thiaselenane"),
        ("C1CSCC[Te]1", "1,4-thiatellurane"),
        ("C1C[Se]CC[Se]1", "1,4-diselenane"),
        ("C1C[Te]CC[Te]1", "1,4-ditellurane"),
    ],
)
def test_two_heteroatom_saturated_ring_names_se_te(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # 1,3-related two-heteroatom saturated 5-membered rings with their
        # own retained/systematic name (P-22.2.1) -- PubChem-verified:
        # imidazolidine (CID 449488), 1,3-oxazolidine (CID 536683),
        # 1,3-thiazolidine (CID 10444), 1,3-dioxolane (CID 12586),
        # 1,3-oxathiolane (CID 65092), 1,3-dithiolane (CID 20970).
        ("C1CNCN1", "imidazolidine"),
        ("C1CNCO1", "1,3-oxazolidine"),
        ("C1CSCN1", "1,3-thiazolidine"),
        ("C1COCO1", "1,3-dioxolane"),
        ("C1CSCO1", "1,3-oxathiolane"),
        ("C1CSCS1", "1,3-dithiolane"),
    ],
)
def test_five_membered_1_3_two_heteroatom_ring_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_five_membered_1_3_two_heteroatoms_se_te_pair_raises():
    # 1,3-related Se+Te pair: out of scope, same reason as the 1,4-ring's
    # Se+Te pair (no name registered in PubChem).
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1C[Se]C[Te]1")


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # 1,3-related two-heteroatom saturated 5-membered rings whose pair
        # includes Se and/or Te -- PubChem-verified: 1,3-selenazolidine,
        # 1,3-tellurazolidine, 1,3-oxaselenolane, 1,3-oxatellurolane,
        # 1,3-thiaselenolane, 1,3-thiatellurolane, 1,3-diselenolane,
        # 1,3-ditellurolane (Se+Te itself has no registered name -- see
        # `test_five_membered_1_3_two_heteroatoms_se_te_pair_raises` above).
        ("C1CNC[Se]1", "1,3-selenazolidine"),
        ("C1CNC[Te]1", "1,3-tellurazolidine"),
        ("C1C[Se]CO1", "1,3-oxaselenolane"),
        ("C1C[Te]CO1", "1,3-oxatellurolane"),
        ("C1C[Se]CS1", "1,3-thiaselenolane"),
        ("C1C[Te]CS1", "1,3-thiatellurolane"),
        ("C1C[Se]C[Se]1", "1,3-diselenolane"),
        ("C1C[Te]C[Te]1", "1,3-ditellurolane"),
    ],
)
def test_five_membered_1_3_two_heteroatom_ring_names_se_te(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # 1,2-related two-heteroatom saturated 5-membered rings -- Blue
        # Book Table 2.3 confirms pyrazolidine is a retained name (no
        # locants, unlike the mixed-heteroatom pairs below); Table 2.3
        # itself gives '1,2-oxazolidine (PIN)'/'1,2-thiazolidine (PIN)'
        # with locants, and P-22.2.2.1.3's own worked example gives
        # '1,2-oxathiolane (PIN)' verbatim (all three sourced directly
        # from the primary text, not PubChem -- PubChem's own computed
        # names for the O/S-only pairs omit the locants, which the
        # primary text's general locant-citation rule and its
        # '1,2-oxathiolane' example both contradict).
        ("C1CCNN1", "pyrazolidine"),
        ("C1CCON1", "1,2-oxazolidine"),
        ("C1CCSN1", "1,2-thiazolidine"),
        ("C1CCOO1", "1,2-dioxolane"),
        ("C1CCOS1", "1,2-oxathiolane"),
        ("C1CCSS1", "1,2-dithiolane"),
    ],
)
def test_five_membered_1_2_two_heteroatom_ring_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # 1,2-related two-heteroatom saturated 5-membered rings whose pair
        # includes Se and/or Te. N+Se/N+Te are Blue Book Table 2.3 retained
        # names, PubChem-verified with locants included ('1,2-selenazolidine'
        # CID 22597361, '1,2-tellurazolidine' CID 18381235). The remaining
        # six (O/S/Se/Te-only pairs) are systematic Hantzsch-Widman names
        # whose locants PubChem's own computed name drops (e.g. 'oxaselenolane'
        # for CID 129737263) -- per P-22.2.2.1.7 the locant can only be
        # omitted when there's no ambiguity, and a 1,2- vs 1,3- placement is
        # always ambiguous here, so the '1,2-' prefix is added per the
        # primary text rather than following PubChem's computed name,
        # mirroring the already-established '1,2-oxathiolane (PIN)' case
        # above.
        ("C1CC[Se]N1", "1,2-selenazolidine"),
        ("C1CC[Te]N1", "1,2-tellurazolidine"),
        ("C1CCO[Se]1", "1,2-oxaselenolane"),
        ("C1CCO[Te]1", "1,2-oxatellurolane"),
        ("C1CCS[Se]1", "1,2-thiaselenolane"),
        ("C1CCS[Te]1", "1,2-thiatellurolane"),
        ("C1CC[Se][Se]1", "1,2-diselenolane"),
        ("C1CC[Te][Te]1", "1,2-ditellurolane"),
    ],
)
def test_five_membered_1_2_two_heteroatom_ring_names_se_te(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_unsupported_heteroatom_element_raises():
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("C1CCCCP1")


def test_imidazole_ring_carbon_substituent_raises():
    # A substituent on imidazole's ring carbon (any locant other than 1,
    # the N-H position itself) is deliberately rejected rather than named:
    # imidazole's N-H can migrate to the *other* nitrogen (a real
    # prototropic tautomer), which reshuffles which carbon counts as
    # adjacent to N1 -- PubChem's own name for this exact SMILES
    # ("5-methyl-1H-imidazole", CID 13195) doesn't match this input's own
    # literal, structurally-fixed N-H position (locant 4), confirming
    # PubChem silently renormalizes to its own canonical tautomer before
    # naming. Without the Blue Book's own tie-breaking rule for this case,
    # the locant can't be trusted yet.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("Cc1cnc[nH]1")


def test_pyrazole_ring_carbon_substituent_raises():
    # Same tautomer ambiguity as imidazole above, one ring family over.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("Cc1cc[nH]n1")


def test_three_heteroatom_mancude_ring_raises():
    # 1,2,4-triazole: three ring heteroatoms, out of scope here.
    with pytest.raises(UnsupportedStructure):
        smiles_to_iupac("c1nc[nH]n1")
