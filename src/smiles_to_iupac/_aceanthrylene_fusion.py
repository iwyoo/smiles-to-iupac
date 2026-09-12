"""Fusion-locant-letter naming for a plain benzo ring ortho-fused onto
aceanthrylene (P-25.3.1.3, https://iupac.qmul.ac.uk/BlueBook/PDF/P2.pdf)
-- same mechanism as `_anthracene_fusion.py`/`_phenanthrene_fusion.py`/
`_pyrene_fusion.py`/`_chrysene_fusion.py`/`_triphenylene_fusion.py`/
`_fluoranthene_fusion.py` (all share `_fusion_locant_letter.py`'s
matching engine), this time for aceanthrylene (`_peri_fused_aromatic.py`'s
own retained-name reference SMILES) as the base component -- the third
peri-fused base (one interior atom, like fluoranthene), but unlike every
prior base, aceanthrylene has *no* nontrivial automorphism at all (its
only substructure self-match is the identity), so none of its eligible
sides collapse into each other by symmetry -- each one is checked and
named independently.

Aceanthrylene's own numbering (1,2,2a,3,4,5,5a,6,6a,7,8,9,10,10a,10b,10c,
with 10c the sole interior atom) was reconstructed here from first
principles (SSSR + bond-ring-membership to find the periphery cycle and
the one interior atom, then anchored against four real
monomethylaceanthrylene PubChem structures -- CID 154217147/23539808/
3025911/23539809 for the 3-/5-/6-/7-methyl isomers respectively -- via
substructure match, which fixed the numbering direction unambiguously
since there is no symmetry to leave it in doubt) gives P-25.3.1.3's
continuous peripheral lettering as: a(1,2), b(2,2a), c(2a,3), d(3,4),
e(4,5), f(5,5a), g(5a,6), h(6,6a), i(6a,7), j(7,8), k(8,9), l(9,10),
m(10,10a), n(10a,10b), o(10b,1). Of the six bonds not touching a
ring-fusion atom (2a, 5a, 6a, 10a, 10b; the interior atom 10c never
touches the periphery at all) -- a, d, e, j, k, l -- letter 'a' turns out
to be the exact same compound (InChI-confirmed) as
`_fluoranthene_fusion.py`'s own letter 'a'/'f' (PubChem CID 9146,
'benzo[a]fluoranthene') and is excluded here to avoid a second route to
it, the same reason `_phenanthrene_fusion.py`/`_chrysene_fusion.py`/
`_triphenylene_fusion.py` each exclude one of their own letters. The
other five are each a distinct, real registered compound, independently
verified by InChI match against PubChem: 'benzo[d]aceanthrylene' (CID
155608), 'benzo[e]aceanthrylene' (CID 114910), 'benzo[j]aceanthrylene'
(CID 104987 -- PubChem's own computed IUPAC name agrees exactly, and is
the only one of the five for which PubChem's algorithm produces a fusion
name rather than falling back to von Baeyer nomenclature),
'benzo[k]aceanthrylene' (CID 146307), and 'benzo[l]aceanthrylene' (CID
105092). (PubChem's synonym lists carry these under the older, elided
CAS style 'benz[x]aceanthrylene'; P-25.3.1.2.2 explicitly says the 2013
recommendations no longer elide the attached-component prefix's final
'o' before a vowel, so 'benzo[x]aceanthrylene' -- matching CID 104987's
own computed IUPAC name -- is what this module produces.)

Scope, deliberately narrow, matching the other fusion modules: exactly
one plain, unsubstituted benzo ring ortho-fused onto aceanthrylene at one
of these five distinct positions. Explicitly out of scope (raise
`UnsupportedStructure`):
- Fusion at any bond touching a ring-fusion atom (a peri-fused system of
  its own, different citation mechanism).
- Any substituent, any heteroatom anywhere, or more than one extra ring.
"""

from rdkit import Chem

from ._fusion_component_registry import make_fusion_component_functions

_ACEANTHRYLENE_REF = Chem.MolFromSmiles("C1=CC=C2C3=C4C(=CC=CC4=CC2=C1)C=C3")
_LETTER_BY_PAIR = {
    frozenset({14, 15}): "a",
    frozenset({7, 8}): "d",
    frozenset({8, 9}): "e",
    frozenset({0, 13}): "j",
    frozenset({0, 1}): "k",
    frozenset({1, 2}): "l",
}
_EXCLUDED_LETTERS = frozenset({"a"})


_NONAROMATIC_REF_ATOMS = frozenset({14, 15})

has_aceanthrylene_fusion_name, name_aceanthrylene_fusion = make_fusion_component_functions(
    "aceanthrylene",
    _ACEANTHRYLENE_REF,
    _LETTER_BY_PAIR,
    _EXCLUDED_LETTERS,
    num_atoms=20,
    num_rings=5,
    nonaromatic_ref_atoms=_NONAROMATIC_REF_ATOMS,
)
