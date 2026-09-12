"""Fusion-locant-letter naming for a plain benzo ring ortho-fused onto
acephenanthrylene (P-25.3.1.3, https://iupac.qmul.ac.uk/BlueBook/PDF/P2.pdf)
-- same mechanism as `_anthracene_fusion.py`/`_phenanthrene_fusion.py`/
`_pyrene_fusion.py`/`_chrysene_fusion.py`/`_triphenylene_fusion.py`/
`_fluoranthene_fusion.py`/`_aceanthrylene_fusion.py` (all share
`_fusion_locant_letter.py`'s matching engine), this time for
acephenanthrylene (`_peri_fused_aromatic.py`'s own retained-name
reference SMILES) as the base component -- the fourth peri-fused base
(one interior atom, like fluoranthene/aceanthrylene) and, like
aceanthrylene, has *no* nontrivial automorphism at all, so each eligible
side is checked and named independently.

Acephenanthrylene's own numbering (1,2,3,3a,4,5,5a,6,6a,7,8,9,10,10a,
10b,10c, with 10c the sole interior atom) was reconstructed here from
first principles (SSSR + bond-ring-membership to find the periphery
cycle and the one interior atom), then anchored against three real
monomethylacephenanthrylene PubChem structures -- CID 14932034/14932039/
14932040 for the 1-/7-/10-methyl isomers respectively (the only three
methylacephenanthrylene isomers PubChem indexes) -- via substructure
match, which fixed the numbering direction unambiguously since there is
no symmetry to leave it in doubt. This gives P-25.3.1.3's continuous
peripheral lettering as: a(1,2), b(2,3), c(3,3a), d(3a,4), e(4,5),
f(5,5a), g(5a,6), h(6,6a), i(6a,7), j(7,8), k(8,9), l(9,10), m(10,10a),
n(10a,10b), o(10b,1). Of the six bonds not touching a ring-fusion atom
(3a, 5a, 6a, 10a, 10b; the interior atom 10c never touches the periphery
at all) -- a, b, e, j, k, l -- letter 'b' turns out to be the exact same
compound (InChI-confirmed) as `_aceanthrylene_fusion.py`'s own letter
'e' (PubChem CID 114910, 'benzo[e]aceanthrylene') and letter 'e' turns
out to be the exact same compound (InChI-confirmed) as
`_fluoranthene_fusion.py`'s own letter 'b'/'e' (PubChem CID 9153,
'benzo[b]fluoranthene') -- both excluded here to avoid a second route to
them, the same reason `_aceanthrylene_fusion.py` itself excludes its own
letter 'a'. The other four are each a distinct, real registered
compound, independently verified by InChI match against PubChem:
'benzo[a]acephenanthrylene' (CID 11482160), 'benzo[j]acephenanthrylene'
(CID 159608), 'benzo[k]acephenanthrylene' (CID 159607), and
'benzo[l]acephenanthrylene' (CID 10848420 -- indexed under the plain,
unlettered synonym 'benzoacephenanthrylene' only, confirmed distinct
from every other letter and every other base's own CIDs by InChI).

Scope, deliberately narrow, matching the other fusion modules: exactly
one plain, unsubstituted benzo ring ortho-fused onto acephenanthrylene
at one of these four distinct positions. Explicitly out of scope (raise
`UnsupportedStructure`):
- Fusion at any bond touching a ring-fusion atom (a peri-fused system of
  its own, different citation mechanism).
- Any substituent, any heteroatom anywhere, or more than one extra ring.
"""

from rdkit import Chem

from ._common import UnsupportedStructure
from ._fusion_locant_letter import find_fusion_letter

_ACEPHENANTHRYLENE_REF = Chem.MolFromSmiles("C1=CC=C2C(=C1)C=C3C=CC4=C3C2=CC=C4")
_LETTER_BY_PAIR = {
    frozenset({13, 14}): "a",
    frozenset({14, 15}): "b",
    frozenset({8, 9}): "e",
    frozenset({0, 5}): "j",
    frozenset({0, 1}): "k",
    frozenset({1, 2}): "l",
}
_EXCLUDED_LETTERS = frozenset({"b", "e"})


_NONAROMATIC_REF_ATOMS = frozenset({8, 9})


def _find_core(mol):
    return find_fusion_letter(
        mol,
        _ACEPHENANTHRYLENE_REF,
        _LETTER_BY_PAIR,
        _EXCLUDED_LETTERS,
        num_atoms=20,
        num_rings=5,
        nonaromatic_ref_atoms=_NONAROMATIC_REF_ATOMS,
    )


def has_acephenanthrylene_fusion_name(mol) -> bool:
    return _find_core(mol) is not None


def name_acephenanthrylene_fusion(mol) -> str:
    letter = _find_core(mol)
    if letter is None:
        raise UnsupportedStructure(
            "this pentacyclic system is not a supported benzo-fused "
            "acephenanthrylene shape (see P-25.3.1.3)"
        )
    return f"benzo[{letter}]acephenanthrylene"
