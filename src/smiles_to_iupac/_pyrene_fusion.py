"""Fusion-locant-letter naming for a plain benzo ring ortho-fused onto
pyrene (P-25.3.1.3, https://iupac.qmul.ac.uk/BlueBook/PDF/P2.pdf) --
mirroring `_anthracene_fusion.py`/`_phenanthrene_fusion.py`'s identical
mechanism (all three share `_fusion_locant_letter.py`'s matching engine),
this time for pyrene (`_peri_fused_aromatic.py`'s own retained-name
reference SMILES) as the base component -- the first *peri*-fused (not
simply catacondensed) base for this general algorithm.

Pyrene's own numbering (P-25.3.3, 1,2,3,3a,4,5,5a,6,7,8,8a,9,10,10a,10b,
10c) gives P-25.3.1.3's continuous peripheral lettering as: a(1,2),
b(2,3), c(3,3a), d(3a,4), e(4,5), f(5,5a), g(5a,6), h(6,7), i(7,8),
j(8,8a), k(8a,9), l(9,10), m(10,10a), n(10a,1). Of the six bonds not
touching a ring-fusion atom (3a, 5a, 8a, 10a; the two further-interior
atoms 10b/10c never touch the periphery at all) -- a, b, e, h, i, l --
pyrene's own high (D2h) symmetry collapses these to just *two* distinct
shapes (confirmed by trying every substructure-match automorphism, same
as the other two fusion modules): {a, b, h, i} are all the same compound
('benzo[a]pyrene', PubChem CID 2336 -- the well-known environmental
carcinogen, letter 'a' the lowest of the four and so the one this
module's own tie-break picks), and {e, l} are the other
('benzo[e]pyrene', PubChem CID 9128). Unlike the anthracene/phenanthrene
cases, neither shape collides with any other already-recognized retained
name, so both are supported here.

Scope, deliberately narrow, matching `_anthracene_fusion.py`: exactly one
plain, unsubstituted benzo ring ortho-fused onto pyrene at one of these
two distinct positions. Explicitly out of scope (raise
`UnsupportedStructure`):
- Fusion at any bond touching a ring-fusion atom (a peri-fused system of
  its own, different citation mechanism).
- Any substituent, any heteroatom anywhere, or more than one extra ring.
"""

from rdkit import Chem

from ._common import UnsupportedStructure
from ._fusion_locant_letter import find_fusion_letter

_PYRENE_REF = Chem.MolFromSmiles("c1cc2ccc3cccc4ccc(c1)c2c34")
_LETTER_BY_PAIR = {
    frozenset({0, 1}): "a",
    frozenset({0, 13}): "b",
    frozenset({10, 11}): "e",
    frozenset({7, 8}): "h",
    frozenset({6, 7}): "i",
    frozenset({3, 4}): "l",
}
_EXCLUDED_LETTERS = frozenset()


def _find_core(mol):
    return find_fusion_letter(mol, _PYRENE_REF, _LETTER_BY_PAIR, _EXCLUDED_LETTERS, num_atoms=20, num_rings=5)


def has_pyrene_fusion_name(mol) -> bool:
    return _find_core(mol) is not None


def name_pyrene_fusion(mol) -> str:
    letter = _find_core(mol)
    if letter is None:
        raise UnsupportedStructure(
            "this pentacyclic system is not a supported benzo-fused "
            "pyrene shape (see P-25.3.1.3)"
        )
    return f"benzo[{letter}]pyrene"
