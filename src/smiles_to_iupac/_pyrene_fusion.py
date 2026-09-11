"""Fusion-locant-letter naming for a plain benzo ring ortho-fused onto
pyrene (P-25.3.1.3, https://iupac.qmul.ac.uk/BlueBook/PDF/P2.pdf) --
mirroring `_anthracene_fusion.py`/`_phenanthrene_fusion.py`'s identical
mechanism, this time for pyrene (`_peri_fused_aromatic.py`'s own
retained-name reference SMILES) as the base component -- the first
*peri*-fused (not simply catacondensed) base for this general algorithm.

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

_PYRENE_REF = Chem.MolFromSmiles("c1cc2ccc3cccc4ccc(c1)c2c34")
_LETTER_BY_PAIR = {
    frozenset({0, 1}): "a",
    frozenset({0, 13}): "b",
    frozenset({10, 11}): "e",
    frozenset({7, 8}): "h",
    frozenset({6, 7}): "i",
    frozenset({3, 4}): "l",
}


def _find_letter(mol):
    matches = mol.GetSubstructMatches(_PYRENE_REF, useChirality=False, uniquify=False)
    best = None
    for match in matches:
        core = set(match)
        extra = [a.GetIdx() for a in mol.GetAtoms() if a.GetIdx() not in core]
        if len(extra) != 4:
            continue
        if any(
            mol.GetAtomWithIdx(a).GetAtomicNum() != 6 or not mol.GetAtomWithIdx(a).GetIsAromatic()
            for a in extra
        ):
            continue
        target_to_ref_idx = {match[i]: i for i in range(len(match))}
        fusion_ref_idxs = set()
        bad = False
        for a in extra:
            for n in mol.GetAtomWithIdx(a).GetNeighbors():
                if n.GetIdx() in target_to_ref_idx:
                    fusion_ref_idxs.add(target_to_ref_idx[n.GetIdx()])
                elif n.GetIdx() not in extra:
                    bad = True
        if bad or len(fusion_ref_idxs) != 2:
            continue
        letter = _LETTER_BY_PAIR.get(frozenset(fusion_ref_idxs))
        if letter is not None and (best is None or letter < best):
            best = letter
    return best


def _find_core(mol):
    if mol.GetNumAtoms() != 20:
        return None
    if mol.GetRingInfo().NumRings() != 5:
        return None
    if any(not atom.GetIsAromatic() for atom in mol.GetAtoms()):
        return None
    for atom in mol.GetAtoms():
        if atom.GetFormalCharge() != 0 or atom.GetIsotope() != 0:
            return None
        if atom.GetAtomicNum() != 6:
            return None
    if len(Chem.GetMolFrags(mol)) > 1:
        return None
    return _find_letter(mol)


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
