"""Fusion-locant-letter naming for a plain benzo ring ortho-fused onto
phenanthrene, at the one periphery bond that doesn't touch a ring-fusion
carbon and doesn't collide with an already-recognized retained name
(P-25.3.1.3, https://iupac.qmul.ac.uk/BlueBook/PDF/P2.pdf) -- mirroring
`_anthracene_fusion.py`'s identical mechanism for anthracene, the other
tricyclic all-carbon base component.

Phenanthrene's own numbering (P-25.3.3, 1,2,3,4,4a,4b,5,6,7,8,8a,9,10,10a)
gives P-25.3.1.3's continuous peripheral lettering as: a(1,2), b(2,3),
c(3,4), d(4,4a), e(4a,4b), f(4b,5), g(5,6), h(6,7), i(7,8), j(8,8a),
k(8a,9), l(9,10), m(10,10a), n(10a,1). Of the seven bonds not touching a
ring-fusion atom (4a, 4b, 8a, 10a) -- a, b, c, g, h, i, l -- phenanthrene's
own mirror symmetry (ring A <-> ring C) makes g/h/i the same compound as
c/b/a respectively (already found automatically by trying every
substructure-match automorphism, same as `_anthracene_fusion.py`), so
only three *distinct* shapes are possible: 'a' (chrysene, CID 9171 --
already a retained name, not yet implemented anywhere in this project),
'b' (benzo[a]anthracene, CID 5954 -- the exact same compound
`_anthracene_fusion.py` already names via anthracene as the base
component instead), and 'l' (triphenylene, CID 9170 -- already a
retained name, `_branched_fused_aromatic.py`). All three are excluded
here so this module claims only the one genuinely new name it adds:
'benzo[c]phenanthrene' (PubChem CID 9136).

Scope, deliberately narrow, matching `_anthracene_fusion.py`: exactly one
plain, unsubstituted benzo ring ortho-fused onto phenanthrene at the 'c'
bond (3,4). Explicitly out of scope (raise `UnsupportedStructure`):
- Fusion at any bond touching a ring-fusion atom (a peri-fused system,
  different citation mechanism), or the 'a'/'b'/'l' shapes above (all
  three either not yet implemented as a retained name, or already
  produced via a different base component).
- Any substituent, any heteroatom anywhere, or more than one extra ring.
"""

from rdkit import Chem

from ._common import UnsupportedStructure

_PHENANTHRENE_REF = Chem.MolFromSmiles("c1ccc2ccc3ccccc3c2c1")
_LETTER_BY_PAIR = {
    frozenset({1, 2}): "a",
    frozenset({0, 1}): "b",
    frozenset({0, 13}): "c",
    frozenset({9, 10}): "g",
    frozenset({8, 9}): "h",
    frozenset({7, 8}): "i",
    frozenset({4, 5}): "l",
}
_EXCLUDED_LETTERS = {"a", "b", "l"}


def _find_letter(mol):
    matches = mol.GetSubstructMatches(_PHENANTHRENE_REF, useChirality=False, uniquify=False)
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
    if best in _EXCLUDED_LETTERS:
        return None
    return best


def _find_core(mol):
    if mol.GetNumAtoms() != 18:
        return None
    if mol.GetRingInfo().NumRings() != 4:
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


def has_phenanthrene_fusion_name(mol) -> bool:
    return _find_core(mol) is not None


def name_phenanthrene_fusion(mol) -> str:
    letter = _find_core(mol)
    if letter is None:
        raise UnsupportedStructure(
            "this tetracyclic system is not a supported benzo-fused "
            "phenanthrene shape (see P-25.3.1.3)"
        )
    return f"benzo[{letter}]phenanthrene"
