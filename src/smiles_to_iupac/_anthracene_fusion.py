"""Fusion-locant-letter naming for a plain benzo ring ortho-fused onto
anthracene, at one of the six periphery bonds that don't touch a ring-
fusion carbon (P-25.3.1.3, https://iupac.qmul.ac.uk/BlueBook/PDF/P2.pdf),
mirroring `_polycyclic_component_fusion.py`'s identical mechanism for
the indole/1-benzofuran base components -- this module is the general
algorithm's first *all-carbon* base component (indole/1-benzofuran are
both heteroatom-containing, so their fixed atom-role table never
collides with an ordinary retained-name PAH).

Anthracene's own numbering (P-25.3.3, 1,2,3,4,4a,5,6,7,8,8a,9,9a,10,10a)
gives P-25.3.1.3's continuous peripheral lettering (a for the '1,2' side,
b for '2,3', etc.) as: a(1,2), b(2,3), c(3,4), d(4,4a), e(4a,10),
f(10,10a), g(10a,5), h(5,6), i(6,7), j(7,8), k(8,8a), l(8a,9), m(9,9a),
n(9a,1). Only the six bonds not touching a ring-fusion atom (4a, 8a, 9a,
10a) support a plain ortho-fusion here: a, b, c, h, i, j.

Because anthracene has two independent symmetries (the long-axis mirror,
1<->4/2<->3/5<->8/6<->7, and the ring-swap rotation, 1<->8/2<->7/3<->6/
4<->5), a substituent at any of {a, c, h, j} gives the *same* compound
(the lowest-alphabet one, 'a', wins per P-25.3.1.3's own tie-break), and
a substituent at either of {b, i} gives the same (linear) compound too.
The linear shape is 'tetracene', already a retained name `_aromatic.py`
recognizes directly -- this module explicitly excludes it (raises
`UnsupportedStructure` via `_find_letter` returning None for it) so
`core.py` can fall through to that check instead, keeping this module
scoped to the one genuinely new name it adds: 'benzo[a]anthracene'
(PubChem CID 5954, `C1=CC=C2C(=C1)C=CC3=CC4=CC=CC=C4C=C32` -- also the
IUPAC Blue Book's own worked example for P-25.3.1.3's fusion-locant-
letter mechanism).

Scope, deliberately narrow, matching `_polycyclic_component_fusion.py`:
exactly one plain, unsubstituted benzo ring ortho-fused onto anthracene
at one of the six candidate periphery bonds, excluding the linear
('tetracene') shape. Explicitly out of scope (raise `UnsupportedStructure`):
- Fusion at any bond touching a ring-fusion atom (a peri-fused system,
  different citation mechanism), or the linear/tetracene shape.
- Any substituent, any heteroatom anywhere, or more than one extra ring.
"""

from rdkit import Chem

from ._common import UnsupportedStructure

_ANTHRACENE_REF = Chem.MolFromSmiles("c1cccc2cc3ccccc3cc12")
_LETTER_BY_PAIR = {
    frozenset({10, 9}): "a",
    frozenset({9, 8}): "b",
    frozenset({8, 7}): "c",
    frozenset({3, 2}): "h",
    frozenset({2, 1}): "i",
    frozenset({1, 0}): "j",
}
_EXCLUDED_LETTERS = {"b"}


def _find_letter(mol):
    matches = mol.GetSubstructMatches(_ANTHRACENE_REF, useChirality=False, uniquify=False)
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


def has_anthracene_fusion_name(mol) -> bool:
    return _find_core(mol) is not None


def name_anthracene_fusion(mol) -> str:
    letter = _find_core(mol)
    if letter is None:
        raise UnsupportedStructure(
            "this tetracyclic system is not a supported benzo-fused "
            "anthracene shape (see P-25.3.1.3)"
        )
    return f"benzo[{letter}]anthracene"
