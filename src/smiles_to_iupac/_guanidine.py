"""Naming of guanidine (HN=C(NH2)2), per the IUPAC 2013 Recommendations
("the Blue Book"):

- Chapter P-6 (https://iupac.qmul.ac.uk/BlueBook/PDF/P6.pdf): 'guanidine'
  is a retained name that is itself the preferred IUPAC name -- not a
  systematic construction -- the same way 'urea' (`_urea.py`) and
  'thiourea' (`_thiourea.py`) are retained rather than derived. Confirmed
  via PubChem structure match: `NC(=N)N` -> "guanidine". Like urea,
  guanidine has no parent-hydride chain to number and no chain-length
  logic of its own.
- N-substituted guanidines do have a confirmed worked example in the Blue
  Book's own text (`tmp/bluebook/P6.txt` line 797:
  'N'-nitro-N-nitroso-N-propylguanidine (PIN)'), but guanidine's three
  nitrogens (one imino =NH, two amino -NH2) make its substitution/locant
  rules a distinct, more involved follow-up from `_urea.py`'s own
  two-nitrogen case -- deferred entirely for this first pass.

Scope, deliberately narrow (a first pass at this parent, mirroring how
`_urea.py` started with unsubstituted urea alone): only the exact
unsubstituted shape (HN=C(NH2)2, no other atoms at all) is supported.
N-substituted guanidines and any ring-fused guanidine are both out of
scope entirely.
"""

from ._common import UnsupportedStructure


def _guanidine_carbon(mol):
    """The guanidine carbon, or None if the molecule isn't shaped like
    plain, unsubstituted guanidine (HN=C(NH2)2, and nothing else)."""
    atoms = list(mol.GetAtoms())
    if len(atoms) != 4:
        return None
    for atom in atoms:
        if atom.GetAtomicNum() != 6:
            continue
        if atom.GetDegree() != 3 or atom.GetFormalCharge() != 0 or atom.GetIsAromatic():
            continue
        nitrogens = [n for n in atom.GetNeighbors() if n.GetAtomicNum() == 7]
        if len(nitrogens) != 3:
            continue
        imino = [
            n
            for n in nitrogens
            if mol.GetBondBetweenAtoms(atom.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 2.0
        ]
        amino = [
            n
            for n in nitrogens
            if mol.GetBondBetweenAtoms(atom.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 1.0
        ]
        if len(imino) != 1 or len(amino) != 2:
            continue
        (imino_n,) = imino
        if imino_n.GetDegree() != 1 or imino_n.GetTotalNumHs() != 1 or imino_n.GetFormalCharge() != 0:
            continue
        if any(
            n.GetDegree() != 1
            or n.GetTotalNumHs() != 2
            or n.GetFormalCharge() != 0
            or n.GetIsotope() != 0
            for n in amino
        ):
            continue
        if imino_n.GetIsotope() != 0:
            continue
        return atom.GetIdx()
    return None


def has_guanidine_shape(mol) -> bool:
    return _guanidine_carbon(mol) is not None


def name_guanidine(mol) -> str:
    if _guanidine_carbon(mol) is None:
        raise UnsupportedStructure(
            "no plain, unsubstituted guanidine (HN=C(NH2)2) shape found; "
            "this module only handles unsubstituted guanidine itself"
        )
    return "guanidine"
