"""Naming of urea (H2N-C(=O)-NH2), per the IUPAC 2013 Recommendations
("the Blue Book"):

- Chapter P-6 (https://iupac.qmul.ac.uk/BlueBook/PDF/P6.pdf): 'urea' is a
  retained name that is itself the preferred IUPAC name -- not a
  systematic construction -- the same way 'carbamic acid' (`_carbamate.py`)
  is retained rather than derived. Confirmed via PubChem structure match:
  `NC(=O)N` -> "urea".
- Unlike every other module in this project, urea has no parent-hydride
  chain to number and no chain-length logic of its own: the name is
  always the single fixed word 'urea'.

Scope, deliberately narrow (a first pass at this parent, mirroring how
`_carbamate.py` originally deferred free carbamic acid itself): only the
exact unsubstituted shape (H2N-C(=O)-NH2, no other atoms at all) is
supported. N-substituted ureas (e.g. 'methylurea', '1,3-dimethylurea') are
deferred -- PubChem's own generated names for those use numeric locants
rather than the classical N/N' convention, and whether that matches the
Blue Book's own PIN locant style needs separate confirmation from the
source text before it's implemented. A ring-fused urea (e.g. hydantoin)
and thiourea (the sulfur analogue) are both out of scope entirely.
"""

from ._common import UnsupportedStructure


def _urea_carbon(mol):
    """The urea carbonyl carbon, or None if the molecule isn't shaped like
    plain, unsubstituted urea (H2N-C(=O)-NH2, and nothing else)."""
    atoms = list(mol.GetAtoms())
    if len(atoms) != 4:
        return None
    for atom in atoms:
        if atom.GetAtomicNum() != 6:
            continue
        if atom.GetDegree() != 3 or atom.GetFormalCharge() != 0 or atom.GetIsAromatic():
            continue
        neighbors = atom.GetNeighbors()
        oxygens = [n for n in neighbors if n.GetAtomicNum() == 8]
        nitrogens = [n for n in neighbors if n.GetAtomicNum() == 7]
        if len(oxygens) != 1 or len(nitrogens) != 2:
            continue
        (oxygen,) = oxygens
        if oxygen.GetDegree() != 1 or mol.GetBondBetweenAtoms(atom.GetIdx(), oxygen.GetIdx()).GetBondTypeAsDouble() != 2.0:
            continue
        if any(
            n.GetDegree() != 1
            or n.GetTotalNumHs() != 2
            or n.GetFormalCharge() != 0
            or n.GetIsotope() != 0
            or mol.GetBondBetweenAtoms(atom.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() != 1.0
            for n in nitrogens
        ):
            continue
        return atom.GetIdx()
    return None


def has_urea_shape(mol) -> bool:
    return _urea_carbon(mol) is not None


def name_urea(mol) -> str:
    if _urea_carbon(mol) is None:
        raise UnsupportedStructure(
            "no plain, unsubstituted urea (H2N-C(=O)-NH2) shape found; "
            "this module only handles unsubstituted urea itself"
        )
    return "urea"
