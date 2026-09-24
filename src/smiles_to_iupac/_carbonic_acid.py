"""Naming of carbonic acid (HO-CO-OH) and its full/partial diesters
(R-O-CO-O-R'/R-O-CO-OH), per the IUPAC 2013 Recommendations ("the Blue
Book"):

- P-65.2.1 (Chapter P-6, https://iupac.qmul.ac.uk/BlueBook/PDF/P6.pdf,
  ~5775-5806): `carbonic acid (PIN)` (HO-CO-OH) is a simple retained
  mononuclear-acid name, no locants, no substituent variability --
  confirmed via PubChem CID 767 (`C(=O)(O)O` -> `"carbonic acid"`).
  Distinct from `_salt.py`'s existing carbonate-anion-in-a-salt coverage
  (e.g. "sodium carbonate") -- this module covers the free acid and its
  organic esters instead.
- P-67.1.3.2's ester citation style (already established by
  `_nitrate_ester.py`/`_sulfate.py`/`_phosphate.py`) applies unchanged:
  the substituent (alkyl/aryl) group(s) are cited as separate word(s),
  in alphanumeric order if more than one, followed by the acid's anion
  name -- reusing `_phosphate.py`'s exported `format_ester_words`
  verbatim, no new word-assembly logic needed. Confirmed via real
  PubChem structures: `COC(=O)OC` -> `"dimethyl carbonate"` (CID 12021),
  `CCOC(=O)OC` -> `"ethyl methyl carbonate"` (CID 522046).
- Partial ester (one -OH, one -O-R remaining): mirrors `_sulfate.py`'s
  own partial-ester citation, inserting the word "hydrogen" between the
  R-group word and "carbonate" -- confirmed worked example via PubChem,
  `COC(=O)O` -> `"methyl hydrogen carbonate"` (CID 78579).

Scope: a single, acyclic (not in any ring -- a cyclic carbonate like
ethylene carbonate is `_ketone.py`'s own Hantzsch-Widman-style ring-
ketone territory, e.g. '1,3-dioxolan-2-one', already handled there)
carbon bonded to exactly one `=O` and exactly two more single-bonded
oxygens, each either a plain `-OH` (free acid) or a plain `-O-R` (ester,
R named via `name_branch`, same alkyl/aryl/halogenated scope
`_sulfate.py`'s own R already covers), with no other atom on the central
carbon.

Explicitly out of scope (raise `UnsupportedStructure`): cyanic acid and
di-/tri-/tetra-/polycarbonic acids (P-65.2.2/.2.3, separate retained-name
families), any chalcogen-replacement (thio-/seleno-/telluro-carbonic
acid, P-65.2.1.2), carbamic acid (a different, N-containing shape), a
salt of a partial ester (P-67.1.3.2's own salt citation, not researched
for this acid), any ring/unsaturated R substituent beyond `name_branch`'s
own existing scope.
"""

from rdkit import Chem

from ._common import HALOGEN_PREFIXES, UnsupportedStructure, adjacency, halogen_substituents
from ._phosphate import format_ester_words
from ._substituents import name_branch

_CARBON = 6
_OXYGEN = 8


def _carbonic_acid_carbons(mol):
    """Carbon atoms shaped like carbonic acid or one of its esters: one
    C=O double bond, and two more single-bonded oxygens, each either an
    ester oxygen (degree 2 -- C plus one R carbon) or a plain uncharged
    hydroxyl (degree 1)."""
    ring_info = mol.GetRingInfo()
    matches = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != _CARBON or atom.GetDegree() != 3:
            continue
        if atom.GetFormalCharge() != 0 or atom.GetIsotope() != 0:
            continue
        if ring_info.NumAtomRings(atom.GetIdx()):
            # A cyclic carbonate (e.g. ethylene carbonate, 1,3-dioxolan-
            # 2-one) is `_ketone.py`'s own ring-ketone territory, not this
            # module's acyclic ester citation style -- see that module's
            # own Hantzsch-Widman-style ring naming.
            continue
        neighbors = atom.GetNeighbors()
        oxygens = [n for n in neighbors if n.GetAtomicNum() == _OXYGEN]
        if len(oxygens) != 3:
            continue
        double_bonded = [
            o for o in oxygens if mol.GetBondBetweenAtoms(atom.GetIdx(), o.GetIdx()).GetBondTypeAsDouble() == 2.0
        ]
        if len(double_bonded) != 1 or double_bonded[0].GetDegree() != 1 or double_bonded[0].GetFormalCharge() != 0:
            continue
        single_bonded = [o for o in oxygens if o.GetIdx() != double_bonded[0].GetIdx()]
        if any(
            o.GetFormalCharge() != 0
            or o.GetIsotope() != 0
            or o.GetDegree() not in (1, 2)
            or mol.GetBondBetweenAtoms(atom.GetIdx(), o.GetIdx()).GetBondTypeAsDouble() != 1.0
            for o in single_bonded
        ):
            continue
        matches.append(atom)
    return matches


def has_carbonic_acid_shape(mol) -> bool:
    return bool(_carbonic_acid_carbons(mol))


def name_carbonic_acid(mol) -> str:
    carbons = _carbonic_acid_carbons(mol)
    if len(carbons) != 1:
        raise UnsupportedStructure("more than one carbonic-acid-shaped carbon is not supported yet")
    (carbon,) = carbons

    group_atom_idxs = {carbon.GetIdx()} | {n.GetIdx() for n in carbon.GetNeighbors() if n.GetAtomicNum() == _OXYGEN}
    for atom in mol.GetAtoms():
        if atom.GetIdx() in group_atom_idxs:
            continue
        if atom.GetFormalCharge() != 0 or atom.GetIsotope() != 0:
            raise UnsupportedStructure("charged or isotopically modified atoms are not supported yet")
        if atom.GetAtomicNum() == _CARBON:
            continue
        if atom.GetAtomicNum() not in (1, *HALOGEN_PREFIXES):
            raise UnsupportedStructure(
                "heteroatoms other than the carbonic acid's own carbon/"
                "oxygens and a halogen substituent are not supported yet"
            )

    group_oxygens = {n.GetIdx() for n in carbon.GetNeighbors() if n.GetAtomicNum() == _OXYGEN}
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == _OXYGEN and atom.GetIdx() not in group_oxygens:
            raise UnsupportedStructure(
                "an oxygen atom not part of the carbonic acid's own "
                "O=C(OR)(OR') group is out of scope for this module"
            )

    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")

    single_bonded_oxygens = [
        idx
        for idx in group_oxygens
        if mol.GetBondBetweenAtoms(carbon.GetIdx(), idx).GetBondTypeAsDouble() == 1.0
    ]
    ester_oxygens = [idx for idx in single_bonded_oxygens if mol.GetAtomWithIdx(idx).GetDegree() == 2]
    hydroxyl_oxygens = [idx for idx in single_bonded_oxygens if mol.GetAtomWithIdx(idx).GetDegree() == 1]

    if not ester_oxygens:
        return "carbonic acid"

    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    aromatic_atoms = {atom.GetIdx() for atom in mol.GetAtoms() if atom.GetIsAromatic()}

    ester_names = []
    for ester_oxygen_idx in ester_oxygens:
        (root,) = [n for n in graph[ester_oxygen_idx] if n != carbon.GetIdx()]
        name, _ = name_branch(graph, root, ester_oxygen_idx, halogens, aromatic_atoms, mol=mol)
        ester_names.append(name)

    if hydroxyl_oxygens:
        return f"{format_ester_words(ester_names)} hydrogen carbonate"
    return format_ester_words(ester_names) + " carbonate"
