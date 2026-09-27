"""One-atom bridge citation (P-25.4.2.1.4) generalized to any
skeleton-dict-matched alicyclic parent hydride (currently the seven
1989 IUPAC steroid parents in `_steroid_parent_hydrides.py`), not just
`_bridged_aromatic.py`'s all-carbon mancude hex-lattice bases.

Structurally this is a different bridge shape from `_bridged_aromatic.py`'s
1,4-type transannular bridge: here the bridge atom's two neighbors are
*already* directly bonded to each other in the bare parent skeleton (e.g.
steroid ring-A locants 5 and 6, adjacent ring carbons), so the bridge adds
a fused three-membered ring (an oxirane, for the O case) across an existing
ring bond, rather than spanning two non-adjacent positions of a larger
ring. No added/subtracted hydrogen prefix is needed either way, since the
bare steroid skeleton is already fully saturated (unlike the aromatic
case, where bridging forces two ring atoms out of the mancude system and
so needs 'dihydro').

- P-25.4.2.1.4 (Chapter P-2, https://iupac.qmul.ac.uk/BlueBook/PDF/P2.pdf):
  the preselected bridge prefix for a divalent -O- bridge is 'epoxy' (not
  'epoxidano'); the vocabulary and locant-citation convention are shared
  verbatim with `_bridged_aromatic.py` (see its own docstring for the
  primary-source citations of 'sulfano'/-S- and 'azano'/-NH-, structurally
  identical in shape to 'epoxy' here, degree-2, no further substitution).
- Verified against a real registered structure: PubChem CID 281912, whose
  synonym list includes the plain (non-stereo) CAS-style name
  '5,6-epoxycholestane' alongside the fully stereo-specified
  '5.alpha.,6.alpha.-Epoxycholestane' -- confirming the bare bridge-prefix
  name (no stereo descriptor) is itself a real, correct name, not merely
  an intermediate form. Stripping CID 281912's bridge oxygen and matching
  the remainder against `_PLAIN_CANONICAL_TO_NAME` reproduces 'cholestane'
  exactly, and the bridge oxygen's two ring neighbors recover as locants
  5 and 6 via `_locant_map` -- independently confirming both the skeleton
  match and the locant pair, not assumed from the CAS name alone.

Explicitly out of scope (falls through to whatever else `core.py` tries,
raising `UnsupportedStructure` if nothing else claims it):
- Stereo descriptor citation for the bridge (e.g. the 'alpha'/'beta' face
  distinguishing which of the two possible epoxide diastereomers a real
  molecule is) -- this module matches and names by constitution only,
  same policy as `_steroid_parent_hydrides.py`'s own `steroid_suffix_name`
  (see that module's docstring). A real input's stereochemistry, if any,
  is simply not reflected in the bridge locants' citation.
- A bridge atom of any other identity than O/S/N (the vocabulary
  `_bridged_aromatic._BRIDGE_PREFIXES` already covers) or with any
  substituent of its own.
- More than one bridge, or a bridge whose two neighbors are *not* already
  directly bonded (a genuine 1,4-type transannular span across a larger
  ring, rather than an ortho-fused three-membered ring) -- a structurally
  different case not attempted here.
- Any skeleton other than the seven bare steroid parent hydrides (no
  alkaloid skeleton dict exists yet, see #1065/#1144's tracked follow-on
  work) or any substituent on the ring system beyond the bridge itself.
"""

from rdkit import Chem

from ._bridged_aromatic import _BRIDGE_PREFIXES
from ._parent_hydride_stripping import strip_substituents
from ._steroid_parent_hydrides import _PLAIN_CANONICAL_TO_NAME, _locant_map

_BRIDGE_ATOM_H_COUNT = {7: 1, 8: 0, 16: 0}


def _bridge_candidates(mol):
    candidates = []
    for atom in mol.GetAtoms():
        if atom.GetIsAromatic() or atom.GetDegree() != 2:
            continue
        atomic_num = atom.GetAtomicNum()
        expected_h = _BRIDGE_ATOM_H_COUNT.get(atomic_num)
        if expected_h is None or atom.GetTotalNumHs() != expected_h:
            continue
        if atom.GetFormalCharge() != 0 or atom.GetIsotope() != 0:
            continue
        a, b = atom.GetNeighbors()
        if mol.GetBondBetweenAtoms(a.GetIdx(), b.GetIdx()) is None:
            continue
        candidates.append(atom)
    return candidates


def find_bridged_steroid_core(mol):
    """Return (bridge_atom_idx, bridge_prefix, skeleton_name, stripped_mol,
    old_to_new) if `mol` is exactly one of the seven bare steroid parent
    skeletons plus one O/S/N one-atom bridge across an already-adjacent
    ring bond, else None."""
    candidates = _bridge_candidates(mol)
    if len(candidates) != 1:
        return None
    bridge_atom = candidates[0]
    bridge_prefix = _BRIDGE_PREFIXES[bridge_atom.GetAtomicNum()]

    stripped, old_to_new = strip_substituents(mol, {bridge_atom.GetIdx()})
    if stripped is None:
        return None
    Chem.RemoveStereochemistry(stripped)
    name = _PLAIN_CANONICAL_TO_NAME.get(Chem.MolToSmiles(stripped))
    if name is None:
        return None
    return bridge_atom.GetIdx(), bridge_prefix, name, stripped, old_to_new


def name_bridged_steroid(mol, core) -> str:
    bridge_idx, bridge_prefix, name, stripped, old_to_new = core
    bridge_atom = mol.GetAtomWithIdx(bridge_idx)
    a, b = (n.GetIdx() for n in bridge_atom.GetNeighbors())
    locant_of_atom = {atom: locant for locant, atom in _locant_map(name, stripped).items()}
    locant_a, locant_b = locant_of_atom[old_to_new[a]], locant_of_atom[old_to_new[b]]
    lower, upper = sorted((locant_a, locant_b))
    return f"{lower},{upper}-{bridge_prefix}{name}"


def has_bridged_steroid_name(mol) -> bool:
    return find_bridged_steroid_core(mol) is not None


def name_bridged_steroid_parent(mol) -> str:
    return name_bridged_steroid(mol, find_bridged_steroid_core(mol))
