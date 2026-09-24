"""Naming of the double-bond-junction two-ring assembly (two identical
monocyclic saturated all-carbon rings joined by a C=C double bond between
one ring atom of each, e.g. bi(cyclopentylidene)), per the IUPAC 2013
Recommendations ("the Blue Book"):

- P-28.2.2 (Chapter P-2, https://iupac.qmul.ac.uk/BlueBook/PDF/P2.pdf):
  "When two cyclic systems are linked by a double bond, method (2)
  [P-28.2.1's substituent-group-name method: the prefix 'bi' before the
  name of the corresponding substituent group, enclosed in parentheses]
  is the only recommended method." Confirmed PIN worked example
  `1,1'-bi(cyclopentylidene)` (`tmp/bluebook/P2.txt` ~7826) -- not the
  older CAS-style Δ-locant form.
- The substituent-group name for a monocyclic ring with a single divalent
  radical center ("-ylidene") is P-71.2.1.1/P-71.2.2.1's own rule, already
  proven by `_radical.py`'s `_name_ring_radical`: for an unsubstituted
  monocyclic saturated ring of size N, the name is always
  "cyclo<N>ylidene" (e.g. "cyclopentylidene") with no locant, since a
  monocyclic ring's radical position is symmetric ("any position", per
  the specific method P-71.2.1.1) -- this module computes that same
  "cyclo" + `_numerals.alkyl_name` + "idene" string directly rather than
  importing `_radical.py`'s private helper, since `_name_ring_radical`
  itself only accepts a *fully* unsubstituted ring (it explicitly rejects
  any substituent, including a halogen) and can't be reused unchanged
  once a halogen substituent (in scope here, see below) is present.
- P-14.3.2 (numbering): each ring is numbered independently, starting at
  its own double-bond-junction atom (locant 1), with one ring's locants
  left unprimed and the other's primed -- the identical primed/unprimed
  lowest-locant-set numbering `_ring_assembly.py` already uses for its own
  single-bond-junction (P-28.2.1) case, reused here via its exported
  `_numberings_from_attachment` building block (the double-bond-junction
  atom plays the same "locant 1" role its single-bond attachment atom
  does there).
- P-35.2.1 (Chapter P-3): halogen substituents hang off a ring atom the
  same way as in every other ring module.

Scope: two disjoint, identical-size monocyclic saturated all-carbon rings
(any ring size), each contributing exactly one ring atom to a C=C double
bond between them (no other inter-ring bond), each optionally bearing
halogen substituents (positions may differ between the two rings, exactly
as `_ring_assembly.py`'s own "identical parent, independently
substituted" biphenyl case already allows). Explicitly out of scope
(raise `UnsupportedStructure` via the generic fallback in `core.py`,
since `find_ring_assembly_ylidene_core` below simply returns None for any
of these): a von Baeyer bicyclic/polycyclic or monospiro parent ring on
either side (`2,2'-bi(bicyclo[2.2.1]heptanylidene)`-style, a follow-up
step), an aromatic ring, two different ring sizes, more than 2 rings, and
indicated hydrogen (P-28.2.3 -- not reachable by any all-carbon saturated
ring anyway).
"""

from ._common import (
    UnsupportedStructure,
    adjacency,
    group_substituents,
    halogen_substituents,
    ring_cycle,
    substituent_locant_set_and_citation,
    validate_atoms_and_bonds,
)
from ._numerals import alkyl_name
from ._ring_assembly import _numberings_from_attachment
from ._substituents import format_substituent_prefixes, name_branch


def find_ring_assembly_ylidene_core(mol):
    """Return (ring0_atoms, ring1_atoms, junction0, junction1) if `mol` is
    exactly two disjoint, identical-size monocyclic saturated all-carbon
    rings joined by one C=C double bond (one ring atom from each), else
    None."""
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    if len(atom_rings) != 2:
        return None
    if len(atom_rings[0]) != len(atom_rings[1]):
        return None
    for ring in atom_rings:
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetIsAromatic() or atom.GetAtomicNum() != 6:
                return None

    ring0, ring1 = set(atom_rings[0]), set(atom_rings[1])
    if ring0 & ring1:
        return None

    connecting = [
        bond
        for bond in mol.GetBonds()
        if (bond.GetBeginAtomIdx() in ring0 and bond.GetEndAtomIdx() in ring1)
        or (bond.GetBeginAtomIdx() in ring1 and bond.GetEndAtomIdx() in ring0)
    ]
    if len(connecting) != 1:
        return None
    bond = connecting[0]
    if bond.GetIsAromatic() or bond.GetBondTypeAsDouble() != 2.0:
        return None

    x, y = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
    junction0, junction1 = (x, y) if x in ring0 else (y, x)
    for atom in (mol.GetAtomWithIdx(junction0), mol.GetAtomWithIdx(junction1)):
        if atom.GetDegree() != 3:
            return None

    return atom_rings[0], atom_rings[1], junction0, junction1


def _ylidene_name(ring_size: int) -> str:
    return "cyclo" + alkyl_name(ring_size) + "idene"


def _candidate_key(locants, ring_atoms, ylidene_name, graph, halogens, mol=None):
    substituents = {}
    for atom, position in locants.items():
        branch_roots = [n for n in graph[atom] if n not in ring_atoms]
        if branch_roots:
            substituents[position] = [name_branch(graph, root, atom, halogens, mol=mol) for root in branch_roots]

    grouped = group_substituents(substituents)
    locant_set, _, citation_locants = substituent_locant_set_and_citation(grouped)
    prefix = format_substituent_prefixes(grouped) if grouped else ""
    if prefix:
        prefix += "-"
    name = f"{prefix}1,1'-bi({ylidene_name})"
    return locant_set, citation_locants, name


def name_ring_assembly_ylidene(mol, core) -> str:
    validate_atoms_and_bonds(mol)

    ring0_atoms, ring1_atoms, junction0, junction1 = core
    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    ring_atoms = set(ring0_atoms) | set(ring1_atoms)
    ylidene_name = _ylidene_name(len(ring0_atoms))

    best_key = None
    best_name = None
    for (ring_a, attach_a), (ring_b, attach_b) in (
        ((ring0_atoms, junction0), (ring1_atoms, junction1)),
        ((ring1_atoms, junction1), (ring0_atoms, junction0)),
    ):
        for locants_a in _numberings_from_attachment(graph, ring_a, attach_a, prime=False):
            for locants_b in _numberings_from_attachment(graph, ring_b, attach_b, prime=True):
                locants = {**locants_a, **locants_b}
                key = _candidate_key(locants, ring_atoms, ylidene_name, graph, halogens, mol=mol)
                if best_key is None or key < best_key:
                    best_key, best_name = key, key[-1]

    return best_name
