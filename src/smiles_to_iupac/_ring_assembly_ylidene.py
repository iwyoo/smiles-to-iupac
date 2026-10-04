"""Naming of the double-bond-junction two-ring assembly (two identical
rings or ring systems joined by a C=C double bond between one atom of
each, e.g. bi(cyclopentylidene), bi(bicyclo[2.2.1]heptanylidene)), per the
IUPAC 2013 Recommendations ("the Blue Book"):

- P-28.2.2 (Chapter P-2, https://iupac.qmul.ac.uk/BlueBook/PDF/P2.pdf):
  "When two cyclic systems are linked by a double bond, method (2)
  [P-28.2.1's substituent-group-name method: the prefix 'bi' before the
  name of the corresponding substituent group, enclosed in parentheses]
  is the only recommended method." Confirmed PIN worked examples
  `1,1'-bi(cyclopentylidene)` and `2,2'-bi(bicyclo[2.2.1]heptanylidene)`
  (`tmp/bluebook/P2.txt` ~7826) -- not the older CAS-style Δ-locant form.
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
- A von Baeyer **bicyclic** ring's own "-ylidene" name instead needs a
  locant (its own numbering isn't symmetric the way a monocyclic ring's
  is, e.g. `bicyclo[2.2.1]heptan-2-ylidene` vs `...heptan-7-ylidene`
  aren't the same position) -- but per the worked example above, that
  locant is cited *outside* the parenthesized name, at the front of the
  whole assembly, exactly where the monocyclic case's own always-"1"
  locant already sits (`2,2'-bi(bicyclo[2.2.1]heptanylidene)`, not
  `bi(bicyclo[2.2.1]heptan-2-ylidene)`), so the parenthesized name itself
  stays locant-free either way (`"bicyclo[2.2.1]heptanylidene"`, elided
  same as the monocyclic "-idene" case). This module finds that locant by
  severing the C=C bond and isolating each side into its own single-ring-
  system RDKit `Mol` (`_isolate_ring_system`), then reuses
  `_bicyclic.py`'s existing `find_bicyclic_core`/`iter_bicyclic_numberings`/
  `bicyclic_parent_name` (already proven end to end by `_polycyclic_suffix
  .name_von_baeyer_suffix`, which every other von Baeyer suffix module in
  this codebase already goes through) to pick, among every von Baeyer-
  valid numbering, whichever gives the junction atom its lowest possible
  locant.
- P-14.3.2 (numbering): each ring is numbered independently, starting at
  its own double-bond-junction atom (locant 1), with one ring's locants
  left unprimed and the other's primed -- the identical primed/unprimed
  lowest-locant-set numbering `_ring_assembly.py` already uses for its own
  single-bond-junction (P-28.2.1) case, reused here via its exported
  `_numberings_from_attachment` building block (the double-bond-junction
  atom plays the same "locant 1" role its single-bond attachment atom
  does there) for the monocyclic case; the bicyclic case's own primed/
  unprimed pair is trivial since both sides are structurally identical
  and each independently minimizes to the same locant value.
- P-35.2.1 (Chapter P-3): halogen substituents hang off a ring atom the
  same way as in every other ring module -- for the monocyclic case only;
  see Scope below for why the bicyclic case doesn't support one yet.

Scope: two disjoint, identical rings or ring systems, each contributing
exactly one atom to a C=C double bond between them (no other inter-ring
bond) -- either (a) monocyclic saturated all-carbon rings of the same
size, optionally bearing halogen substituents (positions may differ
between the two rings, exactly as `_ring_assembly.py`'s own "identical
parent, independently substituted" biphenyl case already allows), or (b)
von Baeyer bicyclic all-carbon ring systems with the same bridge-length
pattern (e.g. two `bicyclo[2.2.1]heptane`s), **unsubstituted** -- a
halogen substituent on a bicyclic side is deferred to a follow-up (it
needs each side's own name built with substituent locants primed
alongside the suffix locant, which isn't needed at all when there are no
substituents to prime, see the module's own history for why this was
scoped out narrowly rather than solved in the same pass). Explicitly out
of scope (raise `UnsupportedStructure` via the generic fallback in
`core.py`, since `find_ring_assembly_ylidene_core` below simply returns
None for any of these): a polycyclic (3+ ring) or monospiro parent ring
on either side (each has its own numbering-candidate source module and
isn't wired in here yet), an aromatic ring, two different ring systems,
more than 2 rings, a halogen substituent on a bicyclic side, and
indicated hydrogen (P-28.2.3 -- not reachable by any all-carbon saturated
ring anyway).
"""

from rdkit import Chem

from ._bicyclic import bicyclic_parent_name, find_bicyclic_core, iter_bicyclic_numberings
from ._common import (
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


def _isolate_ring_system(mol, cross_bond, keep_atom):
    """Sever `cross_bond` (the C=C double bond joining the two ring
    systems) and return the isolated, sanitized `Mol` fragment containing
    `keep_atom`, plus `keep_atom`'s own index within that fragment.
    `find_bicyclic_core`/`iter_bicyclic_numberings` derive ring count from
    the whole passed-in molecule's own graph (`adjacency(mol)`), so they
    can't be run on the two-ring assembly directly (its combined cyclomatic
    number is double a genuine single bicyclic core's) -- each side has to
    be extracted as its own standalone molecule first."""
    rw = Chem.RWMol(mol)
    a, b = cross_bond.GetBeginAtomIdx(), cross_bond.GetEndAtomIdx()
    rw.RemoveBond(a, b)
    for idx in (a, b):
        atom = rw.GetAtomWithIdx(idx)
        atom.SetNoImplicit(False)
        atom.SetNumExplicitHs(0)
    orig_idx_frags = Chem.GetMolFrags(rw, sanitizeFrags=False)
    mol_frags = Chem.GetMolFrags(rw, asMols=True, sanitizeFrags=False)
    for orig_idxs, frag in zip(orig_idx_frags, mol_frags):
        if keep_atom in orig_idxs:
            Chem.SanitizeMol(frag)
            return frag, orig_idxs.index(keep_atom)
    raise AssertionError("keep_atom must belong to one of the two severed fragments")


def find_ring_assembly_ylidene_core(mol):
    """Return (ring0_atoms, ring1_atoms, junction0, junction1, kind) if
    `mol` is exactly two disjoint, identical rings/ring systems joined by
    one C=C double bond (one atom from each) -- `kind` is `("monocyclic",
    size)` for two identical-size monocyclic saturated all-carbon rings,
    or `("bicyclic", bicyclic_parent_name)` for two von Baeyer bicyclic
    all-carbon ring systems with the same bridge-length pattern -- else
    None. `ring0_atoms`/`ring1_atoms` are each side's *entire* atom set
    (both rings, for the bicyclic case). Ring-carbon-only/non-aromatic is
    checked per ring below (not globally here), since a halogen substituent
    hanging off a ring atom is in scope for the monocyclic case."""
    double_bonds = [bond for bond in mol.GetBonds() if not bond.GetIsAromatic() and bond.GetBondTypeAsDouble() == 2.0]
    if len(double_bonds) != 1:
        return None
    bond = double_bonds[0]
    x, y = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
    for atom in (mol.GetAtomWithIdx(x), mol.GetAtomWithIdx(y)):
        if atom.GetDegree() != 3:
            return None
    # The two ring systems are connected to each other by `bond` itself, so
    # `Chem.GetMolFrags(mol)` sees one connected molecule, not two -- side
    # membership has to be found by graph traversal with that one edge
    # excluded instead.
    graph = adjacency(mol)
    side0_atoms = {x}
    frontier = [x]
    while frontier:
        current = frontier.pop()
        for neighbor in graph[current]:
            if neighbor == y and current == x:
                continue
            if neighbor not in side0_atoms:
                side0_atoms.add(neighbor)
                frontier.append(neighbor)
    side1_atoms = set(range(mol.GetNumAtoms())) - side0_atoms
    if y not in side1_atoms:
        return None
    junction0, junction1 = x, y

    ring_info = mol.GetRingInfo()
    rings0 = [r for r in ring_info.AtomRings() if set(r) <= side0_atoms]
    rings1 = [r for r in ring_info.AtomRings() if set(r) <= side1_atoms]

    if len(rings0) == 1 and len(rings1) == 1 and len(rings0[0]) == len(rings1[0]):
        if junction0 not in rings0[0] or junction1 not in rings1[0]:
            return None
        for ring in (rings0[0], rings1[0]):
            ring_set = set(ring)
            for idx in ring:
                atom = mol.GetAtomWithIdx(idx)
                if atom.GetIsAromatic() or atom.GetAtomicNum() != 6:
                    return None
            for bond2 in mol.GetBonds():
                a, b = bond2.GetBeginAtomIdx(), bond2.GetEndAtomIdx()
                if a in ring_set and b in ring_set and bond2.GetBondTypeAsDouble() != 1.0:
                    return None
        return set(rings0[0]), set(rings1[0]), junction0, junction1, ("monocyclic", len(rings0[0]))

    if len(rings0) == 2 and len(rings1) == 2:
        for idx in side0_atoms | side1_atoms:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetIsAromatic() or atom.GetAtomicNum() != 6:
                return None
        frag0, _ = _isolate_ring_system(mol, bond, junction0)
        frag1, _ = _isolate_ring_system(mol, bond, junction1)
        core0 = find_bicyclic_core(frag0)
        core1 = find_bicyclic_core(frag1)
        if core0 is None or core1 is None:
            return None
        parent0, parent1 = bicyclic_parent_name(core0), bicyclic_parent_name(core1)
        if parent0 != parent1:
            return None
        return side0_atoms, side1_atoms, junction0, junction1, ("bicyclic", parent0)

    return None


def _ylidene_name(ring_size: int) -> str:
    return "cyclo" + alkyl_name(ring_size) + "idene"


def _bicyclic_ylidene_name(parent_name: str) -> str:
    return parent_name[:-1] + "ylidene"


def _bicyclic_junction_locant(mol, cross_bond, junction) -> int:
    frag, local_idx = _isolate_ring_system(mol, cross_bond, junction)
    core = find_bicyclic_core(frag)
    return min(full_order.index(local_idx) + 1 for full_order in iter_bicyclic_numberings(core))


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


def _name_monocyclic_ring_assembly_ylidene(mol, ring0_atoms, ring1_atoms, junction0, junction1, ring_size) -> str:
    validate_atoms_and_bonds(mol)

    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    ring_atoms = ring0_atoms | ring1_atoms
    ylidene_name = _ylidene_name(ring_size)

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


def _name_bicyclic_ring_assembly_ylidene(mol, junction0, junction1, parent_name) -> str:
    bond = next(
        b
        for b in mol.GetBonds()
        if not b.GetIsAromatic() and b.GetBondTypeAsDouble() == 2.0
    )
    locant0 = _bicyclic_junction_locant(mol, bond, junction0)
    locant1 = _bicyclic_junction_locant(mol, bond, junction1)
    lower, higher = sorted((locant0, locant1))
    return f"{lower},{higher}'-bi({_bicyclic_ylidene_name(parent_name)})"


def name_ring_assembly_ylidene(mol, core) -> str:
    ring0_atoms, ring1_atoms, junction0, junction1, (kind, detail) = core
    if kind == "monocyclic":
        return _name_monocyclic_ring_assembly_ylidene(mol, ring0_atoms, ring1_atoms, junction0, junction1, detail)
    return _name_bicyclic_ring_assembly_ylidene(mol, junction0, junction1, detail)
