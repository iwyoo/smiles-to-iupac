"""Naming of simple Group 13 (aluminium, gallium, indium, thallium), Group
14 (germanium, tin, lead), and Group 15 (arsenic, antimony, bismuth)
mononuclear organometallics (a single metal atom bearing 1-3 Group 13/15 /
1-4 Group 14 alkyl/phenyl substituents), per the IUPAC 2013 Recommendations
("the Blue Book"):

- P-69.1 (Chapter P-6a, https://iupac.qmul.ac.uk/BlueBook/PDF/P6a.pdf):
  organometallic compounds of Groups 13-16 are named substitutively, the
  same way `_borane.py`/`_phosphane.py` already name boron/phosphorus --
  the metal atom itself is the mononuclear parent hydride (P-68's own
  element table), with substituents cited as prefixes via
  `format_mononuclear_prefixes`. Confirmed worked examples
  (`tmp/bluebook/P6a.txt` lines 8602-8613): `Al(CH2-CH3)3` ->
  'triethylalumane', `Pb(CH2-CH3)4` -> 'tetraethylplumbane',
  `BrSb(CH=CH2)2` -> 'bromodi(ethenyl)stibane', `HIn(CH3)2` ->
  'dimethylindigane'. Thallium's own stem name, 'thallane', is confirmed
  at `tmp/bluebook/P6a.txt` line 5459 (preselected name). The Group 15
  stems themselves are confirmed separately: `ethylarsane (PIN)` (line
  7810), `trimethylbismuthane (PIN)` (line 8010); `stibane` is confirmed
  by the `bromodi(ethenyl)stibane` example above.
- An unsaturated acyclic substituent (e.g. the `bromodi(ethenyl)stibane`
  example's vinyl groups) is named the same way `name_branch` already
  names one anywhere else in this project (confirmed by real PubChem
  structures: CID 23271262 "tris(ethenyl)arsane", CID 81998
  "tributyl(ethenyl)stannane") -- only a multiple bond *directly to the
  metal itself* (a fundamentally different, ylidene-shaped bonding
  pattern this module was never built for) is rejected.
- This module is structurally identical to `_borane.py` with the boron
  atomic number swapped for one of Al/Ga/In/Tl/Ge/Sn/Pb/As/Sb/Bi -- same
  validation shape (plain/halogenated-phenyl support, halogen-on-metal
  support), same `format_mononuclear_prefixes` assembly, differing only
  in the stems dict (atomic number -> parent-hydride name) and the
  maximum substituent count (3 for the trivalent Group 13/15 elements, 4
  for the tetravalent Group 14 ones) -- generalized here as ONE shared
  mechanism parameterized by both, rather than a separate near-duplicate
  module per group (this milestone's own generalization check, applied a
  third time: first to unify Al/Ga/In/Tl into one module instead of four,
  then to share that same module's logic with Ge/Sn/Pb, now with As/Sb/Bi
  too -- each addition purely a stems dict + a valence number, no new
  logic).

Scope: exactly the same as `_borane.py`'s (see that module's own
docstring for the full derivation of every rule below, all unchanged with
the metal swapped in): a plain/branched alkyl substituent, a plain phenyl
substituent, a halogen-substituted-phenyl substituent, a halogen bonded
directly to the metal, up to `max_substituents` substituents total (3 for
Group 13/15, 4 for Group 14 -- never more, these are all trivalent/
tetravalent respectively, like boron/carbon), mixed freely except a
halogenated-phenyl group may not mix with a differently-named substituent.

Explicitly out of scope (raise `UnsupportedStructure`): more than one
metal atom, a multiple bond directly to the metal, a non-phenyl aromatic
ring, charged/isotopically modified atoms, coexistence with any other
heteroatom including another Group 13/14/15/16 element, plus silicon/
carbon (organosilicon naming is a separate, already-established area of
this project, `_silane_chain.py`, not part of this mononuclear-
organometallic mechanism) and Group 16 elements (a different scope, a
separate milestone).
"""

from rdkit import Chem

from ._common import (
    HALOGEN_PREFIXES,
    UnsupportedStructure,
    adjacency,
    carbon_adjacency,
    halogen_substituents,
    non_single_bonds,
    plain_phenyl_substituent_atoms,
)
from ._substituents import format_mononuclear_prefixes, halogenated_phenyl_substituent, name_branch

GROUP_13_STEMS = {
    13: "alumane",
    31: "gallane",
    49: "indigane",
    81: "thallane",
}

GROUP_14_STEMS = {
    32: "germane",
    50: "stannane",
    82: "plumbane",
}

GROUP_15_STEMS = {
    33: "arsane",
    51: "stibane",
    83: "bismuthane",
}


def _branch_atoms(graph, root, boundary):
    """Every atom reachable from `root` without crossing into `boundary`
    (here, just the metal atom) -- mirrors `_aromatic.py`'s identical
    local helper, used the same way: to isolate one substituent's own
    branch so an unsaturated bond can be attributed to the branch that
    contains it."""
    seen = {root}
    stack = [root]
    while stack:
        current = stack.pop()
        for neighbor in graph[current]:
            if neighbor in boundary or neighbor in seen:
                continue
            seen.add(neighbor)
            stack.append(neighbor)
    return seen


def has_group13_hydride_shape(mol) -> bool:
    return any(atom.GetAtomicNum() in GROUP_13_STEMS for atom in mol.GetAtoms())


def has_group14_hydride_shape(mol) -> bool:
    return any(atom.GetAtomicNum() in GROUP_14_STEMS for atom in mol.GetAtoms())


def has_group15_hydride_shape(mol) -> bool:
    return any(atom.GetAtomicNum() in GROUP_15_STEMS for atom in mol.GetAtoms())


def _validate_and_collect_substituents(mol, metal, stems, max_substituents):
    metal_atomic_num = metal.GetAtomicNum()
    if metal.GetFormalCharge() != 0 or metal.GetIsotope() != 0:
        raise UnsupportedStructure("charged or isotopically modified atoms are not supported yet")
    if metal.GetDegree() > max_substituents:
        raise UnsupportedStructure(
            f"a {stems[metal_atomic_num]}-forming atom with more than "
            f"{max_substituents} substituents is not supported"
        )

    graph = adjacency(mol)
    roots = set(graph[metal.GetIdx()])
    phenyl_atoms = plain_phenyl_substituent_atoms(mol, graph, roots)

    halogens = halogen_substituents(mol)
    aromatic_atoms = {atom.GetIdx() for atom in mol.GetAtoms() if atom.GetIsAromatic()}
    halophenyl = {}
    for root in roots - phenyl_atoms:
        result = halogenated_phenyl_substituent(graph, aromatic_atoms, root, metal.GetIdx(), halogens)
        if result is not None:
            halophenyl[root] = result
    halophenyl_ring_atoms = {a for _, ring_atoms, _ in halophenyl.values() for a in ring_atoms}
    halophenyl_halogen_atoms = {a for _, _, halogen_atoms in halophenyl.values() for a in halogen_atoms}

    for atom in mol.GetAtoms():
        idx = atom.GetIdx()
        is_halogen = atom.GetAtomicNum() in HALOGEN_PREFIXES
        if atom.GetAtomicNum() not in (metal_atomic_num, 6) and not is_halogen:
            raise UnsupportedStructure(
                f"heteroatoms other than the {stems[metal_atomic_num]}'s "
                "own metal atom are not supported yet (see P-69.1)"
            )
        if is_halogen and idx in halophenyl_halogen_atoms:
            continue
        if is_halogen and (atom.GetDegree() != 1 or atom.GetNeighbors()[0].GetIdx() != metal.GetIdx()):
            raise UnsupportedStructure(
                "a halogen-substituted alkyl chain is out of scope for this "
                "module -- only a halogen bonded directly to the metal, or to "
                "a phenyl ring directly on the metal, is supported so far"
            )
        if atom.GetFormalCharge() != 0 or atom.GetIsotope() != 0:
            raise UnsupportedStructure("charged or isotopically modified atoms are not supported yet")
        if (
            atom.GetAtomicNum() == 6
            and atom.GetIsAromatic()
            and idx not in phenyl_atoms
            and idx not in halophenyl_ring_atoms
        ):
            raise UnsupportedStructure(
                "an aromatic substituent other than a plain or halogen-"
                "substituted phenyl group is out of scope for this module"
            )
    all_ring_atoms = {a for ring in mol.GetRingInfo().AtomRings() for a in ring}
    if all_ring_atoms - phenyl_atoms - halophenyl_ring_atoms:
        raise UnsupportedStructure(
            "a ring other than a plain or halogen-substituted phenyl "
            "substituent directly on the metal is out of scope for this "
            "module"
        )
    metal_idx = metal.GetIdx()
    if any(b[0] == metal_idx or b[1] == metal_idx for b in non_single_bonds(mol)):
        raise UnsupportedStructure("a multiple bond directly to the metal is out of scope for this module")
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")

    unsaturated_bonds = [b for b in non_single_bonds(mol) if b[0] != metal_idx and b[1] != metal_idx]
    carbon_graph = carbon_adjacency(mol) if unsaturated_bonds else None

    substituent_names = []
    for root in roots:
        if root in phenyl_atoms:
            substituent_names.append(("phenyl", False))
            continue
        if root in halophenyl:
            name, _, _ = halophenyl[root]
            substituent_names.append((name, True))
            continue
        root_atomic_num = mol.GetAtomWithIdx(root).GetAtomicNum()
        if root_atomic_num in HALOGEN_PREFIXES:
            substituent_names.append((HALOGEN_PREFIXES[root_atomic_num], False))
            continue
        if unsaturated_bonds:
            branch_atoms = _branch_atoms(graph, root, {metal_idx})
            branch_bonds = [b for b in unsaturated_bonds if b[0] in branch_atoms and b[1] in branch_atoms]
            if branch_bonds:
                # The Blue Book's worked example cites it as `bromodi(ethenyl)stibane (PIN)`.
                name, _ = name_branch(graph, root, metal_idx, {}, mol=mol)
                substituent_names.append((name, True))
                continue
        substituent_names.append(name_branch(graph, root, metal.GetIdx(), {}, mol=mol))

    distinct_names = {name for name, _ in substituent_names}
    if len(distinct_names) > 1 and any(
        is_compound and name[0].isdigit() for name, is_compound in substituent_names
    ):
        raise UnsupportedStructure(
            "a halogen-substituted phenyl group mixed with any differently-"
            "named substituent is out of scope for this module (mirrors "
            "`_borane.py`'s identical restriction)"
        )
    return substituent_names


def _name_mononuclear_hydride(mol, stems, max_substituents, error_label) -> str:
    metal_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() in stems]
    if len(metal_atoms) != 1:
        raise UnsupportedStructure(f"more than one {error_label} metal atom is not supported yet")
    (metal,) = metal_atoms
    stem = stems[metal.GetAtomicNum()]

    substituent_names = _validate_and_collect_substituents(mol, metal, stems, max_substituents)
    if not substituent_names:
        return stem
    return format_mononuclear_prefixes(substituent_names) + stem


def name_group13_hydride(mol) -> str:
    return _name_mononuclear_hydride(mol, GROUP_13_STEMS, 3, "Group 13")


def name_group14_hydride(mol) -> str:
    return _name_mononuclear_hydride(mol, GROUP_14_STEMS, 4, "Group 14")


def name_group15_hydride(mol) -> str:
    return _name_mononuclear_hydride(mol, GROUP_15_STEMS, 3, "Group 15")
