"""Naming of simple Group 13 organometallics (a single aluminium, gallium,
indium, or thallium atom bearing 1-3 alkyl/phenyl substituents), per the
IUPAC 2013 Recommendations ("the Blue Book"):

- P-69.1 (Chapter P-6a, https://iupac.qmul.ac.uk/BlueBook/PDF/P6a.pdf):
  organometallic compounds of Groups 13-16 are named substitutively, the
  same way `_borane.py`/`_phosphane.py` already name boron/phosphorus --
  the metal atom itself is the mononuclear parent hydride (P-68's own
  element table), with substituents cited as prefixes via
  `format_mononuclear_prefixes`. Confirmed worked examples
  (`tmp/bluebook/P6a.txt` lines 8602-8613): `Al(CH2-CH3)3` ->
  'triethylalumane', `Pb(CH2-CH3)4` -> 'tetraethylplumbane' (Group 14,
  out of scope here), `BrSb(CH=CH2)2` -> 'bromodi(ethenyl)stibane'
  (Group 15, out of scope here), `HIn(CH3)2` -> 'dimethylindigane'.
  Thallium's own stem name, 'thallane', is confirmed at
  `tmp/bluebook/P6a.txt` line 5459 (preselected name).
- This module is structurally identical to `_borane.py` with the boron
  atomic number swapped for one of Al/Ga/In/Tl (each a Group 13, trivalent
  element like boron itself) -- same validation shape (max 3
  substituents, plain/halogenated-phenyl support, halogen-on-metal
  support), same `format_mononuclear_prefixes` assembly, generalized here
  as ONE shared mechanism parameterized by atomic number rather than 4
  separate near-duplicate modules (this milestone's own generalization
  check).

Scope: exactly the same as `_borane.py`'s (see that module's own
docstring for the full derivation of every rule below, all unchanged with
the metal swapped in): a plain/branched alkyl substituent, a plain phenyl
substituent, a halogen-substituted-phenyl substituent, a halogen bonded
directly to the metal, up to 3 substituents total (never more -- these
are all trivalent, like boron), mixed freely except a halogenated-phenyl
group may not mix with a differently-named substituent.

Explicitly out of scope (raise `UnsupportedStructure`): any of
`_borane.py`'s own out-of-scope cases (more than one metal atom, an
unsaturated substituent, a non-phenyl aromatic ring, charged/isotopically
modified atoms, coexistence with any other heteroatom including another
Group 13/14/15/16 element), plus Group 14 elements (Ge/Sn/Pb -- a
different valence, a separate milestone step).
"""

from rdkit import Chem

from ._common import (
    HALOGEN_PREFIXES,
    UnsupportedStructure,
    adjacency,
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


def has_group13_hydride_shape(mol) -> bool:
    return any(atom.GetAtomicNum() in GROUP_13_STEMS for atom in mol.GetAtoms())


def _validate_and_collect_substituents(mol, metal):
    metal_atomic_num = metal.GetAtomicNum()
    if metal.GetFormalCharge() != 0 or metal.GetIsotope() != 0:
        raise UnsupportedStructure("charged or isotopically modified atoms are not supported yet")
    if metal.GetDegree() > 3:
        raise UnsupportedStructure(
            f"a {GROUP_13_STEMS[metal_atomic_num]}-forming atom with more than "
            "three substituents is not supported"
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
                f"heteroatoms other than the {GROUP_13_STEMS[metal_atomic_num]}'s "
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
    non_ring_unsaturation = [
        b
        for b in non_single_bonds(mol)
        if b[0] not in phenyl_atoms
        and b[1] not in phenyl_atoms
        and b[0] not in halophenyl_ring_atoms
        and b[1] not in halophenyl_ring_atoms
    ]
    if non_ring_unsaturation:
        raise UnsupportedStructure("an unsaturated substituent is out of scope for this module (see P-69.1)")
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")

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


def name_group13_hydride(mol) -> str:
    metal_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() in GROUP_13_STEMS]
    if len(metal_atoms) != 1:
        raise UnsupportedStructure("more than one Group 13 metal atom is not supported yet")
    (metal,) = metal_atoms
    stem = GROUP_13_STEMS[metal.GetAtomicNum()]

    substituent_names = _validate_and_collect_substituents(mol, metal)
    if not substituent_names:
        return stem
    return format_mononuclear_prefixes(substituent_names) + stem
