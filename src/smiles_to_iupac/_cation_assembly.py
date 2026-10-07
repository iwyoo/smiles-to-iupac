"""Ring assemblies of two identical cationic ring systems joined by one single bond (P-28.2.1, P-73.5.1.1, P-73.5.1.3):
'1,1'-dimethyl[4,4'-bipyridine]-1,1'-diium', "2,2',5,5'-tetraoxo[3,3'-bipyrrolidine]-1,1'-diium". One unit is numbered
as a cation with the other unit replaced by a probe atom, so the junction takes its locant with the cationic centres
ahead of the prefixes; the same numbering is applied to the primed unit.
"""

import re

from rdkit import Chem

from ._common import (
    UnsupportedStructure,
    adjacency,
    group_substituents,
    halogen_substituents,
)
from ._multiplicative import _fragment_key, _ring_systems
from ._multiplicative_text import PrimedLocant
from ._substituents import format_substituent_prefixes, name_branch

_COUNT = {2: "di", 4: "tetra", 6: "hexa"}


class NotAnAssembly(UnsupportedStructure):
    pass


def _units(mol, systems):
    if len(systems) != 2:
        raise NotAnAssembly("not a pair of ring systems")
    joins = [
        b
        for b in mol.GetBonds()
        if (b.GetBeginAtomIdx() in systems[0] and b.GetEndAtomIdx() in systems[1])
        or (b.GetBeginAtomIdx() in systems[1] and b.GetEndAtomIdx() in systems[0])
    ]
    if len(joins) != 1 or joins[0].GetBondTypeAsDouble() != 1.0 or joins[0].GetIsAromatic():
        raise NotAnAssembly("the ring systems are not joined by one single bond")
    bond = joins[0]
    first, second = (
        (bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
        if bond.GetBeginAtomIdx() in systems[0]
        else (bond.GetEndAtomIdx(), bond.GetBeginAtomIdx())
    )
    units = []
    for start, other in ((first, second), (second, first)):
        atoms, stack = {start}, [start]
        while stack:
            for n in mol.GetAtomWithIdx(stack.pop()).GetNeighbors():
                if n.GetIdx() not in atoms and n.GetIdx() != other:
                    atoms.add(n.GetIdx())
                    stack.append(n.GetIdx())
        units.append((start, other, atoms))
    if _fragment_key(mol, units[0][2], first) != _fragment_key(mol, units[1][2], second):
        raise UnsupportedStructure("the ring systems of an assembly must be identical")
    return units


def _probe(mol, unit, centres):
    start, other, atoms = unit
    editable = Chem.RWMol(mol)
    for index in centres:
        atom = editable.GetAtomWithIdx(index)
        hydrogens = atom.GetTotalNumHs()
        atom.SetFormalCharge(0)
        atom.SetNoImplicit(True)
        atom.SetNumExplicitHs(max(hydrogens - 1, 0))
        if atom.GetIsAromatic():
            atom.SetBoolProp("_ring_cation_centre", True)
    marker = editable.AddAtom(Chem.Atom(53))
    editable.RemoveBond(start, other)
    editable.AddBond(start, marker, Chem.BondType.SINGLE)
    keep = set(atoms) | {marker}
    for index in sorted(set(range(editable.GetNumAtoms())) - keep, reverse=True):
        editable.RemoveAtom(index)
    order = sorted(keep)
    probe = editable.GetMol()
    probe.UpdatePropertyCache(strict=False)
    Chem.FastFindRings(probe)
    return probe, {old: new for new, old in enumerate(order)}


def name_cation_ring_assembly(mol) -> str:
    from ._diester_ring_diyl import _system_of, evaluate_skeleton
    from ._functional_prefixes import functional_names

    centres = [a.GetIdx() for a in mol.GetAtoms() if a.GetFormalCharge()]
    if len(Chem.GetMolFrags(mol)) != 1 or any(mol.GetAtomWithIdx(c).GetFormalCharge() != 1 for c in centres):
        raise NotAnAssembly("an assembly needs one fragment of singly charged centres")
    systems = _ring_systems(mol)
    if any(not any(c in s for s in systems) for c in centres):
        raise NotAnAssembly("a cationic centre outside the ring systems")
    unit = _units(mol, systems)[0]
    unit_centres = [c for c in centres if c in unit[2]]
    if not unit_centres or len(unit_centres) * 2 != len(centres):
        raise UnsupportedStructure("the two ring systems must carry the same cationic centres")
    probe, mapping = _probe(mol, unit, unit_centres)
    centre_atoms = [mapping[c] for c in unit_centres]
    junction = mapping[unit[0]]
    rings, atoms = _system_of(probe, junction)
    graph = adjacency(probe)
    found = evaluate_skeleton(probe, graph, "ring", rings, atoms, centre_atoms, set(), "ium")
    if found is None:
        raise UnsupportedStructure("this cationic ring system has no supported name yet")
    position_of = found[2]
    marker = probe.GetNumAtoms() - 1
    halogens = halogen_substituents(probe)
    seeds = [(a, n) for a in atoms for n in graph[a] if n not in atoms and n != marker]
    named, shown, _ = functional_names(probe, graph, seeds, set(atoms) | {marker}, halogens)
    entries = {}
    for atom, root in seeds:
        entry = named[root] if root in named else name_branch(graph, root, atom, shown, mol=probe)
        for primes in (0, 1):
            entries.setdefault(PrimedLocant(primes, position_of[atom]), []).append(entry)
    prefix = format_substituent_prefixes(group_substituents(entries)) if entries else ""
    parent = found[1]
    iodo_prefix = format_substituent_prefixes(
        group_substituents(
            {
                position_of[a]: [named[r] if r in named else name_branch(graph, r, a, shown, mol=probe)]
                for a, r in [(a, n) for a in atoms for n in graph[a] if n not in atoms]
            }
        )
    )
    tail = parent[len(iodo_prefix):].lstrip("-") if parent.startswith(iodo_prefix) else None
    match = re.fullmatch(r"(?P<stem>[a-z]+)-(?P<locants>\d+(?:,\d+)*)-(?P<multiplier>di|tri|tetra)?ium", tail or "")
    if match is None or match.group("stem").startswith(("cyclo",)):
        raise UnsupportedStructure("the indicated hydrogen or hydro prefixes of this cation are not supported in an assembly yet")
    hydride = match.group("stem") + ("" if match.group("multiplier") else "e")
    join = position_of[junction]
    centre_locants = sorted(
        [PrimedLocant(p, position_of[c]) for c in centre_atoms for p in (0, 1)]
    )
    suffix_count = _COUNT.get(len(centre_locants))
    if suffix_count is None:
        raise UnsupportedStructure("this number of cationic centres in an assembly is not supported yet")
    cited = ",".join(str(loc) for loc in centre_locants)
    return f"{prefix}[{join},{join}'-bi{hydride}]-{cited}-{suffix_count}ium"
