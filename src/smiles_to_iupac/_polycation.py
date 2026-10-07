"""Several identical 'ium' centres on one parent hydride (P-73.1.1.2): the final 'e' of the parent hydride name is kept
before 'di', 'tri' (an 'e' is lost only before a vowel) and the locants are lowest for the cationic centres, ahead of
the prefixes (P-31.1.4.2.4): 1,4-dioxane-1,4-diium, tetramethyldiazene-1,2-diium."""

import re

from rdkit import Chem

from ._common import (
    UnsupportedStructure,
    adjacency,
    group_substituents,
    halogen_substituents,
    substituent_locant_set_and_citation,
)

class _DifferentSystems(UnsupportedStructure):
    pass


class _NotAPair(UnsupportedStructure):
    pass


_PAIR_PARENTS = {7: ("hydrazine", "diazene"), 15: ("diphosphane", None), 16: ("disulfane", None), 8: ("dioxidane", None)}
_MULTIPLIED = {2: "di", 3: "tri", 4: "tetra"}


def _centres(mol):
    return [a for a in mol.GetAtoms() if a.GetFormalCharge()]


def has_polycation_shape(mol) -> bool:
    centres = _centres(mol)
    return (
        bool(centres)
        and sum(a.GetFormalCharge() for a in centres) >= 2
        and len(Chem.GetMolFrags(mol)) == 1
        and all(a.GetFormalCharge() > 0 and not a.GetIsotope() for a in centres)
        and not any(a.GetNumRadicalElectrons() for a in mol.GetAtoms())
        and (len(centres) > 1 or centres[0].GetAtomicNum() == 6)
    )


def name_polycation(mol) -> str:
    from ._multiplicative_cation import name_cation_assembly

    centres = _centres(mol)
    if all(a.GetAtomicNum() == 6 for a in centres):
        try:
            return _name_polycarbenium(mol, centres)
        except _NoSkeletonName:
            if len(centres) > 1:
                return name_cation_assembly(mol)
            raise
    if any(a.GetFormalCharge() != 1 for a in centres):
        raise UnsupportedStructure("a multiply charged heteroatom centre is not supported yet")
    if all(a.IsInRing() for a in centres) and not any(a.GetAtomicNum() == 6 for a in centres):
        try:
            return _name_ring_polycation(mol, centres)
        except _DifferentSystems:
            from ._cation_assembly import NotAnAssembly, name_cation_ring_assembly

            try:
                return name_cation_ring_assembly(mol)
            except NotAnAssembly:
                pass
    elif len(centres) == 2 and not any(a.IsInRing() or a.GetAtomicNum() == 6 for a in centres):
        try:
            return _name_pair_polycation(mol, centres)
        except _NotAPair:
            pass
    return name_cation_assembly(mol)


class _NoSkeletonName(UnsupportedStructure):
    pass


_DIIDE = re.compile(r"(?P<head>.+?)-(?P<locants>\d+(?:,\d+)*)-(?P<multiplier>di|tri|tetra)ide$")
_BIS = {"di": "bis", "tri": "tris", "tetra": "tetrakis"}


def _name_polycarbenium(mol, centres):
    """Several carbon centres that lost hydride ions on one parent hydride are cited as 'bis(ylium)' (P-73.2.2.1.2): the
    carbanion analogue (same bonds, the charge reversed) is named by the anion namer and its 'ide' ending is swapped."""
    from ._anion import name_anion

    analogue = Chem.RWMol(mol)
    for atom in centres:
        analogue.GetAtomWithIdx(atom.GetIdx()).SetFormalCharge(-atom.GetFormalCharge())
    analogue = analogue.GetMol()
    try:
        Chem.SanitizeMol(analogue)
        name = name_anion(analogue)
    except UnsupportedStructure as error:
        raise _NoSkeletonName(str(error)) from error
    match = _DIIDE.match(name)
    if match is None:
        raise _NoSkeletonName("the carbanion analogue has no multiplied 'ide' name")
    return f"{match.group('head')}-{match.group('locants')}-{_BIS[match.group('multiplier')]}(ylium)"


def _name_ring_polycation(mol, centres):
    from ._diester_ring_diyl import _system_of, evaluate_skeleton

    if any(a.GetIsotope() or a.GetAtomicNum() not in _RING_CENTRE_ELEMENTS for a in centres):
        raise UnsupportedStructure("this ring heteroatom is not supported as a cationic centre yet")
    editable = Chem.RWMol(mol)
    indices = [a.GetIdx() for a in centres]
    for index in indices:
        atom = editable.GetAtomWithIdx(index)
        hydrogens = atom.GetTotalNumHs()
        atom.SetFormalCharge(0)
        atom.SetNoImplicit(True)
        atom.SetNumExplicitHs(max(hydrogens - 1, 0))
        if atom.GetIsAromatic():
            atom.SetBoolProp("_ring_cation_centre", True)
    base = editable.GetMol()
    base.UpdatePropertyCache(strict=False)
    Chem.FastFindRings(base)
    rings, atoms = _system_of(base, indices[0])
    if any(i not in atoms for i in indices):
        raise _DifferentSystems("cationic centres in different ring systems are named multiplicatively")
    found = evaluate_skeleton(base, adjacency(base), "ring", rings, atoms, indices, set(), "ium")
    if found is None:
        raise UnsupportedStructure("this cationic ring system has no supported name yet")
    return found[1]


_RING_CENTRE_ELEMENTS = {7, 8, 15, 16, 33, 34, 52}


def _name_pair_polycation(mol, centres):
    """Two adjacent identical heteroatoms (hydrazine, diazene, diphosphane, disulfane, dioxidane) that both carry the
    cationic centre."""
    from ._substituents import format_substituent_prefixes, name_branch

    first, second = centres
    element = first.GetAtomicNum()
    bond = mol.GetBondBetweenAtoms(first.GetIdx(), second.GetIdx())
    if element not in _PAIR_PARENTS or second.GetAtomicNum() != element or bond is None:
        raise _NotAPair("the cationic centres are not an identical adjacent heteroatom pair")
    single, double = _PAIR_PARENTS[element]
    parent = {1.0: single, 2.0: double}.get(bond.GetBondTypeAsDouble())
    if parent is None:
        raise UnsupportedStructure("this bond between the cationic centres has no parent hydride name")
    pair = (first.GetIdx(), second.GetIdx())
    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    options = []
    for numbering in (pair, pair[::-1]):
        entries = {}
        for position, atom in enumerate(numbering, start=1):
            for n in graph[atom]:
                if n not in pair:
                    entries.setdefault(position, []).append(name_branch(graph, n, atom, halogens, mol=mol))
        grouped = group_substituents(entries)
        locant_set, _, citation = substituent_locant_set_and_citation(grouped)
        options.append((locant_set, citation, format_substituent_prefixes(grouped)))
    _, _, prefix = min(options, key=lambda o: (o[0], o[1]))
    valence = 4 if element in (7, 15) else 3
    free = 2 * (valence - int(bond.GetBondTypeAsDouble()))
    used = sum(
        int(mol.GetBondBetweenAtoms(a, n).GetBondTypeAsDouble()) for a in pair for n in graph[a] if n not in pair
    )
    if used == free:
        grouped = group_substituents(
            {position: [name_branch(graph, n, atom, halogens, mol=mol) for n in graph[atom] if n not in pair]
             for position, atom in enumerate(pair, start=1)}
        )
        prefix = format_substituent_prefixes(grouped, omit_locants=True)
    return f"{prefix}{parent}-1,2-diium"
