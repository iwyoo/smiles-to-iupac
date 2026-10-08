"""Several identical 'ium' centres on one parent hydride (P-73.1.1.2): the final 'e' of the parent hydride name is kept
before 'di', 'tri' (an 'e' is lost only before a vowel) and the locants are lowest for the cationic centres, ahead of
the prefixes (P-31.1.4.2.4): 1,4-dioxane-1,4-diium, tetramethyldiazene-1,2-diium."""

import re

from rdkit import Chem

from ._common import (
    CITE_SKELETAL_LAMBDA,
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
        and len(Chem.GetMolFrags(mol)) == 1
        and all(a.GetFormalCharge() > 0 and not a.GetIsotope() for a in centres)
        and not any(
            a.GetNumRadicalElectrons() and not (a.GetAtomicNum() == 7 and _ylium_centre(a) and a.GetFormalCharge() == 1)
            for a in mol.GetAtoms()
        )
        and (
            (sum(a.GetFormalCharge() for a in centres) >= 2 and (len(centres) > 1 or centres[0].GetAtomicNum() == 6))
            or (len(centres) == 1 and _single_ring_heteroatom_cation(centres[0]))
        )
    )


def _ylium_centre(atom) -> bool:
    """A ring carbon or a ring nitrogen that has lost a hydride ion (carbenium, nitrenium): an 'ylium' centre."""
    return atom.GetAtomicNum() == 6 or (
        atom.GetAtomicNum() == 7 and atom.IsInRing() and atom.GetDegree() == 2 and not atom.GetTotalNumHs() and not atom.GetIsAromatic()
    )


def has_ring_nitrenium_shape(mol) -> bool:
    """A cation with a ring nitrogen that lost a hydride ion (RDKit reports its missing valences as radical electrons)."""
    return any(
        a.GetNumRadicalElectrons() and a.GetAtomicNum() == 7 and _ylium_centre(a) and a.GetFormalCharge() == 1
        for a in mol.GetAtoms()
    ) and has_polycation_shape(mol)


def _single_ring_heteroatom_cation(atom) -> bool:
    if not atom.IsInRing() or atom.GetAtomicNum() not in _RING_CENTRE_ELEMENTS:
        return False
    if atom.GetAtomicNum() == 7:
        return atom.GetDegree() == 4 and all(n.IsInRing() for n in atom.GetNeighbors()) and atom.GetIsAromatic() is False
    return _hydron_added(atom)


def _hydron_added(atom) -> bool:
    """A ring chalcogen or pnictogen that gained a hydron: three sigma bonds for sulfur, selenium and tellurium, or five
    when the parent hydride is the lambda4 one (P-73.8.2)."""
    sigma = atom.GetDegree() + atom.GetTotalNumHs()
    z = atom.GetAtomicNum()
    return sigma == _CATION_VALENCE[z] or (z in (16, 34, 52) and sigma == _CATION_VALENCE[z] + 2)


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
    if all(a.IsInRing() for a in centres) and any(_ylium_centre(a) for a in centres) and any(
        not _ylium_centre(a) for a in centres
    ):
        return _name_ring_ium_ylium(mol, centres)
    if all(a.IsInRing() for a in centres) and not any(_ylium_centre(a) for a in centres):
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


def _name_ring_ium_ylium(mol, centres):
    """'ium' centres on ring heteroatoms and 'ylium' centres on ring carbons of one parent hydride (P-73.5.2): the
    suffixes follow the name in that order, the lowest locants go to all the cationic centres whatever their type and
    then to the 'ylium' centres."""
    from ._diester_ring_diyl import _system_of, evaluate_skeleton

    ium = [a.GetIdx() for a in centres if not _ylium_centre(a)]
    ylium = [a.GetIdx() for a in centres if _ylium_centre(a)]
    if any(a.GetFormalCharge() != 1 or a.GetAtomicNum() not in _RING_CENTRE_ELEMENTS | {6} for a in centres):
        raise UnsupportedStructure("this combination of ring centres is not supported yet")
    stage = Chem.RWMol(mol)
    Chem.Kekulize(stage, clearAromaticFlags=True)
    for index in ylium:
        atom = stage.GetAtomWithIdx(index)
        atom.SetFormalCharge(0)
        atom.SetNoImplicit(True)
        atom.SetNumExplicitHs(atom.GetTotalNumHs() + 1)
    stage = stage.GetMol()
    Chem.SanitizeMol(stage)
    editable = Chem.RWMol(stage)
    for index in ium:
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
    rings, atoms = _system_of(base, ium[0])
    if any(i not in atoms for i in ium + ylium):
        raise UnsupportedStructure("cationic centres in different ring systems are named multiplicatively")
    key_centers = [(i, "ium") for i in ium] + [(i, "uide") for i in ylium]
    found = evaluate_skeleton(base, adjacency(base), "ring", rings, atoms, ium, set(), "ium", key_centers=key_centers)
    if found is None:
        raise UnsupportedStructure("this cationic ring system has no supported name yet")
    locants = ",".join(str(x) for x in sorted(found[2][i] for i in ylium))
    word = "ylium" if len(ylium) == 1 else f"{'bis' if len(ylium) == 2 else 'tris'}(ylium)"
    return f"{found[1]}-{locants}-{word}"


def _name_ring_polycation(mol, centres):
    from ._diester_ring_diyl import _system_of, evaluate_skeleton

    if any(a.GetIsotope() or a.GetAtomicNum() not in _RING_CENTRE_ELEMENTS for a in centres):
        raise UnsupportedStructure("this ring heteroatom is not supported as a cationic centre yet")
    if any(a.GetAtomicNum() != 7 and not _hydron_added(a) for a in centres):
        raise UnsupportedStructure("a ring centre that is not a hydron-added heteroatom is a 'ylium' centre")
    indices = [a.GetIdx() for a in centres]
    rings, atoms = _system_of(mol, indices[0])
    if any(i not in atoms for i in indices):
        raise _DifferentSystems("cationic centres in different ring systems are named multiplicatively")
    ring_sigma = {i: sum(1 for n in mol.GetAtomWithIdx(i).GetNeighbors() if n.GetIdx() in atoms) for i in indices}
    lambda_centres = [i for i in indices if ring_sigma[i] > _STANDARD_VALENCE[mol.GetAtomWithIdx(i).GetAtomicNum()]]
    if lambda_centres and (len(lambda_centres) != len(indices) or len(indices) != 1):
        raise UnsupportedStructure("lambda ylium centres beside other cationic centres are not supported yet")
    editable = Chem.RWMol(mol)
    for index in indices:
        if index in lambda_centres and mol.GetAtomWithIdx(index).GetAtomicNum() == 7:
            continue
        atom = editable.GetAtomWithIdx(index)
        hydrogens = atom.GetTotalNumHs()
        atom.SetFormalCharge(0)
        atom.SetNoImplicit(True)
        atom.SetNumExplicitHs(max(hydrogens - 1, 0))
        if mol.GetAtomWithIdx(index).GetDegree() + hydrogens == _CATION_VALENCE[mol.GetAtomWithIdx(index).GetAtomicNum()]:
            atom.SetBoolProp("_ring_cation_centre", True)
    base = editable.GetMol()
    base.UpdatePropertyCache(strict=False)
    Chem.FastFindRings(base)
    token = CITE_SKELETAL_LAMBDA.set(True)
    try:
        found = evaluate_skeleton(
            base, adjacency(base), "ring", rings, atoms, indices, set(), "ylium" if lambda_centres else "ium"
        )
    finally:
        CITE_SKELETAL_LAMBDA.reset(token)
    if found is None:
        raise UnsupportedStructure("this cationic ring system has no supported name yet")
    if lambda_centres and "λ" not in found[1]:
        from ._anion_center import _insert_lambda

        centre = lambda_centres[0]
        bonding = ring_sigma[centre] + base.GetAtomWithIdx(centre).GetTotalNumHs() + 1
        return _insert_lambda(found[1], {found[2][centre]: bonding})
    return found[1]


_RING_CENTRE_ELEMENTS = {7, 8, 15, 16, 33, 34, 52}
_STANDARD_VALENCE = {7: 3, 8: 2, 15: 3, 16: 2, 33: 3, 34: 2, 52: 2}
_CATION_VALENCE = {7: 4, 8: 3, 15: 4, 16: 3, 33: 4, 34: 3, 52: 3}  # sigma bonds of the hydron-added atom (N: with pi)


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
