"""Amides, acid halides and aci-nitro compounds of mononuclear inorganic acids (P-41, P-61.3.2.2, P-61.5.3, P-67.1.2.5,
P-67.1.2.6): a nitroso or nitro group on an amine nitrogen makes an amide of nitrous or nitric acid ('methylnitrous
amide', 'methyl(nitro)nitramide'), a halogen on P, As or Sb an acid halide ('methylphosphinous chloride') and the group
=N(O)OH or >N(O)OH an azinic acid ('ethylideneazinic acid'). Carbon groups are cited as prefixes without locants."""

from ._common import HALOGEN_PREFIXES, UnsupportedStructure, adjacency, group_substituents, halogen_substituents
from ._phosphonic_acid import CENTER_STEMS
from ._substituents import format_mononuclear_prefixes, format_substituent_prefixes, name_branch

_HALIDE = {9: "fluoride", 17: "chloride", 35: "bromide", 53: "iodide"}
_MULTIPLIER = {1: "", 2: "di"}


def _bond(mol, a, b):
    return mol.GetBondBetweenAtoms(a, b).GetBondTypeAsDouble()


def _is_nitroso(mol, atom, parent):
    others = [n for n in atom.GetNeighbors() if n.GetIdx() != parent]
    return (
        atom.GetAtomicNum() == 7
        and not atom.GetFormalCharge()
        and len(others) == 1
        and others[0].GetAtomicNum() == 8
        and others[0].GetDegree() == 1
        and _bond(mol, atom.GetIdx(), others[0].GetIdx()) == 2.0
    )


def _is_nitro(mol, atom, parent):
    others = [n for n in atom.GetNeighbors() if n.GetIdx() != parent]
    return (
        atom.GetAtomicNum() == 7
        and atom.GetFormalCharge() == 1
        and len(others) == 2
        and all(o.GetAtomicNum() == 8 and o.GetDegree() == 1 for o in others)
        and sorted(o.GetFormalCharge() for o in others) == [-1, 0]
    )


def _prefix_entries(mol, centre, excluded, graph):
    halogens = halogen_substituents(mol)
    aromatic = frozenset(a.GetIdx() for a in mol.GetAtoms() if a.GetIsAromatic())
    return [
        name_branch(graph, n.GetIdx(), centre.GetIdx(), halogens, aromatic, mol=mol, unsaturated=True)
        for n in centre.GetNeighbors()
        if n.GetIdx() not in excluded
    ]


def _terminal_nitrogen_ok(mol, nitrogen, centre):
    """The far nitrogen of a hydrazide carries only hydrogen, carbon groups or one ylidene carbon."""
    others = [n for n in nitrogen.GetNeighbors() if n.GetIdx() != centre.GetIdx()]
    if any(n.GetAtomicNum() != 6 for n in others):
        return False
    return sum(_bond(mol, nitrogen.GetIdx(), n.GetIdx()) for n in others) + nitrogen.GetTotalNumHs() == 2


def _nitrogen_amide(mol):
    for centre in mol.GetAtoms():
        if centre.GetAtomicNum() != 7 or centre.GetFormalCharge() or centre.IsInRing() or centre.GetDegree() < 1:
            continue
        nitro = [n for n in centre.GetNeighbors() if _is_nitro(mol, n, centre.GetIdx())]
        nitroso = [n for n in centre.GetNeighbors() if _is_nitroso(mol, n, centre.GetIdx())]
        if not nitro and not nitroso:
            continue
        taken = {n.GetIdx() for n in (*nitro, *nitroso)}
        others = [n for n in centre.GetNeighbors() if n.GetIdx() not in taken]
        amino = [n for n in others if n.GetAtomicNum() == 7 and not n.GetFormalCharge() and not n.IsInRing()]
        if len(amino) == 1 and len(nitro) + len(nitroso) == 1 and _terminal_nitrogen_ok(mol, amino[0], centre):
            carbons = [n for n in others if n.GetIdx() != amino[0].GetIdx()]
            if all(n.GetAtomicNum() == 6 and _bond(mol, centre.GetIdx(), n.GetIdx()) == 1.0 for n in carbons):
                covered = {centre.GetIdx(), amino[0].GetIdx(), *taken, *(o.GetIdx() for g in (*nitro, *nitroso) for o in g.GetNeighbors())}
                if all(a.GetIdx() in covered or a.GetAtomicNum() == 6 or a.GetAtomicNum() in HALOGEN_PREFIXES for a in mol.GetAtoms()):
                    return centre, nitro, nitroso, carbons, amino
        if any(n.GetAtomicNum() != 6 or _bond(mol, centre.GetIdx(), n.GetIdx()) != 1.0 for n in others):
            continue
        covered = {centre.GetIdx()}
        for group in (*nitro, *nitroso):
            covered |= {group.GetIdx(), *(o.GetIdx() for o in group.GetNeighbors())}
        carbons_ok = all(
            a.GetIdx() in covered or a.GetAtomicNum() == 6 or a.GetAtomicNum() in HALOGEN_PREFIXES for a in mol.GetAtoms()
        )
        if carbons_ok:
            return centre, nitro, nitroso, others, []
    return None


def _name_nitrogen_amide(mol, found):
    centre, nitro, nitroso, others, hydrazide = found
    graph = adjacency(mol)
    if hydrazide:
        parent = "nitric hydrazide" if nitro else "nitrous hydrazide"
        far = hydrazide[0]
        entries = {}
        for letter, atom, skip in (("N", centre, {far.GetIdx(), *(n.GetIdx() for n in (*nitro, *nitroso))}), ("N'", far, {centre.GetIdx()})):
            for n in atom.GetNeighbors():
                if n.GetIdx() in skip:
                    continue
                entries.setdefault(letter, []).append(
                    name_branch(graph, n.GetIdx(), atom.GetIdx(), halogen_substituents(mol), mol=mol, unsaturated=True)
                )
        if not entries:
            return parent
        return format_substituent_prefixes(group_substituents(entries)) + parent
    if nitro:
        parent, spare_nitro, spare_nitroso = "nitramide", nitro[1:], nitroso
    else:
        parent, spare_nitro, spare_nitroso = "nitrous amide", [], nitroso[1:]
    excluded = {n.GetIdx() for n in (*nitro, *nitroso)}
    entries = _prefix_entries(mol, centre, excluded, graph)
    entries += [("nitro", False)] * len(spare_nitro) + [("nitroso", False)] * len(spare_nitroso)
    prefix = format_mononuclear_prefixes(entries) if entries else ""
    return prefix + parent


def _halide_centre(mol):
    centres = [a for a in mol.GetAtoms() if a.GetAtomicNum() in CENTER_STEMS]
    if len(centres) != 1:
        return None
    centre = centres[0]
    if centre.GetFormalCharge() or centre.IsInRing() or centre.GetIsotope():
        return None
    halides = [n for n in centre.GetNeighbors() if n.GetAtomicNum() in HALOGEN_PREFIXES and n.GetDegree() == 1]
    carbons = [n for n in centre.GetNeighbors() if n.GetAtomicNum() == 6]
    if (
        not halides
        or not carbons
        or len(halides) + len(carbons) != centre.GetDegree()
        or len({h.GetAtomicNum() for h in halides}) != 1
        or any(_bond(mol, centre.GetIdx(), n.GetIdx()) != 1.0 for n in centre.GetNeighbors())
    ):
        return None
    hydrogens = centre.GetTotalNumHs()
    if len(halides) + len(carbons) + hydrogens != 3 or len(halides) > 2:
        return None
    covered = {centre.GetIdx(), *(h.GetIdx() for h in halides)}
    if not all(a.GetIdx() in covered or a.GetAtomicNum() == 6 or a.GetAtomicNum() in HALOGEN_PREFIXES for a in mol.GetAtoms()):
        return None
    return centre, halides, carbons


def _name_halide(mol, found):
    centre, halides, carbons = found
    graph = adjacency(mol)
    excluded = {h.GetIdx() for h in halides}
    entries = _prefix_entries(mol, centre, excluded, graph)
    prefix = format_mononuclear_prefixes(entries) if entries else ""
    ending = "inous" if len(halides) == 1 else "onous"
    halide = _MULTIPLIER[len(halides)] + _HALIDE[halides[0].GetAtomicNum()]
    return f"{prefix}{CENTER_STEMS[centre.GetAtomicNum()]}{ending} {halide}"


def _azinic_centre(mol):
    for centre in mol.GetAtoms():
        if centre.GetAtomicNum() != 7 or centre.GetFormalCharge() != 1 or centre.IsInRing() or centre.GetDegree() not in (3, 4):
            continue
        oxygens = [n for n in centre.GetNeighbors() if n.GetAtomicNum() == 8 and n.GetDegree() == 1]
        oxide = [o for o in oxygens if o.GetFormalCharge() == -1]
        hydroxyl = [o for o in oxygens if not o.GetFormalCharge() and o.GetTotalNumHs() == 1]
        if len(oxide) != 1 or len(hydroxyl) != 1:
            continue
        rest = [n for n in centre.GetNeighbors() if n.GetAtomicNum() != 8 or n.GetDegree() != 1]
        orders = sorted(_bond(mol, centre.GetIdx(), n.GetIdx()) for n in rest)
        if any(n.GetAtomicNum() != 6 for n in rest) or orders not in ([2.0], [1.0, 1.0]):
            continue
        covered = {centre.GetIdx(), oxide[0].GetIdx(), hydroxyl[0].GetIdx()}
        if all(a.GetIdx() in covered or a.GetAtomicNum() == 6 or a.GetAtomicNum() in HALOGEN_PREFIXES for a in mol.GetAtoms()):
            return centre, {oxide[0].GetIdx(), hydroxyl[0].GetIdx()}
    return None


def _name_azinic(mol, found):
    centre, excluded = found
    graph = adjacency(mol)
    entries = _prefix_entries(mol, centre, excluded, graph)
    prefix = format_mononuclear_prefixes(entries) if len(entries) > 1 else (entries[0][0] if entries else "")
    if len(entries) == 1 and entries[0][1]:
        prefix = f"({prefix})"
    return f"{prefix}azinic acid"


def has_inorganic_acid_derivative_shape(mol) -> bool:
    return (
        _nitrogen_amide(mol) is not None or _halide_centre(mol) is not None or _azinic_centre(mol) is not None
    )


def name_inorganic_acid_derivative(mol) -> str:
    from rdkit import Chem

    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")
    if any(a.GetIsotope() or a.GetNumRadicalElectrons() for a in mol.GetAtoms()):
        raise UnsupportedStructure("isotopes and radicals are not supported on an acid derivative of this kind")
    found = _nitrogen_amide(mol)
    if found is not None:
        return _name_nitrogen_amide(mol, found)
    found = _halide_centre(mol)
    if found is not None:
        return _name_halide(mol, found)
    found = _azinic_centre(mol)
    if found is not None:
        return _name_azinic(mol, found)
    raise UnsupportedStructure("no mononuclear inorganic acid derivative")
