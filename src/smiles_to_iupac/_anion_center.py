"""Anionic centers on a mononuclear parent hydride (P-72.2.2.1, P-72.3, P-72.8):
'ide' when the center has lost hydrons (bonding number below the standard one)
and 'uide' when it has gained hydrides (bonding number above it), with the
lambda convention for a non-standard neutral parent. Substituents are cited as
prefixes of the parent hydride (P-16.5.1.3).
"""

import re

from rdkit import Chem

from ._common import UnsupportedStructure, adjacency, halogen_substituents
from ._substituents import format_mononuclear_prefixes, name_branch

_PARENTS = {
    5: ("boran", 3),
    7: ("azan", 3),
    8: ("oxidan", 2),
    9: ("fluoran", 1),
    13: ("aluman", 3),
    14: ("silan", 4),
    15: ("phosphan", 3),
    16: ("sulfan", 2),
    17: ("chloran", 1),
    31: ("gallan", 3),
    32: ("german", 4),
    33: ("arsan", 3),
    34: ("selan", 2),
    35: ("broman", 1),
    49: ("indigan", 3),
    50: ("stannan", 4),
    51: ("stiban", 3),
    52: ("tellan", 2),
    53: ("iodan", 1),
    81: ("thallan", 3),
    82: ("plumban", 4),
    83: ("bismuthan", 3),
}
_MULTIPLIED = {1: "", 2: "di", 3: "tri"}


def has_center_anion_shape(mol):
    return any(a.GetFormalCharge() < 0 and a.GetAtomicNum() in _PARENTS for a in mol.GetAtoms())


_SENIORITY = [7, 15, 33, 51, 83, 14, 32, 50, 82, 5, 13, 31, 49, 81, 8, 16, 34, 52, 9, 17, 35, 53]


def name_center_anion(mol):
    from ._anion import marked_neutral

    if len(Chem.GetMolFrags(mol)) != 1:
        raise UnsupportedStructure("a multi-fragment anionic structure is not supported here")
    centers = [a for a in mol.GetAtoms() if a.GetFormalCharge() < 0 and not _is_nitro_oxygen(a)]
    cations = [a for a in mol.GetAtoms() if a.GetFormalCharge() > 0 and not _is_nitro_nitrogen(a)]
    if cations:
        if len(cations) == 1 and cations[0].IsInRing() and cations[0].GetFormalCharge() == 1 and all(a.IsInRing() for a in centers):
            return _name_ring_zwitterion(mol, centers, cations[0])
        raise UnsupportedStructure("cationic centers beside an anionic center are not supported here")
    if any(a.IsInRing() for a in centers):
        if len(centers) > 1 and not all(a.IsInRing() for a in centers):
            raise UnsupportedStructure("ring and chain anionic centers together are not supported yet")
        return _name_skeletal_anion(mol, centers)
    if any(a.GetAtomicNum() not in _PARENTS or a.GetNumRadicalElectrons() for a in centers):
        raise UnsupportedStructure("only anionic mononuclear centers are supported here")
    if len(centers) == 1:
        special = _single_center_special(mol, centers[0])
        if special is not None:
            return special
    neutral = marked_neutral(mol)
    from ._anion_chain import is_chain_anion, name_chain_anion

    if is_chain_anion(neutral):
        return name_chain_anion(neutral)
    candidates = []
    for atom in neutral.GetAtoms():
        if not atom.HasProp("_anion_word"):
            continue
        charge, word = int(atom.GetProp("_anion_charge")), atom.GetProp("_anion_word")
        name = _mononuclear_name(neutral, atom.GetIdx())
        candidates.append(((-charge, word != "uide", _SENIORITY.index(atom.GetAtomicNum()), name), name))
    return min(candidates)[1]


def _mononuclear_name(mol, root):
    atom = mol.GetAtomWithIdx(root)
    stem = _PARENTS[atom.GetAtomicNum()][0]
    word, charge = atom.GetProp("_anion_word"), int(atom.GetProp("_anion_charge"))
    lam = atom.GetProp("_anion_lambda") if atom.HasProp("_anion_lambda") else None
    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    aromatic_atoms = frozenset(a.GetIdx() for a in mol.GetAtoms() if a.GetIsAromatic())
    entries = [
        name_branch(graph, n, root, halogens, aromatic_atoms, mol=mol, unsaturated=True) for n in graph[root]
    ]
    prefixes = format_mononuclear_prefixes(entries) if entries else ""
    parent = f"{prefixes}{'-' if prefixes and lam else ''}{f'λ{lam}-' if lam else ''}{stem}"
    return parent + ("e" if charge > 1 else "") + _MULTIPLIED[charge] + word


def _single_center_special(mol, center):
    z, charge = center.GetAtomicNum(), -center.GetFormalCharge()
    heavy = list(center.GetNeighbors())
    if not heavy and z == 8 and charge == 1:
        return "hydroxide"
    if z == 8 and charge == 1 and len(heavy) == 1 and heavy[0].GetAtomicNum() == 8 and heavy[0].GetDegree() == 1:
        return "hydroperoxide"
    if z == 7 and charge == 1 and len(heavy) == 1 and heavy[0].GetAtomicNum() in (15, 33, 51):
        host = heavy[0]
        if mol.GetBondBetweenAtoms(center.GetIdx(), host.GetIdx()).GetBondTypeAsDouble() == 2.0:
            return _imide_of_hydride(mol, center, host)
    if z == 8 and charge == 1 and len(heavy) == 1 and heavy[0].GetAtomicNum() in _HYDROXY_HYDRIDE_HOSTS:
        named = _hydride_olate(mol, center, heavy[0])
        if named is not None:
            return named
    if z == 8 and charge == 1 and len(heavy) == 1 and _is_plain_amino_nitrogen(heavy[0], center):
        graph = adjacency(mol)
        halogens = halogen_substituents(mol)
        aromatic_atoms = frozenset(a.GetIdx() for a in mol.GetAtoms() if a.GetIsAromatic())
        nitrogen = heavy[0]
        subs = [
            name_branch(graph, n.GetIdx(), nitrogen.GetIdx(), halogens, aromatic_atoms, mol=mol, unsaturated=True)
            for n in nitrogen.GetNeighbors()
            if n.GetIdx() != center.GetIdx()
        ]
        return (format_mononuclear_prefixes(subs) if subs else "") + "aminoxide"
    return None


_HYDROXY_HYDRIDE_HOSTS = (13, 14, 31, 32, 49, 50, 81, 82)


def _hydride_olate(mol, center, host):
    """'dimethylalumanolate', 'dimethylthallanolate' (P-68.1.4.1, P-68.1.5.1): the anion of a hydroxy group on a hydride
    of Group 13 or 14 carrying only organyl groups or halogens; the hydroxy parent has the suffix 'ol' (P-72.2.2.2.2)."""
    if host.IsInRing() or host.GetFormalCharge() or host.GetIsotope() or host.GetDegree() > _PARENTS[host.GetAtomicNum()][1]:
        return None
    if any(a.GetAtomicNum() not in (1, 6, 9, 17, 35, 53, host.GetAtomicNum(), 8) or a.GetIsAromatic() for a in mol.GetAtoms()):
        return None
    if sum(a.GetAtomicNum() == 8 for a in mol.GetAtoms()) != 1:
        return None
    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    aromatic_atoms = frozenset()
    entries = [
        name_branch(graph, n, host.GetIdx(), halogens, aromatic_atoms, mol=mol, unsaturated=True)
        for n in graph[host.GetIdx()]
        if n != center.GetIdx()
    ]
    prefixes = format_mononuclear_prefixes(entries) if entries else ""
    return prefixes + _PARENTS[host.GetAtomicNum()][0] + "olate"


def _is_nitro_oxygen(atom):
    return any(n.GetFormalCharge() > 0 for n in atom.GetNeighbors())


def _is_nitro_nitrogen(atom):
    return atom.GetAtomicNum() == 7 and any(n.GetFormalCharge() < 0 for n in atom.GetNeighbors())


def _is_plain_amino_nitrogen(nitrogen, oxide):
    return (
        nitrogen.GetAtomicNum() == 7
        and not nitrogen.GetFormalCharge()
        and not nitrogen.IsInRing()
        and all(
            n.GetIdx() == oxide.GetIdx() or (n.GetAtomicNum() == 6 and not _is_acyl_carbon(n))
            for n in nitrogen.GetNeighbors()
        )
    )


def _is_acyl_carbon(carbon):
    return any(
        b.GetBondTypeAsDouble() >= 2.0 and b.GetOtherAtom(carbon).GetAtomicNum() in (7, 8, 16)
        for b in carbon.GetBonds()
    )


def is_center_atom(atom):
    return atom.GetAtomicNum() in _PARENTS and not atom.GetNumRadicalElectrons()


def center_kind(center, preferred):
    return _center_kind(center, _standard(center), preferred)


def _bond_count(atom):
    mol = Chem.Mol(atom.GetOwningMol())
    try:
        Chem.Kekulize(mol, clearAromaticFlags=True)
    except Exception:
        pass
    center = mol.GetAtomWithIdx(atom.GetIdx())
    return int(sum(b.GetBondTypeAsDouble() for b in center.GetBonds()) + center.GetTotalNumHs())


def _center_kind(center, standard, preferred):
    """('ide' | 'uide', lambda number or None) of an anionic center (P-72.2.2.1, P-72.3, P-72.8)."""
    charge = -center.GetFormalCharge()
    bonds = _bond_count(center)
    n_ide, n_uide = bonds + charge, bonds - charge
    if n_ide == standard:
        return "ide", None
    if n_uide == standard and n_ide > standard and any(b.GetBondTypeAsDouble() > 1.0 for b in center.GetBonds()):
        # P-72.2.1: a centre with a multiple bond (oxo(phenyl)-λ4-sulfanide) is the hydron-loss anion of a λ-hydride
        return "ide", n_ide
    if n_uide == standard:
        return "uide", None
    if n_uide > standard and n_ide > standard:
        if preferred == "uide":
            return "uide", n_uide
        return "ide", n_ide
    if n_ide > standard:
        return "ide", n_ide
    if n_uide > standard:
        return "uide", n_uide
    raise UnsupportedStructure("this bonding number is not an ide or uide center")


def _name_skeletal_anion(mol, centers):
    from ._diester_ring_diyl import _system_of, evaluate_skeleton

    if any(a.GetAtomicNum() not in _PARENTS and a.GetAtomicNum() != 6 for a in centers):
        raise UnsupportedStructure("this skeletal anionic atom is not supported yet")
    kinds = {a.GetIdx(): _center_kind(a, _standard(a), preferred="ide") for a in centers}
    words = {kind for kind, _ in kinds.values()}
    lams = {i: lam for i, (_, lam) in kinds.items()}
    mixed = len(words) == 2
    editable = Chem.RWMol(mol)
    attach, uide_extra = [], []
    for center in centers:
        atom = editable.GetAtomWithIdx(center.GetIdx())
        weight = -center.GetFormalCharge()
        if mixed and kinds[center.GetIdx()][0] == "uide":
            uide_extra += [center.GetIdx()] * weight
        else:
            attach += [center.GetIdx()] * weight
        hydrogens = atom.GetTotalNumHs()
        atom.SetFormalCharge(0)
        atom.SetNoImplicit(True)
        atom.SetNumExplicitHs(hydrogens + (center.GetAtomicNum() == 7 and center.IsInRing() and center.GetDegree() == 2 and weight))
    neutral = editable.GetMol()
    neutral.UpdatePropertyCache(strict=False)
    graph = adjacency(neutral)
    rings, atoms = _system_of(neutral, (attach or uide_extra)[0])
    if any(a not in atoms for a in attach + uide_extra):
        raise UnsupportedStructure("anionic centers in different ring systems are not supported yet")
    word = "ide" if mixed or "ide" in words else "uide"
    key_centers = [(a.GetIdx(), kinds[a.GetIdx()][0]) for a in centers] if mixed else ()
    found = evaluate_skeleton(
        neutral, graph, "ring", rings, atoms, attach or uide_extra, set(), word if attach else "uide",
        key_centers=key_centers,
    )
    if found is None:
        raise UnsupportedStructure("this anionic ring system has no supported name yet")
    name, position_of = found[1], found[2]
    if mixed:
        uide_locants = sorted(position_of[i] for i in uide_extra)
        count = {1: "", 2: "di", 3: "tri"}[len(uide_locants)]
        name = name[: -len("e")] if name.endswith("e") else name
        name = f"{name}-{','.join(str(x) for x in uide_locants)}-{count}uide"
    lam_locants = {position_of[i]: n for i, n in lams.items() if n is not None}
    return _insert_lambda(name, lam_locants) if lam_locants else name


def _name_ring_zwitterion(mol, centers, cation):
    """One ring cation and ring anions of one parent hydride: the suffixes are cumulative, 'ium' before 'ide' (P-74.1.1)."""
    from ._diester_ring_diyl import _system_of, evaluate_skeleton

    if any(a.GetAtomicNum() not in _PARENTS and a.GetAtomicNum() != 6 for a in centers):
        raise UnsupportedStructure("this skeletal anionic atom is not supported yet")
    kinds = [_center_kind(a, _standard(a), preferred="ide") for a in centers]
    if any(word != "ide" or lam for word, lam in kinds):
        raise UnsupportedStructure("only plain ide centers beside a ring cation are supported here")
    editable = Chem.RWMol(mol)
    for center in centers:
        atom = editable.GetAtomWithIdx(center.GetIdx())
        hydrogens = atom.GetTotalNumHs()
        atom.SetFormalCharge(0)
        atom.SetNoImplicit(True)
        atom.SetNumExplicitHs(hydrogens)
    hydron_added = editable.GetAtomWithIdx(cation.GetIdx())
    hydrogens = hydron_added.GetTotalNumHs()
    hydron_added.SetFormalCharge(0)
    hydron_added.SetNoImplicit(True)
    hydron_added.SetNumExplicitHs(max(hydrogens - 1, 0))
    hydron_added.SetBoolProp("_ring_cation_centre", True)
    base = editable.GetMol()
    base.UpdatePropertyCache(strict=False)
    Chem.FastFindRings(base)
    rings, atoms = _system_of(base, cation.GetIdx())
    indices = [a.GetIdx() for a in centers]
    if any(i not in atoms for i in indices):
        raise UnsupportedStructure("the ionic centers lie in different ring systems")
    key_centers = [(cation.GetIdx(), "ium")] + [(i, "ide") for i in indices]
    found = evaluate_skeleton(
        base, adjacency(base), "ring", rings, atoms, [cation.GetIdx()], set(), "ium", key_centers=key_centers
    )
    if found is None:
        raise UnsupportedStructure("this zwitterionic ring system has no supported name yet")
    locants = sorted(found[2][i] for i in indices)
    if not 1 <= len(locants) <= 3:
        raise UnsupportedStructure("this zwitterionic ring system has no supported number of anionic centers")
    count = {1: "", 2: "di", 3: "tri"}[len(locants)]
    return f"{found[1]}-{','.join(str(x) for x in locants)}-{count}ide"


def _standard(atom):
    return 4 if atom.GetAtomicNum() == 6 else _PARENTS[atom.GetAtomicNum()][1]


_RING_HETERO_STEM = (
    r"(?:di|tri|tetra)?(?:thi|selen|tellur|phosphin|phosphol|phosphor|silin|silol|german|arsin|arsol|stibin|stibol|borin|borol|"
    r"alumin|azin|azol|oxin|oxol|thia|phospha|bora|sila|aza|oxa|selena|tellura|arsa|stiba)"
)


def _insert_lambda(name, lam_locants):
    def mark(locants):
        return ",".join(
            f"{loc}\u03bb{lam_locants[int(loc)]}" if loc.isdigit() and int(loc) in lam_locants else loc for loc in locants
        )

    grouped = re.search(rf"(\d+(?:,\d+)*)-(?={_RING_HETERO_STEM})", name)
    if grouped is not None and {int(x) for x in grouped.group(1).split(",")} >= set(lam_locants):
        return name[: grouped.start()] + mark(grouped.group(1).split(",")) + "-" + name[grouped.end():]
    bare = re.search(_RING_HETERO_STEM, name)
    if bare is None:
        raise UnsupportedStructure("the lambda convention is not supported for this ring name yet")
    text = ",".join(f"{loc}\u03bb{lam_locants[loc]}" for loc in sorted(lam_locants))
    return name[: bare.start()] + text + "-" + name[bare.start():]


def center_prefix(mol, root, coming_from):
    """'azanidyl', 'boranuidyl', 'azanediidyl', ... for an anionic mononuclear center cited as a prefix (P-72.6.3)."""
    atom = mol.GetAtomWithIdx(root)
    if atom.IsInRing():
        raise UnsupportedStructure("an anionic ring atom is named with its ring system")
    stem = _PARENTS[atom.GetAtomicNum()][0]
    word, charge = atom.GetProp("_anion_word"), int(atom.GetProp("_anion_charge"))
    lam = atom.GetProp("_anion_lambda") if atom.HasProp("_anion_lambda") else None
    order = mol.GetBondBetweenAtoms(root, coming_from).GetBondTypeAsDouble()
    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    aromatic_atoms = frozenset(a.GetIdx() for a in mol.GetAtoms() if a.GetIsAromatic())
    entries = [
        name_branch(graph, n, root, halogens, aromatic_atoms, mol=mol, unsaturated=True)
        for n in graph[root]
        if n != coming_from
    ]
    prefixes = format_mononuclear_prefixes(entries) if entries else ""
    parent = f"{prefixes}{'-' if prefixes and lam else ''}{f'λ{lam}-' if lam else ''}{stem}"
    body = parent + ("e" if charge > 1 else "") + _MULTIPLIED[charge] + word[:-1]
    return body + {1.0: "yl", 2.0: "ylidene", 3.0: "ylidyne"}[order], bool(entries)


def _imide_of_hydride(mol, center, host):
    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    aromatic_atoms = frozenset(a.GetIdx() for a in mol.GetAtoms() if a.GetIsAromatic())
    entries = [
        name_branch(graph, n.GetIdx(), host.GetIdx(), halogens, aromatic_atoms, mol=mol, unsaturated=True)
        for n in host.GetNeighbors()
        if n.GetIdx() != center.GetIdx()
    ]
    stem, standard = _PARENTS[host.GetAtomicNum()]
    bonds = int(sum(b.GetBondTypeAsDouble() for b in host.GetBonds()) + host.GetTotalNumHs())
    prefixes = format_mononuclear_prefixes(entries) if entries else ""
    lam = f"-λ{bonds}-" if bonds != standard else ""
    return f"{prefixes}{lam}{stem}iminide"
