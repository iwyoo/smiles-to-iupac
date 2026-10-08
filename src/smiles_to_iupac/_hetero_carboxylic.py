"""Acid groups on an acyclic heteroatom parent hydride (P-65.1.2.2.2): 'carboxylic
acid' and its functional-replacement analogues are cited as the suffix of a
heteroacyclic parent such as silane, phosphane, hydrazine or disiloxane
('hydrazinecarboxylic acid', 'disiloxanecarboxylic acid'). The parent hydride
must have a single, unambiguous attachment position, so no locant is cited."""

import re

from rdkit import Chem

from ._acid_groups import acid_group_at
from ._acid_lexicon import carbo_suffix
from ._common import UnsupportedStructure, adjacency, halogen_substituents

_PARENT_ELEMENTS = {5, 7, 13, 14, 15, 32, 33, 50, 51, 82, 83}
_ACID_CENTERS = {6, 16, 34, 52}


def _hydrazine_acid(mol, found):
    """'2-methylhydrazine-1-carboxylic acid', 'hydrazine-1,2-dicarboxylic acid': the two nitrogens are numbered so the
    suffix gets the lowest locants, then the prefixes; None when the acid groups are not on an N-N parent."""
    from ._common import group_substituents
    from ._substituents import format_substituent_prefixes, name_branch

    first = found[0][2]
    partners = [
        n.GetIdx()
        for n in mol.GetAtomWithIdx(first).GetNeighbors()
        if n.GetAtomicNum() == 7 and not n.IsInRing()
    ]
    if len(partners) != 1:
        return None
    if mol.GetBondBetweenAtoms(first, partners[0]).GetBondTypeAsDouble() != 1.0:
        raise UnsupportedStructure("a diazene is not a hydrazine parent")
    pair = (first, partners[0])
    if any(h not in pair for _, _, h in found) or len({g.spec for _, g, _ in found}) != 1:
        raise UnsupportedStructure("this hydrazine acid is not supported yet")
    centers = {c for c, _, _ in found}
    owned = {a for _, g, _ in found for a in g.owned}
    cations = [i for i in pair if mol.GetAtomWithIdx(i).GetFormalCharge() == 1 and mol.GetAtomWithIdx(i).GetTotalDegree() == 4]
    anions = {a for a in owned if mol.GetAtomWithIdx(a).GetFormalCharge()}
    if {a.GetIdx() for a in mol.GetAtoms() if a.GetFormalCharge()} - anions - set(cations) or (
        cations and (len(cations) > 1 or not found[0][1].spec.anion)
    ):
        raise UnsupportedStructure("this charged hydrazine acid is not supported yet")
    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    options = []
    for numbering in (pair, pair[::-1]):
        entries = {}
        for position, nitrogen in enumerate(numbering, start=1):
            for n in graph[nitrogen]:
                if n in pair or n in centers or n in owned:
                    continue
                entries.setdefault(position, []).append(name_branch(graph, n, nitrogen, halogens, mol=mol))
        suffix_locants = sorted(numbering.index(h) + 1 for _, _, h in found)
        prefix_locants = sorted(p for p, names in entries.items() for _ in names)
        ium_locants = [numbering.index(c) + 1 for c in cations]
        options.append((ium_locants, suffix_locants, prefix_locants, format_substituent_prefixes(group_substituents(entries))))
    ium_locants, suffix_locants, _, prefix = min(options)
    spec = found[0][1].spec
    suffix = carbo_suffix(spec, len(found))
    suffix_text = "-" + ",".join(map(str, suffix_locants)) + "-"
    if ium_locants:
        return f"{prefix}hydrazin-{ium_locants[0]}-ium{suffix_text}{suffix}"
    bare = len(found) == 1 and not prefix
    return f"{prefix}hydrazine{'' if bare else suffix_text}{suffix}"


def hydrazine_acyl_prefix(mol, graph, nitrogen, acyl_atom, halogens, ending):
    """'hydrazinesulfonyl', '2-methylhydrazine-1-sulfonyl' (P-66.3.2.1): the acylated nitrogen of a hydrazine takes
    locant 1 as the free valence; None unless it sits in an acyclic N-N pair."""
    from ._common import group_substituents
    from ._substituents import format_substituent_prefixes, name_branch

    first = mol.GetAtomWithIdx(nitrogen)
    partners = [n for n in graph[nitrogen] if n != acyl_atom and mol.GetAtomWithIdx(n).GetAtomicNum() == 7]
    if first.IsInRing() or first.GetFormalCharge() or len(partners) != 1:
        return None
    partner = mol.GetAtomWithIdx(partners[0])
    if partner.IsInRing() or partner.GetFormalCharge() or mol.GetBondBetweenAtoms(nitrogen, partners[0]).GetBondTypeAsDouble() != 1.0:
        return None
    entries = {}
    for position, atom in enumerate((nitrogen, partners[0]), start=1):
        for n in graph[atom]:
            if n not in (acyl_atom, nitrogen, partners[0]):
                entries.setdefault(position, []).append(name_branch(graph, n, atom, halogens, mol=mol))
    prefix = format_substituent_prefixes(group_substituents(entries))
    return f"{prefix}hydrazine{'-1-' if prefix else ''}{ending}"


def name_hetero_parent_acid(mol):
    from .core import smiles_to_iupac

    if len(Chem.GetMolFrags(mol)) != 1:
        raise UnsupportedStructure("several fragments")
    found = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in _ACID_CENTERS or atom.IsInRing():
            continue
        group = acid_group_at(mol, atom.GetIdx(), hetero_attach=True)
        if group is None:
            continue
        attach = [n for n in atom.GetNeighbors() if n.GetIdx() not in group.owned and n.GetAtomicNum() != 1]
        if (
            len(attach) == 1
            and attach[0].GetAtomicNum() in _PARENT_ELEMENTS
            and not attach[0].IsInRing()
            and not _is_acylated(mol, attach[0], atom.GetIdx())
        ):
            found.append((atom.GetIdx(), group, attach[0].GetIdx()))
    if found and all(mol.GetAtomWithIdx(h).GetAtomicNum() == 7 for _, _, h in found):
        hydrazine = _hydrazine_acid(mol, found)
        if hydrazine is not None:
            return hydrazine
    if not found or len({g.spec for _, g, _ in found}) != 1:
        raise UnsupportedStructure("acid groups of one kind on a heteroatom parent are required")
    if found[0][1].spec.anion or any(mol.GetAtomWithIdx(c).GetAtomicNum() != 6 for c, _, _ in found):
        raise UnsupportedStructure("only a hydrazine parent takes anionic or chalcogen acid groups here")
    if any(a.GetFormalCharge() for a in mol.GetAtoms()):
        raise UnsupportedStructure("a charged atom on a heteroatom parent is not supported here")
    spec = found[0][1].spec
    hosts = {h for _, _, h in found}
    editable = Chem.RWMol(mol)
    for _, _, host in found:
        host_atom = editable.GetAtomWithIdx(host)
        host_atom.SetNumExplicitHs(host_atom.GetTotalNumHs() + 1)
        host_atom.SetNoImplicit(True)
    removed = {a for c, g, _ in found for a in (c, *g.owned)}
    for idx in sorted(removed, reverse=True):
        editable.RemoveAtom(idx)
    hydride = editable.GetMol()
    Chem.SanitizeMol(hydride)
    parent = smiles_to_iupac(Chem.MolToSmiles(hydride))
    if len(hosts) == 1 and not (any(ch.isdigit() for ch in parent) or " " in parent):
        host = next(iter(hosts))
        new_host = host - sum(1 for r in removed if r < host)
        if _single_site_kind(hydride, new_host) or _mononuclear(hydride):
            return parent + carbo_suffix(spec, len(found))
    elif len(hosts) == 1 and _mononuclear(hydride):
        raise UnsupportedStructure("the heteroatom parent has several possible attachment positions")
    return _located_parent_acid(mol, found, spec)


def _is_acylated(mol, host, acid_carbon):
    """A host carrying a carbon acyl group (an amide, R-CO-N) makes the group a carbamic or carbonic acid derivative
    (P-65.2), not a carboxylic acid on a hetero parent hydride."""
    return any(
        n.GetIdx() != acid_carbon
        and n.GetAtomicNum() == 6
        and any(b.GetBondTypeAsDouble() == 2.0 and b.GetOtherAtom(n).GetAtomicNum() != 6 for b in n.GetBonds())
        and any(m.GetAtomicNum() == 6 for m in n.GetNeighbors())
        for n in host.GetNeighbors()
    )


def _mononuclear(hydride):
    return sum(1 for a in hydride.GetAtoms() if a.GetAtomicNum() not in (1, 6)) == 1


def _single_site_kind(hydride, host):
    """Whether every hydrogen-bearing non-carbon atom of the hydride is equivalent to `host` (so no locant is needed)."""
    ranks = list(Chem.CanonicalRankAtoms(hydride, breakTies=False))
    return all(
        ranks[a.GetIdx()] == ranks[host]
        for a in hydride.GetAtoms()
        if a.GetAtomicNum() not in (1, 6) and a.GetTotalNumHs()
    )


_IODO = re.compile(r"^(?P<locants>\d+(?:,\d+)*)-(?:di|tri|tetra)?iodo(?P<parent>[a-z]+)$")


def _located_parent_acid(mol, found, spec):
    """'disilane-1,2-dicarboxylic acid', 'disiloxane-1,3-dicarboxylic acid': the acid groups are replaced by iodine,
    which cites after every other prefix and so numbers like the suffix, and the locants are read from that name;
    parents that carry other substituents are not named."""
    from .core import smiles_to_iupac

    editable = Chem.RWMol(mol)
    for _, _, host in found:
        editable.AddBond(host, editable.AddAtom(Chem.Atom(53)), Chem.BondType.SINGLE)
    removed = {a for c, g, _ in found for a in (c, *g.owned)}
    for idx in sorted(removed, reverse=True):
        editable.RemoveAtom(idx)
    surrogate = editable.GetMol()
    Chem.SanitizeMol(surrogate)
    match = _IODO.match(smiles_to_iupac(Chem.MolToSmiles(surrogate)))
    if match is None:
        raise UnsupportedStructure("the heteroatom parent carries other substituents or has an unusual name")
    return f"{match.group('parent')}-{match.group('locants')}-{carbo_suffix(spec, len(found))}"
