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
    pair = (first, partners[0])
    if any(h not in pair for _, _, h in found) or len({g.spec for _, g, _ in found}) != 1:
        raise UnsupportedStructure("this hydrazine acid is not supported yet")
    centers = {c for c, _, _ in found}
    owned = {a for _, g, _ in found for a in g.owned}
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
        options.append((suffix_locants, prefix_locants, format_substituent_prefixes(group_substituents(entries))))
    suffix_locants, _, prefix = min(options)
    spec = found[0][1].spec
    suffix = carbo_suffix(spec, len(found))
    bare = len(found) == 1 and not prefix
    locant_text = "" if bare else "-" + ",".join(map(str, suffix_locants)) + "-"
    return f"{prefix}hydrazine{locant_text}{suffix}"


def name_hetero_parent_acid(mol):
    from .core import smiles_to_iupac

    if len(Chem.GetMolFrags(mol)) != 1:
        raise UnsupportedStructure("several fragments")
    found = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6 or atom.IsInRing():
            continue
        group = acid_group_at(mol, atom.GetIdx(), hetero_attach=True)
        if group is None or group.spec.anion:
            continue
        attach = [n for n in atom.GetNeighbors() if n.GetIdx() not in group.owned and n.GetAtomicNum() != 1]
        if len(attach) == 1 and attach[0].GetAtomicNum() in _PARENT_ELEMENTS and not attach[0].IsInRing():
            found.append((atom.GetIdx(), group, attach[0].GetIdx()))
    if found and all(mol.GetAtomWithIdx(h).GetAtomicNum() == 7 for _, _, h in found):
        hydrazine = _hydrazine_acid(mol, found)
        if hydrazine is not None:
            return hydrazine
    if not found or len({g.spec for _, g, _ in found}) != 1:
        raise UnsupportedStructure("acid groups of one kind on a heteroatom parent are required")
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
