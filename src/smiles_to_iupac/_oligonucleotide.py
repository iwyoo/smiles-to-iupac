"""Oligonucleotides: nucleosides joined by phosphodiester links (P-106.3.4, P-106.3.5).

Cutting out each link phosphorus leaves nucleosides; every unit but the last becomes its acyl group
('2′-deoxyguanylyl', 'P-thio' before the base stem for a sulfur-bearing link) plus '(3′→5′)', read from the end
that gives the lower link locants.
"""

import re

from rdkit import Chem

from ._common import UnsupportedStructure
from ._nucleoside_substituted import _NUCLEOTIDYL_STEM, _best_name

_LINK = Chem.MolFromSmarts("[PX4](=[O,S])([OX2H1,SX2H1])([OX2][#6])[OX2][#6]")
_SUGAR = Chem.MolFromSmarts("[#7]-[CH1]1[CH2,CH1][CH2,CH1][CH1]([CH2,CH3,CH1])O1")
_UNIT = re.compile(
    r"(2′-deoxy)?(adenosine|guanosine|inosine|xanthosine|cytidine|thymidine|uridine)"
)
_CACHE = {}


def _link_atoms(mol):
    links = []
    for p, oxo, acid, o1, c1, o2, c2 in (
        (m[0], m[1], m[2], m[3], m[4], m[5], m[6])
        for m in mol.GetSubstructMatches(_LINK)
    ):
        atom = mol.GetAtomWithIdx(p)
        if atom.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED:
            raise UnsupportedStructure("stereo phosphorus in an oligonucleotide link")
        sulfur = sum(mol.GetAtomWithIdx(i).GetAtomicNum() == 16 for i in (oxo, acid))
        links.append((p, (oxo, acid), o1, o2, sulfur))
    return links


def _unit_position(fragment, oxygen_orig):
    (oxygen,) = [a.GetIdx() for a in fragment.GetAtoms() if a.GetIntProp("orig") == oxygen_orig]
    carbon = next(n.GetIdx() for n in fragment.GetAtomWithIdx(oxygen).GetNeighbors())
    for match in fragment.GetSubstructMatches(_SUGAR):
        ring = {"2": match[2], "3": match[3], "5": match[5]}
        for position, index in ring.items():
            if index == carbon:
                return position
    raise UnsupportedStructure("a link not on a sugar carbon of a nucleoside")


def _compute(mol):
    links = _link_atoms(mol)
    if not links:
        return None
    rw = Chem.RWMol(mol)
    for atom in rw.GetAtoms():
        atom.SetIntProp("orig", atom.GetIdx())
    removed = set()
    for p, terminals, *_ in links:
        removed |= {p, *terminals}
    for idx in sorted(removed, reverse=True):
        rw.RemoveAtom(idx)
    link_oxygens = {i for _, _, o1, o2, _ in links for i in (o1, o2)}
    for atom in rw.GetAtoms():
        if atom.GetIntProp("orig") in link_oxygens:
            atom.SetNoImplicit(False)
    cut = rw.GetMol()
    Chem.SanitizeMol(cut)
    fragments = Chem.GetMolFrags(cut, asMols=True, sanitizeFrags=True)
    if len(fragments) != len(links) + 1:
        return None
    owner = {}
    for f_idx, fragment in enumerate(fragments):
        for atom in fragment.GetAtoms():
            owner[atom.GetIntProp("orig")] = f_idx
    names = []
    for fragment in fragments:
        name = _best_name(fragment)
        if name is None or not _UNIT.fullmatch(name):
            return None
        names.append(name)
    edges = []
    for _, _, o1, o2, sulfur in links:
        a, b = owner[o1], owner[o2]
        if a == b:
            return None
        edges.append(
            (a, _unit_position(fragments[a], o1), b, _unit_position(fragments[b], o2), sulfur)
        )
    degree = {}
    for a, _, b, _, _ in edges:
        degree[a] = degree.get(a, 0) + 1
        degree[b] = degree.get(b, 0) + 1
    ends = [u for u, d in degree.items() if d == 1]
    if len(ends) != 2 or any(d > 2 for d in degree.values()):
        return None
    candidates = []
    for start in ends:
        sequence, here, used = [start], start, set()
        steps = []
        while len(used) < len(edges):
            step = next(
                (
                    i
                    for i, (a, _, b, _, _) in enumerate(edges)
                    if i not in used and here in (a, b)
                ),
                None,
            )
            if step is None:
                return None
            used.add(step)
            a, pa, b, pb, sulfur = edges[step]
            if here == a:
                steps.append((pa, pb, sulfur))
                here = b
            else:
                steps.append((pb, pa, sulfur))
                here = a
            sequence.append(here)
        candidates.append((steps, sequence))
    steps, sequence = min(
        candidates, key=lambda c: ([(int(x), int(y)) for x, y, _ in c[0]], [names[i] for i in c[1]])
    )
    text = []
    for (first, second, sulfur), unit in zip(steps, sequence):
        deoxy, base = _UNIT.fullmatch(names[unit]).groups()
        acyl = (
            (deoxy or "")
            + ("-" if deoxy and sulfur else "")
            + ("P-thio" if sulfur else "")
            + _NUCLEOTIDYL_STEM[base]
        )
        text.append(f"{acyl}-({first}′→{second}′)-")
    return "".join(text) + names[sequence[-1]]


def oligonucleotide_name(mol):
    key = Chem.MolToSmiles(mol)
    if key not in _CACHE:
        try:
            _CACHE[key] = _compute(mol) if mol.HasSubstructMatch(_LINK) else None
        except (UnsupportedStructure, ValueError):
            _CACHE[key] = None
    return _CACHE[key]
