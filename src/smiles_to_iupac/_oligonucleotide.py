"""Oligonucleotides: nucleosides joined by phosphodiester links (P-106.3.4, P-106.3.5).

Cutting out each link phosphorus leaves nucleosides; every unit but the last becomes its acyl group
('2′-deoxyguanylyl', 'P-thio' before the base stem for a sulfur-bearing link) plus '(3′→5′)', read from the end
that gives the lower link locants.
"""

import re

from rdkit import Chem

from ._common import UnsupportedStructure
from ._nucleoside_substituted import _NUCLEOTIDE_STEM, _NUCLEOTIDYL_STEM, _best_name

_STEM_OF_ACID = {_NUCLEOTIDE_STEM[nucleoside]: stem for nucleoside, stem in _NUCLEOTIDYL_STEM.items()}

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


_NUCLEOTIDE_ACID = re.compile(r"(?P<head>.*?)-?\d+′-(?P<stem>[a-z]+ic) acid")


def _nucleotidyl_prefix(fragment, oxygen_orig, terminals):
    """P-106.3.3: the unit as a nucleotide whose phosphate sits on the link oxygen, its '-ic acid' ending turned into
    '-yl' and the phosphate locant left to the (x′→y′) bracket; None when the unit is not named as one nucleotide."""
    editable = Chem.RWMol(fragment)
    (oxygen,) = [a.GetIdx() for a in editable.GetAtoms() if a.GetIntProp("orig") == oxygen_orig]
    phosphorus = editable.AddAtom(Chem.Atom(15))
    editable.AddBond(oxygen, phosphorus, Chem.BondType.SINGLE)
    oxo, acid = terminals
    editable.AddBond(phosphorus, editable.AddAtom(Chem.Atom(oxo)), Chem.BondType.DOUBLE)
    editable.AddBond(phosphorus, editable.AddAtom(Chem.Atom(acid)), Chem.BondType.SINGLE)
    editable.AddBond(phosphorus, editable.AddAtom(Chem.Atom(8)), Chem.BondType.SINGLE)
    nucleotide = editable.GetMol()
    Chem.SanitizeMol(nucleotide)
    match = _NUCLEOTIDE_ACID.fullmatch(_best_name(nucleotide) or "")
    if match is None or match["stem"] not in _STEM_OF_ACID:
        return None
    return match["head"] + _STEM_OF_ACID[match["stem"]]


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
    names = [_best_name(fragment) for fragment in fragments]
    if any(name is None for name in names):
        return None
    edges = []
    for (_, (oxo, acid), o1, o2, sulfur) in links:
        a, b = owner[o1], owner[o2]
        if a == b:
            return None
        edges.append(
            (a, _unit_position(fragments[a], o1), b, _unit_position(fragments[b], o2), sulfur, (o1, o2, mol.GetAtomWithIdx(oxo).GetAtomicNum(), mol.GetAtomWithIdx(acid).GetAtomicNum()))
        )
    degree = {}
    for a, _, b, _, _, _ in edges:
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
                    for i, (a, _, b, _, _, _) in enumerate(edges)
                    if i not in used and here in (a, b)
                ),
                None,
            )
            if step is None:
                return None
            used.add(step)
            a, pa, b, pb, sulfur, atoms = edges[step]
            if here == a:
                steps.append((pa, pb, sulfur, atoms[0], atoms[2:]))
                here = b
            else:
                steps.append((pb, pa, sulfur, atoms[1], atoms[2:]))
                here = a
            sequence.append(here)
        candidates.append((steps, sequence))
    valid = []
    for steps, sequence in candidates:
        text = _chain_text(fragments, names, steps, sequence)
        if text is not None:
            key = ([(int(x), int(y)) for x, y, *_ in steps], [names[i] for i in sequence])
            valid.append((key, text + names[sequence[-1]]))
    return min(valid)[1] if valid else None


def _chain_text(fragments, names, steps, sequence):
    text = []
    for (first, second, sulfur, oxygen, terminals), unit in zip(steps, sequence):
        plain = _UNIT.fullmatch(names[unit])
        if plain:
            deoxy, base = plain.groups()
            acyl = (
                (deoxy or "")
                + ("-" if deoxy and sulfur else "")
                + ("P-thio" if sulfur else "")
                + _NUCLEOTIDYL_STEM[base]
            )
        else:
            acyl = _nucleotidyl_prefix(fragments[unit], oxygen, terminals)
            if acyl is None:
                return None
        text.append(f"{acyl}-({first}′→{second}′)-")
    return "".join(text)


def oligonucleotide_name(mol):
    key = Chem.MolToSmiles(mol)
    if key not in _CACHE:
        try:
            _CACHE[key] = _compute(mol) if mol.HasSubstructMatch(_LINK) else None
        except (UnsupportedStructure, ValueError):
            _CACHE[key] = None
    return _CACHE[key]
