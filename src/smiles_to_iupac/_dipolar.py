"""Dipolar compounds: ylides, imides and oxides of the allyl and propargyl types, azoxy compounds, S-oxides (P-74.2).

Three constructions follow the clauses. A cationic centre bonded to a carbon or nitrogen anion is cited as a cationic
prefix on the anionic parent ('2-(trimethylphosphaniumyl)propan-2-ide', P-74.2.1.1, P-74.2.2.1.3 to .1.7, .2.1.3). Cation
and anion on one chain of identical heteroatoms are the 'ium'/'ide' pair of that parent hydride ('hydrazin-2-ium-1-ide',
P-74.2.2.1.1, .1.2, .1.6, .2.1.1). Azoxy compounds take the class term 'oxide' (P-68.3.1.3.3.1), and a thioaldehyde or
thioketone S-oxide the lambda-convention name 'ylidene-lambda4-sulfanone' (P-74.2.2.1.8).
"""

import re

from rdkit import Chem

from ._common import UnsupportedStructure, adjacency, alpha_sort_key
from ._substituents import FORCED_BRANCH_NAMES, format_mononuclear_prefixes, name_branch

_ONIUM_STEMS = {7: "azanium", 8: "oxidanium", 15: "phosphanium", 16: "sulfanium", 34: "selanium", 52: "telluranium", 33: "arsanium"}
_CHAIN_STEMS = {(7, 2): "hydrazin", (8, 2): "dioxidan", (7, 3): "triaz"}


def _group(mol, graph, root, parent):
    """(name, is_compound) of the group at `root` hanging off `parent`; a multiple bond to the parent is named from a
    carbon stand-in so that an ylidene or ylidyne group is not mistaken for a carbonyl-type group."""
    order = mol.GetBondBetweenAtoms(root, parent).GetBondTypeAsDouble()
    if order == 1.0:
        return name_branch(graph, root, parent, {}, frozenset(), mol=mol)
    surrogate = Chem.RWMol(mol)
    atom = surrogate.GetAtomWithIdx(parent)
    atom.SetAtomicNum(6)
    atom.SetFormalCharge(0)
    standin = surrogate.GetMol()
    standin.UpdatePropertyCache(strict=False)
    return name_branch(adjacency(standin), root, parent, {}, frozenset(), mol=standin, unsaturated=True)


def _side(graph, start, blocked):
    seen, stack = {start}, [start]
    while stack:
        for n in graph[stack.pop()]:
            if n != blocked and n not in seen:
                seen.add(n)
                stack.append(n)
    return seen


def _cation_prefix(mol, graph, x, y):
    """Prefix 'yl' name of the cationic group at `x` that is bonded to the anion at `y` (method 1 of P-73.6)."""
    from .core import smiles_to_iupac

    order = int(mol.GetBondBetweenAtoms(x, y).GetBondTypeAsDouble())
    editable = Chem.RWMol(mol)
    editable.RemoveBond(x, y)
    center = editable.GetAtomWithIdx(x)
    center.SetNoImplicit(True)
    center.SetNumExplicitHs(mol.GetAtomWithIdx(x).GetTotalNumHs() + order)
    keep = _side(graph, x, y)
    if any(y in graph[a] for a in keep if a != x):
        raise UnsupportedStructure("the ionic centres lie in one ring")
    for index in sorted(set(range(mol.GetNumAtoms())) - keep, reverse=True):
        editable.RemoveAtom(index)
    fragment = editable.GetMol()
    Chem.SanitizeMol(fragment)
    try:
        name = smiles_to_iupac(Chem.MolToSmiles(fragment))
        single = re.fullmatch(r"(.+?ylidene)(oxidanium|sulfanium|selanium|telluranium|chloranium|bromanium|iodanium)", name)
        if single and "," not in single.group(1) and single.group(1)[0] not in "([":
            return f"({single.group(1)}){single.group(2)}yl"
        if name.endswith("ium"):
            return name + "yl"
    except UnsupportedStructure:
        pass
    stem = _ONIUM_STEMS.get(mol.GetAtomWithIdx(x).GetAtomicNum())
    if stem is None:
        raise UnsupportedStructure("an unsupported cationic centre")
    inner = adjacency(fragment)
    centre_index = next(a.GetIdx() for a in fragment.GetAtoms() if a.GetFormalCharge() > 0)
    names = [
        (name, compound or bool(re.search(r"[\d(]", name)))
        for name, compound in (_group(fragment, inner, n, centre_index) for n in inner[centre_index])
    ]
    if len(names) == 1 and names[0][1]:
        return f"({names[0][0]}){stem}yl"
    return (format_mononuclear_prefixes(names) if names else "") + stem + "yl"


def _enclose(text):
    return f"[{text}]" if "(" in text else f"({text})"


def _smiles_with_order(mol):
    smiles = Chem.MolToSmiles(mol, canonical=False)
    order = [int(i) for i in mol.GetProp("_smilesAtomOutputOrder").strip("[]").split(",") if i.strip()]
    return smiles, {original: position for position, original in enumerate(order)}


def _carbon_anion_name(mol, graph, x, y, prefix):
    """The anion parent with the cationic group stood in for by a methoxy group whose name is then replaced."""
    from .core import smiles_to_iupac

    editable = Chem.RWMol(mol)
    for index in sorted(_side(graph, x, y), reverse=True):
        editable.RemoveAtom(index)
    anion = y - sum(1 for index in _side(graph, x, y) if index < y)
    oxygen = editable.AddAtom(Chem.Atom(8))
    carbon = editable.AddAtom(Chem.Atom(6))
    editable.AddBond(anion, oxygen, Chem.BondType.SINGLE)
    editable.AddBond(oxygen, carbon, Chem.BondType.SINGLE)
    stand_in = editable.GetMol()
    Chem.SanitizeMol(stand_in)
    smiles, position = _smiles_with_order(stand_in)
    token = FORCED_BRANCH_NAMES.set((stand_in.GetNumAtoms(), {position[oxygen]: (prefix, True)}))
    try:
        return smiles_to_iupac(smiles)
    finally:
        FORCED_BRANCH_NAMES.reset(token)


def _ylide_pair(mol, graph, cation, anion):
    prefix = _cation_prefix(mol, graph, cation.GetIdx(), anion.GetIdx())
    if anion.GetAtomicNum() == 6:
        return _carbon_anion_name(mol, graph, cation.GetIdx(), anion.GetIdx(), prefix)
    if anion.GetAtomicNum() == 7 and anion.GetDegree() == 2:
        from .core import smiles_to_iupac

        (root,) = [n for n in graph[anion.GetIdx()] if n != cation.GetIdx()]
        editable = Chem.RWMol(mol)
        for index in sorted(_side(graph, cation.GetIdx(), anion.GetIdx()), reverse=True):
            editable.RemoveAtom(index)
        remainder = editable.GetMol()
        for atom in remainder.GetAtoms():
            if atom.GetAtomicNum() == 7 and atom.GetFormalCharge() < 0:
                atom.SetNoImplicit(True)
                atom.SetNumExplicitHs(1)
        Chem.SanitizeMol(remainder)
        parent = smiles_to_iupac(Chem.MolToSmiles(remainder))
        if not parent.endswith("aminide") or "-" in parent.replace("aminide", "") or parent.startswith("N"):
            raise UnsupportedStructure("the anionic amine parent is not a plain alkanaminide")
        return f"N-{_enclose(prefix)}{parent}"
    raise UnsupportedStructure("an unsupported anionic centre")


def _remote_pair(mol, graph, cation, anion):
    """Anionic and cationic centres on different parent structures (P-74.1.3): the cationic group is cited as a prefix of
    the parent that holds the anion, here a boranuide or a carbanide joined to it through an acyclic chain."""
    if cation.IsInRing() or anion.IsInRing() or anion.GetAtomicNum() not in (5, 6):
        raise UnsupportedStructure("not a remote pair of an acyclic onium group and a boranuide or carbanide")
    toward = [n for n in graph[cation.GetIdx()] if anion.GetIdx() in _side(graph, n, cation.GetIdx())]
    if len(toward) != 1:
        raise UnsupportedStructure("the anionic parent is not reached through one bond of the cationic group")
    attach = toward[0]
    prefix = _cation_prefix(mol, graph, cation.GetIdx(), attach)
    return _carbon_anion_name(mol, graph, cation.GetIdx(), attach, prefix)


def _chain_pair(mol, graph, cation, anion):
    """Cation and anion on a chain of identical heteroatoms: the 'ium'/'ide' pair of hydrazine, dioxidane or triazene."""
    z = anion.GetAtomicNum()
    chain = [anion.GetIdx(), cation.GetIdx()]
    if z == 7 and cation.GetAtomicNum() == 7:
        tail = [n.GetIdx() for n in cation.GetNeighbors() if n.GetAtomicNum() == 7 and n.GetIdx() != anion.GetIdx()]
        if len(tail) == 1 and mol.GetBondBetweenAtoms(cation.GetIdx(), tail[0]).GetBondTypeAsDouble() == 2.0:
            chain.append(tail[0])
    if cation.GetAtomicNum() != z:
        raise UnsupportedStructure("not a zwitterionic chain of identical heteroatoms")
    stem = _CHAIN_STEMS.get((z, len(chain)))
    if stem is None:
        raise UnsupportedStructure("an unsupported heteroatom chain")
    entries = []
    for locant, atom in enumerate(chain, 1):
        for n in graph[atom]:
            if n in chain:
                continue
            if mol.GetAtomWithIdx(n).GetAtomicNum() == 1:
                continue
            name, compound = _group(mol, graph, n, atom)
            entries.append((locant, name, compound))
    grouped = {}
    for locant, name, compound in entries:
        grouped.setdefault((name, compound), []).append(locant)
    from ._substituents import format_substituent_prefixes

    prefixes = format_substituent_prefixes({n: {"locants": sorted(l), "compound": c} for (n, c), l in grouped.items()})
    en = "-2-en" if len(chain) == 3 else ""
    base = f"{stem}{en}-2-ium-1-ide"
    return f"{prefixes}{base}"


def _azoxy(mol, graph):
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 7 and atom.GetFormalCharge() == 1:
            oxygens = [n for n in atom.GetNeighbors() if n.GetAtomicNum() == 8 and n.GetFormalCharge() == -1 and n.GetDegree() == 1]
            partner = [n for n in atom.GetNeighbors() if n.GetAtomicNum() == 7 and mol.GetBondBetweenAtoms(atom.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 2.0]
            if len(oxygens) == 1 and len(partner) == 1 and atom.GetDegree() == 3:
                return atom, partner[0], oxygens[0]
    return None


def _azoxy_name(mol, graph):
    oxide_nitrogen, other_nitrogen, oxygen = _azoxy(mol, graph)
    if any(b.GetStereo() != Chem.BondStereo.STEREONONE for b in mol.GetBonds()):
        raise UnsupportedStructure("the configuration of the azoxy double bond is not cited yet")
    groups = []
    for nitrogen in (oxide_nitrogen, other_nitrogen):
        (root,) = [n for n in graph[nitrogen.GetIdx()] if n not in (oxygen.GetIdx(), oxide_nitrogen.GetIdx(), other_nitrogen.GetIdx())]
        groups.append(_group(mol, graph, root, nitrogen.GetIdx()))
    if groups[0] == groups[1]:
        text = format_mononuclear_prefixes([groups[0], groups[1]])
        return f"{text}diazene oxide"
    ordered = sorted(enumerate(groups), key=lambda item: alpha_sort_key(item[1][0]))
    first, second = ordered
    oxide_locant = 1 if first[0] == 0 else 2
    left, right = (_enclose(g[0]) if g[1] else g[0] for g in (first[1], second[1]))
    return f"1-{left}-2-{right}diazene {oxide_locant}-oxide"


def _sulfur_oxide(mol, graph):
    sulfurs = [a for a in mol.GetAtoms() if a.GetAtomicNum() == 16 and a.GetFormalCharge() == 1]
    if len(sulfurs) != 1 or mol.GetNumAtoms() < 4:
        return None
    sulfur = sulfurs[0]
    oxygen = [n for n in sulfur.GetNeighbors() if n.GetAtomicNum() == 8 and n.GetFormalCharge() == -1 and n.GetDegree() == 1]
    carbon = [n for n in sulfur.GetNeighbors() if n.GetAtomicNum() == 6 and mol.GetBondBetweenAtoms(sulfur.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 2.0]
    if len(oxygen) != 1 or len(carbon) != 1 or sulfur.GetDegree() != 2:
        return None
    name, _ = _group(mol, graph, carbon[0].GetIdx(), sulfur.GetIdx())
    return f"{name}-λ4-sulfanone"


def dipolar_name(mol):
    if len(Chem.GetMolFrags(mol)) != 1:
        return None
    from ._radical_hetero import isodiazene_name

    isodiazene = isodiazene_name(mol) if not any(a.GetNumRadicalElectrons() for a in mol.GetAtoms()) else None
    if isodiazene is not None:
        return isodiazene
    charged = [a for a in mol.GetAtoms() if a.GetFormalCharge()]
    if len(charged) != 2 or sorted(a.GetFormalCharge() for a in charged) != [-1, 1]:
        return None
    cation = next(a for a in charged if a.GetFormalCharge() > 0)
    anion = next(a for a in charged if a.GetFormalCharge() < 0)
    link = mol.GetBondBetweenAtoms(cation.GetIdx(), anion.GetIdx())
    if any(a.GetIsotope() or a.GetNumRadicalElectrons() for a in mol.GetAtoms()):
        return None
    graph = adjacency(mol)
    if link is None:
        if cation.GetAtomicNum() not in _ONIUM_STEMS:
            return None
        try:
            return _remote_pair(mol, graph, cation, anion)
        except (UnsupportedStructure, ValueError, KeyError):
            return None
    if link.GetBondTypeAsDouble() != 1.0:
        return None
    try:
        if _azoxy(mol, graph) is not None:
            return _azoxy_name(mol, graph)
        sulfur = _sulfur_oxide(mol, graph)
        if sulfur is not None:
            return sulfur
        if cation.GetAtomicNum() == anion.GetAtomicNum() and anion.GetAtomicNum() in (7, 8):
            return _chain_pair(mol, graph, cation, anion)
        if anion.GetAtomicNum() in (6, 7) and cation.GetAtomicNum() in _ONIUM_STEMS:
            return _ylide_pair(mol, graph, cation, anion)
    except (UnsupportedStructure, ValueError, KeyError):
        return None
    return None


def has_dipolar_shape(mol) -> bool:
    return dipolar_name(mol) is not None


def name_dipolar(mol) -> str:
    return dipolar_name(mol)
