"""Radicals on parent hydrides of elements other than carbon (P-71.2.1.1, P-71.2.1.2, P-71.2.2, P-71.2.3).

A mononuclear centre with carbon substituents is named through the hydride that results when its radical electrons are
filled with hydrogens, and the final 'ane' takes the suffix 'yl', 'ylidene' or 'ylidyne' ('silane' gives 'silyl',
'phosphane' gives 'phosphanyl'). An unsubstituted chain of identical atoms keeps the locants of its radical centres
('trisilan-2-yl', 'hydrazine-1,2-diyl').
"""

from rdkit import Chem

from ._common import UnsupportedStructure

_SUFFIX = {1: "yl", 2: "ylidene", 3: "ylidyne"}
_GROUP_14 = {14, 32, 50, 82}
_HYDRIDE_NAMES = {"ammonia": "azane", "water": "oxidane", "hydrogen sulfide": "sulfane", "hydrogen selenide": "selane", "hydrogen telluride": "tellane"}
_RETAINED = {("oxidane", 1): "hydroxyl", ("dioxidane", 1): "hydroperoxyl"}
_MULTIPLIER = {2: "di", 3: "tri", 4: "tetra"}


def _hydride_name(smiles):
    from .core import smiles_to_iupac

    name = smiles_to_iupac(smiles)
    return _HYDRIDE_NAMES.get(name, name)


def _healed(mol, radicals):
    editable = Chem.RWMol(mol)
    for atom in radicals:
        target = editable.GetAtomWithIdx(atom.GetIdx())
        target.SetNoImplicit(True)
        target.SetNumExplicitHs(atom.GetTotalNumHs() + atom.GetNumRadicalElectrons())
        target.SetNumRadicalElectrons(0)
    healed = editable.GetMol()
    Chem.SanitizeMol(healed)
    return healed


def _mononuclear(mol, radical):
    count = radical.GetNumRadicalElectrons()
    if count not in _SUFFIX or radical.GetFormalCharge() or radical.GetIsotope() or radical.IsInRing():
        return None
    others = [a for a in mol.GetAtoms() if a.GetIdx() != radical.GetIdx()]
    if any(a.GetAtomicNum() != 6 or a.GetFormalCharge() or a.GetNumRadicalElectrons() for a in others):
        return None
    from ._hetero_prefixes import MONONUCLEAR_HYDRIDES

    standard = MONONUCLEAR_HYDRIDES.get(radical.GetAtomicNum(), (None, None, None))[2]
    if standard is not None and radical.GetTotalDegree() + count > standard:
        return None
    if any(mol.GetBondBetweenAtoms(radical.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() != 1.0 for n in radical.GetNeighbors()):
        return None
    try:
        name = _hydride_name(Chem.MolToSmiles(_healed(mol, [radical])))
    except (UnsupportedStructure, Chem.rdchem.AtomValenceException):
        return None
    if not name.endswith("ane"):
        return None
    retained = _RETAINED.get((name, count))
    if retained:
        return retained
    suffix = _SUFFIX[count]
    if radical.GetAtomicNum() in _GROUP_14:
        return name[:-3] + suffix
    return name[:-1] + suffix


def _chain(mol, radicals):
    elements = {a.GetAtomicNum() for a in mol.GetAtoms()}
    if len(elements) != 1 or elements == {6} or mol.GetRingInfo().NumRings() or mol.GetNumAtoms() < 2:
        return None
    if any(a.GetFormalCharge() or a.GetIsotope() or any(b.GetBondTypeAsDouble() != 1.0 for b in a.GetBonds()) for a in mol.GetAtoms()):
        return None
    if len(radicals) > 2 or any(a.GetDegree() > 2 for a in mol.GetAtoms()):
        return None
    ends = [a.GetIdx() for a in mol.GetAtoms() if a.GetDegree() == 1]
    order, previous = [ends[0]], None
    while len(order) < mol.GetNumAtoms():
        nxt = [n.GetIdx() for n in mol.GetAtomWithIdx(order[-1]).GetNeighbors() if n.GetIdx() != previous]
        previous = order[-1]
        order.append(nxt[0])
    try:
        base = _hydride_name(Chem.MolToSmiles(_healed(mol, radicals)))
    except (UnsupportedStructure, Chem.rdchem.AtomValenceException):
        return None
    if not base.endswith("ane") and base != "hydrazine":
        return None
    valences = {a.GetIdx(): a.GetNumRadicalElectrons() for a in radicals}
    best = None
    for sequence in (order, order[::-1]):
        key = [i + 1 for i, atom in enumerate(sequence) if atom in valences]
        if best is None or key < best[0]:
            best = (key, [valences[a] for a in sequence if a in valences])
    locants, counts = best
    if len(set(counts)) != 1 or counts[0] not in _SUFFIX:
        return None
    suffix, number = _SUFFIX[counts[0]], len(counts)
    if (base, number) in _RETAINED:
        return _RETAINED[(base, number)]
    if base == "hydrazine" and number == 1:
        return "hydrazin" + suffix
    joined = ",".join(map(str, locants))
    if number == 1:
        if mol.GetNumAtoms() == 2:
            return base[:-1] + suffix
        return f"{base[:-1]}-{joined}-{suffix}"
    return f"{base}-{joined}-{_MULTIPLIER[number]}{suffix}"


def _substituted_chain(mol, radicals):
    """Radical centres of one kind on an unbranched chain of one heteroatom that carries other groups: the chain is the
    parent (P-44.1.2.2), its radical locants come first and every other group is a prefix (P-71.2.3)."""
    from ._common import adjacency, group_substituents, halogen_substituents, substituent_locant_set_and_citation
    from ._hydride_chain import _STEMS
    from ._numerals import multiplying_prefix
    from ._substituents import format_substituent_prefixes, name_branch

    z = radicals[0].GetAtomicNum()
    if z not in _STEMS or any(a.GetAtomicNum() != z or a.GetNumRadicalElectrons() not in _SUFFIX for a in radicals):
        return None
    if any(a.GetFormalCharge() or a.GetIsotope() or a.IsInRing() for a in mol.GetAtoms()) or any(
        b.GetBondTypeAsDouble() != 1.0 for b in mol.GetBonds()
    ):
        return None
    chain_atoms = {radicals[0].GetIdx()}
    stack = [radicals[0]]
    while stack:
        for n in stack.pop().GetNeighbors():
            if n.GetAtomicNum() == z and n.GetIdx() not in chain_atoms:
                chain_atoms.add(n.GetIdx())
                stack.append(n)
    ends = [i for i in chain_atoms if sum(n.GetIdx() in chain_atoms for n in mol.GetAtomWithIdx(i).GetNeighbors()) <= 1]
    if len(chain_atoms) < 2 or len(ends) != 2 or not {a.GetIdx() for a in radicals} <= chain_atoms:
        return None
    if any(sum(n.GetIdx() in chain_atoms for n in mol.GetAtomWithIdx(i).GetNeighbors()) > 2 for i in chain_atoms):
        return None
    order, previous = [ends[0]], None
    while len(order) < len(chain_atoms):
        order.append(next(n.GetIdx() for n in mol.GetAtomWithIdx(order[-1]).GetNeighbors() if n.GetIdx() in chain_atoms and n.GetIdx() != previous))
        previous = order[-2]
    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    valence_of = {a.GetIdx(): a.GetNumRadicalElectrons() for a in radicals}
    best = None
    for candidate in (order, order[::-1]):
        substituents = {}
        for position, atom in enumerate(candidate, start=1):
            for n in graph[atom]:
                if n in chain_atoms:
                    continue
                if mol.GetBondBetweenAtoms(atom, n).GetBondTypeAsDouble() != 1.0:
                    return None
                substituents.setdefault(position, []).append(name_branch(graph, n, atom, halogens, frozenset(), mol=mol, unsaturated=True))
        grouped = group_substituents(substituents)
        locant_set, _, citation = substituent_locant_set_and_citation(grouped)
        by_valence = {v: [i + 1 for i, atom in enumerate(candidate) if valence_of.get(atom) == v] for v in _SUFFIX}
        key = (sorted(sum(by_valence.values(), [])), by_valence[1], by_valence[2], by_valence[3], locant_set, citation)
        if best is None or key < best[0]:
            best = (key, grouped, by_valence)
    _, grouped, by_valence = best
    parent = "hydrazine" if z == 7 and len(order) == 2 else f"{multiplying_prefix(len(order))}{_STEMS[z]}"
    pieces = [
        f"{','.join(map(str, locants))}-{'' if len(locants) == 1 else _MULTIPLIER[len(locants)]}{_SUFFIX[v]}"
        for v, locants in by_valence.items()
        if locants
    ]
    stem = parent if len(by_valence[min(v for v, l in by_valence.items() if l)]) > 1 else parent[:-1]
    return f"{format_substituent_prefixes(grouped, omit_locants=False)}{stem}-{'-'.join(pieces)}"


def _isodiazene(mol):
    """R2N-N: (an N-N compound with a divalent terminal nitrogen) is the parent radical hydrazinylidene (P-68.3.1.3.7)."""
    from ._cited_group import cited_group
    from ._common import adjacency
    from ._substituents import format_mononuclear_prefixes

    radicals = [a for a in mol.GetAtoms() if a.GetNumRadicalElectrons()]
    if len(radicals) != 1 or radicals[0].GetAtomicNum() != 7 or radicals[0].GetNumRadicalElectrons() != 2 or radicals[0].GetDegree() != 1:
        return None
    terminal = radicals[0]
    (root,) = terminal.GetNeighbors()
    if root.GetAtomicNum() != 7 or root.GetFormalCharge() or terminal.GetFormalCharge() or root.IsInRing():
        return None
    graph = adjacency(mol)
    substituents = [n.GetIdx() for n in root.GetNeighbors() if n.GetIdx() != terminal.GetIdx()]
    if any(mol.GetAtomWithIdx(n).GetAtomicNum() != 6 or mol.GetBondBetweenAtoms(root.GetIdx(), n).GetBondTypeAsDouble() != 1.0 for n in substituents):
        return None
    if root.GetTotalNumHs() + len(substituents) != 2:
        return None
    if any(a.GetAtomicNum() != 6 for a in mol.GetAtoms() if a.GetIdx() not in (terminal.GetIdx(), root.GetIdx())):
        return None
    try:
        names = [cited_group(mol, graph, n, root.GetIdx()) for n in substituents]
    except UnsupportedStructure:
        return None
    return (format_mononuclear_prefixes(names) if names else "") + "hydrazinylidene"


def isodiazene_name(mol):
    """The name of an isodiazene given as R2N-N: or as the zwitterion R2N(+)=N(-)."""
    if any(a.GetNumRadicalElectrons() for a in mol.GetAtoms()):
        return _isodiazene(mol)
    charged = [a for a in mol.GetAtoms() if a.GetFormalCharge()]
    if len(charged) != 2 or sorted(a.GetFormalCharge() for a in charged) != [-1, 1]:
        return None
    cation = next(a for a in charged if a.GetFormalCharge() > 0)
    anion = next(a for a in charged if a.GetFormalCharge() < 0)
    bond = mol.GetBondBetweenAtoms(cation.GetIdx(), anion.GetIdx())
    if bond is None or bond.GetBondTypeAsDouble() != 2.0 or cation.GetAtomicNum() != 7 or anion.GetAtomicNum() != 7 or anion.GetDegree() != 1:
        return None
    editable = Chem.RWMol(mol)
    editable.GetBondBetweenAtoms(cation.GetIdx(), anion.GetIdx()).SetBondType(Chem.BondType.SINGLE)
    for atom in (cation, anion):
        target = editable.GetAtomWithIdx(atom.GetIdx())
        target.SetFormalCharge(0)
        target.SetNoImplicit(True)
        target.SetNumExplicitHs(atom.GetTotalNumHs())
    editable.GetAtomWithIdx(anion.GetIdx()).SetNumRadicalElectrons(2)
    neutral = editable.GetMol()
    neutral.UpdatePropertyCache(strict=False)
    return _isodiazene(neutral)


def hetero_radical_name(mol):
    isodiazene = _isodiazene(mol)
    if isodiazene is not None:
        return isodiazene
    if all(a.GetAtomicNum() == 6 for a in mol.GetAtoms()) or len(Chem.GetMolFrags(mol)) > 1:
        return None
    radicals = [a for a in mol.GetAtoms() if a.GetNumRadicalElectrons()]
    if not radicals:
        return None
    if len(radicals) > 1:
        return _chain(mol, radicals) or _substituted_chain(mol, radicals)
    if sum(a.GetNumRadicalElectrons() for a in radicals) > 3:
        return None
    if mol.GetNumAtoms() == 1 or radicals[0].GetAtomicNum() != 6:
        found = _mononuclear(mol, radicals[0])
        if found is not None:
            return found
    return _chain(mol, radicals)
