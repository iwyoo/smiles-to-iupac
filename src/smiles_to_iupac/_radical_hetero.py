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
_MULTIPLIER = {2: "di", 3: "tri"}


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
        return f"{base[:-1]}-{joined}-{suffix}"
    return f"{base}-{joined}-{_MULTIPLIER[number]}{suffix}"


def hetero_radical_name(mol):
    if all(a.GetAtomicNum() == 6 for a in mol.GetAtoms()) or len(Chem.GetMolFrags(mol)) > 1:
        return None
    radicals = [a for a in mol.GetAtoms() if a.GetNumRadicalElectrons()]
    if not radicals or sum(a.GetNumRadicalElectrons() for a in radicals) > 3:
        return None
    if len(radicals) == 1 and (mol.GetNumAtoms() == 1 or radicals[0].GetAtomicNum() != 6):
        found = _mononuclear(mol, radicals[0])
        if found is not None:
            return found
    return _chain(mol, radicals)
