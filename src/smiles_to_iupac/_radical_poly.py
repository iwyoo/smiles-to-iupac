"""Radicals with several centres (P-71.2.3, P-71.2.4, P-71.4, P-71.6, P-71.7): the parent holding every centre carries
the suffixes 'yl', 'ylidene' and 'ylidyne' in that order, and identical monovalent groups on a shared skeleton are
named as an assembly of parent radicals, 'cyclopropane-1,2-diyl' with two 'methyl' groups."""

import re

from rdkit import Chem

from ._common import UnsupportedStructure, adjacency, specified_stereo_elements
from ._hetero_prefixes import is_functional_carbon
from ._numerals import multiplying_prefix

_MAX_ATOMS = 80
_BOND = {1: Chem.BondType.SINGLE, 2: Chem.BondType.DOUBLE, 3: Chem.BondType.TRIPLE}
_CHALCOGEN_WORDS = {16: "sulfanyl", 34: "selanyl", 52: "tellanyl"}


def _plain(mol):
    return (
        len(Chem.GetMolFrags(mol)) == 1
        and mol.GetNumAtoms() <= _MAX_ATOMS
        and not any(a.GetIsotope() for a in mol.GetAtoms())
        and sum(a.GetFormalCharge() for a in mol.GetAtoms()) == 0
        and not specified_stereo_elements(mol)
    )


def _capped(mol, caps, removed=()):
    """`mol` without the atoms `removed`, each (atom, valence) of `caps` ended by a dummy hydrogen bonded by that many
    bonds; (capped mol, [(atom, dummy, valence)]) in the numbering of the capped molecule."""
    editable = Chem.RWMol(mol)
    dummies = []
    for atom, valence in caps:
        target = editable.GetAtomWithIdx(atom)
        target.SetNumRadicalElectrons(0)
        dummies.append((atom, editable.AddAtom(Chem.Atom(1)), valence))
    for atom, dummy, valence in dummies:
        editable.AddBond(atom, dummy, _BOND[valence])
    keep = [i for i in range(mol.GetNumAtoms()) if i not in set(removed)]
    renumber = {old: new for new, old in enumerate(keep + [d for _, d, _ in dummies])}
    for index in sorted(set(removed), reverse=True):
        editable.RemoveAtom(index)
    out = editable.GetMol()
    out.UpdatePropertyCache(strict=False)
    return out, [(renumber[a], renumber[d], v) for a, d, v in dummies]


def _context(mol):
    from ._chain_multiplicative import _make_context

    return _make_context(mol, adjacency(mol))


def _one_acyclic_carbon_run(mol, atoms):
    seen, stack = {atoms[0]}, [atoms[0]]
    while stack:
        for n in mol.GetAtomWithIdx(stack.pop()).GetNeighbors():
            if n.GetAtomicNum() == 6 and not n.IsInRing() and n.GetIdx() not in seen:
                seen.add(n.GetIdx())
                stack.append(n.GetIdx())
    return set(atoms) <= seen


def _component(mol, attachments):
    """Part naming the skeleton that holds every attachment, or None."""
    from ._chain_multiplicative import _component_part, _ring_system_of
    from ._multiplicative_linker import DecompositionRejected

    atoms = [a for a, _, _ in attachments]
    ring_atoms = [a for a in atoms if mol.GetAtomWithIdx(a).IsInRing()]
    if ring_atoms:
        system = _ring_system_of(mol, ring_atoms[0])
        if len(ring_atoms) != len(atoms) or not set(atoms) <= set(system):
            return None
        kind, members = "ring", sorted(system)
    else:
        if any(mol.GetAtomWithIdx(a).GetAtomicNum() != 6 for a in atoms):
            return None
        if not _one_acyclic_carbon_run(mol, atoms) or any(
            a.GetAtomicNum() == 6 and not a.IsInRing() and is_functional_carbon(mol, a.GetIdx()) for a in mol.GetAtoms()
        ):
            return None
        kind, members = "carbon", atoms
    graph = adjacency(mol)
    try:
        return _component_part(mol, graph, kind, members, attachments, _context(mol), None)
    except (DecompositionRejected, UnsupportedStructure):
        return None


def _benzene(text):
    return re.sub(r"-?(\d[\d,]*)-phenylene", r"benzene-\1-diyl", text)


def _all_on_one_parent(mol, centres):
    capped, attachments = _capped(mol, [(a.GetIdx(), a.GetNumRadicalElectrons()) for a in centres])
    if any(capped.GetAtomWithIdx(a).GetIsAromatic() for a, _, _ in attachments):
        return None
    part = _component(capped, attachments)
    return None if part is None else _benzene(part.text)


def _unit(mol, centre):
    """(word, atoms from the centre toward the skeleton, skeleton atom) of a monovalent radical group, or None."""
    z = centre.GetAtomicNum()
    if centre.GetNumRadicalElectrons() != 1 or centre.GetDegree() != 1 or centre.IsInRing():
        return None
    (neighbour,) = centre.GetNeighbors()
    if mol.GetBondBetweenAtoms(centre.GetIdx(), neighbour.GetIdx()).GetBondTypeAsDouble() != 1.0:
        return None
    if z == 6:
        if centre.GetTotalNumHs() != 2 or neighbour.GetAtomicNum() != 6:
            return None
        return "methyl", [centre.GetIdx()], neighbour.GetIdx()
    if z == 7:
        if centre.GetTotalNumHs() != 1 or neighbour.GetAtomicNum() != 6:
            return None
        return "aminyl", [centre.GetIdx()], neighbour.GetIdx()
    if z not in (8, *_CHALCOGEN_WORDS) or centre.GetTotalNumHs():
        return None
    chain = [centre.GetIdx()]
    current = neighbour
    while current.GetAtomicNum() == z:
        if current.GetDegree() != 2 or current.GetTotalNumHs() or current.GetNumRadicalElectrons() or current.IsInRing():
            return None
        chain.append(current.GetIdx())
        onward = [n for n in current.GetNeighbors() if n.GetIdx() != chain[-2]]
        current = onward[0]
    if current.GetAtomicNum() != 6:
        return None
    count = len(chain)
    if z == 8:
        word = {1: "oxyl", 2: "peroxyl"}.get(count) or multiplying_prefix(count) + "oxidanyl"
    else:
        word = (multiplying_prefix(count) if count > 1 else "") + _CHALCOGEN_WORDS[z]
    return word, chain, current.GetIdx()


def _assembly(mol, centres):
    units = [_unit(mol, c) for c in centres]
    if len(centres) < 2 or any(u is None for u in units) or len({u[0] for u in units}) != 1:
        return None
    if any(is_functional_carbon(mol, attach) for _, _, attach in units):
        return None
    removed = {a for _, atoms, _ in units for a in atoms}
    skeleton = Chem.RWMol(mol)
    for index in sorted(removed, reverse=True):
        skeleton.RemoveAtom(index)
    if len(Chem.GetMolFrags(skeleton)) != 1:
        return None
    capped, _ = _capped(mol, [], removed)
    renumber = {old: new for new, old in enumerate(i for i in range(mol.GetNumAtoms()) if i not in removed)}
    editable = Chem.RWMol(capped)
    pairs = []
    for _, _, attach in units:
        dummy = editable.AddAtom(Chem.Atom(1))
        editable.AddBond(renumber[attach], dummy, Chem.BondType.SINGLE)
        pairs.append((renumber[attach], dummy, 1))
    out = editable.GetMol()
    out.UpdatePropertyCache(strict=False)
    if any(out.GetAtomWithIdx(a).GetIsAromatic() and not out.GetAtomWithIdx(a).IsInRing() for a, _, _ in pairs):
        return None
    part = _component(out, pairs)
    if part is None:
        return None
    word, count = units[0][0], len(units)
    central = _benzene(part.text)
    if part.has_locants or part.has_prefix:
        central = f"({central})"
    unit_text = multiplying_prefix(count) + word if word == "methyl" else multiplying_prefix(count, compound=True) + f"({word})"
    return central + unit_text


_ACID_ENDINGS = (
    ("carboxylic acid", "carbonyl"),
    ("oic acid", "oyl"),
    ("sulfinic acid", "sulfinyl"),
    ("sulfonic acid", "sulfonyl"),
)
_NITROGEN_ENDINGS = ("carboxamide", "amide", "imine")
_MULTIPLIER_WORDS = {2: "bis", 3: "tris", 4: "tetrakis"}


def _is_acyl(atom):
    if atom.GetNumRadicalElectrons() != 1 or atom.GetAtomicNum() not in (6, 16):
        return False
    oxygens = [
        n
        for n in atom.GetNeighbors()
        if n.GetAtomicNum() == 8 and n.GetDegree() == 1 and atom.GetOwningMol().GetBondBetweenAtoms(atom.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 2.0
    ]
    return len(oxygens) >= 1


def _healed(mol, centres, kind):
    editable = Chem.RWMol(mol)
    for atom in centres:
        target = editable.GetAtomWithIdx(atom.GetIdx())
        target.SetNumRadicalElectrons(0)
        target.SetNoImplicit(True)
        target.SetNumExplicitHs(atom.GetTotalNumHs())
        if kind == "acyl":
            editable.AddBond(atom.GetIdx(), editable.AddAtom(Chem.Atom(8)), Chem.BondType.SINGLE)
        else:
            target.SetNumExplicitHs(atom.GetTotalNumHs() + atom.GetNumRadicalElectrons())
    healed = editable.GetMol()
    Chem.SanitizeMol(healed)
    return healed


def _group_radical_name(mol, centres):
    """Several acyl radicals, or several imine or amide nitrogens (P-71.3.1, P-71.3.3): the name of the acid, imine or
    amide that results when each centre is filled, with its ending read as the radical suffix."""
    from .core import smiles_to_iupac

    count = len(centres)
    if all(_is_acyl(a) for a in centres):
        kind = "acyl"
    elif all(a.GetAtomicNum() == 7 and a.GetNumRadicalElectrons() == centres[0].GetNumRadicalElectrons() for a in centres):
        kind = "nitrogen"
    else:
        return None
    try:
        name = smiles_to_iupac(Chem.MolToSmiles(_healed(mol, centres, kind)))
    except (UnsupportedStructure, Chem.rdchem.MolSanitizeException):
        return None
    multiplier = multiplying_prefix(count)
    if kind == "acyl":
        for ending, radical in _ACID_ENDINGS:
            if name.endswith(multiplier + ending):
                return name[: -len(ending)] + radical
        return None
    suffix = "yl" if centres[0].GetNumRadicalElectrons() == 1 else "ylidene"
    for ending in _NITROGEN_ENDINGS:
        if name.endswith(multiplier + ending):
            stem = name[: -len(multiplier + ending)]
            return f"{stem}{_MULTIPLIER_WORDS[count]}({ending[:-1]}{suffix})"
    return None


def polyradical_name(mol):
    """Name of a radical with several centres, or None when no supported parent or assembly describes it."""
    centres = [a for a in mol.GetAtoms() if a.GetNumRadicalElectrons()]
    if len(centres) < 2 or not _plain(mol):
        return None
    name = _group_radical_name(mol, centres)
    if name is not None:
        return name
    if all(a.GetAtomicNum() == 6 for a in centres):
        name = _all_on_one_parent(mol, centres)
        if name is not None:
            return name
    return _assembly(mol, centres)
