"""Isotopically modified mononuclear hydrides of Groups 13 to 16 and their chains, and hydrides with substituents
(P-82.2.1, P-82.3, P-82.5.2, P-82.6)."""

import re

from rdkit import Chem

from ._common import InternalInvariantError, UnsupportedStructure, group_substituents, substituent_locant_set_and_citation
from ._substituents import format_substituent_prefixes

_MONONUCLEAR_HYDRIDES = {
    5: "borane", 7: "azane", 8: "oxidane", 13: "alumane", 14: "silane", 15: "phosphane", 16: "sulfane", 31: "gallane",
    32: "germane", 33: "arsane", 34: "selane", 49: "indigane", 50: "stannane", 51: "stibane", 52: "tellane",
    81: "thallane", 82: "plumbane", 83: "bismuthane",
}
_HALANES = {9: "fluorane", 17: "chlorane", 35: "bromane", 53: "iodane"}


def _heavy_chain(mol):
    """The ordered heavy atoms of an unbranched acyclic chain of one hydride-forming element, else None."""
    heavy = [a for a in mol.GetAtoms() if a.GetAtomicNum() != 1]
    if not heavy or len({a.GetAtomicNum() for a in heavy}) != 1:
        return None
    if heavy[0].GetAtomicNum() in _HALANES:
        single = len(heavy) == 1 and len(Chem.GetMolFrags(mol)) == 1 and heavy[0].GetFormalCharge() in (0, -1)
        return [heavy[0].GetIdx()] if single and not heavy[0].GetNumRadicalElectrons() else None
    if heavy[0].GetAtomicNum() not in _MONONUCLEAR_HYDRIDES:
        return None
    if len(Chem.GetMolFrags(mol)) != 1 or any(a.IsInRing() for a in mol.GetAtoms()):
        return None
    if len(heavy) > 1 and any(a.GetFormalCharge() or a.GetNumRadicalElectrons() for a in mol.GetAtoms()):
        return None
    if len(heavy) == 1 and heavy[0].GetFormalCharge() and heavy[0].GetNumRadicalElectrons():
        return None
    if any(b.GetBondTypeAsDouble() != 1.0 for b in mol.GetBonds()):
        return None
    degrees = {a.GetIdx(): sum(1 for n in a.GetNeighbors() if n.GetAtomicNum() != 1) for a in heavy}
    ends = [i for i, d in degrees.items() if d <= 1]
    if any(d > 2 for d in degrees.values()) or (len(heavy) > 1 and len(ends) != 2) or (len(heavy) == 1 and not ends):
        return None
    order, previous = [ends[0]], None
    while len(order) < len(heavy):
        onward = [n.GetIdx() for n in mol.GetAtomWithIdx(order[-1]).GetNeighbors() if n.GetAtomicNum() != 1 and n.GetIdx() != previous]
        previous = order[-1]
        order.append(onward[0])
    return order


_STANDARD_VALENCE = {5: 3, 7: 3, 8: 2, 13: 3, 14: 4, 15: 3, 16: 2, 31: 3, 32: 4, 33: 3, 34: 2, 49: 3, 50: 4, 51: 3, 52: 2, 81: 3, 82: 4, 83: 3, 9: 1, 17: 1, 35: 1, 53: 1}


def _is_radical(mol, atom):
    return not atom.GetFormalCharge() and atom.GetTotalValence() < _STANDARD_VALENCE[atom.GetAtomicNum()]


def _ion_or_radical_name(clean):
    """The name of the unmodified ion or radical of a mononuclear hydride, 'azanium', 'azanide', 'azanyl', 'hydroxyl'."""
    from .core import smiles_to_iupac

    return smiles_to_iupac(Chem.MolToSmiles(clean))


def has_isotope_hydride_shape(mol) -> bool:
    """A parent hydride of one element from Groups 13 to 16, a single atom or an unbranched chain, whose skeletal atoms
    or hydrogens are nuclides (P-82.2.1, P-21.2.2)."""
    return any(a.GetIsotope() for a in mol.GetAtoms()) and _heavy_chain(mol) is not None


def name_isotope_hydride(mol) -> str:
    """'(2H3)azane', '(13N)azane', '(18O)oxidane', '(1,1-2H2)hydrazine': the descriptor of the nuclides, symbols in
    alphabetical order (P-82.3.1), then the hydride name; the numbering gives the modified atoms the lowest locants
    (P-82.5.2) and a position that every atom shares needs none (P-82.6.1.3)."""
    from ._isotope_labels import descriptor, split_isotopes
    from ._numerals import multiplying_prefix

    order = _heavy_chain(mol)
    clean, labels, index = split_isotopes(mol)
    z = mol.GetAtomWithIdx(order[0]).GetAtomicNum()
    chain = [index[i] for i in order]
    best = None
    for candidate in (chain, chain[::-1]):
        positions = {atom: place for place, atom in enumerate(candidate, start=1)}
        key = sorted(positions[a] for a in labels)
        if best is None or key < best[0]:
            best = (key, positions)
    positions = best[1]
    uniform = (
        set(labels) == set(chain)
        and len({repr(sorted(e["H"].items())) + str(e["skeleton"]) for e in labels.values()}) == 1
        and all(sum(e["H"].values()) == clean.GetAtomWithIdx(a).GetTotalNumHs() for a, e in labels.items())
    )
    capacity = {a: clean.GetAtomWithIdx(a).GetTotalNumHs() for a in chain}
    text = descriptor(labels, positions, len(chain) == 1 or uniform, capacity=capacity)
    stem = _MONONUCLEAR_HYDRIDES.get(z) or _HALANES[z]
    centre = clean.GetAtomWithIdx(chain[0])
    if len(chain) == 1 and (centre.GetFormalCharge() or centre.GetNumRadicalElectrons() or _is_radical(clean, centre)):
        return text + _ion_or_radical_name(clean)
    parent = stem if len(chain) == 1 else "hydrazine" if z == 7 and len(chain) == 2 else multiplying_prefix(len(chain)) + stem
    return text + parent


def _substituted_hydride(mol):
    """(chain atoms, substituent roots by chain atom) of an unbranched acyclic chain of one hydride-forming element that
    carries only carbon groups, else None: 'methylhydrazine', 'dimethylphosphane'."""
    if len(Chem.GetMolFrags(mol)) != 1 or any(a.GetFormalCharge() or a.GetNumRadicalElectrons() for a in mol.GetAtoms()):
        return None
    elements = {a.GetAtomicNum() for a in mol.GetAtoms() if a.GetAtomicNum() not in (1, 6)}
    if len(elements) != 1 or (z := next(iter(elements))) not in _MONONUCLEAR_HYDRIDES:
        return None
    chain_atoms = {a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == z}
    degrees = {i: sum(1 for n in mol.GetAtomWithIdx(i).GetNeighbors() if n.GetIdx() in chain_atoms) for i in chain_atoms}
    ends = [i for i, d in degrees.items() if d <= 1]
    if any(d > 2 for d in degrees.values()) or len(ends) != (2 if len(chain_atoms) > 1 else 1):
        return None
    if z in (8, 16, 34, 52) and (len(chain_atoms) == 1 or z == 8):
        return None
    if len(chain_atoms) > 1 and z in (5, 13):
        return None
    order, previous = [ends[0]], None
    while len(order) < len(chain_atoms):
        onward = [n.GetIdx() for n in mol.GetAtomWithIdx(order[-1]).GetNeighbors() if n.GetIdx() in chain_atoms and n.GetIdx() != previous]
        if not onward:
            return None
        previous = order[-1]
        order.append(onward[0])
    roots = {}
    for i in order:
        for n in mol.GetAtomWithIdx(i).GetNeighbors():
            if n.GetIdx() in chain_atoms or n.GetAtomicNum() == 1:
                continue
            bond = mol.GetBondBetweenAtoms(i, n.GetIdx())
            ylidene = z == 7 and bond.GetBondTypeAsDouble() == 2.0
            if n.GetAtomicNum() != 6 or (bond.GetBondTypeAsDouble() != 1.0 and not ylidene) or bond.IsInRing():
                return None
            roots.setdefault(i, []).append(n.GetIdx())
    chain_bonds = (mol.GetBondBetweenAtoms(a, b) for a, b in zip(order, order[1:]))
    if not roots or any(b.GetBondTypeAsDouble() != 1.0 for b in chain_bonds):
        return None
    return order, roots


def _letters(name):
    """The letters of a name: locants, hyphens and enclosing marks differ between the cited and the bare form."""
    return re.sub(r"[^A-Za-z]", "", name)


def has_isotope_substituted_hydride_shape(mol) -> bool:
    """A nuclide on the atoms or hydrogens of the parent hydride of a hydrazine or phosphane-type compound that carries
    carbon groups, no other atom of the molecule being modified."""
    from ._isotope_labels import split_isotopes

    split = split_isotopes(mol)
    if split is None:
        return False
    clean, labels, _ = split
    found = _substituted_hydride(clean)
    if found is None or not all(a in found[0] for a in labels):
        return False
    try:
        return _labelled_hydride_name(mol) is not None
    except UnsupportedStructure:
        return False


def name_isotope_substituted_hydride(mol) -> str:
    """'1-methyl(2-2H1)hydrazine', 'methyl(2H1)phosphane', 'trimethyl(31P)phosphane': the unmodified name with the descriptor of
    the hydride's own atoms and hydrogens in front of the parent hydride, after the prefixes; once the descriptor needs
    locants every locant is cited (P-82.2.1, P-82.6.1.1, P-82.5.2)."""
    name = _labelled_hydride_name(mol)
    if name is None:
        raise InternalInvariantError("the substituted hydride shape was accepted for a structure named in another way")
    return name


def _labelled_hydride_name(mol):
    """The name of `name_isotope_substituted_hydride`, or None when the unmodified structure is named in another way."""
    from .core import smiles_to_iupac
    from ._common import adjacency, halogen_substituents
    from ._isotope_labels import descriptor, split_isotopes
    from ._numerals import multiplying_prefix
    from ._substituents import name_branch

    clean, labels, _ = split_isotopes(mol)
    chain, roots = _substituted_hydride(clean)
    graph = adjacency(clean)
    halogens = halogen_substituents(clean)
    z = clean.GetAtomWithIdx(chain[0]).GetAtomicNum()
    stem = _MONONUCLEAR_HYDRIDES[z]
    parent = stem if len(chain) == 1 else "hydrazine" if z == 7 and len(chain) == 2 else multiplying_prefix(len(chain)) + stem
    named = {a: [name_branch(graph, root, a, halogens, mol=clean) for root in group] for a, group in roots.items()}
    best = None
    for candidate in (chain, chain[::-1]):
        positions = {atom: place for place, atom in enumerate(candidate, start=1)}
        grouped = group_substituents({positions[a]: entries for a, entries in named.items()})
        locant_set, _, citation = substituent_locant_set_and_citation(grouped)
        key = (locant_set, citation, sorted(positions[a] for a in labels))
        if best is None or key < best[0]:
            best = (key, positions)
    positions = best[1]
    substituents = {positions[a]: entries for a, entries in named.items()}
    prefixes = format_substituent_prefixes(group_substituents(substituents), omit_locants=len(chain) == 1)
    capacity = {a: clean.GetAtomWithIdx(a).GetTotalNumHs() for a in chain}
    text = descriptor(labels, positions, len(chain) == 1, capacity=capacity)
    name = prefixes + text + parent
    if _letters(smiles_to_iupac(Chem.MolToSmiles(clean))) != _letters(prefixes + parent):
        return None
    return name
