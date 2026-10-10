"""Acid-anion naming and citation shared by polyesters of one polyol (P-65.6.3.3.3): identical
anions multiply ('di' unsubstituted, 'bis' substituted), differing anions are listed
alphanumerically with locants when needed. Used by `_diester_acyloxy.py` and `_diester_ring_diyl.py`.
"""

import re

from rdkit import Chem
from rdkit.Chem import rdCIPLabeler

from ._common import (
    UnsupportedStructure,
    adjacency,
    alpha_sort_key,
    carbon_adjacency,
    component_subgraph,
    halogen_substituents,
    longest_chains,
)
from ._ester import _name_acyl_part
from ._functional_prefixes import functional_names
from ._numerals import multiplying_prefix

_MULTIPLICATIVE_START = re.compile(r"(?:di|do|tri|tetra|penta|hexa|hepta|octa|nona|dec)")


def find_ester_carbons(mol):
    """Every carbon shaped like an ester acyl carbon (a carbonyl oxygen plus
    a second, carbon-bonded ester oxygen), as (acyl_carbon, carbonyl_oxygen,
    ester_oxygen, alcohol_carbon) tuples."""
    matches = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6:
            continue
        oxygens = [n for n in atom.GetNeighbors() if n.GetAtomicNum() == 8]
        if len(oxygens) < 2:
            continue
        carbonyls = [
            o
            for o in oxygens
            if o.GetDegree() == 1 and mol.GetBondBetweenAtoms(atom.GetIdx(), o.GetIdx()).GetBondTypeAsDouble() == 2.0
        ]
        ester_oxygens = [
            o
            for o in oxygens
            if o.GetDegree() == 2
            and mol.GetBondBetweenAtoms(atom.GetIdx(), o.GetIdx()).GetBondTypeAsDouble() == 1.0
            and any(n.GetAtomicNum() == 6 for n in o.GetNeighbors() if n.GetIdx() != atom.GetIdx())
        ]
        if len(carbonyls) == 1 and len(ester_oxygens) == 1:
            alcohol_carbon = next(n for n in ester_oxygens[0].GetNeighbors() if n.GetIdx() != atom.GetIdx())
            if any(
                n.GetAtomicNum() == 8 and mol.GetBondBetweenAtoms(alcohol_carbon.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 2.0
                for n in alcohol_carbon.GetNeighbors()
            ):
                continue
            matches.append((atom, carbonyls[0], ester_oxygens[0], alcohol_carbon))
    return matches


def cip_labels(mol, side, chain_double_bonds=False):
    """[(atom or (begin, end) of a double bond, CIP code)] for the specified stereo elements within `side`; raises
    if `side` also has an unspecified element or a specified one that is neither tetrahedral nor a ring C=C (nor a
    chain double bond when `chain_double_bonds`)."""
    in_side = []
    for element in Chem.FindPotentialStereo(mol):
        if element.type == Chem.StereoType.Atom_Tetrahedral:
            inside = element.centeredOn in side
        else:
            bond = mol.GetBondWithIdx(element.centeredOn)
            inside = bond.GetBeginAtomIdx() in side and bond.GetEndAtomIdx() in side
        if inside:
            in_side.append(element)
    specified = [e for e in in_side if e.specified == Chem.StereoSpecified.Specified]
    if not specified:
        return []
    allowed = _is_double_bond if chain_double_bonds else _is_ring_double_bond
    if any(e.type != Chem.StereoType.Atom_Tetrahedral and not allowed(mol, e) for e in specified) or any(
        e.type == Chem.StereoType.Bond_Double and e.specified != Chem.StereoSpecified.Specified for e in in_side
    ):
        raise UnsupportedStructure(
            "stereochemistry beyond fully specified tetrahedral stereocenters and ring double bonds is not supported "
            "on the polyol side"
        )
    rdCIPLabeler.AssignCIPLabels(mol)
    pseudo = _pseudoasymmetric_in_group(mol, side)
    labels = []
    for element in specified:
        if element.type == Chem.StereoType.Atom_Tetrahedral:
            atom = mol.GetAtomWithIdx(element.centeredOn)
            if not atom.HasProp("_CIPCode"):
                raise UnsupportedStructure("could not determine a CIP label for a stereocenter")
            labels.append((element.centeredOn, pseudo.get(element.centeredOn, atom.GetProp("_CIPCode"))))
        else:
            bond = mol.GetBondWithIdx(element.centeredOn)
            if not bond.HasProp("_CIPCode"):
                raise UnsupportedStructure("could not determine a CIP label for a ring double bond")
            labels.append(((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()), bond.GetProp("_CIPCode")))
    return labels


def _is_double_bond(mol, element):
    return element.type == Chem.StereoType.Bond_Double and mol.GetBondWithIdx(element.centeredOn).GetBondType() == Chem.BondType.DOUBLE


def _is_ring_double_bond(mol, element):
    bond = mol.GetBondWithIdx(element.centeredOn)
    return element.type == Chem.StereoType.Bond_Double and bond.IsInRing() and bond.GetBondType() == Chem.BondType.DOUBLE


def _pseudoasymmetric_in_group(mol, side):
    """{atom: 'r'|'s'} for the centres of the substituent group `side` that are pseudoasymmetric in the group on its
    own (free valences as phantom atoms), where the rest of the molecule would otherwise break the symmetry (P-93.5.1)."""
    outside = {n.GetIdx() for a in side for n in mol.GetAtomWithIdx(a).GetNeighbors() if n.GetIdx() not in side}
    if not outside or len(side) == mol.GetNumAtoms():
        return {}
    keep = sorted(set(side) | outside)
    new_of = {old: new for new, old in enumerate(keep)}
    editable = Chem.RWMol(mol)
    for idx in sorted(set(range(mol.GetNumAtoms())) - set(keep), reverse=True):
        editable.RemoveAtom(idx)
    for a in outside:
        for b in outside:
            if a < b and editable.GetBondBetweenAtoms(new_of[a], new_of[b]) is not None:
                editable.RemoveBond(new_of[a], new_of[b])
    for a in outside:
        phantom = editable.GetAtomWithIdx(new_of[a])
        phantom.SetAtomicNum(0)
        phantom.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
        phantom.SetFormalCharge(0)
        phantom.SetNoImplicit(True)
        phantom.SetNumExplicitHs(0)
        phantom.SetIsAromatic(False)
        for bond in phantom.GetBonds():
            bond.SetIsAromatic(False)
    group = editable.GetMol()
    try:
        Chem.SanitizeMol(group)
    except Chem.rdchem.MolSanitizeException:
        return {}
    rdCIPLabeler.AssignCIPLabels(group)
    return {
        old: group.GetAtomWithIdx(new_of[old]).GetProp("_CIPCode")
        for old in side
        if group.GetAtomWithIdx(new_of[old]).HasProp("_CIPCode") and group.GetAtomWithIdx(new_of[old]).GetProp("_CIPCode") in "rs"
    }


def _acyl_component(graph, acyl_idx, ester_oxygen_idx):
    component = set()
    stack = [acyl_idx]
    while stack:
        node = stack.pop()
        if node in component:
            continue
        component.add(node)
        stack.extend(n for n in graph[node] if n != ester_oxygen_idx and n not in component)
    return component


def _specified_stereo_in(mol, atoms):
    for element in Chem.FindPotentialStereo(mol):
        if element.specified != Chem.StereoSpecified.Specified:
            continue
        if element.type == Chem.StereoType.Atom_Tetrahedral:
            if element.centeredOn in atoms:
                return True
        else:
            bond = mol.GetBondWithIdx(element.centeredOn)
            if bond.GetBeginAtomIdx() in atoms and bond.GetEndAtomIdx() in atoms:
                return True
    return False


def _fast_anion(mol, graph, carbon_graph, halogens, match, component):
    acyl, carbonyl, ester_o, _ = match
    acyl_idx = acyl.GetIdx()
    ring_neighbors = [n for n in graph[acyl_idx] if n in component and mol.GetRingInfo().NumAtomRings(n) > 0]
    if ring_neighbors:
        from ._diester_ring_diyl import ring_carboxylate_name

        name, atoms = ring_carboxylate_name(mol, graph, acyl_idx, ring_neighbors[0])
        return name, bool(component - atoms - {acyl_idx, carbonyl.GetIdx()})
    blocked = {acyl_idx, carbonyl.GetIdx(), ester_o.GetIdx()}
    seeds = [(acyl_idx, n) for n in graph[acyl_idx] if n not in blocked]
    named, shown, covered = functional_names(mol, graph, seeds, blocked, halogens, chain_seeds=True)
    chain_graph = {a: [n for n in ns if n not in covered] for a, ns in carbon_graph.items() if a not in covered}
    chains = longest_chains(component_subgraph(chain_graph, acyl_idx))
    substituted = len([a for a in component if a != carbonyl.GetIdx()]) > len(chains[0])
    name = _name_acyl_part(
        mol,
        acyl,
        carbonyl.GetIdx(),
        ester_o.GetIdx(),
        None,
        extra_names={a: shown[a] for a in named},
        extra_excluded_atoms=(frozenset(range(mol.GetNumAtoms())) - component) | covered,
    )
    return name, substituted


def _is_acid_group_carbon(mol, atom_idx):
    atom = mol.GetAtomWithIdx(atom_idx)
    if atom.GetAtomicNum() != 6:
        return False
    carbonyl = hydroxyl = False
    for n in atom.GetNeighbors():
        if n.GetAtomicNum() != 8:
            continue
        order = mol.GetBondBetweenAtoms(atom_idx, n.GetIdx()).GetBondTypeAsDouble()
        carbonyl |= order == 2.0
        hydroxyl |= order == 1.0 and n.GetTotalNumHs() == 1
    return carbonyl and hydroxyl


def _general_anion(mol, graph, match, component):
    """Anion name from the parent acid's own name (P-65.6.3.2.1: '-ic acid' becomes '-ate')."""
    from .core import smiles_to_iupac

    acyl, _, ester_o, _ = match
    if any(_is_acid_group_carbon(mol, a) for a in component if a != acyl.GetIdx()):
        raise UnsupportedStructure("a second free carboxylic acid group on an ester acyl group is not supported yet")
    rw = Chem.RWMol(mol)
    new_o = rw.AddAtom(Chem.Atom(8))
    rw.AddBond(acyl.GetIdx(), new_o, Chem.BondType.SINGLE)
    for idx in sorted(set(range(mol.GetNumAtoms())) - component, reverse=True):
        rw.RemoveAtom(idx)
    sub = rw.GetMol()
    Chem.SanitizeMol(sub)
    acid_name = smiles_to_iupac(Chem.MolToSmiles(sub))
    if not acid_name.endswith("ic acid"):
        raise UnsupportedStructure("the acyl group's parent acid has no 'ic acid' name to convert to an anion")
    name = acid_name[: -len("ic acid")] + "ate"
    carbonyl_idx = match[1].GetIdx()
    ring_atoms_in = {a for a in component if mol.GetRingInfo().NumAtomRings(a) > 0}
    if ring_atoms_in:
        substituted = bool(component - ring_atoms_in - {acyl.GetIdx(), carbonyl_idx})
    else:
        chains = longest_chains(component_subgraph(carbon_adjacency(mol), acyl.GetIdx()))
        substituted = len(component - {carbonyl_idx}) > len(chains[0])
    return name, substituted


_INORGANIC_ANIONS = {
    (16, 2): "hydrogen sulfate",
    (16, 1): "hydrogen sulfite",
    (15, 1): "dihydrogen phosphate",
    (15, 0): "dihydrogen phosphite",
}


def inorganic_anion_name(mol, centre):
    """Anion of a sulfur or phosphorus acid esterified once, its remaining hydroxyl groups cited as 'hydrogen' (P-67.1.3.2)."""
    doubled = sum(
        1 for n in centre.GetNeighbors() if mol.GetBondBetweenAtoms(centre.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 2.0
    )
    return _INORGANIC_ANIONS[(centre.GetAtomicNum(), doubled)]


def acid_anions(mol, matches):
    """[(anion_name, is_substituted), ...] aligned with `matches`."""
    graph = adjacency(mol)
    carbon_graph = carbon_adjacency(mol)
    halogens = halogen_substituents(mol)
    ester_atoms = set()
    for acyl, carbonyl, ester_o, _ in matches:
        ester_atoms |= {acyl.GetIdx(), carbonyl.GetIdx(), ester_o.GetIdx()}

    results = []
    for match in matches:
        acyl, carbonyl, ester_o, _ = match
        if acyl.GetAtomicNum() != 6:
            results.append((inorganic_anion_name(mol, acyl), True))
            continue
        component = _acyl_component(graph, acyl.GetIdx(), ester_o.GetIdx())
        if (component & ester_atoms) != {acyl.GetIdx(), carbonyl.GetIdx()}:
            raise UnsupportedStructure("an acyl group connected to another ester group is not supported yet")
        if not _specified_stereo_in(mol, component):
            try:
                results.append(_fast_anion(mol, graph, carbon_graph, halogens, match, component))
                continue
            except UnsupportedStructure:
                pass
        results.append(_general_anion(mol, graph, match, component))
    return results


def _enclose(name):
    if "[" in name:
        return f"{{{name}}}"
    if "(" in name:
        return f"[{name}]"
    return f"({name})"


def multiplied_anion(name, substituted, count):
    if count == 1:
        return name
    if substituted:
        return multiplying_prefix(count, compound=True) + _enclose(name)
    if any(ch.isdigit() for ch in name) or name.startswith("cyclo") or _MULTIPLICATIVE_START.match(name):
        return multiplying_prefix(count) + _enclose(name)
    return multiplying_prefix(count) + name


def cite_anions(anions, locants_by_match, cite_locants):
    """Anion words of the ester name; locant sets are cited only when `cite_locants`."""
    groups = {}
    for (name, substituted), locant in zip(anions, locants_by_match):
        entry = groups.setdefault(name, {"substituted": substituted, "locants": []})
        entry["locants"].append(locant)
    if len(groups) == 1:
        (name, entry), = groups.items()
        return multiplied_anion(name, entry["substituted"], len(entry["locants"]))
    parts = []
    for name in sorted(groups, key=alpha_sort_key):
        entry = groups[name]
        word = multiplied_anion(name, entry["substituted"], len(entry["locants"]))
        if cite_locants:
            if word[0].isdigit() or word[0] in "([{":
                word = _enclose(word)
            loc_str = ",".join(str(loc) for loc in sorted(entry["locants"]))
            word = f"{loc_str}-{word}"
        parts.append(word)
    return " ".join(parts)


def anion_locant_key(anions, locants_by_match):
    groups = {}
    for (name, _), locant in zip(anions, locants_by_match):
        groups.setdefault(name, []).append(locant)
    return tuple(tuple(sorted(groups[name])) for name in sorted(groups, key=alpha_sort_key))


__all__ = [
    "acid_anions",
    "anion_locant_key",
    "cip_labels",
    "cite_anions",
    "find_ester_carbons",
    "inorganic_anion_name",
    "multiplied_anion",
]
