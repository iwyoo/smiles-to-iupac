"""Group 3-12 organometallic naming (P-69.2.1-P-69.2.3, `tmp/bluebook/P6a.txt`
lines 8614-8710): ligands in alphanumerical order, then the metal
('trichlorido(methyl)titanium'); charged complexes, M-M bonds and
mu-bridging atoms included. Not SMILES-expressible, hence out of this
repo's scope: eta-n hapto ligands (other than the separate [CH-] Cp fragment).
"""

import re

from rdkit import Chem

from ._common import (
    HALOGEN_PREFIXES,
    UnsupportedStructure,
    adjacency,
    alpha_sort_key,
    non_single_bonds,
    plain_phenyl_substituent_atoms,
)
from ._metal_pair import _STEMS as _CLASS2_STEMS, _brackets, _metal_group_name, substituted_aryl_name
from ._numerals import multiplying_prefix
from ._substituents import name_branch

_METAL_NAMES = {
    21: "scandium", 22: "titanium", 23: "vanadium", 24: "chromium",
    25: "manganese", 26: "iron", 27: "cobalt", 28: "nickel", 29: "copper",
    30: "zinc", 39: "yttrium", 40: "zirconium", 41: "niobium",
    42: "molybdenum", 43: "technetium", 44: "ruthenium", 45: "rhodium",
    46: "palladium", 47: "silver", 48: "cadmium", 72: "hafnium",
    73: "tantalum", 74: "tungsten", 75: "rhenium", 76: "osmium",
    77: "iridium", 78: "platinum", 79: "gold", 80: "mercury",
}
_S_BLOCK_NAMES = {
    3: "lithium", 11: "sodium", 12: "magnesium", 19: "potassium", 20: "calcium",
    37: "rubidium", 38: "strontium", 55: "caesium", 56: "barium",
}
_METAL_NAMES.update(_S_BLOCK_NAMES)
_HALIDO = {9: "fluorido", 17: "chlorido", 35: "bromido", 53: "iodido"}
_NEUTRAL_VALENCE = {7: 3, 8: 2, 15: 3, 16: 2, 33: 3}


def _ate(z):
    if z not in _ATE_NAMES:
        raise UnsupportedStructure("an anionic complex of this metal is not supported yet")
    return _ATE_NAMES[z]


def _is_metal(atom) -> bool:
    z = atom.GetAtomicNum()
    if z in _S_BLOCK_NAMES:
        return sum(1 for n in atom.GetNeighbors() if n.GetAtomicNum() in (7, 8, 15, 16)) >= 2
    return z in _METAL_NAMES


def has_coordination_shape(mol) -> bool:
    return any(_is_metal(atom) for atom in mol.GetAtoms())


def _is_simple(name: str) -> bool:
    return not any(ch.isdigit() for ch in name) and "(" not in name


def _component(graph, start, metal_idx):
    seen, stack = set(), [start]
    while stack:
        a = stack.pop()
        if a in seen:
            continue
        seen.add(a)
        stack.extend(n for n in graph[a] if n != metal_idx)
    return seen


def _carbonyl(mol, donor, atoms):
    if len(atoms) != 2:
        return False
    other = next(i for i in atoms if i != donor.GetIdx())
    bond = mol.GetBondBetweenAtoms(donor.GetIdx(), other)
    return mol.GetAtomWithIdx(other).GetAtomicNum() == 8 and bond.GetBondTypeAsDouble() >= 2


def _monodentate_anion(mol, metal, donor, atoms):
    if donor.GetAtomicNum() in HALOGEN_PREFIXES or donor.GetAtomicNum() == 6:
        return None
    from ._anion_ligands import monodentate_anion

    return monodentate_anion(mol, metal, donor, atoms)


def _chelate_label(mol, metal, donors_in, atoms):
    symbols = []
    for donor in donors_in:
        z = donor.GetAtomicNum()
        own = [n for n in donor.GetNeighbors() if n.GetIdx() != metal.GetIdx()]
        valence = sum(mol.GetBondBetweenAtoms(donor.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() for n in own)
        if z not in _NEUTRAL_VALENCE or donor.GetFormalCharge() != 0 or valence + donor.GetNumExplicitHs() != _NEUTRAL_VALENCE[z]:
            raise UnsupportedStructure("only neutral heteroatom donors are supported in a chelating ligand")
        symbols.append(donor.GetSymbol())
    if any(mol.GetAtomWithIdx(i).GetFormalCharge() != 0 for i in atoms):
        raise UnsupportedStructure("charged ligands are not supported yet")
    name = _neutral_ligand_name(mol, metal, donors_in[0], atoms)
    seen: dict[str, int] = {}
    cited = []
    for symbol in sorted(symbols):
        cited.append(symbol + "'" * seen.get(symbol, 0))
        seen[symbol] = seen.get(symbol, 0) + 1
    return f"{name}-\u03ba{len(symbols)}{','.join(cited)}"


def _neutral_ligand_name(mol, metal, donor, atoms):
    atomic_num = donor.GetAtomicNum()
    if atomic_num not in _NEUTRAL_VALENCE:
        raise UnsupportedStructure("this donor atom is not supported as a ligand yet")
    own = [n for n in donor.GetNeighbors() if n.GetIdx() != metal.GetIdx()]
    valence = sum(mol.GetBondBetweenAtoms(donor.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() for n in own)
    if valence + donor.GetNumExplicitHs() != _NEUTRAL_VALENCE[atomic_num] or donor.GetFormalCharge() != 0:
        raise UnsupportedStructure("an anionic or multiply bonded heteroatom ligand is not supported yet")
    if atomic_num == 8 and not own and donor.GetNumExplicitHs() == 2:
        return "aqua"
    if atomic_num == 7 and not own and donor.GetNumExplicitHs() == 3:
        return "ammine"

    rw = Chem.RWMol(mol)
    keep = sorted(atoms)
    for idx in sorted(set(range(mol.GetNumAtoms())) - set(keep), reverse=True):
        rw.RemoveAtom(idx)
    ligand = rw.GetMol()
    for a in ligand.GetAtoms():
        a.SetNoImplicit(False)
    Chem.SanitizeMol(ligand)
    from .core import smiles_to_iupac

    return smiles_to_iupac(Chem.MolToSmiles(ligand))


_ATE_NAMES = {
    21: "scandate", 22: "titanate", 23: "vanadate", 24: "chromate",
    25: "manganate", 26: "ferrate", 27: "cobaltate", 28: "nickelate",
    29: "cuprate", 30: "zincate", 39: "yttrate", 40: "zirconate",
    41: "niobate", 42: "molybdate", 43: "technetate", 44: "ruthenate",
    45: "rhodate", 46: "palladate", 47: "argentate", 48: "cadmate",
    72: "hafnate", 73: "tantalate", 74: "tungstate", 75: "rhenate",
    76: "osmate", 77: "iridate", 78: "platinate", 79: "aurate",
    80: "mercurate",
}
_CATIONS = {3: "lithium", 11: "sodium", 19: "potassium"}
_ANIONS = {9: "fluoride", 17: "chloride", 35: "bromide", 53: "iodide"}


def _net_charge(mol) -> int:
    carbonyl_o = {
        n.GetIdx()
        for a in mol.GetAtoms()
        if a.GetAtomicNum() == 6 and any(_is_metal(m) for m in a.GetNeighbors())
        for n in a.GetNeighbors()
        if n.GetAtomicNum() == 8 and n.GetDegree() == 1 and mol.GetBondBetweenAtoms(a.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() >= 2.0
    }
    return sum(a.GetFormalCharge() for a in mol.GetAtoms() if a.GetIdx() not in carbonyl_o)


def _charge_text(charge: int) -> str:
    return f"({abs(charge)}{'+' if charge > 0 else '-'})"


def _counter_ion_words(frags):
    counts: dict[str, int] = {}
    for frag in frags:
        atom = frag.GetAtomWithIdx(0)
        if frag.GetNumAtoms() != 1:
            raise UnsupportedStructure("polyatomic counter-ions are not supported yet")
        table, charge = (_CATIONS, 1) if atom.GetFormalCharge() > 0 else (_ANIONS, -1)
        if atom.GetFormalCharge() != charge or atom.GetAtomicNum() not in table:
            raise UnsupportedStructure("this counter-ion is not supported yet")
        word = table[atom.GetAtomicNum()]
        counts[word] = counts.get(word, 0) + 1
    return " ".join((multiplying_prefix(n) if n > 1 else "") + w for w, n in sorted(counts.items()))


def _simple_anion_label(mol, donor, atoms):
    z, own = donor.GetAtomicNum(), [n for n in donor.GetNeighbors() if _is_metal(n)]
    if len(atoms) == 1 and z in (8, 16):
        bond_orders = [mol.GetBondBetweenAtoms(donor.GetIdx(), m.GetIdx()).GetBondTypeAsDouble() for m in own]
        hydrogens = donor.GetNumExplicitHs() + donor.GetNumImplicitHs()
        if hydrogens == 1 and donor.GetFormalCharge() == 0 and bond_orders == [1.0]:
            return "hydroxido" if z == 8 else "sulfanido"
        if hydrogens == 0 and donor.GetFormalCharge() in (0, -1):
            return "oxido" if z == 8 else "sulfido"
    if len(atoms) == 1 and z == 7 and donor.GetNumExplicitHs() == 2 and donor.GetDegree() == 1:
        return "azanido"
    if len(atoms) == 2:
        other = next(mol.GetAtomWithIdx(i) for i in atoms if i != donor.GetIdx())
        bond = mol.GetBondBetweenAtoms(donor.GetIdx(), other.GetIdx()).GetBondTypeAsDouble()
        if z == 6 and other.GetAtomicNum() == 7 and bond == 3.0 and donor.GetFormalCharge() == 0:
            return "cyanido"
        if z == 7 and other.GetAtomicNum() == 8 and bond == 2.0:
            return "nitrosyl"
    return None


def _ligand_fragment_name(mol, donor, atoms, cap_atomic_num, cap_bond=Chem.BondType.SINGLE):
    rw = Chem.RWMol(mol)
    for idx in sorted(set(range(mol.GetNumAtoms())) - set(atoms), reverse=True):
        rw.RemoveAtom(idx)
    new_donor = sorted(atoms).index(donor.GetIdx())
    cap = rw.AddAtom(Chem.Atom(cap_atomic_num))
    rw.AddBond(new_donor, cap, cap_bond)
    ligand = rw.GetMol()
    for a in ligand.GetAtoms():
        a.SetFormalCharge(0)
        a.SetNoImplicit(False)
    Chem.SanitizeMol(ligand)
    from .core import smiles_to_iupac

    return smiles_to_iupac(Chem.MolToSmiles(ligand))


def _alkoxido_name(mol, metal, donor, atoms):
    alcohol = _ligand_fragment_name(mol, donor, atoms, 1)
    if not alcohol.endswith("ol") or any(ch.isdigit() for ch in alcohol[-4:]):
        raise UnsupportedStructure("this alkoxide ligand is not supported yet")
    return alcohol[:-2] + "olato"


def _is_acyl(mol, donor):
    return any(
        n.GetAtomicNum() == 8 and mol.GetBondBetweenAtoms(donor.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 2.0
        for n in donor.GetNeighbors()
    )


def _acyl_name(mol, metal, donor, atoms):
    if any(mol.GetAtomWithIdx(i).GetFormalCharge() != 0 for i in atoms):
        raise UnsupportedStructure("charged ligands are not supported yet")
    name = _ligand_fragment_name(mol, donor, atoms, 17)
    if not name.endswith(" chloride"):
        raise UnsupportedStructure("this acyl ligand is not supported yet")
    return {"ethanoyl": "acetyl", "benzenecarbonyl": "benzoyl", "methanoyl": "formyl"}.get(name[: -len(" chloride")], name[: -len(" chloride")])


def collect_ligands(mol, metal, graph, skip=frozenset()):
    """Ligand labels on `metal` -> counts, plus the label classes used for
    enclosure; neighbors in `skip` (e.g. ring atoms) are not ligands."""
    counts: dict[str, int] = {}
    simple_labels: set[str] = set()
    organic: set[str] = set()
    neutral: set[str] = set()
    donors: dict[str, str] = {}
    if metal.GetTotalNumHs():
        counts["hydrido"] = metal.GetTotalNumHs()
        simple_labels.add("hydrido")

    seen: set[int] = set()
    for donor in metal.GetNeighbors():
        if donor.GetIdx() in skip:
            continue
        if donor.GetIdx() in seen:
            continue
        atoms = _component(graph, donor.GetIdx(), metal.GetIdx())
        seen |= atoms
        if any(_is_metal(mol.GetAtomWithIdx(i)) for i in atoms):
            raise UnsupportedStructure("bridging ligands and metal-metal bonds are not supported here")
        donors_in = [n for n in metal.GetNeighbors() if n.GetIdx() in atoms]
        if len(donors_in) > 1:
            if all(d.GetAtomicNum() == 6 for d in donors_in):
                from ._hapto import hapto_label

                label = hapto_label(mol, metal, donors_in, atoms)
                donors[label] = "\u03b7"
            else:
                from ._anion_ligands import bidentate_anion

                label = bidentate_anion(mol, metal, donors_in, atoms) or _chelate_label(mol, metal, donors_in, atoms)
                donors[label] = "\u03ba"
            organic.add(label)
            neutral.add(label)
            counts[label] = counts.get(label, 0) + 1
            continue
        atomic_num = donor.GetAtomicNum()
        if atomic_num in HALOGEN_PREFIXES and len(atoms) == 1:
            label = _HALIDO[atomic_num]
            simple_labels.add(label)
        elif _simple_anion_label(mol, donor, atoms) is not None:
            label = _simple_anion_label(mol, donor, atoms)
            simple_labels.add(label)
        elif (anion := _monodentate_anion(mol, metal, donor, atoms)) is not None:
            label = anion
            organic.add(label)
        elif atomic_num == 8 and len(atoms) > 1 and donor.GetDegree() == 2 and donor.GetNumExplicitHs() == 0 and donor.GetFormalCharge() in (0, -1):
            label = _alkoxido_name(mol, metal, donor, atoms)
            organic.add(label)
        elif atomic_num == 6 and _is_acyl(mol, donor):
            label = _acyl_name(mol, metal, donor, atoms)
            organic.add(label)
        elif atomic_num in _CLASS2_STEMS:
            if any(mol.GetAtomWithIdx(i).GetFormalCharge() != 0 for i in atoms):
                raise UnsupportedStructure("charged ligands are not supported yet")
            seen: set[int] = set()
            label = _brackets(_metal_group_name(mol, graph, donor.GetIdx(), metal.GetIdx(), seen))
            if seen != atoms:
                raise UnsupportedStructure("this metal-group ligand is not supported yet")
            organic.add(label)
        elif atomic_num == 6 and _carbonyl(mol, donor, atoms):
            label = "carbonyl"
            simple_labels.add(label)
        elif atomic_num == 6:
            if any(mol.GetAtomWithIdx(i).GetFormalCharge() != 0 for i in atoms):
                raise UnsupportedStructure("charged ligands are not supported yet")
            if donor.GetIsAromatic():
                seen = set()
                label = _brackets(substituted_aryl_name(mol, graph, donor.GetIdx(), metal.GetIdx(), seen)[0])
                if seen != atoms:
                    raise UnsupportedStructure("this aryl ligand is not supported yet")
            else:
                for i in atoms:
                    if mol.GetAtomWithIdx(i).GetAtomicNum() != 6 or mol.GetAtomWithIdx(i).IsInRing():
                        raise UnsupportedStructure("only acyclic alkyl and aryl carbon ligands are supported here")
                if any(b[0] in atoms or b[1] in atoms for b in non_single_bonds(mol)):
                    raise UnsupportedStructure("an unsaturated ligand is out of scope here")
                label = name_branch(graph, donor.GetIdx(), metal.GetIdx(), {}, mol=mol)[0]
            organic.add(label)
        else:
            if any(mol.GetAtomWithIdx(i).GetFormalCharge() != 0 for i in atoms):
                raise UnsupportedStructure("charged ligands are not supported yet")
            label = _neutral_ligand_name(mol, metal, donor, atoms)
            if label in ("aqua", "ammine"):
                simple_labels.add(label)
            else:
                organic.add(label)
                neutral.add(label)
        if "\u03ba" in label or (label.endswith(("azanido", "phosphanido")) and label not in ("azanido", "phosphanido")):
            neutral.add(label)
        counts[label] = counts.get(label, 0) + 1
        donors[label] = donor.GetSymbol()
    if "hydrido" in counts:
        donors["hydrido"] = "H"
    return counts, simple_labels, organic, neutral, donors


_CP_LABEL = "\u03b75-cyclopenta-2,4-dien-1-yl"
_CP_ANION = Chem.CanonSmiles("[CH-]1C=CC=C1")
_CP_RADICAL = Chem.CanonSmiles("[CH]1C=CC=C1")


def _cp_charge(frag):
    smiles = Chem.MolToSmiles(frag)
    if smiles == _CP_ANION:
        return -1
    if smiles == _CP_RADICAL:
        return 0
    return None


def name_coordination(mol) -> str:
    frags = Chem.GetMolFrags(mol, asMols=True)
    if len(frags) == 1:
        return _name_complex(mol)
    complexes = [f for f in frags if has_coordination_shape(f)]
    if len(complexes) != 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")
    (complex_mol,) = complexes
    others = [f for f in frags if f is not complex_mol]
    cp_charges = [_cp_charge(f) for f in others]
    cp_count = sum(c is not None for c in cp_charges)
    others = [f for f, c in zip(others, cp_charges) if c is None]
    metal = next(a for a in complex_mol.GetAtoms() if _is_metal(a))
    charge = _net_charge(complex_mol) + sum(c for c in cp_charges if c is not None)
    name = _name_complex(complex_mol, {_CP_LABEL: cp_count} if cp_count else None, charge)
    if not others:
        return name
    words = _counter_ion_words(others)
    if charge == 0:
        raise UnsupportedStructure("a neutral complex with counter-ions is not supported")
    return f"{name} {words}" if charge > 0 else f"{words} {name}"


def _branch_atoms(graph, start, blocked):
    seen, stack = set(), [start]
    while stack:
        x = stack.pop()
        if x in seen or x in blocked:
            continue
        seen.add(x)
        stack.extend(graph[x])
    return seen


def _bridge_label(mol, graph, atom):
    z = atom.GetAtomicNum()
    metal_ids = {n.GetIdx() for n in atom.GetNeighbors() if _is_metal(n)}
    others = [n for n in atom.GetNeighbors() if not _is_metal(n)]
    hydrogens = atom.GetNumExplicitHs() + atom.GetNumImplicitHs()
    if z == 1 and not others:
        return "hydrido"
    if z in _HALIDO and not others:
        return _HALIDO[z]
    if z in (8, 16) and not others:
        if hydrogens == 1 and z == 8:
            return "hydroxido"
        if hydrogens == 0:
            return "oxido" if z == 8 else "sulfido"
    if z == 6 and len(others) == 1 and others[0].GetAtomicNum() == 8 and others[0].GetDegree() == 1:
        return "carbonyl"
    if z in (8, 16) and len(others) == 1 and others[0].GetAtomicNum() == 6:
        atoms = {atom.GetIdx()} | _branch_atoms(graph, others[0].GetIdx(), metal_ids | {atom.GetIdx()})
        from ._anion_ligands import _neutral_name

        name = _neutral_name(mol, atoms)
        suffix = "ol" if z == 8 else "thiol"
        if name.endswith(suffix):
            return name + "ato"
    if z in (7, 15) and others and hydrogens + len(others) == 2:
        from ._anion_ligands import _prefixed

        metal_idx = next(iter(metal_ids))
        return _prefixed(mol, graph, atom.GetIdx(), metal_idx, "azanido" if z == 7 else "phosphanido", set())
    raise UnsupportedStructure("this bridging ligand is not supported yet")


_CLUSTER_WORDS = {(3, 3): "triangulo", (4, 6): "tetrahedro", (4, 4): "quadro"}


def _name_polynuclear(mol, metals) -> str:
    if any(a.GetIsotope() != 0 for a in mol.GetAtoms()):
        raise UnsupportedStructure("isotopically modified atoms are not supported yet")
    graph = adjacency(mol)
    ids = [m.GetIdx() for m in metals]
    metal_set = set(ids)
    bonds = [
        (a, b)
        for i, a in enumerate(ids)
        for b in ids[i + 1:]
        if mol.GetBondBetweenAtoms(a, b) is not None
    ]
    bridge_atoms = [
        a for a in mol.GetAtoms()
        if not _is_metal(a) and len({n.GetIdx() for n in a.GetNeighbors()} & metal_set) >= 2
    ]
    parent = {i: i for i in ids}

    def find(x):
        while parent[x] != x:
            x = parent[x]
        return x

    for a, b in bonds:
        parent[find(a)] = find(b)
    for atom in bridge_atoms:
        linked = [n.GetIdx() for n in atom.GetNeighbors() if n.GetIdx() in metal_set]
        for other in linked[1:]:
            parent[find(other)] = find(linked[0])
    if len({find(i) for i in ids}) != 1:
        raise UnsupportedStructure("the metal atoms are not all connected to each other")

    bridge_info = {}
    for atom in bridge_atoms:
        label = _bridge_label(mol, graph, atom)
        mu = len({n.GetIdx() for n in atom.GetNeighbors()} & metal_set)
        bridge_info[(label, mu)] = bridge_info.get((label, mu), 0) + 1
    skip = metal_set | {a.GetIdx() for a in bridge_atoms}

    raw = {i: collect_ligands(mol, mol.GetAtomWithIdx(i), graph, skip=skip) for i in ids}
    order = sorted(ids, key=lambda i: (_METAL_NAMES[mol.GetAtomWithIdx(i).GetAtomicNum()], -sum(raw[i][0].values()), i))
    number = {atom: k for k, atom in enumerate(order, start=1)}

    counts, simple_labels, organic, neutral, donors, per_label = {}, set(), set(), set(), {}, {}
    for i in order:
        c, sl, org, neu, don = raw[i]
        for label, n in c.items():
            shown = label.replace("-κ", f"-{number[i]}κ", 1) if "κ" in label else label
            counts[shown] = counts.get(shown, 0) + n
            per_label.setdefault(shown, {})[number[i]] = n
            if label in sl:
                simple_labels.add(shown)
            if label in org:
                organic.add(shown)
            if label in neu or "κ" in label:
                neutral.add(shown)
            donors[shown] = don[label]
    tags = {}
    for label, by_metal in per_label.items():
        if "κ" in label or donors[label] in ("κ", "η"):
            continue
        parts = [f"{k}κ{n if n > 1 else ''}{donors[label]}" for k, n in sorted(by_metal.items())]
        tags[label] = "-" + ",".join(parts) + "-"
    for label in counts:
        if "κ" in label and label not in tags:
            tags[label] = "-"
    simple_labels |= {label for label, _ in bridge_info}
    ligands = _format_ligands(counts, simple_labels, organic, neutral, tags, bridge_info)

    charge = _net_charge(mol)
    names = [_METAL_NAMES[mol.GetAtomWithIdx(i).GetAtomicNum()] for i in order]
    distinct = sorted(set(names))
    if charge < 0:
        if len(distinct) != 1:
            raise UnsupportedStructure("a heteronuclear anionic complex is not supported yet")
        metal_part = multiplying_prefix(len(order)) + _ate(mol.GetAtomWithIdx(order[0]).GetAtomicNum())
    else:
        metal_part = "".join(
            (multiplying_prefix(names.count(n)) if names.count(n) > 1 else "") + n for n in distinct
        )
    metal_part += _charge_text(charge) if charge else ""
    cluster = _CLUSTER_WORDS.get((len(order), len(bonds)), "")
    if cluster:
        ligands += cluster + "-"
    pair_counts: dict[str, int] = {}
    for a, b in bonds:
        sa, sb = sorted((mol.GetAtomWithIdx(a).GetSymbol(), mol.GetAtomWithIdx(b).GetSymbol()))
        pair_counts[f"{sa}—{sb}"] = pair_counts.get(f"{sa}—{sb}", 0) + 1
    descriptor = ""
    if pair_counts:
        parts = [f"{n} {text}" if n > 1 else text for text, n in sorted(pair_counts.items())]
        descriptor = "(" + ", ".join(parts) + ")"
    return f"{ligands}{metal_part}{descriptor}"


def _name_complex(mol, extra=None, charge=None) -> str:
    metals = [a for a in mol.GetAtoms() if _is_metal(a)]
    if len(metals) >= 2:
        return _name_polynuclear(mol, metals)
    (metal,) = metals
    if any(a.GetIsotope() != 0 for a in mol.GetAtoms()):
        raise UnsupportedStructure("isotopically modified atoms are not supported yet")
    if charge is None:
        charge = _net_charge(mol)

    graph = adjacency(mol)
    counts, simple_labels, organic, neutral, _ = collect_ligands(mol, metal, graph)
    for label, n in (extra or {}).items():
        counts[label] = counts.get(label, 0) + n
        organic.add(label)


    out = _format_ligands(counts, simple_labels, organic, neutral)
    if charge < 0:
        return out + _ate(metal.GetAtomicNum()) + _charge_text(charge)
    metal_name = _METAL_NAMES[metal.GetAtomicNum()]
    return out + metal_name + (_charge_text(charge) if charge else "")


def _sort_key(label: str) -> str:
    stripped = re.sub(r"^[\[(]*(?:[\d,]+-\u03b7\)-)?(?:\u03b7\d+-)?[\d,\-]*", "", label)
    return alpha_sort_key(stripped)


def _format_ligands(counts, simple_labels, organic, neutral, tags=None, bridges=None) -> str:
    tags = tags or {}
    entries = [(label, n, False, 0) for label, n in counts.items()] + [
        (label, n, True, mu) for (label, mu), n in (bridges or {}).items()
    ]
    entries.sort(key=lambda e: (_sort_key(e[0]), 0 if e[2] else 1, e[3]))
    out = []
    for position, (label, n, is_bridge, mu) in enumerate(entries):
        simple = label in simple_labels or (label in organic and label not in neutral and _is_simple(label))
        if simple or label.startswith("["):
            wrapped = label
        elif "{" in label:
            wrapped = f"({label})"
        elif "[" in label:
            wrapped = "{" + label + "}"
        elif "(" in label:
            wrapped = "[" + label + "]"
        else:
            wrapped = f"({label})"
        text = (multiplying_prefix(n, compound=not simple) if n > 1 else "") + wrapped
        if is_bridge:
            out.append((multiplying_prefix(n) + "-" if n > 1 else "") + "\u03bc" + (str(mu) if mu > 2 else "") + "-" + label + "-")
            continue
        if label in organic and position > 0 and simple and not text.startswith("("):
            text = f"({text})"
        out.append(text + tags.get(label, ""))
    return "".join(out)
