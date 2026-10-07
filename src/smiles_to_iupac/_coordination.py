"""Group 3-12 organometallic naming (P-69.2.1-P-69.2.3, the Blue Book): ligands in alphanumerical order, then the metal
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
    3: "lithium", 4: "beryllium", 11: "sodium", 12: "magnesium", 19: "potassium", 20: "calcium",
    37: "rubidium", 38: "strontium", 55: "caesium", 56: "barium",
}
_METAL_NAMES.update(_S_BLOCK_NAMES)
_F_BLOCK = {
    57: ("lanthanum", "lanthanate"), 58: ("cerium", "cerate"), 59: ("praseodymium", "praseodymate"),
    60: ("neodymium", "neodymate"), 61: ("promethium", "promethate"), 62: ("samarium", "samarate"),
    63: ("europium", "europate"), 64: ("gadolinium", "gadolinate"), 65: ("terbium", "terbate"),
    66: ("dysprosium", "dysprosate"), 67: ("holmium", "holmate"), 68: ("erbium", "erbate"),
    69: ("thulium", "thulate"), 70: ("ytterbium", "ytterbate"), 71: ("lutetium", "lutetate"),
    89: ("actinium", "actinate"), 90: ("thorium", "thorate"), 91: ("protactinium", "protactinate"),
    92: ("uranium", "uranate"), 93: ("neptunium", "neptunate"), 94: ("plutonium", "plutonate"),
    95: ("americium", "americate"), 96: ("curium", "curate"),
}
_METAL_NAMES.update({z: names[0] for z, names in _F_BLOCK.items()})
_HALIDO = {9: "fluorido", 17: "chlorido", 35: "bromido", 53: "iodido"}
_NEUTRAL_VALENCE = {7: 3, 8: 2, 15: 3, 16: 2, 33: 3}


def _ate(z):
    if z not in _ATE_NAMES:
        raise UnsupportedStructure("an anionic complex of this metal is not supported yet")
    return _ATE_NAMES[z]


def _is_metal(atom) -> bool:
    z = atom.GetAtomicNum()
    if z in _S_BLOCK_NAMES:
        donors = sum(1 for n in atom.GetNeighbors() if n.GetAtomicNum() in (7, 8, 15, 16))
        covalent = sum(1 for n in atom.GetNeighbors() if n.GetAtomicNum() in (1, 6))
        return donors >= 2 or covalent >= 1 or (atom.GetTotalNumHs() > 0 and atom.GetDegree() > 0)
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


def _covalent_root(mol, metal, donors_in):
    """The single donor bonded by an ordinary covalent bond (an aryl/alkyl
    carbon or an anionic N/P) when the other donors are lone-pair donors."""
    from ._metallacycle import _is_sigma

    roots = []
    for d in donors_in:
        i = d.GetIdx()
        bond = mol.GetBondBetweenAtoms(metal.GetIdx(), i)
        if bond.GetBondType() == Chem.BondType.DATIVE:
            continue
        z = d.GetAtomicNum()
        if z == 6 or (z in (7, 15) and _is_sigma(mol, metal.GetIdx(), i)):
            roots.append(i)
    return roots[0] if len(roots) == 1 else None


def _root_chelate_label(mol, metal, donors_in, atoms):
    """Chelating ligand with one covalent root: the root's substituent group
    or anion name carries the kappa terms ('2-(pyridin-2-yl-kN)phenyl-kC1',
    IR-10.2.3.3)."""
    from ._prefix_groups import PrefixNamer

    root = _covalent_root(mol, metal, donors_in)
    if root is None:
        return None
    sigma = {d.GetIdx() for d in donors_in} - {root}
    if any(mol.GetAtomWithIdx(i).GetFormalCharge() for i in atoms):
        raise UnsupportedStructure("charged atoms in a chelating ligand are not supported yet")
    for i in sigma:
        if mol.GetAtomWithIdx(i).GetAtomicNum() not in _NEUTRAL_VALENCE:
            raise UnsupportedStructure("this donor atom is not supported as a ligand yet")
    graph = {a: [n for n in adjacency(mol)[a] if n in atoms] for a in atoms}
    namer = PrefixNamer(mol, graph, sigma)
    atom = mol.GetAtomWithIdx(root)
    if atom.GetAtomicNum() == 6:
        name, _ = namer.name(root, metal.GetIdx())
        return f"{name}-\u03baC1"
    from ._substituents import format_mononuclear_prefixes

    entries = [namer.name(q, root) for q in graph[root]]
    word = "azanido" if atom.GetAtomicNum() == 7 else "phosphanido"
    return f"{format_mononuclear_prefixes(entries)}{word}-\u03ba{atom.GetSymbol()}"


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
    anionic = {7: "azanido", 8: "hydroxido"}
    if (
        atomic_num in anionic
        and not own
        and not donor.GetFormalCharge()
        and donor.GetTotalNumHs() == _NEUTRAL_VALENCE[atomic_num] - 1
    ):
        return anionic[atomic_num]
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
        a.SetNumRadicalElectrons(0)
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
_ATE_NAMES.update({z: names[1] for z, names in _F_BLOCK.items()})
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
        a.SetNumRadicalElectrons(0)
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


def _split_pi_sigma(mol, donors_in, metal=None):
    """(pi donors, sigma donors) when exactly one connected group of two or
    more donors (with at least two carbons) is accompanied only by single
    donors; None for an ordinary chelate or several pi groups."""
    ids = {d.GetIdx() for d in donors_in}
    if metal is not None:
        multiple = {
            i for i in ids
            if mol.GetAtomWithIdx(i).GetAtomicNum() == 6
            and mol.GetBondBetweenAtoms(metal.GetIdx(), i).GetBondTypeAsDouble() >= 2.0
        }
        if multiple and len(ids - multiple) >= 2:
            return ids - multiple, multiple
    if all(d.GetAtomicNum() == 6 for d in donors_in):
        return ids, set()
    groups, left = [], set(ids)
    while left:
        start = left.pop()
        group, stack = {start}, [start]
        while stack:
            for n in mol.GetAtomWithIdx(stack.pop()).GetNeighbors():
                if n.GetIdx() in left:
                    left.discard(n.GetIdx())
                    group.add(n.GetIdx())
                    stack.append(n.GetIdx())
        groups.append(group)
    big = [g for g in groups if len(g) > 1]
    if len(big) != 1 or sum(mol.GetAtomWithIdx(i).GetAtomicNum() == 6 for i in big[0]) < 2:
        return None
    return big[0], ids - big[0]


def _carbon_ligand_name(mol, graph, metal, donor, atoms):
    if donor.GetIsAromatic():
        seen = set()
        label = _brackets(substituted_aryl_name(mol, graph, donor.GetIdx(), metal.GetIdx(), seen)[0])
        if seen != atoms:
            raise UnsupportedStructure("this aryl ligand is not supported yet")
        return label
    for i in atoms:
        if mol.GetAtomWithIdx(i).GetAtomicNum() != 6 or mol.GetAtomWithIdx(i).IsInRing():
            raise UnsupportedStructure("only acyclic alkyl and aryl carbon ligands are supported here")
    if any(b[0] in atoms or b[1] in atoms for b in non_single_bonds(mol)):
        raise UnsupportedStructure("an unsaturated ligand is out of scope here")
    return name_branch(graph, donor.GetIdx(), metal.GetIdx(), {}, mol=mol)[0]


def _ylidene_name(mol, metal, donor, atoms, graph):
    if any(mol.GetAtomWithIdx(i).GetFormalCharge() != 0 or mol.GetAtomWithIdx(i).GetAtomicNum() != 6 for i in atoms):
        raise UnsupportedStructure("only a hydrocarbon alkylidene or alkylidyne ligand is supported here")
    for a, b, *_ in non_single_bonds(mol):
        if (a in atoms or b in atoms) and metal.GetIdx() not in (a, b):
            raise UnsupportedStructure("an unsaturated carbene substituent is out of scope here")
    return name_branch(graph, donor.GetIdx(), metal.GetIdx(), {}, mol=mol)[0]


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
            split = _split_pi_sigma(mol, donors_in, metal)
            if split:
                from ._hapto import hapto_label

                label = hapto_label(mol, metal, [mol.GetAtomWithIdx(i) for i in split[0]], atoms, split[1])
                donors[label] = "\u03b7"
            else:
                from ._anion_ligands import bidentate_anion

                label = (
                    bidentate_anion(mol, metal, donors_in, atoms)
                    or _root_chelate_label(mol, metal, donors_in, atoms)
                    or _chelate_label(mol, metal, donors_in, atoms)
                )
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
        elif atomic_num == 6 and mol.GetBondBetweenAtoms(metal.GetIdx(), donor.GetIdx()).GetBondTypeAsDouble() >= 2.0:
            label = _ylidene_name(mol, metal, donor, atoms, graph)
            organic.add(label)
        elif atomic_num == 6:
            if any(mol.GetAtomWithIdx(i).GetFormalCharge() != 0 for i in atoms):
                raise UnsupportedStructure("charged ligands are not supported yet")
            try:
                label = _carbon_ligand_name(mol, graph, metal, donor, atoms)
            except UnsupportedStructure:
                from ._prefix_groups import PrefixNamer

                ligand_graph = {a: [n for n in graph[a] if n in atoms] for a in atoms}
                label = PrefixNamer(mol, ligand_graph).name(donor.GetIdx(), metal.GetIdx())[0]
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
    if z == 6 and not others and hydrogens < 4:
        return ("carbido", "methanetriyl", "methanediyl", "methyl")[hydrogens]
    if z in (7, 15) and others and hydrogens + len(others) == 2:
        from ._anion_ligands import _prefixed

        metal_idx = next(iter(metal_ids))
        return _prefixed(mol, graph, atom.GetIdx(), metal_idx, "azanido" if z == 7 else "phosphanido", set())
    raise UnsupportedStructure("this bridging ligand is not supported yet")


_TABLE_VI_GROUP_ORDER = [18, 17, 16, 0, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1]
_TABLE_VI_GROUP_RANK = {g: i for i, g in enumerate(_TABLE_VI_GROUP_ORDER)}
_F_BLOCK_Z = set(range(57, 72)) | set(range(89, 104))
_GROUP_OF = {
    **{z: 1 for z in (3, 11, 19, 37, 55, 87)}, **{z: 2 for z in (4, 12, 20, 38, 56, 88)},
    **{z: 3 for z in (21, 39)}, **{z: 3 for z in _F_BLOCK_Z},
    **{z: 4 for z in (22, 40, 72)}, **{z: 5 for z in (23, 41, 73)}, **{z: 6 for z in (24, 42, 74)},
    **{z: 7 for z in (25, 43, 75)}, **{z: 8 for z in (26, 44, 76)}, **{z: 9 for z in (27, 45, 77)},
    **{z: 10 for z in (28, 46, 78)}, **{z: 11 for z in (29, 47, 79)}, **{z: 12 for z in (30, 48, 80)},
}


def _table_vi_rank(z):
    """Position in the Red Book Table VI traversal (later = more electropositive)."""
    return (_TABLE_VI_GROUP_RANK[_GROUP_OF.get(z, 0)], z)


def _hapto_bridges(mol, graph, metal_set):
    """Ligand fragments that reach two or more metals without any single atom
    bridging them (mu-eta ligands): [(atoms, {metal: donor atoms})]."""
    seen, found = set(), []
    for start in range(mol.GetNumAtoms()):
        if start in seen or start in metal_set:
            continue
        comp, stack = set(), [start]
        while stack:
            x = stack.pop()
            if x in comp or x in metal_set:
                continue
            comp.add(x)
            stack.extend(graph[x])
        seen |= comp
        reach = {m: {a for a in comp if m in graph[a]} for m in metal_set}
        reach = {m: d for m, d in reach.items() if d}
        multi = [a for a in comp if len([m for m in reach if a in reach[m]]) > 1]
        donors = set().union(*reach.values()) if reach else set()

        def pi_linked(a):
            return mol.GetAtomWithIdx(a).GetAtomicNum() == 6 and any(
                n in donors and mol.GetBondBetweenAtoms(a, n).GetBondTypeAsDouble() >= 2.0 for n in graph[a]
            )

        if len(reach) >= 2 and len(comp) > 1 and all(pi_linked(a) for a in multi):
            found.append((comp, reach))
    return found


_CLUSTER_WORDS = {(3, 3): "triangulo", (4, 6): "tetrahedro", (4, 4): "quadro"}


def _central_atom_order(mol, graph, ids, metal_set, raw, bridge_counts, bridge_sites):
    """Central-atom numbering of IR-9.2.5.6: Table VI (later element first),
    higher coordination number, then the ligands in alphabetical order (more
    ligating atoms of the first uneven ligand first), then neighbouring
    central atoms' locants."""
    labels = sorted(
        {label for i in ids for label in raw[i][0]} | {label for i in ids for label in bridge_counts[i]},
        key=lambda lab: (_sort_key(lab), lab),
    )

    def profile(i):
        terminal = raw[i][0]
        return tuple(-(terminal.get(lab, 0) + bridge_counts[i].get(lab, 0)) for lab in labels)

    def coordination_number(i):
        counts, _, _, _, donors = raw[i]
        total = 0
        for label, n in counts.items():
            if donors.get(label) == "\u03b7":
                total += n
            else:
                match = re.search(r"\u03ba(\d+)", label)
                total += n * (int(match.group(1)) if match else 1)
        return total + sum(bridge_sites.get(i, {}).values())

    def base(i):
        z = mol.GetAtomWithIdx(i).GetAtomicNum()
        return (tuple(-x for x in _table_vi_rank(z)), -coordination_number(i), profile(i))

    neighbours = {i: [n for n in graph[i] if n in metal_set] for i in ids}
    ordered = sorted(ids, key=lambda i: (base(i), i))
    for _ in range(len(ids)):
        classes = {}
        for i in ids:
            classes.setdefault(base(i), []).append(i)
        position = {}
        for k, key in enumerate(sorted(classes)):
            for i in classes[key]:
                position[i] = k
        refined = sorted(
            ids,
            key=lambda i: (base(i), tuple(sorted(position[n] for n in neighbours[i])), i),
        )
        if refined == ordered:
            break
        ordered = refined
    return ordered


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
    hapto_bridges = _hapto_bridges(mol, graph, metal_set)
    hapto_atoms = set().union(*(comp for comp, _ in hapto_bridges)) if hapto_bridges else set()
    bridge_atoms = [a for a in bridge_atoms if a.GetIdx() not in hapto_atoms]
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
    for _, reach in hapto_bridges:
        linked = list(reach)
        for other in linked[1:]:
            parent[find(other)] = find(linked[0])
    if len({find(i) for i in ids}) != 1:
        raise UnsupportedStructure("the metal atoms are not all connected to each other")

    bridge_info = {}
    bridge_counts: dict[int, dict[str, int]] = {i: {} for i in ids}
    bridge_sites: dict[int, dict[str, int]] = {i: {} for i in ids}
    for atom in bridge_atoms:
        label = _bridge_label(mol, graph, atom)
        linked = {n.GetIdx() for n in atom.GetNeighbors()} & metal_set
        bridge_info[(label, len(linked))] = bridge_info.get((label, len(linked)), 0) + 1
        for m in linked:
            bridge_counts[m][label] = bridge_counts[m].get(label, 0) + 1
            bridge_sites[m][label] = bridge_sites[m].get(label, 0) + 1
    skip = metal_set | {a.GetIdx() for a in bridge_atoms}
    hapto_results = []
    for comp, reach in hapto_bridges:
        from ._hapto import bridge_result

        skip |= comp
        res = bridge_result(mol, set().union(*reach.values()), comp)
        hapto_results.append((res, reach))
        for m, donors_m in reach.items():
            bridge_counts[m][res.stem] = bridge_counts[m].get(res.stem, 0) + len(donors_m)
            bridge_sites[m][res.stem] = bridge_sites[m].get(res.stem, 0) + 1

    raw = {i: collect_ligands(mol, mol.GetAtomWithIdx(i), graph, skip=skip) for i in ids}
    order = _central_atom_order(mol, graph, ids, metal_set, raw, bridge_counts, bridge_sites)
    number = {atom: k for k, atom in enumerate(order, start=1)}
    for res, reach in hapto_results:
        from ._hapto_ext import lkey, render

        pieces = sorted(
            (min((lkey(res.labels[d]) for d in donors)), number[m], render(res, donors)) for m, donors in reach.items()
        )
        if res.prefix:
            raise UnsupportedStructure("a substituted bridging hapto ligand is not supported here")
        label = "[" + ":".join(f"{num}({text})" for _, num, text in pieces) + "]" + res.stem
        bridge_info[(label, len(reach))] = bridge_info.get((label, len(reach)), 0) + 1

    from ._complex_stereo import complex_stereo_descriptor

    stereo_parts = [
        f"{number[m]}({body})" for m in order if (body := complex_stereo_descriptor(mol, m, metal_set))
    ]
    stereo = f"[{','.join(stereo_parts)}]-" if stereo_parts else ""

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
    distinct = list(dict.fromkeys(names))
    listing = "".join((multiplying_prefix(names.count(n)) if names.count(n) > 1 else "") + n for n in distinct)
    if charge < 0:
        if len(distinct) != 1:
            metal_part = f"({listing})ate"
        else:
            metal_part = multiplying_prefix(len(order)) + _ate(mol.GetAtomWithIdx(order[0]).GetAtomicNum())
    else:
        metal_part = listing
    metal_part += _charge_text(charge) if charge else ""
    cluster = _CLUSTER_WORDS.get((len(order), len(bonds)), "")
    if cluster:
        ligands += cluster + "-"
    pair_counts: dict[tuple[int, int], int] = {}
    for a, b in bonds:
        key = tuple(sorted((number[a], number[b])))
        pair_counts[key] = pair_counts.get(key, 0) + 1
    descriptor = ""
    if pair_counts:
        symbol = {number[i]: mol.GetAtomWithIdx(i).GetSymbol() for i in order}
        by_text: dict[str, int] = {}
        for (x, y), n in sorted(pair_counts.items()):
            text = f"{symbol[x]}—{symbol[y]}"
            by_text[text] = by_text.get(text, 0) + n
        parts = [f"{n} {text}" if n > 1 else text for text, n in by_text.items()]
        descriptor = "(" + ", ".join(parts) + ")"
    return f"{stereo}{ligands}{metal_part}{descriptor}"


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


    from ._complex_stereo import complex_stereo_prefix

    stereo = complex_stereo_prefix(mol, metal.GetIdx())
    out = _format_ligands(counts, simple_labels, organic, neutral)
    if charge < 0:
        return stereo + out + _ate(metal.GetAtomicNum()) + _charge_text(charge)
    metal_name = _METAL_NAMES[metal.GetAtomicNum()]
    return stereo + out + metal_name + (_charge_text(charge) if charge else "")


def _is_enclosed(label: str) -> bool:
    """True when the first enclosing mark closes at the very end."""
    pairs = {"[": "]", "(": ")", "{": "}"}
    if label[:1] not in pairs:
        return False
    depth = 0
    for i, ch in enumerate(label):
        if ch in "([{":
            depth += 1
        elif ch in ")]}":
            depth -= 1
            if depth == 0:
                return i == len(label) - 1
    return False


def _sort_key(label: str) -> str:
    """Alphanumerical key of a ligand name: enclosing marks, locants,
    stereodescriptors and eta terms are ignored (P-14.5)."""
    if re.match(r"\[\d+\(", label):
        label = label[label.index("]") + 1:]
    label = label.lstrip("[(")
    label = re.sub(r"^[\dEZRSrsez,' ]+\)-", "", label)
    label = re.sub(r"^(?:[\d\u2013,a-z]+-\u03b7\)-)?(?:\u03b7\d+-)?", "", label)
    return alpha_sort_key(label)


def _format_ligands(counts, simple_labels, organic, neutral, tags=None, bridges=None) -> str:
    tags = tags or {}
    hapto_bridge = any(label.startswith("[") for label, _ in (bridges or {}))
    entries = [(label, n, False, 0) for label, n in counts.items()] + [
        (label, n, True, mu) for (label, mu), n in (bridges or {}).items()
    ]
    entries.sort(key=lambda e: (_sort_key(e[0]), 0 if e[2] else 1, e[3]))
    out = []
    for position, (label, n, is_bridge, mu) in enumerate(entries):
        simple = label in simple_labels or (label in organic and label not in neutral and _is_simple(label))
        if simple or _is_enclosed(label):
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
            head = multiplying_prefix(n) + "-" if n > 1 else ""
            mu_text = "\u03bc" + (str(mu) if mu > 2 else "") + "-"
            out.append(head + ("{" + mu_text + label + "}-" if label.startswith("[") else mu_text + label + "-"))
            continue
        if label in organic and position > 0 and simple and not text.startswith("("):
            text = f"({text})"
        tag = tags.get(label, "")
        if hapto_bridge and "\u03ba" in tag and not text.startswith(("(", "[", "{")):
            out.append(f"({text}{tag.rstrip('-')})")
            continue
        out.append(text + tag)
    return "".join(out)
