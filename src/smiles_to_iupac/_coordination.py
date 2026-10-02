"""Mononuclear Group 3-12 organometallic naming (P-69.2.1-P-69.2.3,
`tmp/bluebook/P6a.txt` lines 8614-8710): ligands in alphanumerical order,
then the metal ('trichlorido(methyl)titanium'); charged complexes take
'(n+)'/'-ate(n-)' plus simple counter-ions. Chelating, hapto and bridging
ligands are out of scope.
"""

from rdkit import Chem

from ._common import (
    HALOGEN_PREFIXES,
    UnsupportedStructure,
    adjacency,
    non_single_bonds,
    plain_phenyl_substituent_atoms,
)
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
_HALIDO = {9: "fluorido", 17: "chlorido", 35: "bromido", 53: "iodido"}
_NEUTRAL_VALENCE = {7: 3, 8: 2, 15: 3, 16: 2, 33: 3}


def has_coordination_shape(mol) -> bool:
    return any(atom.GetAtomicNum() in _METAL_NAMES for atom in mol.GetAtoms())


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
        if any(mol.GetAtomWithIdx(i).GetAtomicNum() in _METAL_NAMES for i in atoms):
            raise UnsupportedStructure("bridging ligands and metal-metal bonds are not supported here")
        if sum(1 for n in metal.GetNeighbors() if n.GetIdx() in atoms) > 1:
            raise UnsupportedStructure("chelating, hapto and bridging ligands are not supported yet")
        atomic_num = donor.GetAtomicNum()
        if atomic_num in HALOGEN_PREFIXES and len(atoms) == 1:
            label = _HALIDO[atomic_num]
            simple_labels.add(label)
        elif atomic_num == 6 and _carbonyl(mol, donor, atoms):
            label = "carbonyl"
            simple_labels.add(label)
        elif atomic_num == 6:
            if any(mol.GetAtomWithIdx(i).GetFormalCharge() != 0 for i in atoms):
                raise UnsupportedStructure("charged ligands are not supported yet")
            phenyl = plain_phenyl_substituent_atoms(mol, graph, {donor.GetIdx()})
            if phenyl:
                if atoms != phenyl:
                    raise UnsupportedStructure("a substituted phenyl ligand is out of scope here")
                label = "phenyl"
            else:
                for i in atoms:
                    if mol.GetAtomWithIdx(i).GetAtomicNum() != 6 or mol.GetAtomWithIdx(i).IsInRing():
                        raise UnsupportedStructure("only acyclic alkyl and phenyl carbon ligands are supported here")
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
    metal = next(a for a in complex_mol.GetAtoms() if a.GetAtomicNum() in _METAL_NAMES)
    charge = metal.GetFormalCharge() + sum(c for c in cp_charges if c is not None)
    name = _name_complex(complex_mol, {_CP_LABEL: cp_count} if cp_count else None, charge)
    if not others:
        return name
    words = _counter_ion_words(others)
    if charge == 0:
        raise UnsupportedStructure("a neutral complex with counter-ions is not supported")
    return f"{name} {words}" if charge > 0 else f"{words} {name}"


def _name_dinuclear(mol, metals) -> str:
    first, second = metals
    if (
        first.GetAtomicNum() != second.GetAtomicNum()
        or mol.GetBondBetweenAtoms(first.GetIdx(), second.GetIdx()) is None
        or first.IsInRing()
        or any(m.GetFormalCharge() != 0 for m in metals)
        or any(a.GetIsotope() != 0 for a in mol.GetAtoms())
    ):
        raise UnsupportedStructure("only a neutral, directly bonded homodinuclear complex is supported here")
    graph = adjacency(mol)
    per_metal = [collect_ligands(mol, m, graph, skip={o.GetIdx()}) for m, o in ((first, second), (second, first))]
    per_metal.sort(key=lambda r: -sum(r[0].values()))
    counts: dict[str, int] = {}
    simple_labels, organic, neutral, donors = set(), set(), set(), {}
    for c, sl, org, neu, don in per_metal:
        for label, n in c.items():
            counts[label] = counts.get(label, 0) + n
        simple_labels |= sl
        organic |= org
        neutral |= neu
        donors.update(don)
    tags = {}
    for label in counts:
        parts = []
        for i, (c, *_rest) in enumerate(per_metal, start=1):
            n = c.get(label, 0)
            if n:
                parts.append(f"{i}\u03ba{n if n > 1 else ''}{donors[label]}")
        tags[label] = "-" + ",".join(parts) + "-"
    ligands = _format_ligands(counts, simple_labels, organic, neutral, tags)
    metal_name = _METAL_NAMES[first.GetAtomicNum()]
    symbol = first.GetSymbol()
    return f"{ligands}di{metal_name}({symbol}\u2014{symbol})"


def _name_complex(mol, extra=None, charge=None) -> str:
    metals = [a for a in mol.GetAtoms() if a.GetAtomicNum() in _METAL_NAMES]
    if len(metals) == 2:
        return _name_dinuclear(mol, metals)
    if len(metals) != 1:
        raise UnsupportedStructure("more than one transition-metal atom is not supported yet")
    (metal,) = metals
    if metal.IsInRing():
        raise UnsupportedStructure("a ring metal atom is not supported here (see P-69.4)")
    if any(a.GetIsotope() != 0 for a in mol.GetAtoms()):
        raise UnsupportedStructure("isotopically modified atoms are not supported yet")
    if charge is None:
        charge = metal.GetFormalCharge()

    graph = adjacency(mol)
    counts, simple_labels, organic, neutral, _ = collect_ligands(mol, metal, graph)
    for label, n in (extra or {}).items():
        counts[label] = counts.get(label, 0) + n
        organic.add(label)


    out = _format_ligands(counts, simple_labels, organic, neutral)
    if charge < 0:
        return out + _ATE_NAMES[metal.GetAtomicNum()] + _charge_text(charge)
    metal_name = _METAL_NAMES[metal.GetAtomicNum()]
    return out + metal_name + (_charge_text(charge) if charge else "")


def _format_ligands(counts, simple_labels, organic, neutral, tags=None) -> str:
    tags = tags or {}
    out = []
    for position, label in enumerate(sorted(counts, key=lambda s: s.lstrip("(").replace("\u03b75-", "").lower())):
        n = counts[label]
        simple = label in simple_labels or (label in organic and label not in neutral and _is_simple(label))
        if n > 1:
            text = multiplying_prefix(n, compound=not simple) + (label if simple else f"({label})")
        else:
            text = label if simple else f"({label})"
        if label in organic and position > 0 and simple and not text.startswith("("):
            text = f"({text})"
        out.append(text + tags.get(label, ""))
    return "".join(out)
