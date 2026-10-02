"""Naming of mononuclear neutral Group 3-12 organometallics with
monodentate sigma-bound ligands (P-69.2.3, Chapter P-6a,
`tmp/bluebook/P6a.txt` lines 8680-8710): ligand names in alphanumerical
order, then the metal name, e.g. 'trichlorido(methyl)titanium'.

Anionic inorganic ligands take the '-ido' form (hydrido, chlorido, ...);
carbon ligands use their substituent-group names, which the Blue Book
lists as an accepted alternative to the '-ide' (methanido) form. A
multiplied simple name takes 'di'/'tri'; a name with a locant takes
'bis'/'tris'. Organic ligands after the first are enclosed in
parentheses, as in the worked examples.

Out of scope (raise `UnsupportedStructure`): charged atoms, more than one
metal, ring metals, ligands other than hydrogen/halogen/alkyl/phenyl
(carbonyl, aqua, ammine, phosphane, hapto and bridging ligands), and
oxidation-state/charge descriptors.
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


def has_coordination_shape(mol) -> bool:
    return any(atom.GetAtomicNum() in _METAL_NAMES for atom in mol.GetAtoms())


def _is_simple(name: str) -> bool:
    return not any(ch.isdigit() for ch in name) and "(" not in name


def name_coordination(mol) -> str:
    metals = [a for a in mol.GetAtoms() if a.GetAtomicNum() in _METAL_NAMES]
    if len(metals) != 1:
        raise UnsupportedStructure("more than one transition-metal atom is not supported yet")
    (metal,) = metals
    if metal.IsInRing():
        raise UnsupportedStructure("a ring metal atom is not supported here (see P-69.4)")
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")
    if any(a.GetFormalCharge() != 0 or a.GetIsotope() != 0 for a in mol.GetAtoms()):
        raise UnsupportedStructure("charged or isotopically modified atoms are not supported yet")

    graph = adjacency(mol)
    neighbors = list(metal.GetNeighbors())
    halogens = [n for n in neighbors if n.GetAtomicNum() in HALOGEN_PREFIXES]
    carbons = [n for n in neighbors if n.GetAtomicNum() == 6]
    if len(halogens) + len(carbons) != len(neighbors):
        raise UnsupportedStructure("ligands other than hydrogen, halogen, alkyl and phenyl are not supported yet")

    roots = {n.GetIdx() for n in carbons}
    phenyl_atoms = plain_phenyl_substituent_atoms(mol, graph, roots)
    halogen_idx = {n.GetIdx() for n in halogens}
    for atom in mol.GetAtoms():
        idx = atom.GetIdx()
        if atom.GetAtomicNum() not in (6, metal.GetAtomicNum()) and idx not in halogen_idx:
            raise UnsupportedStructure("heteroatoms other than the metal and its halido ligands are not supported yet")
        if atom.GetAtomicNum() == 6 and atom.GetIsAromatic() and idx not in phenyl_atoms:
            raise UnsupportedStructure("an aromatic ligand other than plain phenyl is out of scope here")
    ring_atoms = {a for ring in mol.GetRingInfo().AtomRings() for a in ring}
    if ring_atoms - phenyl_atoms:
        raise UnsupportedStructure("a ring ligand other than plain phenyl is out of scope here")
    if any(b[0] not in phenyl_atoms and b[1] not in phenyl_atoms for b in non_single_bonds(mol)):
        raise UnsupportedStructure("an unsaturated ligand is out of scope here")

    counts: dict[str, int] = {}
    organic: set[str] = set()
    if metal.GetTotalNumHs():
        counts["hydrido"] = metal.GetTotalNumHs()
    for n in halogens:
        label = _HALIDO[n.GetAtomicNum()]
        counts[label] = counts.get(label, 0) + 1
    for root in roots:
        label = "phenyl" if root in phenyl_atoms else name_branch(graph, root, metal.GetIdx(), {}, mol=mol)[0]
        counts[label] = counts.get(label, 0) + 1
        organic.add(label)

    out = []
    for position, label in enumerate(sorted(counts, key=lambda s: s.lstrip("(").lower())):
        n = counts[label]
        simple = _is_simple(label)
        if n > 1:
            text = multiplying_prefix(n, compound=not simple) + (label if simple else f"({label})")
        else:
            text = label if simple else f"({label})"
        if label in organic and position > 0 and simple:
            text = f"({text})"
        out.append(text)
    return "".join(out) + _METAL_NAMES[metal.GetAtomicNum()]
