"""Naming of the bare metallacyclic parent hydride (P-69.4, Chapter P-6a,
`tmp/bluebook/P6a.txt` lines 8873-8888): a monocyclic saturated all-carbon
ring in which exactly one ring carbon is replaced by a single Group 4-12
transition metal, named via the nondetachable skeletal-replacement ('a')
prefix, e.g. '1-titanacyclobutane'. The replacement locant is always cited
(never omitted, unlike the analogous single-heteroatom Hantzsch-Widman
case in `_hetero_monocyclic.py`) -- P-69.4's own worked examples always
show it: '1-platinacyclopenta-2,4-diene', '1-sila-2-ferracyclopentane',
'6-titanabicyclo[3.2.0]heptane', '1-iridabenzene'.

The 'a'-prefix stems below are drawn directly from that primary text
('titana', 'ferra', 'irida', 'platina') plus the same stems this project's
own `_metallocene.py` already uses for its "-ocene" names (ferrocene,
ruthenocene, osmocene, nickelocene, chromocene, cobaltocene, vanadocene --
stripping their shared '-ocene' suffix gives 'ferra'/'ruthena'/'osma'/
'nickela'/'chroma'/'cobalta'/'vanada', the identical element stem this
section's own 'a' prefix is built from).

Scope, deliberately narrow (the bare-parent-hydride case only, per the
issue): exactly one monocyclic, saturated, all-carbon-plus-one-metal ring
(3+ members), the metal bearing no substituent beyond its two ring bonds
and no charge/isotope, and no atom anywhere in the molecule outside the
ring (i.e. no ligand on the metal, no ring substituent). Any ligand/
substituent on the metal (P-69.2's coordination-nomenclature ligand
naming), any second ring (fused/bridged/spiro), any ring unsaturation, two
or more metal atoms, or a metal outside the ten-element table above is out
of scope for this pass.
"""

from ._numerals import alkane_name

_METAL_A_PREFIXES = {
    "Ti": "titana",
    "V": "vanada",
    "Cr": "chroma",
    "Fe": "ferra",
    "Co": "cobalta",
    "Ni": "nickela",
    "Ru": "ruthena",
    "Os": "osma",
    "Ir": "irida",
    "Pt": "platina",
}


def _find_ring_metal(mol):
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() != 1:
        return None
    ring_atoms = ring_info.AtomRings()[0]
    if len(ring_atoms) < 3 or mol.GetNumAtoms() != len(ring_atoms):
        return None

    metals = [idx for idx in ring_atoms if mol.GetAtomWithIdx(idx).GetSymbol() in _METAL_A_PREFIXES]
    if len(metals) != 1:
        return None
    (metal_idx,) = metals
    metal_atom = mol.GetAtomWithIdx(metal_idx)
    if metal_atom.GetDegree() != 2 or metal_atom.GetFormalCharge() != 0 or metal_atom.GetIsotope() != 0:
        return None

    for idx in ring_atoms:
        if idx == metal_idx:
            continue
        atom = mol.GetAtomWithIdx(idx)
        if atom.GetAtomicNum() != 6 or atom.GetIsAromatic() or atom.GetFormalCharge() != 0 or atom.GetIsotope() != 0:
            return None

    for bond in mol.GetBonds():
        if bond.GetBondTypeAsDouble() != 1.0:
            return None

    return metal_idx, len(ring_atoms)


def has_metallacycle_shape(mol) -> bool:
    return _find_ring_metal(mol) is not None


def name_metallacycle(mol) -> str:
    metal_idx, ring_size = _find_ring_metal(mol)
    prefix = _METAL_A_PREFIXES[mol.GetAtomWithIdx(metal_idx).GetSymbol()]
    return f"1-{prefix}cyclo{alkane_name(ring_size)}"
