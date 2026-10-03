"""Shared pieces of the metallacycle namers: metal ligands cited as ring
prefixes, ionic suffixes (-ium/-ide) and ring stereodescriptors (P-93).
"""

import re

from rdkit import Chem
from rdkit.Chem import rdCIPLabeler

from ._common import UnsupportedStructure
from ._numerals import multiplying_prefix

_HALO_FOR_LIGAND = {"fluorido": "fluoro", "chlorido": "chloro", "bromido": "bromo", "iodido": "iodo"}


def _is_plain(name: str) -> bool:
    return not any(ch.isdigit() for ch in name) and "(" not in name


def ring_ligand_entries(counts, simple_labels, neutral):
    """(name, compound, copies, donors) per ligand on the ring metal; a
    chelate is cited once without its kappa tag and its metal locant is
    repeated per donor ('9,9-[methylenebis(dimethylphosphane)]')."""
    entries = []
    for label, n in counts.items():
        compound = label in neutral or not (label in simple_labels or _is_plain(label))
        name = _HALO_FOR_LIGAND.get(label, label)
        donors = 1
        match = re.search(r"-κ(\d+)", label)
        if match:
            name, donors = label[: match.start()], int(match.group(1))
        entries.append((name, compound, n, donors))
    return entries


def add_ligands(grouped, entries, metal_locant):
    for name, compound, n, donors in entries:
        entry = grouped.setdefault(name, {"locants": [], "compound": compound})
        entry["locants"] += [metal_locant] * (n * donors)
        if donors > 1:
            entry["multiplier"] = multiplying_prefix(n, compound=True) if n > 1 else ""


def ionic_stem(stem, metal_locant, charge):
    """Cation/anion suffix on the metal ('platinacyclobutan-1-ium')."""
    if not charge:
        return stem
    if abs(charge) > 1:
        raise UnsupportedStructure("a metallacycle carrying more than one charge on the metal is not supported")
    word = "ium" if charge > 0 else "ide"
    return f"{stem[:-1] if stem.endswith('e') else stem}-{metal_locant}-{word}"


def check_charge(mol, metal_idx, ring_atoms, net_charge):
    metal_charge = mol.GetAtomWithIdx(metal_idx).GetFormalCharge()
    if any(mol.GetAtomWithIdx(i).GetFormalCharge() for i in ring_atoms if i != metal_idx):
        raise UnsupportedStructure("a charge on a ring atom other than the metal is not supported")
    if net_charge != metal_charge:
        raise UnsupportedStructure("a charge on a ligand of a metallacycle is not supported here")
    return metal_charge


_NONMETAL_RING_ATOMS = {5, 6, 7, 8, 14, 15, 16, 32, 33, 34}


def _branch_atoms(mol, ring):
    """Atoms of the substituent groups hanging off non-metal ring atoms."""
    found, stack = set(), [a for a in ring if mol.GetAtomWithIdx(a).GetAtomicNum() in _NONMETAL_RING_ATOMS]
    seen = set(stack)
    while stack:
        for n in mol.GetAtomWithIdx(stack.pop()).GetNeighbors():
            i = n.GetIdx()
            if i in ring or i in seen or n.GetAtomicNum() not in _NONMETAL_RING_ATOMS:
                continue
            seen.add(i)
            found.add(i)
            stack.append(i)
    return found


def ring_stereo(mol, ring_atoms, locant):
    """Sorted (locant, 'R'/'S') for the specified ring stereocentres, or [].

    A carbon stereocentre bonded to the ring is left to the substituent
    namer; any other one, or a double bond, is rejected so no descriptor is
    silently dropped."""
    ring = set(ring_atoms)
    branches = _branch_atoms(mol, ring)
    probe = Chem.Mol(mol)
    rdCIPLabeler.AssignCIPLabels(probe)
    labels = []
    for atom in probe.GetAtoms():
        if atom.GetChiralTag() == Chem.ChiralType.CHI_UNSPECIFIED:
            continue
        if not atom.HasProp("_CIPCode"):
            raise UnsupportedStructure("a stereocentre outside the ring skeleton is not supported here")
        if atom.GetIdx() not in ring:
            if atom.GetIdx() in branches:
                continue
            raise UnsupportedStructure("a stereocentre outside the ring skeleton is not supported here")
        labels.append((locant[atom.GetIdx()], atom.GetProp("_CIPCode")))
    if any(b.GetStereo() != Chem.BondStereo.STEREONONE for b in mol.GetBonds()):
        raise UnsupportedStructure("double-bond stereochemistry is not supported here")
    return sorted(labels, key=lambda x: (int(re.match(r'\d+', str(x[0])).group()), str(x[0])))


def stereo_prefix(labels):
    return f"({','.join(f'{l}{c}' for l, c in labels)})-" if labels else ""
