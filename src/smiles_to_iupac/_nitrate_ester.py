"""Naming of nitrate esters (R-O-NO2, a single R group), per the IUPAC
2013 Recommendations ("the Blue Book"):

- P-67.1.3.2 (Chapter P-6a, https://iupac.qmul.ac.uk/BlueBook/PDF/P6a.pdf):
  esters of mononuclear noncarbon oxoacids are named by citing the
  substituent (alkyl/aryl) group as a separate word, followed by the
  acid's anion name -- mirrors `_sulfate.py`/`_sulfite.py`'s own
  citation style, differing only in the anion word ("nitrate") and the
  central atom's shape.
- Word-assembly logic is shared verbatim with `_phosphate.py` via its
  exported `format_ester_words` -- exercised here only in its simplest
  (single word, no multiplying prefix) form, since nitrogen's valence
  leaves room for exactly one ester oxygen (unlike sulfate/sulfite's
  two).
- Structurally distinct from `_nitro.py`'s own nitrogen shape: a nitro
  group's nitrogen is bonded to one carbon plus two oxygens (a direct
  N-C bond); a nitrate ester's nitrogen is bonded to zero carbons and
  three oxygens instead (one ester O-R, plus the same two-representation
  terminal-oxygen pair `_nitro_nitrogen_atoms` already handles -- a
  charged `(=O)(-O^-)` pair or the neutral `(=O)(=O)` tautomer some
  SMILES writers use) -- not directly reusable, though the same
  dual-representation check applies.

A chalcogen analogue with S, Se or Te in the ester position is 'S-methyl thionitrate' (P-67.1.3.2).

Scope: a single nitrogen atom shaped like a nitrate ester -- one N-O-R
single bond, two terminal oxygens (either the charged or neutral
representation), R named via `name_branch`. The ester ranks above the
alcohol, amine, nitrile, aldehyde, ketone, amide, ether and ring groups of R,
which stay prefixes (P-41); R holding an acid or another ester raises `UnsupportedStructure`. Several nitrate
groups on one organyl group are named like the acetates of the same polyol (P-65.6.3.2.2).
"""

from rdkit import Chem

from ._common import UnsupportedStructure, adjacency, halogen_substituents
from ._substituents import centre_nuclides
from ._phosphate import check_ester_organyl, cited_ester_groups, format_ester_words, name_polyester_through_acetates

_NITROGEN = 7
_OXYGEN = 8
_ESTER_ATOMS = {8: ("", "O"), 16: ("thio", "S"), 34: ("seleno", "Se"), 52: ("telluro", "Te")}


def _nitrate_ester_nitrogen_atoms(mol):
    """Nitrogen atoms shaped like a fully-esterified nitrate: bonded to
    zero carbons, one single-bonded ester oxygen (degree 2 -- N plus one
    R carbon), and two more oxygens matching either the charged
    `[N+](=O)[O-]` pair or the neutral `N(=O)=O` tautomer, mirroring
    `_nitro_nitrogen_atoms`'s own dual-representation check."""
    matches = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != _NITROGEN or atom.GetDegree() != 3:
            continue
        neighbors = atom.GetNeighbors()
        if any(n.GetAtomicNum() == 6 for n in neighbors):
            continue
        oxygens = [n for n in neighbors if n.GetAtomicNum() in _ESTER_ATOMS]
        if len(oxygens) != 3:
            continue
        ester_os = [
            o
            for o in oxygens
            if o.GetDegree() == 2 and mol.GetBondBetweenAtoms(atom.GetIdx(), o.GetIdx()).GetBondTypeAsDouble() == 1.0
        ]
        if len(ester_os) != 1:
            continue
        terminal_os = [o for o in oxygens if o.GetIdx() != ester_os[0].GetIdx()]
        if any(o.GetDegree() != 1 or o.GetAtomicNum() != _OXYGEN for o in terminal_os):
            continue
        bond_orders = sorted(
            mol.GetBondBetweenAtoms(atom.GetIdx(), o.GetIdx()).GetBondTypeAsDouble() for o in terminal_os
        )
        if bond_orders == [2.0, 2.0]:
            if atom.GetFormalCharge() != 0:
                continue
        elif bond_orders == [1.0, 2.0]:
            single_o = terminal_os[
                0
                if mol.GetBondBetweenAtoms(atom.GetIdx(), terminal_os[0].GetIdx()).GetBondTypeAsDouble() == 1.0
                else 1
            ]
            if atom.GetFormalCharge() != 1 or single_o.GetFormalCharge() != -1:
                continue
        else:
            continue
        matches.append(atom)
    return matches


def has_nitrate_ester_shape(mol) -> bool:
    return bool(_nitrate_ester_nitrogen_atoms(mol))


def name_nitrate_ester(mol) -> str:
    nitrogen_atoms = _nitrate_ester_nitrogen_atoms(mol)
    if len(nitrogen_atoms) > 1:
        group_atoms = {
            n.GetIdx() for nitrogen in nitrogen_atoms for n in [nitrogen, *nitrogen.GetNeighbors()] if n.GetAtomicNum() in (7, 8)
        }
        if any(n.GetAtomicNum() != _OXYGEN for nitrogen in nitrogen_atoms for n in nitrogen.GetNeighbors()):
            raise UnsupportedStructure("a chalcogen analogue of several nitrate groups is not supported yet")
        check_ester_organyl(mol, group_atoms)
        return name_polyester_through_acetates(mol, [n.GetIdx() for n in nitrogen_atoms], "nitrate", "nitrooxy")
    (nitrogen,) = nitrogen_atoms

    group_atom_idxs = {nitrogen.GetIdx()} | {n.GetIdx() for n in nitrogen.GetNeighbors() if n.GetAtomicNum() in _ESTER_ATOMS}
    check_ester_organyl(mol, group_atom_idxs)
    group_oxygens = {n.GetIdx() for n in nitrogen.GetNeighbors() if n.GetAtomicNum() in _ESTER_ATOMS}

    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")

    (ester_oxygen_idx,) = [
        idx
        for idx in group_oxygens
        if mol.GetAtomWithIdx(idx).GetDegree() == 2
        and mol.GetBondBetweenAtoms(nitrogen.GetIdx(), idx).GetBondTypeAsDouble() == 1.0
    ]

    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    aromatic_atoms = {atom.GetIdx() for atom in mol.GetAtoms() if atom.GetIsAromatic()}

    (root,) = [n for n in graph[ester_oxygen_idx] if n != nitrogen.GetIdx()]
    ((name, _),) = cited_ester_groups(mol, graph, [(root, ester_oxygen_idx)], halogens, aromatic_atoms)
    word, symbol = _ESTER_ATOMS[mol.GetAtomWithIdx(ester_oxygen_idx).GetAtomicNum()]
    nuclide = centre_nuclides(nitrogen.GetIdx())
    if word:
        return f"{symbol}-{format_ester_words([name])} {word}{nuclide}nitrate"
    return f"{format_ester_words([name])} {nuclide}nitrate"
