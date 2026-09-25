"""Naming of open-chain (Fischer-projection) monosaccharides -- currently
just plain aldoses of 3-6 carbons, named via their retained stem name --
per the IUPAC 2013 Recommendations ("the Blue Book"), Chapter P-10
(https://iupac.qmul.ac.uk/BlueBook/P10.html):

- P-102.3.2/P-102.3.3: the 'D'/'L' stereodescriptor is assigned from the
  CIP configuration of the highest-numbered chirality center (the carbon
  nearest the terminal -CH2OH, i.e. C-(n-1) for an n-carbon aldose
  counted from the aldehyde) -- confirmed empirically (see P-102 epic
  #1039 M1 step 1) against 19 real PubChem aldose structures that this
  always maps CIP 'R' to 'D' and 'S' to 'L' for this specific
  substituent pattern (every chain carbon bears exactly one -OH, a
  terminal aldehyde at one end and -CH2OH at the other).
- P-102.5.1/P-102.5.2.2 (Table 10.2): each configuration of the
  remaining chirality centers (C-2 through C-(n-1) inclusive, counted
  from the aldehyde -- the configurational atom is always the pattern's
  own last entry) maps to one retained stem name (glyceraldehyde,
  erythrose/threose, ribose/arabinose/xylose/lyxose,
  allose/altrose/glucose/mannose/gulose/idose/galactose/talose). The
  L-series pattern is always the exact CIP mirror (every R flipped to S
  and vice versa) of its D-series counterpart -- confirmed against 4
  more real L-series structures (#1039 M1 step 1).

Scope, deliberately narrow (first pass at carbohydrate nomenclature, WS1/
M1 step 1 of #1039): a plain, unbranched, acyclic aldose backbone only --
a terminal aldehyde (C-1), 1 to 4 -CH(OH)- chirality-bearing carbons, and
a terminal -CH2OH, 3 to 6 carbons total, with every stereocenter's
configuration specified and no substituent anywhere beyond each chain
carbon's own single -OH (or, at C-1/the terminal carbon, none beyond
what the aldehyde/-CH2OH shape itself requires). Ketoses (P-102.5.2,
#1041), aldoses beyond 6 carbons (P-102.5.1.1.2, #1042), any cyclic/ring
form (P-102.3.4, this project's M2/#85), deoxy/amino sugars, glycosides,
and any other substituent are all out of scope here -- a molecule
matching any of those still falls through to the existing generic
acyclic-aldehyde/polyol naming unchanged, same as it does today.
"""

from rdkit.Chem import rdCIPLabeler

from ._common import adjacency, ordered_chain

# D-series CIP pattern (C-2 .. C-(n-1), in that order, counted from the
# aldehyde) -> retained stem name. A pattern's length is (carbon count -
# 2): 1 for the triose (where the sole entry is also the configurational
# atom itself), 2 for tetroses, 3 for pentoses, 4 for hexoses. Confirmed
# via `rdCIPLabeler` against real PubChem aldehydo-structures (see module
# docstring; #1039 M1 step 1).
_D_ALDOSE_PATTERNS = {
    ("R",): "glyceraldehyde",
    ("R", "R"): "erythrose",
    ("S", "R"): "threose",
    ("R", "R", "R"): "ribose",
    ("S", "R", "R"): "arabinose",
    ("R", "S", "R"): "xylose",
    ("S", "S", "R"): "lyxose",
    ("R", "S", "R", "R"): "glucose",
    ("S", "S", "R", "R"): "mannose",
    ("R", "S", "S", "R"): "galactose",
    ("R", "R", "R", "R"): "allose",
    ("S", "R", "R", "R"): "altrose",
    ("R", "R", "S", "R"): "gulose",
    ("S", "R", "S", "R"): "idose",
    ("S", "S", "S", "R"): "talose",
}

_FLIP_CIP = {"R": "S", "S": "R"}


def _open_chain_aldose_backbone(mol):
    """The aldose backbone's carbon chain, aldehyde carbon first and
    terminal -CH2OH carbon last, or None if `mol` doesn't have this exact
    shape (see module docstring): acyclic, exactly one terminal aldehyde,
    an unbranched carbon chain of 3-6 atoms from it, every chain carbon
    bearing exactly one single-bonded oxygen substituent (an -OH, or the
    aldehyde's own =O at C-1), and no heavy atom anywhere in the molecule
    outside that chain and those oxygens."""
    if mol.GetRingInfo().NumRings() > 0:
        return None
    aldehyde_carbons = [
        atom.GetIdx()
        for atom in mol.GetAtoms()
        if atom.GetSymbol() == "C"
        and atom.GetDegree() == 2
        and any(
            bond.GetBondTypeAsDouble() == 2.0 and bond.GetOtherAtom(atom).GetSymbol() == "O"
            for bond in atom.GetBonds()
        )
    ]
    if len(aldehyde_carbons) != 1:
        return None
    (aldehyde_carbon,) = aldehyde_carbons

    graph = adjacency(mol)
    oxygens = {atom.GetIdx() for atom in mol.GetAtoms() if atom.GetSymbol() == "O"}
    chain = ordered_chain(graph, aldehyde_carbon, -1, oxygens)
    if chain is None or not 3 <= len(chain) <= 6:
        return None
    if mol.GetNumAtoms() != 2 * len(chain):
        return None

    for atom_idx in chain[1:-1]:
        chain_oxygens = [n for n in graph[atom_idx] if n in oxygens]
        if len(chain_oxygens) != 1:
            return None
        (oh_oxygen,) = chain_oxygens
        if mol.GetBondBetweenAtoms(atom_idx, oh_oxygen).GetBondTypeAsDouble() != 1.0:
            return None

    terminal = chain[-1]
    terminal_oxygens = [n for n in graph[terminal] if n in oxygens]
    if len(terminal_oxygens) != 1:
        return None
    (terminal_oxygen,) = terminal_oxygens
    if mol.GetBondBetweenAtoms(terminal, terminal_oxygen).GetBondTypeAsDouble() != 1.0:
        return None
    if mol.GetAtomWithIdx(terminal).GetTotalNumHs() != 2:
        return None
    return chain


def _aldose_cip_pattern(mol, chain):
    """`(cip, ...)` for `chain[1:-1]` (C-2 through C-(n-1)), or None if
    any of those stereocenters is left unspecified."""
    rdCIPLabeler.AssignCIPLabels(mol)
    pattern = []
    for atom_idx in chain[1:-1]:
        atom = mol.GetAtomWithIdx(atom_idx)
        if not atom.HasProp("_CIPCode"):
            return None
        pattern.append(atom.GetProp("_CIPCode"))
    return tuple(pattern)


def has_open_chain_aldose_shape(mol) -> bool:
    chain = _open_chain_aldose_backbone(mol)
    if chain is None:
        return False
    pattern = _aldose_cip_pattern(mol, chain)
    if pattern is None:
        return False
    if pattern in _D_ALDOSE_PATTERNS:
        return True
    flipped = tuple(_FLIP_CIP[c] for c in pattern)
    return flipped in _D_ALDOSE_PATTERNS


def name_open_chain_aldose(mol) -> str:
    chain = _open_chain_aldose_backbone(mol)
    pattern = _aldose_cip_pattern(mol, chain)
    if pattern in _D_ALDOSE_PATTERNS:
        return f"D-{_D_ALDOSE_PATTERNS[pattern]}"
    flipped = tuple(_FLIP_CIP[c] for c in pattern)
    return f"L-{_D_ALDOSE_PATTERNS[flipped]}"
