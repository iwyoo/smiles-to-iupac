"""Naming of open-chain (Fischer-projection) monosaccharides -- currently
plain aldoses of 3-6 carbons and plain 2-ketoses of 4-6 carbons, each
named via their retained stem name -- per the IUPAC 2013 Recommendations
("the Blue Book"), Chapter P-10 (https://iupac.qmul.ac.uk/BlueBook/P10.html):

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
- P-102.5.2.1/P-102.5.2.2 (Table 10.3): a *2*-ketose (the carbonyl at
  C-2, one carbon in from the -CH2OH end that becomes C-1) follows the
  identical D/L rule at its own highest-numbered chirality center, and
  its remaining chirality centers (C-3 through C-(n-1)) map to one of a
  much shorter list of retained names (erythrulose, ribulose/xylulose,
  fructose/sorbose/psicose/tagatose) -- confirmed against 13 more real
  PubChem structures (#1039 M1 step 2). A ketose whose carbonyl sits at
  C-3 or higher has no retained name at all (P-102.5.2.3's own examples
  compose a configurational-prefix name instead), out of scope here.

Scope, deliberately narrow (first pass at carbohydrate nomenclature, WS1
of #1039): a plain, unbranched, acyclic aldose (M1 step 1) or 2-ketose
(M1 step 2) backbone only -- a terminal aldehyde or a C-2 carbonyl
flanked by a terminal -CH2OH, 1 to 4 more -CH(OH)- chirality-bearing
carbons, and a second terminal -CH2OH, 3 to 6 (aldose) or 4 to 6
(2-ketose) carbons total, with every stereocenter's configuration
specified and no substituent anywhere beyond each chain carbon's own
single -OH (or, at the terminal carbons, none beyond what the aldehyde/
-CH2OH/carbonyl shape itself requires). A ketose with its carbonyl at
C-3 or higher (#1041's own later scope), any aldose/ketose beyond 6
carbons (P-102.5.1.1.2, #1042), any cyclic/ring form (P-102.3.4, this
project's M2/#85), deoxy/amino sugars, glycosides, and any other
substituent are all out of scope here -- a molecule matching any of
those still falls through to the existing generic acyclic-aldehyde/
ketone/polyol naming unchanged, same as it does today.
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


def _cip_pattern(mol, atom_indices):
    """`(cip, ...)` for `atom_indices`, in order, or None if any of them
    is left unspecified. Shared by the aldose and ketose backbones below
    -- each just picks a different slice of its own chain."""
    rdCIPLabeler.AssignCIPLabels(mol)
    pattern = []
    for atom_idx in atom_indices:
        atom = mol.GetAtomWithIdx(atom_idx)
        if not atom.HasProp("_CIPCode"):
            return None
        pattern.append(atom.GetProp("_CIPCode"))
    return tuple(pattern)


def has_open_chain_aldose_shape(mol) -> bool:
    chain = _open_chain_aldose_backbone(mol)
    if chain is None:
        return False
    pattern = _cip_pattern(mol, chain[1:-1])
    if pattern is None:
        return False
    if pattern in _D_ALDOSE_PATTERNS:
        return True
    flipped = tuple(_FLIP_CIP[c] for c in pattern)
    return flipped in _D_ALDOSE_PATTERNS


def name_open_chain_aldose(mol) -> str:
    chain = _open_chain_aldose_backbone(mol)
    pattern = _cip_pattern(mol, chain[1:-1])
    if pattern in _D_ALDOSE_PATTERNS:
        return f"D-{_D_ALDOSE_PATTERNS[pattern]}"
    flipped = tuple(_FLIP_CIP[c] for c in pattern)
    return f"L-{_D_ALDOSE_PATTERNS[flipped]}"


# D-series CIP pattern (C-3 .. C-(n-1), in that order, counted from the
# terminal -CH2OH at C-1) -> retained 2-ketose ('-ulose') stem name. Only
# the 2-ketose retained names (Table 10.3) are covered here -- a ketose
# whose carbonyl sits at C-3 or higher has no single-word retained name
# at all (P-102.5.2.3's own worked examples compose a configurational-
# prefix name instead), out of scope for this step. Confirmed via
# `rdCIPLabeler` against real PubChem open-chain (keto-) structures (see
# module docstring; #1039 M1 step 2).
_D_2_KETOSE_PATTERNS = {
    ("R",): "erythrulose",
    ("R", "R"): "ribulose",
    ("S", "R"): "xylulose",
    ("S", "R", "R"): "fructose",
    ("R", "S", "R"): "sorbose",
    ("R", "R", "R"): "psicose",
    ("S", "S", "R"): "tagatose",
}


def _is_terminal_ch2oh(mol, graph, oxygens, atom_idx):
    atom = mol.GetAtomWithIdx(atom_idx)
    if atom.GetDegree() != 2 or atom.GetTotalNumHs() != 2:
        return False
    atom_oxygens = [n for n in graph[atom_idx] if n in oxygens]
    if len(atom_oxygens) != 1:
        return False
    (oh_oxygen,) = atom_oxygens
    return mol.GetBondBetweenAtoms(atom_idx, oh_oxygen).GetBondTypeAsDouble() == 1.0


def _open_chain_2_ketose_backbone(mol):
    """`[C1, C2, C3, ..., Cn]` (C1/Cn both -CH2OH, C2 the carbonyl), or
    None if `mol` doesn't have this exact shape: acyclic, exactly one
    non-terminal carbonyl carbon with a terminal -CH2OH on one side (C1,
    making this a *2*-ketose specifically, P-102.5.2.1) and an unbranched
    chain of 2-4 more carbons on the other, every one of those bearing
    exactly one single-bonded -OH except the far terminal, which is
    itself a -CH2OH, 4 to 6 carbons total, no heavy atom anywhere in the
    molecule outside that chain and its oxygens."""
    if mol.GetRingInfo().NumRings() > 0:
        return None
    carbonyl_carbons = [
        atom.GetIdx()
        for atom in mol.GetAtoms()
        if atom.GetSymbol() == "C"
        and atom.GetDegree() == 3
        and any(
            bond.GetBondTypeAsDouble() == 2.0 and bond.GetOtherAtom(atom).GetSymbol() == "O"
            for bond in atom.GetBonds()
        )
    ]
    if len(carbonyl_carbons) != 1:
        return None
    (carbonyl,) = carbonyl_carbons

    graph = adjacency(mol)
    oxygens = {atom.GetIdx() for atom in mol.GetAtoms() if atom.GetSymbol() == "O"}
    carbon_neighbors = [n for n in graph[carbonyl] if n not in oxygens]
    if len(carbon_neighbors) != 2:
        return None
    terminal_candidates = [n for n in carbon_neighbors if _is_terminal_ch2oh(mol, graph, oxygens, n)]
    if len(terminal_candidates) != 1:
        return None
    (c1,) = terminal_candidates
    (chain_start,) = [n for n in carbon_neighbors if n != c1]

    rest = ordered_chain(graph, chain_start, carbonyl, oxygens)
    if rest is None or not 2 <= len(rest) <= 4:
        return None
    chain = [c1, carbonyl] + rest
    if mol.GetNumAtoms() != 2 * len(chain):
        return None

    for atom_idx in rest[:-1]:
        chain_oxygens = [n for n in graph[atom_idx] if n in oxygens]
        if len(chain_oxygens) != 1:
            return None
        (oh_oxygen,) = chain_oxygens
        if mol.GetBondBetweenAtoms(atom_idx, oh_oxygen).GetBondTypeAsDouble() != 1.0:
            return None
    if not _is_terminal_ch2oh(mol, graph, oxygens, rest[-1]):
        return None
    return chain


def has_open_chain_2_ketose_shape(mol) -> bool:
    chain = _open_chain_2_ketose_backbone(mol)
    if chain is None:
        return False
    pattern = _cip_pattern(mol, chain[2:-1])
    if pattern is None:
        return False
    if pattern in _D_2_KETOSE_PATTERNS:
        return True
    flipped = tuple(_FLIP_CIP[c] for c in pattern)
    return flipped in _D_2_KETOSE_PATTERNS


def name_open_chain_2_ketose(mol) -> str:
    chain = _open_chain_2_ketose_backbone(mol)
    pattern = _cip_pattern(mol, chain[2:-1])
    if pattern in _D_2_KETOSE_PATTERNS:
        return f"D-{_D_2_KETOSE_PATTERNS[pattern]}"
    flipped = tuple(_FLIP_CIP[c] for c in pattern)
    return f"L-{_D_2_KETOSE_PATTERNS[flipped]}"
