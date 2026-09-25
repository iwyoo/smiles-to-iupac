"""Naming of open-chain (Fischer-projection) monosaccharides -- currently
plain aldoses of 3-7 carbons and plain 2-ketoses of 4-6 carbons, each
named via their retained stem name (3-6 carbons) or, for the 7-carbon
heptoses, a two-segment configurational-prefix name -- per the IUPAC 2013
Recommendations ("the Blue Book"), Chapter P-10
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
- P-102.5.1.1.2: a heptose (7 carbons, 5 chirality centers) has no
  single retained name -- its C-2..C-5 centers (adjacent to the
  aldehyde) and its lone C-6 center (adjacent to the terminal -CH2OH)
  each get their own independent D/L descriptor and configurational-
  prefix word, with the C-6 group's prefix cited first. The C-2..C-5
  group reuses the hexose stem-name table above stripped to its prefix
  form (`glucose` -> `gluco`, etc.); the lone C-6 group uses its own
  table (`glyceraldehyde`'s prefix `glycero`, which doesn't follow that
  same strip rule) -- confirmed against 3 real PubChem heptose
  structures spanning both matching and mixed D/L series (#1042).
- P-102.5.2.1/P-102.5.2.2 (Table 10.3): a *2*-ketose (the carbonyl at
  C-2, one carbon in from the -CH2OH end that becomes C-1) follows the
  identical D/L rule at its own highest-numbered chirality center, and
  its remaining chirality centers (C-3 through C-(n-1)) map to one of a
  much shorter list of retained names (erythrulose, ribulose/xylulose,
  fructose/sorbose/psicose/tagatose) -- confirmed against 13 more real
  PubChem structures (#1039 M1 step 2). A ketose whose carbonyl sits at
  C-3 or higher has no retained name at all (P-102.5.2.3's own examples
  compose a configurational-prefix name instead), out of scope here.
- P-102.3.4.1/P-102.3.4.2: a D-aldohexopyranose (the 6-membered cyclic
  hemiacetal form -- C-1 through C-5 in the ring, the ring oxygen
  bridging C-1 and C-5, C-6 an exocyclic -CH2OH on C-5) reuses the
  open-chain D-hexose machinery above unchanged: the CIP pattern at
  C-2..C-5 (in that ring order) is looked up in the exact same
  `_D_ALDOSE_PATTERNS` table for its D/L descriptor and stem name, and
  the newly-formed anomeric center at C-1 gets its own alpha/beta
  descriptor from its own CIP label alone -- confirmed empirically
  (against 14 real PubChem D-hexopyranose structures spanning 7 stem
  names, both anomers, #1039 M2 step 1) that C-1's CIP is always 'S' for
  alpha and 'R' for beta, regardless of stem name (C-1's own CIP
  priority order never depends on the configuration further round the
  ring: the ring-oxygen branch always outranks the exocyclic-hydroxyl
  branch on reaching a carbon vs. a lone hydrogen one atom out, a purely
  local comparison). A D-aldohexofuranose (the 5-membered ring form --
  C-1 through C-4 in the ring, C-5/C-6 an exocyclic -CH(OH)-CH2OH tail on
  C-4) reuses the identical mechanism, with the ring-closing carbon one
  position earlier (C-4 instead of C-5): the same alpha='S'/beta='R'
  anomeric mapping holds unchanged, but the CIP-flip position shifts to
  C-3 instead of C-4 (confirmed against 4 real PubChem D-hexofuranose
  structures, #1039 M2 step 2) -- see `_ring_stem_pattern`'s own
  `flip_index` parameter. L-series pyranoses/furanoses and ketofuranoses
  are out of scope here (later M2 steps).

Scope, deliberately narrow (first pass at carbohydrate nomenclature,
#1039): a plain, unbranched, acyclic aldose (WS1/M1 steps 1/3) or
2-ketose (WS1/M1 step 2) backbone -- a terminal aldehyde or a C-2
carbonyl flanked by a terminal -CH2OH, 1 to 5 more -CH(OH)-
chirality-bearing carbons, and a second terminal -CH2OH, 3 to 7 (aldose)
or 4 to 6 (2-ketose) carbons total; or a D-aldohexopyranose or
D-aldohexofuranose cyclic form (WS2/M2 steps 1/2) -- see above. Every
stereocenter's configuration must be specified, with no substituent
anywhere beyond each ring/chain carbon's own single -OH (or, at the
terminal/anomeric carbons, none beyond what the aldehyde/-CH2OH/
carbonyl/ring-hemiacetal shape itself requires). A ketose with its
carbonyl at C-3 or higher (#1041's own later scope), any aldose/ketose
beyond 7 carbons (octoses/nonoses/decoses each have their own
group-count shape, later M1 steps), an L-series cyclic form or a cyclic
ketose (later M2 steps), deoxy/amino sugars, glycosides, and any other
substituent are all out of scope here -- a molecule
matching any of those still falls through to whatever it names today
(generic acyclic-aldehyde/ketone/polyol naming for an open chain, or
`_ketone.py`'s hetero-ring-ketone path for most cyclic shapes -- see that
function's own docstring for why a plain pyranose collides with it),
unchanged.
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

# The heptose's lone C-6 group (P-102.5.1.1.2) -- 'glycero', the
# configurational-prefix form of glyceraldehyde, doesn't follow the
# hexose table's own strip-to-prefix rule so it needs its own table.
_D_GLYCERO_PATTERNS = {
    ("R",): "glycero",
}


def _segment_descriptor(pattern, patterns_table):
    """`(D/L, stem-name)` for one heptose group's own CIP pattern against
    `patterns_table` -- the stem word itself never changes between the D
    and L series, only which descriptor letter is cited (P-102.3.3
    applied to a single group rather than the whole molecule)."""
    if pattern in patterns_table:
        return "D", patterns_table[pattern]
    flipped = tuple(_FLIP_CIP[c] for c in pattern)
    return "L", patterns_table[flipped]


def _open_chain_aldose_backbone(mol):
    """The aldose backbone's carbon chain, aldehyde carbon first and
    terminal -CH2OH carbon last, or None if `mol` doesn't have this exact
    shape (see module docstring): acyclic, exactly one terminal aldehyde,
    an unbranched carbon chain of 3-7 atoms from it, every chain carbon
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
    if chain is None or not 3 <= len(chain) <= 7:
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
    if len(pattern) == 5:
        # A heptose's near/far group split (P-102.5.1.1.2) always
        # resolves: the hexose and glycero tables' D-series entries plus
        # their CIP mirrors cover every possible 4- and 1-length pattern.
        return True
    if pattern in _D_ALDOSE_PATTERNS:
        return True
    flipped = tuple(_FLIP_CIP[c] for c in pattern)
    return flipped in _D_ALDOSE_PATTERNS


def name_open_chain_aldose(mol) -> str:
    chain = _open_chain_aldose_backbone(mol)
    pattern = _cip_pattern(mol, chain[1:-1])
    if len(pattern) == 5:
        near, far = pattern[:4], pattern[4:]
        near_dl, near_stem = _segment_descriptor(near, _D_ALDOSE_PATTERNS)
        far_dl, far_prefix = _segment_descriptor(far, _D_GLYCERO_PATTERNS)
        return f"{far_dl}-{far_prefix}-{near_dl}-{near_stem[:-2]}-heptose"
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


# C-1's own CIP label (no D/L-series or stem-name dependence -- see module
# docstring) -> anomeric descriptor.
_ANOMERIC_DESCRIPTORS = {"S": "alpha", "R": "beta"}


def _exocyclic_oxygen(graph, oxygens, ring_set, atom_idx):
    """The single exocyclic oxygen substituent on ring atom `atom_idx`, or
    None if it doesn't have exactly one. Shared by the pyranose/furanose
    ring-shape detectors below."""
    candidates = [n for n in graph[atom_idx] if n in oxygens and n not in ring_set]
    return candidates[0] if len(candidates) == 1 else None


def _exocyclic_carbon(graph, oxygens, ring_set, atom_idx):
    """The single exocyclic carbon substituent on ring atom `atom_idx`, or
    None if it doesn't have exactly one. Shared by the pyranose/furanose
    ring-shape detectors below."""
    candidates = [n for n in graph[atom_idx] if n not in ring_set and n not in oxygens]
    return candidates[0] if len(candidates) == 1 else None


def _pyranose_ring_order(mol):
    """`(C-1, C-2, C-3, C-4, C-5, C-6)` for a plain D-aldohexopyranose
    ring, or None if `mol` doesn't have this exact shape (see module
    docstring): a single saturated, non-aromatic 6-membered ring with
    exactly one ring oxygen; one ring carbon adjacent to it (C-1, the
    anomeric carbon) bearing one exocyclic -OH; the ring oxygen's other
    neighbor (C-5) bearing one exocyclic plain -CH2OH (C-6); and the
    remaining 3 ring carbons (C-2..C-4) each bearing exactly one exocyclic
    -OH -- no other substituent, and no heavy atom anywhere in the
    molecule outside the ring and these oxygens."""
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() != 1:
        return None
    ring = ring_info.AtomRings()[0]
    if len(ring) != 6 or any(mol.GetAtomWithIdx(a).GetIsAromatic() for a in ring):
        return None
    ring_set = set(ring)
    ring_oxygens = [a for a in ring if mol.GetAtomWithIdx(a).GetSymbol() == "O"]
    if len(ring_oxygens) != 1:
        return None
    (ring_oxygen,) = ring_oxygens
    if any(mol.GetAtomWithIdx(a).GetSymbol() != "C" for a in ring if a != ring_oxygen):
        return None

    graph = adjacency(mol)
    oxygens = {atom.GetIdx() for atom in mol.GetAtoms() if atom.GetSymbol() == "O"}
    ring_neighbors = [n for n in graph[ring_oxygen] if n in ring_set]
    if len(ring_neighbors) != 2:
        return None

    def exocyclic_oxygen(atom_idx):
        return _exocyclic_oxygen(graph, oxygens, ring_set, atom_idx)

    def exocyclic_carbon(atom_idx):
        return _exocyclic_carbon(graph, oxygens, ring_set, atom_idx)

    c1 = c5 = None
    for candidate in ring_neighbors:
        if mol.GetAtomWithIdx(candidate).GetDegree() != 3:
            return None
        if exocyclic_oxygen(candidate) is not None:
            c1 = candidate
        elif exocyclic_carbon(candidate) is not None:
            c5 = candidate
    if c1 is None or c5 is None or c1 == c5:
        return None

    anomeric_oxygen = exocyclic_oxygen(c1)
    if mol.GetBondBetweenAtoms(c1, anomeric_oxygen).GetBondTypeAsDouble() != 1.0:
        return None
    if mol.GetAtomWithIdx(anomeric_oxygen).GetTotalNumHs() != 1:
        return None

    c6 = exocyclic_carbon(c5)
    if not _is_terminal_ch2oh(mol, graph, oxygens, c6):
        return None

    for middle_carbon in ring_set - {ring_oxygen, c1, c5}:
        if mol.GetAtomWithIdx(middle_carbon).GetDegree() != 3:
            return None
        if exocyclic_oxygen(middle_carbon) is None:
            return None

    if mol.GetNumAtoms() != 12:
        return None

    rotated = ring[ring.index(ring_oxygen) :] + ring[: ring.index(ring_oxygen)]
    order = rotated[1:]
    if order[0] == c5:
        order = tuple(reversed(order))
    return order + (c6,)


def _furanose_ring_order(mol):
    """`(C-1, C-2, C-3, C-4, C-5, C-6)` for a plain D-aldohexofuranose
    ring, or None if `mol` doesn't have this exact shape (P-102.3.4.1): a
    single saturated, non-aromatic 5-membered ring with exactly one ring
    oxygen; one ring carbon adjacent to it (C-1, the anomeric carbon)
    bearing one exocyclic -OH; the ring oxygen's other neighbor (C-4)
    bearing one exocyclic -CH(OH)- (C-5) that itself bears one plain
    -CH2OH (C-6); and the remaining 2 ring carbons (C-2, C-3 -- one fewer
    than pyranose's 3, since the ring itself is one atom smaller) each
    bearing exactly one exocyclic -OH -- no other substituent, and no
    heavy atom anywhere in the molecule outside the ring, C-5/C-6, and
    these oxygens."""
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() != 1:
        return None
    ring = ring_info.AtomRings()[0]
    if len(ring) != 5 or any(mol.GetAtomWithIdx(a).GetIsAromatic() for a in ring):
        return None
    ring_set = set(ring)
    ring_oxygens = [a for a in ring if mol.GetAtomWithIdx(a).GetSymbol() == "O"]
    if len(ring_oxygens) != 1:
        return None
    (ring_oxygen,) = ring_oxygens
    if any(mol.GetAtomWithIdx(a).GetSymbol() != "C" for a in ring if a != ring_oxygen):
        return None

    graph = adjacency(mol)
    oxygens = {atom.GetIdx() for atom in mol.GetAtoms() if atom.GetSymbol() == "O"}
    ring_neighbors = [n for n in graph[ring_oxygen] if n in ring_set]
    if len(ring_neighbors) != 2:
        return None

    def exocyclic_oxygen(atom_idx):
        return _exocyclic_oxygen(graph, oxygens, ring_set, atom_idx)

    def exocyclic_carbon(atom_idx):
        return _exocyclic_carbon(graph, oxygens, ring_set, atom_idx)

    c1 = c4 = None
    for candidate in ring_neighbors:
        if mol.GetAtomWithIdx(candidate).GetDegree() != 3:
            return None
        if exocyclic_oxygen(candidate) is not None:
            c1 = candidate
        elif exocyclic_carbon(candidate) is not None:
            c4 = candidate
    if c1 is None or c4 is None or c1 == c4:
        return None

    anomeric_oxygen = exocyclic_oxygen(c1)
    if mol.GetBondBetweenAtoms(c1, anomeric_oxygen).GetBondTypeAsDouble() != 1.0:
        return None
    if mol.GetAtomWithIdx(anomeric_oxygen).GetTotalNumHs() != 1:
        return None

    for middle_carbon in ring_set - {ring_oxygen, c1, c4}:
        if mol.GetAtomWithIdx(middle_carbon).GetDegree() != 3:
            return None
        if exocyclic_oxygen(middle_carbon) is None:
            return None

    c5 = exocyclic_carbon(c4)
    if c5 is None or mol.GetAtomWithIdx(c5).GetDegree() != 3:
        return None
    c5_oxygen = _exocyclic_oxygen(graph, oxygens, ring_set | {c5}, c5)
    if c5_oxygen is None or mol.GetBondBetweenAtoms(c5, c5_oxygen).GetBondTypeAsDouble() != 1.0:
        return None
    if mol.GetAtomWithIdx(c5_oxygen).GetTotalNumHs() != 1:
        return None
    c6_candidates = [n for n in graph[c5] if n not in ring_set and n != c5_oxygen]
    if len(c6_candidates) != 1:
        return None
    (c6,) = c6_candidates
    if not _is_terminal_ch2oh(mol, graph, oxygens, c6):
        return None

    if mol.GetNumAtoms() != 12:
        return None

    rotated = ring[ring.index(ring_oxygen) :] + ring[: ring.index(ring_oxygen)]
    order = rotated[1:]
    if order[0] == c4:
        order = tuple(reversed(order))
    return order + (c5, c6)


def _ring_stem_pattern(mol, order, flip_index):
    """The open-chain-equivalent CIP pattern at C-2..C-5, for lookup in
    `_D_ALDOSE_PATTERNS`, or None if any is unspecified. The entry at
    `flip_index` must be flipped relative to what `_cip_pattern` reads
    directly off the ring atom: its priority order depends on which
    neighboring branch (towards the anomeric C-1 side vs. towards the
    ring-closing carbon side) ranks higher, and the ring closure changes
    that ranking without changing the spatial configuration at all -- in
    the open-chain form the ring-closing carbon's own position is a plain
    -CH(OH)-, but in the ring it's bonded to the ring oxygen (which
    itself continues on to the anomeric C-1), a higher-priority
    substituent than the open-chain form's own -OH ever was, which
    reorders that carbon's ring-ward neighbor's two branches and flips
    its CIP label even though nothing physically moved. `flip_index` is
    2 for a pyranose (the ring-closing carbon is C-5, so C-4's label
    flips, #1039 M2 step 1) and 1 for a furanose (the ring-closing carbon
    is C-4, so C-3's label flips, #1039 M2 step 2) -- confirmed
    empirically for each against real D-hexose structures spanning all 8
    stem names (pyranose) or 3 stem names (furanose); every other
    position is unaffected in both ring sizes."""
    pattern = _cip_pattern(mol, order[1:5])
    if pattern is None:
        return None
    return pattern[:flip_index] + (_FLIP_CIP[pattern[flip_index]],) + pattern[flip_index + 1 :]


def has_cyclic_aldopyranose_shape(mol) -> bool:
    order = _pyranose_ring_order(mol)
    if order is None:
        return False
    anomeric_pattern = _cip_pattern(mol, [order[0]])
    if anomeric_pattern is None or anomeric_pattern[0] not in _ANOMERIC_DESCRIPTORS:
        return False
    pattern = _ring_stem_pattern(mol, order, flip_index=2)
    return pattern is not None and pattern in _D_ALDOSE_PATTERNS


def name_cyclic_aldopyranose(mol) -> str:
    order = _pyranose_ring_order(mol)
    anomeric_descriptor = _ANOMERIC_DESCRIPTORS[_cip_pattern(mol, [order[0]])[0]]
    stem = _D_ALDOSE_PATTERNS[_ring_stem_pattern(mol, order, flip_index=2)]
    return f"{anomeric_descriptor}-D-{stem[:-2]}pyranose"


def has_cyclic_aldofuranose_shape(mol) -> bool:
    order = _furanose_ring_order(mol)
    if order is None:
        return False
    anomeric_pattern = _cip_pattern(mol, [order[0]])
    if anomeric_pattern is None or anomeric_pattern[0] not in _ANOMERIC_DESCRIPTORS:
        return False
    pattern = _ring_stem_pattern(mol, order, flip_index=1)
    return pattern is not None and pattern in _D_ALDOSE_PATTERNS


def name_cyclic_aldofuranose(mol) -> str:
    order = _furanose_ring_order(mol)
    anomeric_descriptor = _ANOMERIC_DESCRIPTORS[_cip_pattern(mol, [order[0]])[0]]
    stem = _D_ALDOSE_PATTERNS[_ring_stem_pattern(mol, order, flip_index=1)]
    return f"{anomeric_descriptor}-D-{stem[:-2]}furanose"
