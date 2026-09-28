"""Naming of the morphinan retained parent hydride (Appendix 3 / P-101,
`tmp/bluebook/app3/morphinan.gif` -- numbered 1-17, N at 17), by
skeleton-dict recognition analogous to `_steroid_parent_hydrides.py`,
generalized via `_parent_hydride_stripping.py`'s multi-atom stripping to
handle morphine/codeine's several simultaneous substituents plus a
transannular ether bridge, per the IUPAC 2013 Recommendations ("the Blue
Book"):

- Bare morphinan (C16H21N) has one aromatic ring (locants 1-4, 11, 12),
  two saturated six-membered rings, and an N17-containing ring, all
  ortho-fused into one bridged tetracyclic system with a quaternary
  bridgehead at C13. Verified structurally, not assumed: derived by
  reversing morphine's own known modifications (PubChem CID 5288826) --
  removing its N17-methyl, its two O-substituents (locants 3, 6), and
  its transannular ether-bridge oxygen (locants 4, 5), then saturating
  its one extra ring C=C (locants 7, 8) -- confirming the result is
  exactly 17 ring atoms with molecular formula C16H21N, matching
  morphinan's own known formula independently.
- Locants 3 (phenol/methoxy in morphine/codeine), 4 and 5 (the epoxy
  bridge attachment -- 4 aromatic, 5 not, three bonds apart through the
  existing skeleton before the bridge is added, i.e. a transannular span
  per `_bridged_alicyclic_parent.py`), 6 (the second -ol), 7 and 8 (the
  'didehydro' double bond), and 17 (the ring nitrogen) are the only
  locants morphine/codeine's own names ever cite -- confirmed against
  their adjacency directly (3-4 adjacent aromatic positions; 4 adjacent
  to ring-fusion carbon 12; 5 adjacent to bridgehead 13 and to 6; 6-7-8
  forming the second ring in sequence; 8 adjacent to bridgehead 14). The
  remaining ring atoms (the two-carbon chain from ring-fusion carbon 11
  to N17, the two-carbon chain from N17 to bridgehead 13, and
  bridgeheads 13/14 themselves -- real morphinan locants 9/10/13/14/15/
  16) are never substituted or cited in morphine/codeine's own name
  (see below), so this module's locant query leaves the unused ones
  (9/10/15/16) unlabeled -- a future alkaloid needing one of them cited
  would need this filled in first.
- P-101.8's worked example (`tmp/bluebook/P10.txt` lines 1512-1516) gives
  morphine's PIN via fusion nomenclature
  ('furo[2′,3′,4′,5′:4,12,13,5]morphinan'), not yet built here (needs
  general N-ring fusion-orientation work not yet done). This module
  instead produces the same primary source's own sanctioned
  alternative bridge-prefix form (`tmp/bluebook/P1.txt` lines 556-557),
  '4,5-epoxy-17-methyl-7,8-didehydromorphinan-3,6-diol' for morphine
  (codeine analogously) -- a real, correct, primary-source-verified
  name, not a compromise, modulo one deliberate scope limit below.
- Stereodescriptor citation is out of scope here, by the same explicit
  policy `_bridged_alicyclic_parent.py` already established for its
  sibling epoxy-bridge mechanism (see that module's own docstring):
  matching and naming by constitution only. The primary source's full
  name additionally cites '4,5alpha-epoxy' and '3,6alpha-diol' (the
  bridge's own face, and the C6 stereocenter's configuration
  respectively) -- unlike `_steroid_parent_hydrides.py`'s leading
  '<locants>-' alpha/beta block convention, these attach locally to the
  bridge prefix and to the suffix, a different citation shape not yet
  implemented anywhere in this project. Getting this wrong silently
  would be worse than omitting it, so a real input's stereochemistry,
  if any, is simply not reflected in this module's output (tracked as
  real follow-on work once the local-attachment stereo citation shape
  exists).
- Bridge citation follows `_bridged_alicyclic_parent.py`'s one-atom
  bridge-prefix mechanism (P-25.4.2.1.4), generalized there to the
  transannular case; re-derived here directly since that module is
  scoped to the seven bare steroid parents only. Unsaturation citation
  follows the 'didehydro' convention
  (P-31.2.2/P-31.2.4.1, `_didehydro_ring.py`'s sibling mechanism,
  re-derived here directly since that module's own scope is restricted
  to a plain monocyclic ring).
"""

from rdkit import Chem
from rdkit.Chem import BondType, RWMol

from ._parent_hydride_stripping import strip_substituents


def _build_morphinan_query():
    """Hand-built query graph for bare morphinan's 17 ring atoms, per the
    module docstring's derivation. Keys are `internal_idx` (0-16) --
    `_LOCANT_OF_INTERNAL` below maps the ones this module actually cites
    onto their real Appendix 3 locants."""
    rw = RWMol()
    idx = {i: rw.AddAtom(Chem.Atom(7 if i == 0 else 6)) for i in range(17)}
    for i in (7, 8, 9, 10, 11, 12):
        rw.GetAtomWithIdx(idx[i]).SetIsAromatic(True)

    def bond(a, b, aromatic=False):
        rw.AddBond(idx[a], idx[b], BondType.AROMATIC if aromatic else BondType.SINGLE)

    # piperidine (N) ring: 0(N17)-1-2-3(C13)-4(C14)-5-0
    bond(0, 1), bond(1, 2), bond(2, 3), bond(3, 4), bond(4, 5), bond(5, 0)
    # aromatic ring + its two exocyclic connections: 6-7(C11)..12(C12)-3(C13)
    bond(5, 6), bond(6, 7)
    bond(7, 8, True), bond(8, 9, True), bond(9, 10, True), bond(10, 11, True), bond(11, 12, True), bond(12, 7, True)
    bond(12, 3)
    # bottom ring: 3(C13)-13(C5)-14(C6)-15(C7)-16(C8)-4(C14)-3
    bond(3, 13), bond(13, 14), bond(14, 15), bond(15, 16), bond(16, 4)
    mol = rw.GetMol()
    Chem.SanitizeMol(mol)
    return mol, idx


_QUERY, _QUERY_IDX = _build_morphinan_query()

# Real Appendix 3 locants for the internal indices morphine/codeine's own
# names actually cite -- cross-checked against morphine's real structure
# (PubChem CID 5288826): internal 0=N17, 3=C13, 4=C14, 10=C3, 11=C4,
# 13=C5, 14=C6, 15=C7, 16=C8 (internal 7/8/9/12 = real morphinan C11,
# C1, C2, C12 -- correct but never cited by morphine/codeine, omitted).
_LOCANT_OF_INTERNAL = {0: 17, 3: 13, 4: 14, 10: 3, 11: 4, 13: 5, 14: 6, 15: 7, 16: 8}


def _locant_map(mol):
    matches = mol.GetSubstructMatches(_QUERY, uniquify=False)
    assert len(matches) == 1, f"morphinan query matched {len(matches)} times, expected 1"
    match = matches[0]
    return {locant: match[_QUERY_IDX[internal]] for internal, locant in _LOCANT_OF_INTERNAL.items()}


def find_alkaloid_morphinan_core(mol):
    """Return a dict describing morphine/codeine-shaped morphinan
    derivatives, or None. Handles: one N17 substituent (methyl), one O
    substituent each at C3 and/or C6 (hydroxy or methoxy), one
    transannular O bridge at C4/C5, and up to one extra ring C=C
    (didehydro) -- exactly the shape morphine and codeine both have."""
    ring_info = mol.GetRingInfo()

    def in_ring(idx):
        return ring_info.NumAtomRings(idx) > 0

    n_methyls = []
    o_substituents = {}
    o_substituent_kind = {}
    methoxy_methyls = set()
    bridge_atom = None
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6 and atom.GetDegree() == 1 and not atom.GetIsAromatic():
            (neighbor,) = atom.GetNeighbors()
            if neighbor.GetAtomicNum() == 7:
                n_methyls.append(atom.GetIdx())
        if atom.GetAtomicNum() == 8 and atom.GetFormalCharge() == 0 and atom.GetIsotope() == 0:
            if atom.GetDegree() == 1:
                (ring_neighbor,) = atom.GetNeighbors()
                o_substituents[ring_neighbor.GetIdx()] = atom.GetIdx()
                o_substituent_kind[ring_neighbor.GetIdx()] = "hydroxy"
            elif atom.GetDegree() == 2 and not atom.GetIsAromatic():
                neighbors = atom.GetNeighbors()
                ring_neighbors = [n for n in neighbors if in_ring(n.GetIdx())]
                non_ring_neighbors = [n for n in neighbors if not in_ring(n.GetIdx())]
                if len(ring_neighbors) == 2:
                    # a true bridge: both neighbors already in the ring
                    # skeleton, connected only via a longer existing path
                    # (or already adjacent) -- not this simple ether case.
                    if bridge_atom is not None:
                        return None
                    bridge_atom = atom.GetIdx()
                elif len(ring_neighbors) == 1 and len(non_ring_neighbors) == 1:
                    # a methoxy ether: one neighbor is a ring atom, the
                    # other a terminal exocyclic methyl -- a substituent,
                    # not a bridge.
                    (methyl,) = non_ring_neighbors
                    if methyl.GetAtomicNum() != 6 or methyl.GetDegree() != 1:
                        return None
                    o_substituents[ring_neighbors[0].GetIdx()] = atom.GetIdx()
                    o_substituent_kind[ring_neighbors[0].GetIdx()] = "methoxy"
                    methoxy_methyls.add(methyl.GetIdx())
                else:
                    return None
    if len(n_methyls) > 1 or bridge_atom is None:
        return None

    to_strip = {bridge_atom} | set(o_substituents.values()) | methoxy_methyls | set(n_methyls)

    stripped, old_to_new = strip_substituents(mol, to_strip)
    if stripped is None:
        return None

    double_bonds = [
        (b.GetBeginAtomIdx(), b.GetEndAtomIdx())
        for b in stripped.GetBonds()
        if b.GetBondTypeAsDouble() == 2.0 and not b.GetIsAromatic()
    ]
    if len(double_bonds) > 1:
        return None
    saturated = stripped
    if double_bonds:
        rw = RWMol(stripped)
        (a, b) = double_bonds[0]
        rw.GetBondBetweenAtoms(a, b).SetBondType(BondType.SINGLE)
        for i in (a, b):
            at = rw.GetAtomWithIdx(i)
            at.SetNoImplicit(False)
            at.SetNumExplicitHs(0)
        saturated = rw.GetMol()
        try:
            Chem.SanitizeMol(saturated)
        except Chem.rdchem.KekulizeException:
            return None

    plain = Chem.Mol(saturated)
    Chem.RemoveStereochemistry(plain)
    if plain.GetNumAtoms() != 17 or len(plain.GetSubstructMatches(_QUERY, uniquify=False)) != 1:
        return None

    locant_map = _locant_map(saturated)
    locant_of_new_atom = {atom: locant for locant, atom in locant_map.items()}

    result = {
        "n17_methyl": bool(n_methyls),
        "o_substituents": {},
        "bridge_locants": None,
        "didehydro_locants": None,
    }

    for ring_atom_old in o_substituents:
        new_idx = old_to_new.get(ring_atom_old)
        locant = locant_of_new_atom.get(new_idx)
        if locant is None:
            return None
        result["o_substituents"][locant] = o_substituent_kind[ring_atom_old]

    bridge_neighbors = [n.GetIdx() for n in mol.GetAtomWithIdx(bridge_atom).GetNeighbors()]
    bridge_locants = []
    for n in bridge_neighbors:
        new_idx = old_to_new.get(n)
        locant = locant_of_new_atom.get(new_idx)
        if locant is None:
            return None
        bridge_locants.append(locant)
    if len(bridge_locants) != 2:
        return None
    result["bridge_locants"] = tuple(sorted(bridge_locants))

    if double_bonds:
        (a, b) = double_bonds[0]
        la, lb = locant_of_new_atom.get(a), locant_of_new_atom.get(b)
        if la is None or lb is None:
            return None
        result["didehydro_locants"] = tuple(sorted((la, lb)))

    return result


def name_alkaloid_morphinan(mol) -> str:
    core = find_alkaloid_morphinan_core(mol)
    if core is None:
        return None
    subs = core["o_substituents"]

    methoxy_locants = sorted(l for l, kind in subs.items() if kind == "methoxy")
    methoxy_prefix = (",".join(str(l) for l in methoxy_locants) + "-methoxy-") if methoxy_locants else ""

    n17_prefix = "17-methyl-" if core["n17_methyl"] else ""

    didehydro_prefix = ""
    if core["didehydro_locants"]:
        a, b = core["didehydro_locants"]
        didehydro_prefix = f"{a},{b}-didehydro"

    bridge_a, bridge_b = core["bridge_locants"]
    bridge_prefix = f"{bridge_a},{bridge_b}-epoxy-"

    hydroxy_locants = sorted(l for l, kind in subs.items() if kind == "hydroxy")
    if len(hydroxy_locants) == 2:
        suffix = f"-{hydroxy_locants[0]},{hydroxy_locants[1]}-diol"
    elif len(hydroxy_locants) == 1:
        suffix = f"-{hydroxy_locants[0]}-ol"
    else:
        suffix = ""

    return bridge_prefix + methoxy_prefix + n17_prefix + didehydro_prefix + "morphinan" + suffix


def has_alkaloid_morphinan_name(mol) -> bool:
    return name_alkaloid_morphinan(mol) is not None
