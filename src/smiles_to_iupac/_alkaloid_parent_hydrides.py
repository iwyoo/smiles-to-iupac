"""Morphinan retained parent hydride (Appendix 3/P-101, `morphinan.gif`,
locants 1-17, N=17), skeleton-dict matched like `_steroid_parent_hydrides.py`.
Cites the transannular O-bridge as P-25.4.2.1.5's `furo[2',3',4',5':...]`
heterocyclic-bridge form when the bridge closes a furan-shaped ring
(morphine/codeine's real shape), falling back to the plain `epoxy` prefix
otherwise. Stereodescriptors (alpha/beta on the bridge/suffix) are out of
scope: no local-attachment citation shape exists yet in this project.
"""

from rdkit import Chem
from rdkit.Chem import BondType, RWMol, rdmolops

from ._parent_hydride_stripping import strip_substituents

_FURAN_PRIMES = ("2′", "3′", "4′", "5′")


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
_LOCANT_OF_INTERNAL = {0: 17, 3: 13, 4: 14, 10: 3, 11: 4, 12: 12, 13: 5, 14: 6, 15: 7, 16: 8}


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
        "furan_bridge_path": None,
    }

    for ring_atom_old in o_substituents:
        new_idx = old_to_new.get(ring_atom_old)
        locant = locant_of_new_atom.get(new_idx)
        if locant is None:
            return None
        result["o_substituents"][locant] = o_substituent_kind[ring_atom_old]

    bridge_neighbors = [n.GetIdx() for n in mol.GetAtomWithIdx(bridge_atom).GetNeighbors()]
    bridge_locants = []
    bridge_new_idx = []
    for n in bridge_neighbors:
        new_idx = old_to_new.get(n)
        locant = locant_of_new_atom.get(new_idx)
        if locant is None:
            return None
        bridge_locants.append(locant)
        bridge_new_idx.append(new_idx)
    if len(bridge_locants) != 2:
        return None
    result["bridge_locants"] = tuple(sorted(bridge_locants))

    # P-25.4.2.1.5: a transannular bridge whose own path plus the bridge
    # atom closes a 5-membered ring (furan's own shape) is cited by furan's
    # conventional 2',3',4',5' locants against that path, not the plain
    # bridge prefix -- pick whichever walk direction gives lower locants.
    path = list(rdmolops.GetShortestPath(saturated, bridge_new_idx[0], bridge_new_idx[1]))
    if len(path) == 4:
        path_locants = [locant_of_new_atom.get(i) for i in path]
        if all(l is not None for l in path_locants):
            reversed_locants = list(reversed(path_locants))
            result["furan_bridge_path"] = min(path_locants, reversed_locants)

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

    hydroxy_locants = sorted(l for l, kind in subs.items() if kind == "hydroxy")
    if len(hydroxy_locants) == 2:
        suffix = f"-{hydroxy_locants[0]},{hydroxy_locants[1]}-diol"
    elif len(hydroxy_locants) == 1:
        suffix = f"-{hydroxy_locants[0]}-ol"
    else:
        suffix = ""

    if core["furan_bridge_path"]:
        path_locants = ",".join(str(l) for l in core["furan_bridge_path"])
        bridge_prefix = f"furo[{','.join(_FURAN_PRIMES)}:{path_locants}]"
        return methoxy_prefix + n17_prefix + didehydro_prefix + bridge_prefix + "morphinan" + suffix

    bridge_a, bridge_b = core["bridge_locants"]
    bridge_prefix = f"{bridge_a},{bridge_b}-epoxy-"
    return bridge_prefix + methoxy_prefix + n17_prefix + didehydro_prefix + "morphinan" + suffix


def has_alkaloid_morphinan_name(mol) -> bool:
    return name_alkaloid_morphinan(mol) is not None
