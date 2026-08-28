"""Naming of a two-superatom phane assembly with *different* component
rings -- one naphthalene "superatom" and one benzene "superatom", joined
by two bridges of any (possibly unequal) length -- per the IUPAC 2013
Recommendations ("the Blue Book"):

- P-26 (Chapter P-2, https://iupac.qmul.ac.uk/BlueBook/PDF/P2.pdf; full
  text confirmed accessible 2026-08-29 via
  https://iupac.qmul.ac.uk/BlueBook/P2.html after earlier sessions'
  repeated access failures -- see tasks/phane-naming.md's own log). This
  is the "different ring kinds, asymmetric attachment points" shape
  `_cyclophane.py` explicitly left for this broader case.
- Confirmed via a direct Blue Book worked example (P-26.4.1.4):
  '1(2,7)-naphthalena-4(1,4)-benzenacycloheptaphane (PIN)' -- no
  PubChem-listed compound exists for this exotic macrocycle (checked,
  including a from-scratch structure search), so this is settled by the
  primary-source worked example itself, the same standard already used
  elsewhere in this project for a primary-source quotation without an
  independent structural cross-check (see e.g. `_phosphanone.py`).
- P-44.2: a polycyclic ring system (naphthalene) is senior to a
  monocyclic one (benzene), so naphthalene always receives the lower
  superatom locant, '1' -- confirmed directly by the worked example's own
  annotation ("the senior amplificant, naphthalene, must be assigned the
  superatom locant '1'").
- P-26.4.1.2: naphthalene's own two attachment atoms are labeled using
  naphthalene's *own* standard numbering, choosing whichever of its four
  automorphism-equivalent numbering assignments gives the lowest locant
  set for the two real attachment atoms (confirmed pattern: '(2,7)', a
  reviewed, not independently confirmed extension for any other pair of
  non-fusion positions).
- P-26.4.1.1 (lowest overall superatom locant set): with naphthalene
  fixed at superatom locant 1, the macrocycle numbering direction is
  chosen so the *shorter* bridge is traversed first, giving benzene the
  lower of its two possible locants (2 + shorter-bridge-length) --
  confirmed directly: bridges of length 2 and 3 give benzene locant 4,
  not 5, matching the worked example. Equal-length bridges (a genuine
  symmetry tie) fall back to P-26.4.1.4's adjacency rule below.
- P-26.3.2.2 / P-26.4.1.4 (attachment-locant citation order): once the
  macrocycle direction is fixed, each superatom's own two attachment
  locants are cited with whichever one sits on the *lower-numbered*
  (forward) macrocycle side written first -- confirmed directly:
  naphthalene's citation '(2,7)' (not '(7,2)') places locant 2 (on the
  shorter, forward bridge) first; benzene's own para/meta/ortho pattern
  is symmetric under reflection, so its citation is always the plain
  '1,4'/'1,3'/'1,2' regardless of bridge direction.
- The 'naphthalena'/'benzena' amplification-prefix spelling elides the
  parent name's own final 'e' before the vowel-starting suffix '-a'
  (P-16.3.3), confirmed directly by the worked example's own spelling.

Scope, deliberately narrow: exactly one naphthalene ring (unsubstituted
except at its own two bridge-attachment atoms, which must not be
ring-fusion atoms 4a/8a) and exactly one benzene ring (same restriction),
joined by exactly two plain unbranched saturated -CH2-...-CH2- bridges of
any length (equal or unequal). Explicitly out of scope (raise
`UnsupportedStructure`, i.e. the detector below simply returns None):
- More than two rings, rings sharing atoms (fused, not phane), or a ring
  bridged back to itself.
- Attachment at a naphthalene ring-fusion atom (4a/8a) or any ortho (1,2)
  naphthalene attachment pattern spanning the two rings in a way this
  module's automorphism search can't resolve unambiguously.
- Any ring other than naphthalene/benzene, any substituent, or a bridge
  that isn't a plain unbranched saturated chain.
- Two ambiguous naphthalene automorphisms both achieving the same lowest
  attachment locant set with different atom-to-locant assignments (a
  genuine numbering ambiguity this module doesn't attempt to resolve
  further) -- not encountered by the confirmed worked example, but
  possible for an unverified attachment pattern.
"""

from rdkit import Chem

from ._common import UnsupportedStructure, adjacency
from ._numerals import numerical_term
from ._polyspiro import _ring_cyclic_order

_NAPHTHALENE_REF = Chem.MolFromSmiles("c1ccc2ccccc2c1")
# Atom index (in the reference SMILES above) -> canonical naphthalene
# locant; verified directly via RDKit bond inspection (see this module's
# own development notes / tasks/naphthalene-benzene-phane-naming.md).
_NAPHTHALENE_ROLE = {9: 1, 0: 2, 1: 3, 2: 4, 3: "4a", 4: 5, 5: 6, 6: 7, 7: 8, 8: "8a"}
_LOCAL_LOCANTS = {"ortho": "1,2", "meta": "1,3", "para": "1,4"}


def _find_attachments(mol, graph, ring_atoms):
    """Classify every atom of `ring_atoms` as a plain periphery CH, a
    ring-fusion atom (all-in-ring neighbors, e.g. naphthalene's own 4a/
    8a -- always allowed, never counted as an attachment), or a single-
    substituent attachment atom. Returns the list of (ring_atom,
    exocyclic_neighbor) attachment pairs, or None if there aren't exactly
    two or any atom doesn't fit one of these three shapes."""
    atts = []
    for a in ring_atoms:
        atom = mol.GetAtomWithIdx(a)
        if atom.GetFormalCharge() != 0 or atom.GetIsotope() != 0:
            return None
        in_ring = [nb for nb in graph[a] if nb in ring_atoms]
        exo = [nb for nb in graph[a] if nb not in ring_atoms]
        if len(in_ring) == 3 and not exo:
            continue
        if len(in_ring) == 2 and not exo:
            if atom.GetTotalNumHs() != 1:
                return None
            continue
        if len(in_ring) == 2 and len(exo) == 1:
            if atom.GetTotalNumHs() != 0:
                return None
            atts.append((a, exo[0]))
            continue
        return None
    if len(atts) != 2:
        return None
    return atts


def _naphthalene_attachment_locants(mol, att_a, att_b, forward_atom):
    """Try every automorphism-consistent mapping of the naphthalene
    reference onto `mol` and return ((lo, hi), {att_a: label, att_b:
    label}) for the lowest achievable locant set, or None if either
    attachment resolves to a ring-fusion position (4a/8a) in every
    matching automorphism. Among numbering choices tied on that lowest
    set, prefer the one that assigns the lower individual locant to
    `forward_atom` (P-26.4.1.4's citation-order rule) -- naphthalene's
    own symmetry lets more than one automorphism reach the same {lo, hi}
    set while swapping which physical attachment atom is 'lo' vs 'hi'."""
    matches = mol.GetSubstructMatches(_NAPHTHALENE_REF, uniquify=False, useChirality=False)
    candidates = []
    for match in matches:
        target_to_ref = {match[i]: i for i in range(len(match))}
        if att_a not in target_to_ref or att_b not in target_to_ref:
            continue
        label_a = _NAPHTHALENE_ROLE[target_to_ref[att_a]]
        label_b = _NAPHTHALENE_ROLE[target_to_ref[att_b]]
        if isinstance(label_a, str) or isinstance(label_b, str):
            continue
        candidates.append({att_a: label_a, att_b: label_b})
    if not candidates:
        return None
    best_pair = min(tuple(sorted(labels.values())) for labels in candidates)
    tied = [labels for labels in candidates if tuple(sorted(labels.values())) == best_pair]
    best_labels = min(tied, key=lambda labels: labels[forward_atom])
    return best_pair, best_labels


def _walk_bridge(mol, graph, start_ring_atom, start_bridge_atom, other_ring_atoms):
    chain = [start_bridge_atom]
    previous, current = start_ring_atom, start_bridge_atom
    while True:
        atom = mol.GetAtomWithIdx(current)
        if (
            atom.GetAtomicNum() != 6
            or atom.GetIsAromatic()
            or atom.GetFormalCharge() != 0
            or atom.GetIsotope() != 0
            or atom.GetTotalNumHs() != 2
        ):
            return None
        neighbors = [nb for nb in graph[current] if nb != previous]
        if len(neighbors) != 1:
            return None
        next_atom = neighbors[0]
        if next_atom in other_ring_atoms:
            return len(chain), next_atom
        if len(chain) > 20:
            return None
        chain.append(next_atom)
        previous, current = current, next_atom


def _find_naphthalene_benzene_phane(mol):
    naph_matches = mol.GetSubstructMatches(_NAPHTHALENE_REF, uniquify=True, useChirality=False)
    if len(naph_matches) != 1:
        return None
    naph_atoms = set(naph_matches[0])

    ring_info = mol.GetRingInfo()
    benzene_rings = [
        set(ring)
        for ring in ring_info.AtomRings()
        if len(ring) == 6
        and not (set(ring) & naph_atoms)
        and all(mol.GetAtomWithIdx(a).GetAtomicNum() == 6 and mol.GetAtomWithIdx(a).GetIsAromatic() for a in ring)
    ]
    if len(benzene_rings) != 1:
        return None
    benzene_atoms = benzene_rings[0]
    if mol.GetNumAtoms() < len(naph_atoms) + len(benzene_atoms):
        return None

    graph = adjacency(mol)
    naph_atts = _find_attachments(mol, graph, naph_atoms)
    benzene_atts = _find_attachments(mol, graph, benzene_atoms)
    if naph_atts is None or benzene_atts is None:
        return None

    b_cyc = _ring_cyclic_order(graph, list(benzene_atoms), benzene_atts[0][0])
    b_dist = b_cyc.index(benzene_atts[1][0])
    b_dist = min(b_dist, len(b_cyc) - b_dist)
    benzene_pattern = {1: "ortho", 2: "meta", 3: "para"}.get(b_dist)
    if benzene_pattern is None:
        return None

    (na, na_exo), (nb, nb_exo) = naph_atts

    bridge_na = _walk_bridge(mol, graph, na, na_exo, benzene_atoms)
    bridge_nb = _walk_bridge(mol, graph, nb, nb_exo, benzene_atoms)
    if bridge_na is None or bridge_nb is None:
        return None
    len_na, benzene_end_from_na = bridge_na
    len_nb, benzene_end_from_nb = bridge_nb
    if benzene_end_from_na == benzene_end_from_nb:
        return None

    # P-26.4.1.1: the shorter bridge goes forward (gives benzene the lower
    # locant). This is decided purely from bridge length -- independent of
    # naphthalene's own numbering choice, which is resolved next.
    if len_nb < len_na:
        forward_atom, forward_len = nb, len_nb
        backward_atom, backward_len = na, len_na
    else:
        forward_atom, forward_len = na, len_na
        backward_atom, backward_len = nb, len_nb

    result = _naphthalene_attachment_locants(mol, na, nb, forward_atom)
    if result is None:
        return None
    _, naph_labels = result

    macrocycle_size = 2 + forward_len + backward_len
    benzene_locant = 2 + forward_len
    naph_citation = f"{naph_labels[forward_atom]},{naph_labels[backward_atom]}"
    benzene_local = _LOCAL_LOCANTS[benzene_pattern]

    return (
        f"1({naph_citation})-naphthalena-{benzene_locant}({benzene_local})-"
        f"benzenacyclo{numerical_term(macrocycle_size)}phane"
    )


def has_naphthalene_benzene_phane_name(mol) -> bool:
    return _find_naphthalene_benzene_phane(mol) is not None


def name_naphthalene_benzene_phane(mol) -> str:
    name = _find_naphthalene_benzene_phane(mol)
    if name is None:
        raise UnsupportedStructure(
            "this is not a supported naphthalene+benzene two-bridge phane shape (see P-26)"
        )
    return name
