"""Naming of a two-superatom phane assembly with *different* component
rings -- one "senior" amplificant superatom (naphthalene or, #1031,
pyridine) and one benzene superatom, joined by two bridges of any
(possibly unequal) length -- per the IUPAC 2013 Recommendations ("the
Blue Book"):

- P-26 (Chapter P-2, https://iupac.qmul.ac.uk/BlueBook/PDF/P2.pdf). This
  is the "different ring kinds, asymmetric attachment points" shape
  `_cyclophane.py` explicitly left for this broader case.
- Confirmed via a direct Blue Book worked example (P-26.4.1.4):
  '1(2,7)-naphthalena-4(1,4)-benzenacycloheptaphane (PIN)' -- no
  PubChem-listed compound exists for this exotic macrocycle (checked,
  including a from-scratch structure search), so this is settled by the
  primary-source worked example itself, the same standard already used
  elsewhere in this project for a primary-source quotation without an
  independent structural cross-check (see e.g. `_phosphanone.py`).
- P-44.2: a polycyclic ring system (naphthalene) or a nitrogen-containing
  heteromonocycle (pyridine) is senior to an all-carbon monocycle
  (benzene), so the non-benzene amplificant always receives the lower
  superatom locant, '1' -- naphthalene's own precedence confirmed
  directly by the worked example's own annotation ("the senior
  amplificant, naphthalene, must be assigned the superatom locant '1'");
  pyridine's precedence over benzene follows the same general
  seniority-of-rings ordering and is consistent with every P-26 worked
  example in `tmp/bluebook/P2.txt` where a listed heterocyclic
  amplificant (pyrimidine, pyridine, furan, quinoline, piperidine) always
  receives a lower superatom locant than a co-occurring carbocyclic one
  (naphthalene, phenanthrene, benzene) -- no literal 2-ring pyridine+
  benzene worked example exists in the source, so this specific pairing's
  expected name is derived by hand from the rules below rather than
  quoted directly (#1031).
- P-26.4.1.2: the senior amplificant's own two attachment atoms are
  labeled using its *own* standard numbering, choosing whichever
  automorphism-equivalent numbering assignment gives the lowest locant
  set for the two real attachment atoms (confirmed pattern for
  naphthalene: '(2,7)'; pyridine's own numbering fixes its nitrogen at
  locant 1 with no fusion atoms to exclude, so any two of its five carbon
  positions are eligible attachment points).
- P-26.4.1.1 (lowest overall superatom locant set): with the senior
  amplificant fixed at superatom locant 1, the macrocycle numbering
  direction is chosen so the *shorter* bridge is traversed first, giving
  benzene the lower of its two possible locants (2 + shorter-bridge-
  length) -- confirmed directly for naphthalene: bridges of length 2 and
  3 give benzene locant 4, not 5, matching the worked example. Equal-
  length bridges (a genuine symmetry tie) fall back to P-26.4.1.4's
  adjacency rule below.
- P-26.3.2.2 / P-26.4.1.4 (attachment-locant citation order): once the
  macrocycle direction is fixed, each superatom's own two attachment
  locants are cited with whichever one sits on the *lower-numbered*
  (forward) macrocycle side written first -- confirmed directly:
  naphthalene's citation '(2,7)' (not '(7,2)') places locant 2 (on the
  shorter, forward bridge) first; benzene's own para/meta/ortho pattern
  is symmetric under reflection, so its citation is always the plain
  '1,4'/'1,3'/'1,2' regardless of bridge direction.
- The 'naphthalena'/'pyridina'/'benzena' amplification-prefix spelling
  elides the parent name's own final 'e' before the vowel-starting suffix
  '-a' (P-16.3.3), confirmed directly by the worked example's own
  spelling ('naphthalena') and by every 'pyridina'/'tripyridina' spelling
  in `tmp/bluebook/P2.txt`.

Scope, deliberately narrow: exactly one senior-amplificant ring
(naphthalene or pyridine, unsubstituted except at its own two bridge-
attachment atoms, which for naphthalene must not be ring-fusion atoms
4a/8a) and exactly one benzene ring (same restriction), joined by exactly
two plain unbranched saturated -CH2-...-CH2- bridges of any length (equal
or unequal). Explicitly out of scope (raise `UnsupportedStructure`, i.e.
the detector below simply returns None):
- More than two rings, rings sharing atoms (fused, not phane), or a ring
  bridged back to itself.
- Attachment at a naphthalene ring-fusion atom (4a/8a) or any ortho (1,2)
  naphthalene attachment pattern spanning the two rings in a way this
  module's automorphism search can't resolve unambiguously.
- Any ring other than {naphthalene, pyridine}/benzene, any substituent, or
  a bridge that isn't a plain unbranched saturated chain. A third
  amplificant kind, or more than two rings/bridges total, is separate,
  unverified follow-up work (see epic #1026 WS2) -- not attempted here.
- Two ambiguous automorphisms of the senior amplificant both achieving the
  same lowest attachment locant set with different atom-to-locant
  assignments (a genuine numbering ambiguity this module doesn't attempt
  to resolve further) -- not encountered by the confirmed worked example,
  but possible for an unverified attachment pattern.
"""

from rdkit import Chem

from ._common import UnsupportedStructure, adjacency
from ._numerals import numerical_term
from ._polyspiro import _ring_cyclic_order


class _AmplificantKind:
    """A registered "senior" amplificant ring kind (see module docstring):
    `reference` is a plain, unsubstituted RDKit mol for substructure
    matching; `role` maps each of `reference`'s own atom indices to its
    standard locant (an int, or a string like '4a' for a ring-fusion
    position that's never an eligible attachment point); `prefix` is the
    P-26.2.2.1 amplification-prefix spelling."""

    def __init__(self, reference_smiles, role, prefix):
        self.reference = Chem.MolFromSmiles(reference_smiles)
        self.role = role
        self.prefix = prefix


# Atom index (in the reference SMILES) -> canonical naphthalene locant;
# verified directly via RDKit bond inspection.
_NAPHTHALENE_KIND = _AmplificantKind(
    "c1ccc2ccccc2c1",
    {9: 1, 0: 2, 1: 3, 2: 4, 3: "4a", 4: 5, 5: 6, 6: 7, 7: 8, 8: "8a"},
    "naphthalena",
)
# Atom index (in the reference SMILES) -> canonical pyridine locant (N is
# always locant 1); verified directly via RDKit bond inspection. No fusion
# atoms, so every role value is a plain int.
_PYRIDINE_KIND = _AmplificantKind(
    "c1ccncc1",
    {3: 1, 2: 2, 1: 3, 0: 4, 5: 5, 4: 6},
    "pyridina",
)
_SENIOR_KINDS = (_NAPHTHALENE_KIND, _PYRIDINE_KIND)
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
            # A plain periphery ring atom carries 1 H if it's carbon (an
            # aromatic CH, benzene/naphthalene), or 0 if it's pyridine's
            # own nitrogen (its lone pair, not an H, completes the
            # aromatic sextet) -- the only two ring-atom kinds this
            # module's registered `_SENIOR_KINDS`/benzene ever produce at
            # an unsubstituted periphery position.
            expected_h = 1 if atom.GetAtomicNum() == 6 else 0
            if atom.GetTotalNumHs() != expected_h:
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


def _senior_attachment_locants(kind, mol, att_a, att_b, forward_atom):
    """Try every automorphism-consistent mapping of `kind`'s reference onto
    `mol` and return ((lo, hi), {att_a: label, att_b: label}) for the
    lowest achievable locant set, or None if either attachment resolves to
    a ring-fusion position (e.g. naphthalene's 4a/8a) in every matching
    automorphism. Among numbering choices tied on that lowest set, prefer
    the one that assigns the lower individual locant to `forward_atom`
    (P-26.4.1.4's citation-order rule) -- a symmetric amplificant like
    naphthalene or pyridine lets more than one automorphism reach the same
    {lo, hi} set while swapping which physical attachment atom is 'lo' vs
    'hi'."""
    matches = mol.GetSubstructMatches(kind.reference, uniquify=False, useChirality=False)
    candidates = []
    for match in matches:
        target_to_ref = {match[i]: i for i in range(len(match))}
        if att_a not in target_to_ref or att_b not in target_to_ref:
            continue
        label_a = kind.role[target_to_ref[att_a]]
        label_b = kind.role[target_to_ref[att_b]]
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


def _find_senior_kind(mol):
    """Return (kind, senior_atoms) for whichever registered `_SENIOR_KINDS`
    entry has exactly one match in `mol`, or None if none does (or more
    than one kind matches -- an ambiguous/ill-formed input, out of scope)."""
    found = None
    for kind in _SENIOR_KINDS:
        matches = mol.GetSubstructMatches(kind.reference, uniquify=True, useChirality=False)
        if len(matches) == 1:
            if found is not None:
                return None
            found = (kind, set(matches[0]))
    return found


def _find_naphthalene_benzene_phane(mol):
    found = _find_senior_kind(mol)
    if found is None:
        return None
    kind, senior_atoms = found

    ring_info = mol.GetRingInfo()
    benzene_rings = [
        set(ring)
        for ring in ring_info.AtomRings()
        if len(ring) == 6
        and not (set(ring) & senior_atoms)
        and all(mol.GetAtomWithIdx(a).GetAtomicNum() == 6 and mol.GetAtomWithIdx(a).GetIsAromatic() for a in ring)
    ]
    if len(benzene_rings) != 1:
        return None
    benzene_atoms = benzene_rings[0]
    if mol.GetNumAtoms() < len(senior_atoms) + len(benzene_atoms):
        return None

    graph = adjacency(mol)
    senior_atts = _find_attachments(mol, graph, senior_atoms)
    benzene_atts = _find_attachments(mol, graph, benzene_atoms)
    if senior_atts is None or benzene_atts is None:
        return None

    b_cyc = _ring_cyclic_order(graph, list(benzene_atoms), benzene_atts[0][0])
    b_dist = b_cyc.index(benzene_atts[1][0])
    b_dist = min(b_dist, len(b_cyc) - b_dist)
    benzene_pattern = {1: "ortho", 2: "meta", 3: "para"}.get(b_dist)
    if benzene_pattern is None:
        return None

    (sa, sa_exo), (sb, sb_exo) = senior_atts

    bridge_sa = _walk_bridge(mol, graph, sa, sa_exo, benzene_atoms)
    bridge_sb = _walk_bridge(mol, graph, sb, sb_exo, benzene_atoms)
    if bridge_sa is None or bridge_sb is None:
        return None
    len_sa, benzene_end_from_sa = bridge_sa
    len_sb, benzene_end_from_sb = bridge_sb
    if benzene_end_from_sa == benzene_end_from_sb:
        return None

    # P-26.4.1.1: the shorter bridge goes forward (gives benzene the lower
    # locant) -- decided purely from bridge length, independent of the
    # senior amplificant's own numbering. When the bridges tie, P-26.4.1.1
    # can't decide either direction, so P-26.4.1.4 breaks the tie instead:
    # whichever direction puts the *lower* attachment locant on the
    # forward (first-cited) side wins -- naphthalene's own extra symmetry
    # happens to make both directions agree in every case tested so far,
    # but pyridine's lower symmetry (#1031) exposed that this must be
    # checked explicitly, not assumed.
    if len_sb < len_sa:
        forward_atom, backward_atom = sb, sa
    elif len_sa < len_sb:
        forward_atom, backward_atom = sa, sb
    else:
        candidates = []
        for candidate_forward, candidate_backward in ((sa, sb), (sb, sa)):
            candidate_result = _senior_attachment_locants(kind, mol, sa, sb, candidate_forward)
            if candidate_result is not None:
                _, candidate_labels = candidate_result
                candidates.append((candidate_labels[candidate_forward], candidate_forward, candidate_backward))
        if not candidates:
            return None
        _, forward_atom, backward_atom = min(candidates)
    forward_len = len_sa if forward_atom == sa else len_sb
    backward_len = len_sb if forward_atom == sa else len_sa

    result = _senior_attachment_locants(kind, mol, sa, sb, forward_atom)
    if result is None:
        return None
    _, senior_labels = result

    macrocycle_size = 2 + forward_len + backward_len
    benzene_locant = 2 + forward_len
    senior_citation = f"{senior_labels[forward_atom]},{senior_labels[backward_atom]}"
    benzene_local = _LOCAL_LOCANTS[benzene_pattern]

    return (
        f"1({senior_citation})-{kind.prefix}-{benzene_locant}({benzene_local})-"
        f"benzenacyclo{numerical_term(macrocycle_size)}phane"
    )


def has_naphthalene_benzene_phane_name(mol) -> bool:
    return _find_naphthalene_benzene_phane(mol) is not None


def name_naphthalene_benzene_phane(mol) -> str:
    name = _find_naphthalene_benzene_phane(mol)
    if name is None:
        raise UnsupportedStructure(
            "this is not a supported {naphthalene,pyridine}+benzene two-bridge phane shape (see P-26)"
        )
    return name
