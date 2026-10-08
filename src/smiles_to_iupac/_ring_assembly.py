"""Naming of the ring assembly of two identical cyclic systems joined by a
single bond, sharing no atom, per the IUPAC 2013 Recommendations ("the Blue
Book"):

- P-28.2.1 (Chapter P-2, https://iupac.qmul.ac.uk/BlueBook/PDF/P2.pdf): two
  identical cyclic parent hydrides joined directly by a single bond are named
  as a ring assembly using the parent hydride's name preceded by 'bi',
  enclosed in parentheses "if necessary" to avoid confusion with a von Baeyer
  name -- confirmed PIN worked examples `1,1'-bi(cyclopropane)` (parens
  needed: 'bicyclopropane' would misread as a von Baeyer bicyclic name),
  `2,2'-bipyridine`, `1,2'-binaphthalene`, `2,3'-bifuran` (no parens needed
  for any of these -- none start with 'cyclo'), the Blue Book. Two benzene rings joined this way instead use the retained
  substituent-group name 'biphenyl' (method (2) of the same rule), giving the
  PIN '1,1'-biphenyl' -- the sole named exception to the plain parent-hydride
  construction. The locants of the two ring-joining atoms are compulsorily
  cited even when the compound is otherwise unsubstituted.
- Numbering: each ring is numbered independently, one with unprimed locants,
  the other with primed locants, and "lowest possible locants must be used
  to denote the positions of attachment" (same P-28.2.1 text
  `_ring_assembly_chain.py` already cites for the 3-6-ring case). For a
  symmetric ring (benzo/cycloalkane) this is trivially satisfied by starting
  each ring's own numbering at its attachment atom (locant 1 always), then
  choosing the numbering direction to minimize substituent locants as a
  tie-break. For a heteroaromatic parent whose own numbering is fixed by its
  role sequence (any `_ring_assembly_chain.py`'s `_NON_NH_ROLE_SEQUENCES`
  parent, e.g. pyridine's nitrogen always at locant 1), the attachment point
  isn't automatically locant 1 -- among the ring's own valid role-sequence
  alignments (`_hetero_ring_alignments`, reused unchanged from that module),
  the one giving the lowest locant to the attachment atom is chosen first,
  confirmed by the worked example `2,2'-bipyridine` (attachment at the lower
  of pyridine's two symmetry-equivalent non-nitrogen alpha positions, 2 not
  6).
- P-35.2.1 (Chapter P-3): halogen substituents hang off a ring atom the same
  way as in every other ring module.

Scope: two *identical* rings connected by exactly one single (non-aromatic)
bond -- 6-membered all-carbon aromatic (benzo), monocyclic saturated
all-carbon (any one ring size), any mancude 5- or 6-membered ring matching a
`_NON_NH_ROLE_SEQUENCES` parent (pyridine, furan, thiophene, selenophene,
tellurophene, pyridazine, pyrimidine, pyrazine -- same exclusions as that
module: no N-H tautomer-ambiguous parent, no locant-prefixed Hantzsch-Widman
parent), 1H-pyrrole (see below), or 1H-imidazole/1H-pyrazole joined through
both rings' own N-H nitrogen only (see below -- a carbon-attached junction
for either of these two stays out of scope), each ring optionally bearing
halogen substituents. Three to six identical rings in an unbranched chain
(terphenyl etc.) are `_ring_assembly_chain.py`'s own job instead -- a
separate composite-locant numbering scheme (P-28.3), not a generalization of
this module's own primed-locant one.

- P-28.2.3 (indicated hydrogen of a two-component ring assembly, same
  the Blue Book): 1H-pyrrole is the one N-H tautomer-
  unambiguous `_ROLE_SEQUENCES` parent (a single ring nitrogen, no
  prototropic choice about which atom is "N1" the way imidazole/pyrazole
  have -- those stay deferred, see #936's own investigation) whose own
  indicated-hydrogen position can be a ring-assembly junction. Per
  P-28.2.3's own text, indicated hydrogen is "ignoring the indicated
  hydrogen atoms of the component rings... maximum number of noncumulative
  double bonds is then added taking into account the junction
  positions... remaining saturated ring positions are designated as
  indicated hydrogen, placed... at the front of the name of the assembly"
  -- for pyrrole specifically this reduces to a simple per-ring check
  (only one candidate saturated position, the N, exists at all): if a
  ring's own junction bond sits at the N itself, that position is already
  occupied and needs no indicated-H citation for that ring (confirmed PIN
  worked example `1,1'-bipyrrole`, "no indicated hydrogen needed"); if the
  junction is elsewhere (through a ring carbon instead), the N keeps its
  own H and that ring's `1`/`1'` locant is cited together with the other
  ring's own (if it also needs one) at the very front of the whole
  assembly name, before any substituent prefix (confirmed via a real
  PubChem structure, CID 260036, canonical SMILES `C1=CNC(=C1)C2=CC=CN2`,
  PIN `1H,1'H-2,2'-bipyrrole`). `_ring_assembly_chain.py`'s own 3-6-ring
  case has the analogous mechanism (P-28.3.1's composite-locant citation
  format) for its own pyrrole support, per its own docstring.
- 1H-imidazole/1H-pyrazole: unlike pyrrole (a single ring nitrogen, no
  tautomer choice), these have two ring nitrogens and a real prototropic-
  tautomer ambiguity about which one is "N1" -- `_hetero_monocyclic.py`'s
  own `_TAUTOMER_AMBIGUOUS_UNLESS_N1` resolves the analogous single-ring
  substituent case by only trusting a substituent that sits directly on
  the N-H position itself. This module trusts the same "safe" condition
  for a ring-assembly junction: when the junction bond itself replaces a
  ring nitrogen's H (`_tautomer_fixed_ring_kind`), that nitrogen *is* N1
  by construction -- no remaining tautomer choice, and (mirroring
  pyrrole's own N-N-attached case) no indicated hydrogen needed for
  either ring, since the one candidate saturated position on each ring is
  already occupied by the junction. Confirmed PIN worked examples
  `1,1'-biimidazole` (PubChem CID 15034216) and `1,1'-bipyrazole` (CID
  21981271). A carbon-attached junction (e.g. `2,2'-biimidazole`, CID
  101463) is deferred -- no confirmed primary-source worked example for
  whether C2 (flanked by both ring nitrogens) is actually safe the same
  way N1 is, or needs the same conservative rejection the single-ring
  case already applies to imidazole's C2-substituted case.
"""

from rdkit import Chem

from ._common import (
    adjacency,
    group_substituents,
    halogen_substituents,
    ring_cycle,
    substituent_locant_set_and_citation,
    validate_atoms_and_bonds,
)
from ._numerals import alkane_name, numerical_term
from ._ring_assembly_chain import (
    _NON_NH_ROLE_SEQUENCES,
    _hetero_ring_alignments,
    _pyrrole_ring_kind,
    _ring_kind,
    hydro_locants,
    hydro_prefix,
    hydro_sort_key,
    saturated_counterpart_kind,
    validate_hetero_ring_assembly_atoms,
)
from ._substituents import format_substituent_prefixes, name_branch

# 1H-pyrrole/1H-imidazole/1H-pyrazole's own P-28.2.3 indicated-hydrogen
# citation (see module docstring) needs the bare parent-hydride name
# ("pyrrole"/"imidazole"/"pyrazole", not the "1H-"-prefixed form) for the
# "bi"-prefixed ring word -- unlike every `_NON_NH_ROLE_SEQUENCES` parent,
# whose own name never carries an indicated-hydrogen prefix to strip in
# the first place.
_INDICATED_HYDROGEN_BARE_NAMES = {
    "1H-pyrrole": "pyrrole",
    "1H-imidazole": "imidazole",
    "1H-pyrazole": "pyrazole",
}


def _bond_between(bond, atoms_a, atoms_b):
    x, y = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
    return (x in atoms_a and y in atoms_b) or (x in atoms_b and y in atoms_a)


def _tautomer_fixed_ring_kind(mol, graph, ring, attach_atom):
    """("1H-imidazole", 5) or ("1H-pyrazole", 5) if `ring` matches one of
    these two N-H tautomer-ambiguous parents' element pattern
    (`_hetero_ring_alignments`, element-only match) and `attach_atom` is
    one of the ring's own two nitrogens, else None. This module only
    trusts a ring-assembly junction landing on a ring nitrogen for these
    two parents (see module docstring): the junction bond itself replaces
    that nitrogen's H, fixing it as N1 with no remaining prototropic-
    tautomer choice -- structurally guaranteed here since an unsubstituted
    imidazole/pyrazole ring's own reflective symmetry lets *either* ring
    nitrogen be numbered role-sequence position 1 (`_hetero_ring_alignments`
    finds a matching alignment for the physically-substituted one
    specifically because of that symmetry, not because both are
    interchangeable at once). A carbon-attached junction is deferred
    (`attach_atom` never matches a role-sequence N-position, so no
    alignment exists and this returns None, same as an outright
    non-match)."""
    if mol.GetAtomWithIdx(attach_atom).GetAtomicNum() != 7:
        return None
    for parent_name in ("1H-imidazole", "1H-pyrazole"):
        if any(True for _ in _hetero_ring_alignments(mol, graph, ring, parent_name)):
            return parent_name, 5
    return None


_REPLACEMENT_ELEMENTS = {8: "oxa", 16: "thia", 34: "selena", 52: "tellura", 7: "aza"}
_REPLACEMENT_ORDER = (8, 16, 34, 52, 7)
_REPLACEMENT_MINIMUM_SIZE = 11


def _replacement_ring_kind(mol, ring, attach_atom):
    """("saturated", n) for a saturated monocycle of more than ten members whose heteroatoms are O, S, Se, Te or NH and
    whose junction atom is carbon: the cycloalkane assembly takes the skeletal replacement prefixes (P-28.4.2)."""
    atoms = [mol.GetAtomWithIdx(i) for i in ring]
    if len(ring) < _REPLACEMENT_MINIMUM_SIZE or mol.GetAtomWithIdx(attach_atom).GetAtomicNum() != 6:
        return None
    if any(a.GetIsAromatic() or a.GetAtomicNum() not in {6, *_REPLACEMENT_ELEMENTS} for a in atoms):
        return None
    if any(a.GetAtomicNum() == 7 and (a.GetFormalCharge() or a.GetDegree() != 2) for a in atoms):
        return None
    ring_set = set(ring)
    if any(
        b.GetBondTypeAsDouble() != 1.0 for b in mol.GetBonds() if b.GetBeginAtomIdx() in ring_set and b.GetEndAtomIdx() in ring_set
    ):
        return None
    return "saturated", len(ring)


def find_ring_assembly_core(mol):
    """Return (ring0_atoms, ring1_atoms, attach0, attach1, kind) if `mol` is
    exactly two disjoint identical-kind rings (`_ring_kind`'s own
    ("aromatic", 6), ("saturated", n), a `_NON_NH_ROLE_SEQUENCES` parent
    name, `_pyrrole_ring_kind`'s own ("1H-pyrrole", 5), or
    `_tautomer_fixed_ring_kind`'s own ("1H-imidazole", 5)/
    ("1H-pyrazole", 5)) joined by one single bond, else None."""
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    if len(atom_rings) != 2:
        return None

    ring0, ring1 = set(atom_rings[0]), set(atom_rings[1])
    if ring0 & ring1:
        return None

    connecting = [bond for bond in mol.GetBonds() if _bond_between(bond, ring0, ring1)]
    if len(connecting) != 1:
        return None
    bond = connecting[0]
    if bond.GetIsAromatic() or bond.GetBondTypeAsDouble() != 1.0:
        return None

    x, y = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
    attach0, attach1 = (x, y) if x in ring0 else (y, x)

    # `_tautomer_fixed_ring_kind` needs the junction atom (unlike every
    # other ring kind here), so it's checked per-ring only after the
    # junction itself is known, not folded into a single kind-per-ring set
    # the way the junction-agnostic kinds are.
    graph = adjacency(mol)
    kind0 = (
        _ring_kind(mol, atom_rings[0])
        or _pyrrole_ring_kind(mol, atom_rings[0])
        or _tautomer_fixed_ring_kind(mol, graph, atom_rings[0], attach0)
        or _replacement_ring_kind(mol, atom_rings[0], attach0)
    )
    kind1 = (
        _ring_kind(mol, atom_rings[1])
        or _pyrrole_ring_kind(mol, atom_rings[1])
        or _tautomer_fixed_ring_kind(mol, graph, atom_rings[1], attach1)
        or _replacement_ring_kind(mol, atom_rings[1], attach1)
    )
    hydro_ring = None
    if kind0 is None and kind1 is not None and kind1[0] in _NON_NH_ROLE_SEQUENCES:
        kind0, hydro_ring = saturated_counterpart_kind(mol, graph, atom_rings[0]), 0
    elif kind1 is None and kind0 is not None and kind0[0] in _NON_NH_ROLE_SEQUENCES:
        kind1, hydro_ring = saturated_counterpart_kind(mol, graph, atom_rings[1]), 1
    if kind0 is None or kind0 != kind1:
        return None

    return atom_rings[0], atom_rings[1], attach0, attach1, kind0, hydro_ring


def _numberings_from_attachment(graph, ring_atoms, attach, prime):
    cycle = ring_cycle(graph, list(ring_atoms))
    start = cycle.index(attach)
    rotated = cycle[start:] + cycle[:start]
    suffix = "'" if prime else ""
    for seq in (rotated, [rotated[0]] + list(reversed(rotated[1:]))):
        yield {atom: f"{position}{suffix}" for position, atom in enumerate(seq, start=1)}


def _hetero_numberings_from_attachment(mol, graph, ring_atoms, parent_name, prime):
    suffix = "'" if prime else ""
    for alignment in _hetero_ring_alignments(mol, graph, ring_atoms, parent_name):
        yield {atom: f"{position}{suffix}" for atom, position in alignment.items()}


def _carbon_skeleton(mol):
    """`mol` with every ring atom that a skeletal replacement prefix will name turned into carbon."""
    editable = Chem.RWMol(mol)
    for atom in editable.GetAtoms():
        if atom.GetAtomicNum() in _REPLACEMENT_ELEMENTS and atom.IsInRing():
            atom.SetAtomicNum(6)
            atom.SetNoImplicit(False)
            atom.SetNumExplicitHs(0)
    skeleton = editable.GetMol()
    skeleton.UpdatePropertyCache(strict=False)
    return skeleton


def _replacement_text(mol, locants):
    """(sort key, 'a' prefix text) of the skeletal replacement prefixes of an assembly: low locants to the heteroatoms
    as a set, then in the order O > S > Se > Te > N (P-28.4.2)."""
    by_element = {}
    for atom, position in locants.items():
        z = mol.GetAtomWithIdx(atom).GetAtomicNum()
        if z in _REPLACEMENT_ELEMENTS and mol.GetAtomWithIdx(atom).IsInRing():
            by_element.setdefault(z, []).append(position)
    if not by_element:
        return ((), ()), ""
    order = lambda position: hydro_sort_key(position)
    as_set = tuple(sorted(order(p) for ps in by_element.values() for p in ps))
    by_seniority = tuple(tuple(sorted(order(p) for p in by_element.get(z, ()))) for z in _REPLACEMENT_ORDER)
    pieces = []
    for z in _REPLACEMENT_ORDER:
        if z in by_element:
            cited = sorted(by_element[z], key=order)
            multiplier = numerical_term(len(cited)) if len(cited) > 1 else ""
            pieces.append(f"{','.join(cited)}-{multiplier}{_REPLACEMENT_ELEMENTS[z]}")
    return (as_set, by_seniority), "-".join(pieces) + "-"


def _candidate_key(
    locants,
    ring_atoms,
    graph,
    halogens,
    attach_a,
    attach_b,
    ring_word,
    mol=None,
    indicated_hydrogen_prefix="",
    hydro=(),
):
    substituents = {}
    for atom, position in locants.items():
        branch_roots = [n for n in graph[atom] if n not in ring_atoms]
        if branch_roots:
            substituents[position] = [name_branch(graph, root, atom, halogens, mol=mol) for root in branch_roots]

    grouped = group_substituents(substituents)
    locant_set, _, citation_locants = substituent_locant_set_and_citation(grouped)
    prefix = format_substituent_prefixes(grouped) if grouped else ""
    if prefix:
        prefix += "-"
    # P-28.2.1's "lowest possible locants must be used to denote the
    # positions of attachment" governs the choice among candidates before
    # substituent locants do -- for benzo/cycloalkane this pair is always
    # ("1", "1'") since `_numberings_from_attachment` starts at the
    # attachment atom, so including it here doesn't change that branch's
    # own behavior; for a heteroaromatic parent whose numbering is fixed by
    # role sequence instead, it's the deciding factor (confirmed by
    # `2,2'-bipyridine`, not `6,2'-bipyridine` or `6,6'-bipyridine`).
    attach_pair = tuple(sorted((locants[attach_a], locants[attach_b])))
    replaced = _replacement_text(mol, locants) if ring_word.startswith("bi(cyclo") and len(ring_atoms) // 2 >= _REPLACEMENT_MINIMUM_SIZE else ((), "")
    # P-28.2.3: indicated hydrogen, if any, is "placed... at the front of
    # the name of the assembly" -- ahead of the substituent prefix too,
    # not folded next to the ring word the way a single ring's own "nH-"
    # sits (see module docstring).
    name = f"{indicated_hydrogen_prefix}{prefix}{replaced[1]}{hydro_prefix(hydro)}{attach_pair[0]},{attach_pair[1]}-{ring_word}"
    return attach_pair, replaced[0], tuple(hydro_sort_key(p) for p in hydro), locant_set, citation_locants, name


def name_ring_assembly(mol, core) -> str:
    ring0_atoms, ring1_atoms, attach0, attach1, ring_kind, hydro_ring = core
    parent_name, ring_size = ring_kind

    ring_atoms = set(ring0_atoms) | set(ring1_atoms)

    if parent_name in _NON_NH_ROLE_SEQUENCES or parent_name in _INDICATED_HYDROGEN_BARE_NAMES:
        validate_hetero_ring_assembly_atoms(mol, ring_atoms, parent_name)
        ring_word = "bi" + _INDICATED_HYDROGEN_BARE_NAMES.get(parent_name, parent_name)
        numberings = lambda graph, ring_atoms_i, attach_i, prime: _hetero_numberings_from_attachment(
            mol, graph, ring_atoms_i, parent_name, prime
        )
    else:
        validate_atoms_and_bonds(_carbon_skeleton(mol))
        if parent_name == "aromatic":
            ring_word = "biphenyl"
        else:
            ring_word = f"bi(cyclo{alkane_name(ring_size)})"
        numberings = _numberings_from_attachment

    graph = adjacency(mol)
    halogens = halogen_substituents(mol)

    best_key = None
    best_name = None
    for (ring_a, attach_a), (ring_b, attach_b) in (
        ((ring0_atoms, attach0), (ring1_atoms, attach1)),
        ((ring1_atoms, attach1), (ring0_atoms, attach0)),
    ):
        indicated_hydrogen_prefix = ""
        if parent_name in _INDICATED_HYDROGEN_BARE_NAMES:
            # The ring's own N is always role-sequence position 1 (see
            # `_pyrrole_ring_kind`'s docstring), regardless of which
            # numbering direction wins below, so whether indicated
            # hydrogen is needed depends only on whether the junction atom
            # itself is that N -- not on any per-candidate numbering
            # choice.
            needed = [
                f"1{suffix}H" for attach, suffix in ((attach_a, ""), (attach_b, "'")) if mol.GetAtomWithIdx(attach).GetAtomicNum() != 7
            ]
            if needed:
                indicated_hydrogen_prefix = ",".join(needed) + "-"
        for locants_a in numberings(graph, ring_a, attach_a, False):
            for locants_b in numberings(graph, ring_b, attach_b, True):
                locants = {**locants_a, **locants_b}
                hydro = ()
                if hydro_ring is not None:
                    saturated = (ring0_atoms, ring1_atoms)[hydro_ring]
                    hydro = hydro_locants(locants[atom] for atom in saturated)
                key = _candidate_key(
                    locants,
                    ring_atoms,
                    graph,
                    halogens,
                    attach_a,
                    attach_b,
                    ring_word,
                    mol=mol,
                    indicated_hydrogen_prefix=indicated_hydrogen_prefix,
                    hydro=hydro,
                )
                if best_key is None or key < best_key:
                    best_key, best_name = key, key[-1]

    return best_name
