"""Naming of a molecule combining exactly one sulfonic acid (-SO3H) with
one or more thiol (-SH) groups on the same acyclic saturated carbon chain,
per the IUPAC 2013 Recommendations ("the Blue Book"):

- P-41/P-43 (`_seniority.py`): sulfonic acids (Table 4.4 class 4) far
  outrank hydroxy compounds/thiols (class 49), so a coexisting thiol is
  demoted to the 'sulfanyl' substituent prefix instead of its own
  '-thiol' suffix -- e.g. 'SCCS(=O)(=O)O' -> '2-sulfanylethane-1-sulfonic
  acid'. `_sulfonic_acid.py`'s own module docstring already flags this
  exact gap ("no other heteroatom anywhere in the molecule except the
  sulfonic acid group's own three oxygens"); this is the first module
  built directly on `_seniority.py` (see that module's docstring for why
  the five older pairwise modules aren't refactored to use it too).
- Otherwise mirrors `_sulfonic_acid.py`'s acyclic-chain path exactly: the
  P-14.3.4.2(a)/(b) locant-omission rules, and the -SO3H locant minimized
  before substituent-prefix locants (P-45.2). 'sulfanyl' is injected into
  the same {atom_idx -> name} map `_sulfonic_acid.py`/`_carboxylic_acid_
  amine.py` already use for halogens/a demoted amine, so it is formatted,
  alphabetized, and multiplied ('2,4-disulfanyl', etc.) by the same
  general substituent-prefix machinery, with no new logic needed there.

Scope, deliberately narrow (first concrete pairwise case built on
`_seniority.py`): a single sulfonic acid plus one or more thiols, all on
one acyclic saturated chain, with halogen substituents allowed. One
narrow *aromatic*-ring exception: `_name_phenyl_chain_sulfonic_acid_thiol`
names a sulfonic acid/thiol chain hanging off a single plain,
unsubstituted benzene ring (e.g. '3-phenyl-3-sulfanylpropane-1-sulfonic
acid', PubChem CID 57312026), mirroring `_sulfonic_acid.py`'s/
`_thiol.py`'s identical benzene-ring-substituent path.
Explicitly out of scope (raise `UnsupportedStructure`): any chain
unsaturation (ene/yne), any ring other than the single benzene-
substituent exception above, more than one sulfonic acid, a coexisting
hydroxyl/ether/other heteroatom, a specified stereocenter, a thiol
directly on the benzene ring (thiophenol-type), and any sulfonic
acid/thiol not captured by a single longest chain.
"""

from rdkit import Chem

from ._coexisting_groups import name_via_senior_acyclic
from ._common import (
    HALOGEN_PREFIXES,
    UnsupportedStructure,
    adjacency,
    halogen_substituents,
    is_plain_benzene_ring,
    longest_branched_chain_through,
    lowest_locant_set,
    non_single_bonds,
    ring_chain_attachment,
    specified_stereocenters,
)
from ._numerals import alkane_name
from ._seniority import senior_class
from ._substituents import alpha_sort_key, format_substituent_prefixes, name_branch
from ._sulfonic_acid import _name_acyclic_sulfonic_acid

_ALLOWED_ATOMIC_NUMS = {6, 8, 16, *HALOGEN_PREFIXES}

# Decided once, at import time, by consulting the shared rank table rather
# than hardcoding the conclusion in prose the way the older pairwise
# modules do -- this module's whole shape assumes sulfonic acid wins, so a
# future rank-table edit that somehow reversed this would need this module
# rewritten anyway, but the assertion makes that discrepancy loud instead
# of silently producing a wrong name.
assert senior_class("sulfonic_acid", "alcohol") == "sulfonic_acid"


def _sulfonic_sulfur_atoms(mol):
    """Sulfur atoms shaped like a sulfonic acid group: bonded to exactly
    one carbon, two double-bonded (terminal) oxygens, and one
    single-bonded hydroxyl oxygen (terminal, one H). Mirrors
    `_sulfonic_acid.py`'s identical helper."""
    matches = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 16 or atom.GetDegree() != 4:
            continue
        neighbors = atom.GetNeighbors()
        carbons = [n for n in neighbors if n.GetAtomicNum() == 6]
        oxygens = [n for n in neighbors if n.GetAtomicNum() == 8]
        if len(carbons) != 1 or len(oxygens) != 3:
            continue
        double_os = [
            o
            for o in oxygens
            if mol.GetBondBetweenAtoms(atom.GetIdx(), o.GetIdx()).GetBondTypeAsDouble() == 2.0
        ]
        hydroxyl_os = [
            o
            for o in oxygens
            if mol.GetBondBetweenAtoms(atom.GetIdx(), o.GetIdx()).GetBondTypeAsDouble() == 1.0
        ]
        if len(double_os) != 2 or len(hydroxyl_os) != 1:
            continue
        if any(o.GetDegree() != 1 for o in double_os):
            continue
        (hydroxyl_o,) = hydroxyl_os
        if hydroxyl_o.GetDegree() != 1 or hydroxyl_o.GetTotalNumHs() != 1:
            continue
        matches.append(atom)
    return matches


def _thiol_sulfur_atoms(mol):
    """Sulfur atoms shaped like a plain thiol group (-SH): degree 1, one
    hydrogen, single-bonded to a carbon. Mirrors `_thiol.py`'s identical
    shape check."""
    matches = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 16 or atom.GetDegree() != 1 or atom.GetTotalNumHs() != 1:
            continue
        (bond,) = atom.GetBonds()
        if bond.GetBondTypeAsDouble() != 1.0:
            continue
        (neighbor,) = atom.GetNeighbors()
        if neighbor.GetAtomicNum() == 6:
            matches.append(atom)
    return matches


def has_sulfonic_acid_thiol_shape(mol) -> bool:
    return bool(_sulfonic_sulfur_atoms(mol)) and bool(_thiol_sulfur_atoms(mol))


def _validate_and_collect(mol, aromatic_ring_atoms=frozenset()):
    sulfonic_atoms = _sulfonic_sulfur_atoms(mol)
    if len(sulfonic_atoms) != 1:
        raise UnsupportedStructure(
            "exactly one sulfonic acid is required; this module only "
            "handles a single sulfonic acid combined with one or more "
            "thiols"
        )
    (sulfonic_sulfur,) = sulfonic_atoms
    (so3h_carbon,) = (n for n in sulfonic_sulfur.GetNeighbors() if n.GetAtomicNum() == 6)
    sulfonic_oxygens = {n.GetIdx() for n in sulfonic_sulfur.GetNeighbors() if n.GetAtomicNum() == 8}

    thiols = _thiol_sulfur_atoms(mol)
    if not thiols:
        raise UnsupportedStructure(
            "no thiol found; this module only handles a sulfonic acid "
            "combined with at least one thiol (see _sulfonic_acid.py for "
            "a plain sulfonic acid)"
        )
    thiol_idxs = {s.GetIdx() for s in thiols}

    accounted_sulfur_idxs = {sulfonic_sulfur.GetIdx()} | thiol_idxs

    has_carbon = False
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atomic_num not in _ALLOWED_ATOMIC_NUMS:
            raise UnsupportedStructure(
                "heteroatoms other than the sulfonic acid group's own "
                "oxygens (P-65.3.1), a plain thiol (-SH, P-41 Table 4.4), "
                "and halogen substituents (P-35.2.1) are not supported yet"
            )
        if atom.GetFormalCharge() != 0 or atom.GetIsotope() != 0:
            raise UnsupportedStructure("charged or isotopically modified atoms are not supported yet")
        if atomic_num == 6:
            has_carbon = True
            if atom.GetIsAromatic() and atom.GetIdx() not in aromatic_ring_atoms:
                raise UnsupportedStructure(
                    "aromatic rings are out of scope for this module (see "
                    "the separate aromatic-ring module)"
                )
        elif atomic_num == 8:
            if atom.GetIdx() not in sulfonic_oxygens:
                raise UnsupportedStructure(
                    "an oxygen that isn't part of the single sulfonic acid "
                    "group is out of scope for this module (e.g. a "
                    "coexisting hydroxyl, ether, or carbonyl)"
                )
        elif atomic_num == 16:
            if atom.GetIdx() not in accounted_sulfur_idxs:
                raise UnsupportedStructure(
                    "a sulfur atom not shaped like the sulfonic acid group "
                    "or a plain thiol is out of scope for this module"
                )
        else:
            if atom.GetDegree() != 1:
                raise UnsupportedStructure(
                    "a halogen atom must be a monovalent substituent (P-35.2.1)"
                )
    if not has_carbon:
        raise UnsupportedStructure(
            "a structure with no carbon atom has no hydrocarbon parent "
            "hydride to substitute"
        )
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")

    return sulfonic_sulfur.GetIdx(), so3h_carbon.GetIdx(), thiol_idxs


def _group(substituents):
    grouped = {}
    for position, entries in substituents.items():
        for name, is_compound in entries:
            info = grouped.setdefault(name, {"locants": [], "compound": is_compound})
            info["locants"].append(position)
    return grouped


def _name_from_substituents(chain_length, so3h_locant, grouped):
    total_subs = sum(len(info["locants"]) for info in grouped.values())

    if chain_length == 1:
        # P-14.3.4.2(a): a mononuclear parent's locant is always '1' and
        # never cited.
        return format_substituent_prefixes(grouped, omit_locants=True) + alkane_name(1) + "sulfonic acid"
    if chain_length == 2 and total_subs == 0:
        # P-14.3.4.2(b): a homogeneous two-carbon chain with exactly one
        # substituent (the sole -SO3H) in total omits the locant.
        return alkane_name(2) + "sulfonic acid"

    prefix = format_substituent_prefixes(grouped)
    return f"{prefix}{alkane_name(chain_length)}-{so3h_locant}-sulfonic acid"


def _candidate_key(chain_length, so3h_locant, substituents):
    grouped = _group(substituents)
    locant_set = lowest_locant_set(loc for info in grouped.values() for loc in info["locants"])
    citation_locants = tuple(
        loc
        for name in sorted(grouped, key=alpha_sort_key)
        for loc in sorted(grouped[name]["locants"])
    )
    name = _name_from_substituents(chain_length, so3h_locant, grouped)
    return (so3h_locant, locant_set, citation_locants, name), name


def _substituents_for_chain(graph, chain, names, excluded):
    chain_set = set(chain)
    substituents = {}
    for position, atom in enumerate(chain, start=1):
        branch_roots = [n for n in graph[atom] if n not in chain_set and n not in excluded]
        if branch_roots:
            substituents[position] = [name_branch(graph, root, atom, names) for root in branch_roots]
    return substituents


def _name_phenyl_chain_sulfonic_acid_thiol(mol, ring_atoms):
    """Name a sulfonic acid plus one or more thiols, all lying on a single
    unbranched chain hanging off one atom of an otherwise-plain,
    unsubstituted benzene ring -- e.g.
    '3-phenyl-3-sulfanylpropane-1-sulfonic acid' (PubChem CID 57312026).
    The ring is cited as a 'phenyl' substituent prefix (via
    `name_branch`'s aromatic-ring recognition) on the chain, which is the
    parent hydride, mirroring `_sulfonic_acid.py`'s
    `_name_phenyl_chain_sulfonic_acid`/`_thiol.py`'s
    `_name_phenyl_chain_thiol`. Narrower than the acyclic path above: no
    chain unsaturation and no specified stereocenter."""
    sulfonic_sulfur_idx, so3h_carbon, thiol_idxs = _validate_and_collect(
        mol, aromatic_ring_atoms=ring_atoms
    )
    if so3h_carbon in ring_atoms:
        raise UnsupportedStructure(
            "a sulfonic acid directly attached to the benzene ring "
            "(benzenesulfonic acid-style naming) is out of scope for this "
            "module (see the separate aromatic-ring module)"
        )
    if specified_stereocenters(mol):
        raise UnsupportedStructure(
            "a specified stereocenter alongside a benzene-ring-substituent "
            "sulfonic acid/thiol chain is not supported yet"
        )
    non_ring_unsaturation = [
        b
        for b in non_single_bonds(mol)
        if sulfonic_sulfur_idx not in (b[0], b[1])
        and b[0] not in ring_atoms
        and b[1] not in ring_atoms
    ]
    if non_ring_unsaturation:
        raise UnsupportedStructure(
            "chain unsaturation alongside a benzene-ring-substituent "
            "sulfonic acid/thiol chain is not supported yet"
        )

    graph = adjacency(mol)
    attachment = ring_chain_attachment(graph, ring_atoms, set())
    if attachment is None:
        raise UnsupportedStructure(
            "a benzene ring with more than one exocyclic substituent "
            "alongside a chain sulfonic acid/thiol is not supported yet"
        )
    ring_atom, chain_root = attachment
    if chain_root in thiol_idxs:
        raise UnsupportedStructure(
            "a thiol directly on the benzene ring (thiophenol-type) uses "
            "a separate construction, out of scope for this chain-parent "
            "module"
        )
    chain, _ = longest_branched_chain_through(graph, so3h_carbon, ring_atoms, {sulfonic_sulfur_idx} | thiol_idxs)
    chain_set = set(chain)
    for s in thiol_idxs:
        (carbon,) = graph[s]
        if carbon not in chain_set:
            raise UnsupportedStructure(
                "a thiol outside the single unbranched chain hanging off "
                "the benzene ring is not supported yet"
            )

    names = {**halogen_substituents(mol), **{s: "sulfanyl" for s in thiol_idxs}}
    chain_length = len(chain)
    excluded = {sulfonic_sulfur_idx, ring_atom}
    best_key = None
    best_name = None
    for candidate in (chain, list(reversed(chain))):
        position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
        so3h_locant = position_of[so3h_carbon]
        substituents = _substituents_for_chain(graph, candidate, names, excluded)
        ring_entry = name_branch(graph, ring_atom, chain_root, names, ring_atoms)
        substituents.setdefault(position_of[chain_root], []).append(ring_entry)
        key, name = _candidate_key(chain_length, so3h_locant, substituents)
        if best_key is None or key < best_key:
            best_key, best_name = key, name
    return best_name


def name_sulfonic_acid_thiol(mol) -> str:
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() == 1:
        ring_atoms = set(ring_info.AtomRings()[0])
        if is_plain_benzene_ring(mol, ring_atoms):
            return _name_phenyl_chain_sulfonic_acid_thiol(mol, ring_atoms)

    sulfonic_sulfur_idx, so3h_carbon, thiol_idxs = _validate_and_collect(mol)

    if mol.GetRingInfo().NumRings() != 0:
        raise UnsupportedStructure(
            "a sulfonic acid/thiol combination on a ring is out of scope "
            "for this acyclic-only module"
        )
    all_non_single = non_single_bonds(mol)
    if any(a != sulfonic_sulfur_idx and b != sulfonic_sulfur_idx for a, b, _ in all_non_single):
        raise UnsupportedStructure(
            "chain unsaturation (ene/yne) alongside a sulfonic acid/thiol "
            "combination is out of scope for this module"
        )

    graph = adjacency(mol)
    thiol_carbons = {next(iter(graph[s])) for s in thiol_idxs}
    return name_via_senior_acyclic(
        _name_acyclic_sulfonic_acid,
        "sulfonic_acid",
        "alcohol",
        (mol, sulfonic_sulfur_idx, so3h_carbon, []),
        {s: "sulfanyl" for s in thiol_idxs},
        required_atoms=thiol_carbons,
    )
