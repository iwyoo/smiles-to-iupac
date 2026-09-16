"""Naming of a molecule combining exactly one terminal aldehyde (-CHO) with
one or more internal ketone (=O) carbonyls on the same acyclic saturated
carbon chain, per the IUPAC 2013 Recommendations ("the Blue Book"):

- P-41, Table 3.3 (Chapter P-4, https://iupac.qmul.ac.uk/BlueBook/PDF/P4.pdf
  for P-41; Table 3.3 lives in Chapter P-3): when two different
  characteristic groups that could each be cited as a suffix are both
  present, only the more senior one is; the rest are cited as prefixes.
  Table 3.3 ranks 'al' (`_aldehyde.py`) senior to 'one' (`_ketone.py`), so
  a coexisting ketone is demoted to the 'oxo' substituent prefix instead of
  its own '-one' suffix, e.g. 'CH3CH2CH2COCH2CHO' -> '3-oxohexanal'
  (confirmed against this worked example). This module is the first
  concrete case of that general suffix-vs-suffix demotion; every other
  parent-hydride module here still rejects a coexisting carbonyl outright
  (see e.g. `_aldehyde.py`'s own module docstring) until demotion is
  implemented for that specific pair too.
- The 'oxo' prefix is injected into the same {atom_idx -> name} map
  `_aldehyde.py` already uses for a demoted hydroxyl ('hydroxy'), so it is
  formatted, alphabetized, and multiplied ('3,5-dioxo', etc.) by the same
  general substituent-prefix machinery, with no new logic needed there.
- Otherwise mirrors `_aldehyde.py` exactly: the aldehyde carbon is always a
  chain terminus and always becomes C1 with its own locant never cited
  (P-14.3.3); a ketone carbon's locant, by contrast, is always cited
  (P-14.3.3 doesn't apply to a prefix substituent the way it does to the
  principal suffix).

Scope, deliberately narrow (first concrete case of the general
suffix-vs-suffix demotion problem): only a single terminal aldehyde plus
one or more ketones, all on one
acyclic saturated chain, with halogen substituents allowed. Explicitly out
of scope (raise `UnsupportedStructure`): any chain unsaturation (ene/yne),
any ring, more than one aldehyde, a coexisting hydroxyl/ether/other
heteroatom, and any carbonyl not captured by a single longest chain.
"""

from rdkit import Chem

from ._common import (
    HALOGEN_PREFIXES,
    UnsupportedStructure,
    adjacency,
    carbon_adjacency,
    group_substituents,
    halogen_substituents,
    is_plain_benzene_ring,
    longest_branched_chain,
    longest_chains,
    name_from_substituents,
    non_single_bonds,
    ring_chain_attachment,
    substituent_locant_set_and_citation,
)
from ._substituents import format_substituent_prefixes, name_branch, substituents_for_chain

_ALLOWED_ATOMIC_NUMS = {6, 8, *HALOGEN_PREFIXES}


def _carbonyl_oxygens_by_shape(mol):
    """(aldehyde_oxygens, ketone_oxygens): every isolated, doubly-bonded
    carbonyl oxygen's atom index, split by whether its carbon has one
    (aldehyde-shaped) or two (ketone-shaped) carbon neighbors -- permissive
    about everything else (see module docstring's find/name split)."""
    aldehydes, ketones = set(), set()
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 8 or atom.GetDegree() != 1:
            continue
        (bond,) = atom.GetBonds()
        if bond.GetBondTypeAsDouble() != 2.0:
            continue
        (carbon,) = atom.GetNeighbors()
        if carbon.GetAtomicNum() != 6:
            continue
        carbon_neighbors = [n for n in carbon.GetNeighbors() if n.GetAtomicNum() == 6]
        if len(carbon_neighbors) == 1:
            aldehydes.add(atom.GetIdx())
        elif len(carbon_neighbors) == 2:
            ketones.add(atom.GetIdx())
    return aldehydes, ketones


def has_aldehyde_ketone_shape(mol) -> bool:
    aldehydes, ketones = _carbonyl_oxygens_by_shape(mol)
    return bool(aldehydes) and bool(ketones)


def _validate_and_collect(mol, aromatic_ring_atoms=frozenset()):
    """`aromatic_ring_atoms`: atom indices already independently verified
    (by the caller, before this function runs) to form a single plain
    benzene ring with exactly one exocyclic attachment -- exempted from
    the aromatic-atom rejection below so `name_aldehyde_ketone`'s
    benzene-ring-substituent path (see
    `_name_phenyl_chain_aldehyde_ketone`) can reuse this same validation
    for the rest of the molecule. Empty by default, so every other
    caller's behavior is unchanged."""
    has_carbon = False
    aldehydes = set()
    ketones = set()
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atomic_num not in _ALLOWED_ATOMIC_NUMS:
            raise UnsupportedStructure(
                "heteroatoms other than aldehyde/ketone carbonyl oxygens "
                "(P-41, Table 3.3) and halogen substituents (P-35.2.1) are "
                "not supported yet"
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
            if atom.GetDegree() != 1:
                raise UnsupportedStructure(
                    "an oxygen bonded to more than one heavy atom (e.g. an "
                    "ether) is out of scope; a coexisting hydroxyl is out "
                    "of scope for this first pass too (see module docstring)"
                )
            (bond,) = atom.GetBonds()
            if bond.GetBondTypeAsDouble() != 2.0:
                raise UnsupportedStructure(
                    "an oxygen that isn't a doubly-bonded aldehyde/ketone "
                    "carbonyl is out of scope for this module"
                )
            (carbon,) = atom.GetNeighbors()
            if carbon.GetAtomicNum() != 6:
                raise UnsupportedStructure("a carbonyl oxygen must be attached to a carbon atom")
            if carbon.GetIsAromatic():
                raise UnsupportedStructure(
                    "a carbonyl on an aromatic ring is out of scope for this "
                    "module"
                )
            carbon_neighbors = [n for n in carbon.GetNeighbors() if n.GetAtomicNum() == 6]
            if len(carbon_neighbors) == 1:
                aldehydes.add(atom.GetIdx())
            elif len(carbon_neighbors) == 2:
                ketones.add(atom.GetIdx())
            else:
                raise UnsupportedStructure(
                    "a carbonyl carbon with zero or more than two carbon "
                    "neighbors is not a valid aldehyde/ketone carbon"
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
    if len(aldehydes) != 1:
        raise UnsupportedStructure(
            "exactly one aldehyde is required; this module only handles "
            "a single aldehyde combined with one or more ketones"
        )
    if not ketones:
        raise UnsupportedStructure(
            "no ketone found; this module only handles an aldehyde "
            "combined with at least one ketone (see _aldehyde.py for a "
            "plain aldehyde)"
        )
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")
    return aldehydes, ketones


def _name_from_substituents(chain_length, grouped):
    return format_substituent_prefixes(grouped) + name_from_substituents(chain_length, [], [], "al")


def _candidate_key(chain_length, al_locant, substituents):
    grouped = group_substituents(substituents)
    locant_set, _, citation_locants = substituent_locant_set_and_citation(grouped)
    name = _name_from_substituents(chain_length, grouped)
    return (al_locant, locant_set, citation_locants, name), name


def _name_phenyl_chain_aldehyde_ketone(mol, ring_atoms):
    """Name an aldehyde+ketone combination whose -CHO lies entirely on a
    single unbranched chain hanging off one atom of an otherwise-plain,
    unsubstituted benzene ring -- e.g. 3-oxo-4-phenylbutanal. The ring is
    cited as a 'phenyl' substituent prefix (via `name_branch`'s aromatic-
    ring recognition) on the chain, which is the parent hydride, mirroring
    `_aldehyde.py`'s `_name_phenyl_chain_aldehyde`. Narrower than the
    acyclic path above: no chain unsaturation."""
    aldehydes, ketones = _validate_and_collect(mol, aromatic_ring_atoms=ring_atoms)
    carbonyl_oxygens = aldehydes | ketones
    all_non_single = non_single_bonds(mol)
    carbonyl_bonds = [b for b in all_non_single if b[0] in carbonyl_oxygens or b[1] in carbonyl_oxygens]
    ring_bonds = [b for b in all_non_single if b[0] in ring_atoms and b[1] in ring_atoms]
    if len(carbonyl_bonds) + len(ring_bonds) != len(all_non_single):
        raise UnsupportedStructure(
            "chain unsaturation (ene/yne) alongside a benzene-ring-"
            "substituent aldehyde/ketone combination is out of scope for "
            "this module"
        )

    graph = adjacency(mol)
    (aldehyde_o,) = aldehydes
    (aldehyde_carbon,) = graph[aldehyde_o]
    names = {**halogen_substituents(mol), **{o: "oxo" for o in ketones}}
    attachment = ring_chain_attachment(graph, ring_atoms, set())
    if attachment is None:
        raise UnsupportedStructure(
            "a benzene ring with more than one exocyclic substituent "
            "alongside a chain aldehyde/ketone combination is not "
            "supported yet"
        )
    chain, _ = longest_branched_chain(graph, aldehyde_carbon, ring_atoms, aldehydes | set(names), halogens=halogen_substituents(mol))
    if any(graph[o][0] not in chain for o in ketones):
        raise UnsupportedStructure(
            "not every ketone-bearing carbon lies on the chain hanging "
            "off the benzene ring for this benzene-substituent path"
        )
    if len(chain) < 2:
        raise UnsupportedStructure(
            "a -CHO group directly attached to the benzene ring (the "
            "'carbaldehyde' suffix, P-33.3.1.2) uses a separate naming "
            "construction, out of scope for this acyclic-chain-parent "
            "module"
        )

    chain_length = len(chain)
    substituents = substituents_for_chain(graph, chain, names, aldehydes, mol=mol, aromatic_atoms=ring_atoms)
    grouped = group_substituents(substituents)
    return _name_from_substituents(chain_length, grouped)


def name_aldehyde_ketone(mol) -> str:
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() == 1:
        ring_atoms = set(ring_info.AtomRings()[0])
        if is_plain_benzene_ring(mol, ring_atoms):
            return _name_phenyl_chain_aldehyde_ketone(mol, ring_atoms)

    aldehydes, ketones = _validate_and_collect(mol)
    if mol.GetRingInfo().NumRings() != 0:
        raise UnsupportedStructure(
            "an aldehyde/ketone combination on a ring is out of scope for "
            "this acyclic-only module"
        )
    carbonyl_oxygens = aldehydes | ketones
    all_non_single = non_single_bonds(mol)
    carbonyl_bonds = [b for b in all_non_single if b[0] in carbonyl_oxygens or b[1] in carbonyl_oxygens]
    if len(carbonyl_bonds) != len(all_non_single):
        # Every non-single bond in scope here is a carbonyl C=O (excluded
        # from `carbon_adjacency` already); anything else is chain
        # unsaturation, out of scope for this first pass.
        raise UnsupportedStructure(
            "chain unsaturation (ene/yne) combined with an aldehyde/ketone "
            "demotion is out of scope for this module"
        )

    graph = adjacency(mol)
    (aldehyde_o,) = aldehydes
    (aldehyde_carbon,) = graph[aldehyde_o]
    names = {**halogen_substituents(mol), **{o: "oxo" for o in ketones}}
    chains = longest_chains(carbon_adjacency(mol))
    chain_length = len(chains[0])

    eligible = []
    for chain in chains:
        if aldehyde_carbon not in chain:
            continue
        chain_set = set(chain)
        if any(graph[o][0] not in chain_set for o in ketones):
            continue
        eligible.append(chain)
    if not eligible:
        raise UnsupportedStructure(
            "not every aldehyde/ketone-bearing carbon lies on a single "
            "longest carbon chain; a shorter principal chain is not "
            "supported yet"
        )

    best_key = None
    best_name = None
    for chain in eligible:
        for candidate in (chain, list(reversed(chain))):
            if candidate[0] != aldehyde_carbon:
                # The aldehyde carbon must sit at C1 (P-14.3.3, see module
                # docstring); a direction that doesn't start there is
                # never valid.
                continue
            position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
            substituents = substituents_for_chain(graph, candidate, names, aldehydes, mol=mol)
            key, name = _candidate_key(chain_length, position_of[aldehyde_carbon], substituents)
            if best_key is None or key < best_key:
                best_key, best_name = key, name
    return best_name
