"""Naming of thiones (the '-thione' suffix, C=S with two carbon
substituents) on acyclic saturated or unsaturated carbon chains and on
simple monocyclic saturated rings, per the IUPAC 2013 Recommendations
("the Blue Book"):

- P-64.6.1 (Chapter P-6, https://iupac.qmul.ac.uk/BlueBook/PDF/P6.pdf):
  "Chalcogen analogues of ketones, pseudoketones and heterones are named
  by using the following suffixes and prefixes: =S '-thione' ... =Se
  '-selone' ... =Te '-tellone'" -- this module covers only the sulfur
  case ('-thione'), mirroring `_ketone.py`'s own '-one' suffix
  machinery exactly (a thione carbon, like a ketone carbon, always has
  two carbon neighbors and no hydrogens). Confirmed via PubChem PUG REST
  structure match: `CC(=S)C` -> "propane-2-thione" (matches the Blue
  Book's own worked example "propane-2-thione (PIN) (not thioacetone)"),
  `CCC(=S)C`/`CC(=S)CC` -> "butane-2-thione" (matches the Blue Book's own
  "butane-2-thione (PIN)"), `C1CCC(=S)CC1` -> "cyclohexanethione".
- Two or more C=S groups (a dithione) mirrors the Blue Book's own worked
  example "pentane-2,4-dithione (PIN)" directly -- the locant/suffix
  machinery (ported unchanged from `_ketone.py`) already generalizes over
  a list of thione locants.
- Unlike '-one' (vowel-initial, elides the parent hydride's terminal 'e'),
  '-thione' begins with a consonant, so the terminal 'e' is never elided
  (P-16.3.3): 'propane' + 'thione' -> 'propanethione'/'propane-2-thione',
  not 'propanthione'. `_ketone.py`'s own elision check (whether the
  suffix word's first letter is a vowel) already produces this correctly
  unchanged, since 't' isn't a vowel -- no new logic needed for this.
- Selenium/tellurium analogues ('-selone'/'-tellone') are confirmed to
  exist in the Blue Book (same P-64.6.1 rule text, worked example
  "hexane-3-selone (PIN)") but are deferred to a follow-up task, the same
  way `_selenol.py`/`_tellurol.py` were split into separate tasks despite
  being the same mechanism.

- P-91.3/P-92: a
  molecule with one or more *specified* tetrahedral stereocenters -- every
  one on the principal chain/ring itself, no unspecified one alongside
  them, and no C=C/C#N double-bond E/Z element -- gets a
  "(<locant><R/S>,...)-" prefix, ascending locant order, same pattern as
  `_ketone.py` (chain and ring both). A thione's C=S carbon is
  double-bonded to sulfur exactly like a ketone's C=O carbon to oxygen
  (confirmed via RDKit's `Chem.FindPotentialStereo` to never itself be a
  potential stereocenter), so this support is unconditional, unlike the
  single-bonded, lone-pair-bearing chalcogen in `_sulfinic_acid.py`/
  `_seleninic_acid.py`.

Scope, deliberately narrow (mirrors `_ketone.py`'s own scope, minus its
Table 3.3 ketone/hydroxyl seniority-coexistence handling -- that
combination is future work, tracked separately, same as `_ketone.py`'s
own note about it): one or more C=S groups on an acyclic chain or a
single saturated monocyclic ring, with no other heteroatom (in
particular no ketone C=O, hydroxyl -OH, or any other chalcogen) anywhere
in the molecule.

P-31.1.3: a monocyclic ring bearing a thione and exactly one C=C ring
double bond -- e.g. 'cyclohex-2-ene-1-thione', 'cyclohex-3-ene-1-thione',
both confirmed via PubChem PUG REST. Mirrors `_ketone.py`'s identical
extension; deliberately narrow: a ring triple bond, and any other
substituent alongside the ring double bond, are both still explicitly
rejected pending further verification.

Explicitly out of scope (raise `UnsupportedStructure`): any oxygen at
all, a thione carbon with fewer than two carbon neighbors (a
thial/thioaldehyde, a different suffix), an aromatic thione carbon,
polycyclic/spiro rings, unsaturation reaching outside the ring or a ring
triple bond, and a thione on a substituent branch off an
otherwise-unsubstituted *saturated* ring. One narrow *aromatic*-
ring exception: `_name_phenyl_chain_thione` names a single thione whose
chain hangs off a plain, unsubstituted benzene ring (e.g.
'1-phenylpropane-2-thione'), mirroring `_ketone.py`'s identical benzene-
ring-substituent path -- narrower than the acyclic path: exactly one
thione, no chain unsaturation, and no specified stereocenter, and a
thione carbon directly attached to the ring (an aryl thione) stays out
of scope for this module.
"""

from rdkit import Chem

from ._common import (
    ENE_BOND_ORDER,
    HALOGEN_PREFIXES,
    UnsupportedStructure,
    YNE_BOND_ORDER,
    adjacency,
    bond_locants,
    carbon_adjacency,
    group_substituents,
    halogen_substituents,
    is_plain_benzene_ring,
    longest_branched_chain_through,
    longest_chains,
    lowest_locant_set,
    multiplied_word,
    name_from_substituents,
    non_single_bonds,
    ordered_chain,
    ring_bond_locant,
    ring_bond_locants,
    ring_chain_attachment,
    ring_cycle,
    ring_name_from_substituents,
    specified_stereocenters,
    substituent_locant_set_and_citation,
)
from ._numerals import alkyl_name
from ._substituents import branch_atom_locant, format_substituent_prefixes, name_branch, ring_branch_stereo_display

_SULFUR = 16
_ALLOWED_ATOMIC_NUMS = {6, _SULFUR, *HALOGEN_PREFIXES}


def has_thione_shape(mol) -> bool:
    """True iff `mol` has at least one sulfur double-bonded to a carbon
    (a thione, C=S) -- unlike a plain -SH/-S- sulfur (`_thiol.py`'s/
    `_sulfide.py`'s own loose "any sulfur atom" checks), this is precise
    enough to route correctly regardless of check order."""
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != _SULFUR:
            continue
        if atom.GetDegree() == 1:
            (bond,) = atom.GetBonds()
            if bond.GetBondTypeAsDouble() == 2.0:
                return True
    return False


def _validate_and_collect_thiones(mol, aromatic_ring_atoms=frozenset()):
    """`aromatic_ring_atoms`: atom indices already independently verified
    (by the caller, before this function runs) to form a single plain
    benzene ring with exactly one exocyclic attachment -- exempted from
    the aromatic-atom rejection below so `name_thione`'s benzene-ring-
    substituent path (see `_name_phenyl_chain_thione`) can reuse this same
    validation for the rest of the molecule. Empty by default, so every
    other caller's behavior is unchanged."""
    thiones = set()
    has_carbon = False
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atomic_num not in _ALLOWED_ATOMIC_NUMS:
            raise UnsupportedStructure(
                "heteroatoms other than a thione sulfur (P-64.6.1) and "
                "halogen substituents (P-35.2.1) are not supported yet -- "
                "in particular, a coexisting ketone C=O or hydroxyl -OH "
                "needs Table 3.3 seniority-coexistence handling not yet "
                "implemented for thiones"
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
        elif atomic_num == _SULFUR:
            if atom.GetDegree() != 1:
                raise UnsupportedStructure(
                    "a sulfur bonded to more than one heavy atom (e.g. a "
                    "thioether/sulfide) is out of scope; only an isolated "
                    "thione is supported (P-64.6.1)"
                )
            (bond,) = atom.GetBonds()
            (carbon,) = atom.GetNeighbors()
            if carbon.GetAtomicNum() != 6:
                raise UnsupportedStructure("a thione sulfur must be attached to a carbon atom")
            if bond.GetBondTypeAsDouble() != 2.0:
                raise UnsupportedStructure(
                    "a sulfur that isn't a thione (C=S) is out of scope for "
                    "this module (e.g. a thiol, -SH)"
                )
            if carbon.GetIsAromatic():
                raise UnsupportedStructure(
                    "a thione on an aromatic ring is out of scope for this module"
                )
            carbon_neighbors = [n for n in carbon.GetNeighbors() if n.GetAtomicNum() == 6]
            if len(carbon_neighbors) != 2:
                raise UnsupportedStructure(
                    "a thione carbon with fewer than two carbon neighbors "
                    "(a thial/thioaldehyde) is a different suffix, which "
                    "this module does not attempt to disambiguate"
                )
            thiones.add(atom.GetIdx())
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
    if not thiones:
        raise UnsupportedStructure("no thione (C=S) group found; this module only handles thiones")
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")
    return thiones


def _name_from_substituents(chain_length, thione_locants, ene_locants, yne_locants, grouped):
    total_subs = sum(len(info["locants"]) for info in grouped.values())
    return format_substituent_prefixes(grouped) + name_from_substituents(
        chain_length, ene_locants, yne_locants, multiplied_word(len(thione_locants), "thione"), thione_locants, total_subs
    )


def _candidate_key(chain_length, thione_locants, ene_locants, yne_locants, substituents):
    grouped = group_substituents(substituents)
    locant_set, total_count, citation_locants = substituent_locant_set_and_citation(grouped)
    thione_locant_set = lowest_locant_set(thione_locants)
    combined_locant_set = lowest_locant_set(ene_locants + yne_locants)
    ene_locant_set = lowest_locant_set(ene_locants)
    name = _name_from_substituents(chain_length, thione_locants, ene_locants, yne_locants, grouped)
    return (
        (
            thione_locant_set,
            combined_locant_set,
            ene_locant_set,
            -total_count,
            locant_set,
            citation_locants,
            name,
        ),
        name,
    )


def _thione_locants(position_of, thiones, graph):
    locants = []
    for s in thiones:
        (carbon,) = graph[s]
        if carbon not in position_of:
            return None
        locants.append(position_of[carbon])
    return locants


def _substituents_for_chain(graph, chain, halogens, thiones, mol=None):
    chain_set = set(chain)
    substituents = {}
    for position, atom in enumerate(chain, start=1):
        branch_roots = [n for n in graph[atom] if n not in chain_set and n not in thiones]
        if branch_roots:
            substituents[position] = [name_branch(graph, root, atom, halogens, mol=mol) for root in branch_roots]
    return substituents


def _name_acyclic_thione(mol, thiones, bonds, stereo=None):
    """`stereo`: None, or a list of (stereocenter_atom_idx, "R"/"S") from
    `specified_stereocenters` -- if given, only chain candidates that
    include every stereocenter are eligible (P-92: a stereocenter on a
    substituent branch rather than the principal chain is out of scope,
    mirroring `_ketone.py`'s identical treatment), and the winning
    candidate's own locants are used to format a "(<locant><R/S>,...)-"
    prefix onto the name, ascending locant order (P-91.3)."""
    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    chains = longest_chains(carbon_adjacency(mol))
    chain_length = len(chains[0])
    stereo_atoms = [atom for atom, _ in stereo] if stereo is not None else []

    eligible = []
    for chain in chains:
        position_of = {atom: i + 1 for i, atom in enumerate(chain)}
        if _thione_locants(position_of, thiones, graph) is None:
            continue
        if bonds and bond_locants(chain, bonds) is None:
            continue
        chain_set = set(chain)
        if stereo is not None and any(atom not in chain_set for atom in stereo_atoms):
            continue
        eligible.append(chain)
    if not eligible:
        if stereo is not None and any(
            _thione_locants({a: i + 1 for i, a in enumerate(c)}, thiones, graph) is not None
            and (not bonds or bond_locants(c, bonds) is not None)
            for c in chains
        ):
            raise UnsupportedStructure(
                "a stereocenter on a substituent branch rather than the "
                "principal chain is not supported yet (see P-92)"
            )
        raise UnsupportedStructure(
            "not every thione-bearing carbon (and/or multiple bond) lies "
            "on a single longest carbon chain; a shorter principal chain "
            "capturing more C=S groups is not supported yet"
        )

    best_key = None
    best_name = None
    best_position_of = None
    for chain in eligible:
        for candidate in (chain, list(reversed(chain))):
            position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
            thione_locants = _thione_locants(position_of, thiones, graph)
            ene_locants, yne_locants = bond_locants(candidate, bonds) if bonds else ([], [])
            substituents = _substituents_for_chain(graph, candidate, halogens, thiones, mol=mol)
            key, name = _candidate_key(chain_length, thione_locants, ene_locants, yne_locants, substituents)
            if best_key is None or key < best_key:
                best_key, best_name, best_position_of = key, name, position_of

    if stereo is not None:
        labels = sorted((best_position_of[atom], code) for atom, code in stereo)
        prefix = ",".join(f"{locant}{code}" for locant, code in labels)
        return f"({prefix})-{best_name}"
    return best_name


def _substituents_for_ring(graph, ring_order, halogens, thiones, mol=None):
    ring_set = set(ring_order)
    substituents = {}
    for position, atom in enumerate(ring_order, start=1):
        branch_roots = [n for n in graph[atom] if n not in ring_set and n not in thiones]
        if branch_roots:
            substituents[position] = [name_branch(graph, root, atom, halogens, mol=mol) for root in branch_roots]
    return substituents


def _ring_name_from_substituents(ring_size, thione_locants, ene_locants, yne_locants, grouped):
    total_subs = sum(len(info["locants"]) for info in grouped.values())
    prefix = format_substituent_prefixes(grouped)
    return ring_name_from_substituents(
        ring_size,
        ene_locants,
        yne_locants,
        prefix,
        total_subs,
        multiplied_word(len(thione_locants), "thione"),
        thione_locants,
    )


def _ring_candidate_key(ring_size, thione_locants, ene_locants, yne_locants, substituents):
    grouped = group_substituents(substituents)
    locant_set, _, citation_locants = substituent_locant_set_and_citation(grouped)
    thione_locant_set = lowest_locant_set(thione_locants)
    combined_locant_set = lowest_locant_set(ene_locants + yne_locants)
    ene_locant_set = lowest_locant_set(ene_locants)
    name = _ring_name_from_substituents(ring_size, thione_locants, ene_locants, yne_locants, grouped)
    return thione_locant_set, combined_locant_set, ene_locant_set, locant_set, citation_locants, name


def _ring_branch_stereo_display(graph, ring_order, thiones, stereo, halogens, mol=None):
    return ring_branch_stereo_display(graph, ring_order, thiones, stereo, halogens, mol=mol)


def _name_cyclic_thione(mol, thiones, stereo=None, bonds=()):
    """`stereo`: None, or a list of (stereocenter_atom_idx, "R"/"S") from
    `specified_stereocenters` -- if given, every stereocenter must normally
    lie on the ring itself (P-92: a stereocenter on a substituent branch is
    out of scope, mirroring `_ketone.py`'s `_name_cyclic_ketone`), and the
    winning ring numbering's own locants for those atoms are used to
    format a "(<locant><R/S>,...)-" prefix onto the name, ascending
    locant order (P-91.3). The one narrow exception
    (`_ring_branch_stereo_display`, mirroring `_ketone.py`/`_selone.py`/
    `_tellone.py`'s own case): exactly one stereocenter on the ring's sole
    substituent branch instead embeds a bracketed descriptor into that
    substituent's own name, in place of the usual ring-locant prefix."""
    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    ring_info = mol.GetRingInfo()
    ring_atoms = list(ring_info.AtomRings()[0])
    ring_order = ring_cycle(graph, ring_atoms)
    ring_size = len(ring_order)
    branch_stereo = None
    if stereo is not None and any(atom not in ring_order for atom, _ in stereo):
        branch_stereo = _ring_branch_stereo_display(graph, ring_order, thiones, stereo, halogens, mol=mol)
        if branch_stereo is None:
            raise UnsupportedStructure(
                "a stereocenter on a substituent branch rather than the ring "
                "itself is not supported yet (see P-92)"
            )
    if bonds and any(_substituents_for_ring(graph, ring_order, halogens, thiones, mol=mol).values()):
        raise UnsupportedStructure(
            "a substituent alongside both a ring double/triple bond and a "
            "thione is not supported yet (see module docstring)"
        )

    best_key = None
    best_name = None
    best_position_of = None
    for start in range(ring_size):
        rotated = ring_order[start:] + ring_order[:start]
        for candidate in (rotated, list(reversed(rotated))):
            position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
            thione_locants = _thione_locants(position_of, thiones, graph)
            if thione_locants is None:
                raise UnsupportedStructure(
                    "a thione not on the ring itself (e.g. on a "
                    "substituent branch) is not supported yet"
                )
            substituents = _substituents_for_ring(graph, candidate, halogens, thiones, mol=mol)
            if branch_stereo is not None:
                branch_ring_atom, display = branch_stereo
                substituents[position_of[branch_ring_atom]] = [(display, False)]
            ene_locants, yne_locants = ring_bond_locants(position_of, bonds, ring_size)
            key = _ring_candidate_key(ring_size, thione_locants, ene_locants, yne_locants, substituents)
            if best_key is None or key < best_key:
                best_key, best_name, best_position_of = key, key[-1], position_of

    if stereo is not None and branch_stereo is None:
        labels = sorted((best_position_of[atom], r_or_s) for atom, r_or_s in stereo)
        prefix = ",".join(f"{locant}{r_or_s}" for locant, r_or_s in labels)
        return f"({prefix})-{best_name}"
    return best_name


def _name_phenyl_chain_thione(mol, ring_atoms):
    """Name a thione whose C=S lies entirely on a single unbranched chain
    hanging off one atom of an otherwise-plain, unsubstituted benzene ring
    -- e.g. 1-phenylpropane-2-thione. The ring is cited as a 'phenyl'
    substituent prefix (via `name_branch`'s aromatic-ring recognition) on
    the chain, which is the parent hydride, mirroring `_ketone.py`'s
    `_name_phenyl_chain_ketone`. Narrower than the acyclic path above:
    exactly one thione, no chain unsaturation, and no specified
    stereocenter -- each is a separate follow-up rather than being
    combined with the ring case in this first slice. A thione carbon
    directly attached to the ring (an aryl thione) stays out of scope,
    same as the module's existing acyclic-carbonyl check above."""
    thiones = _validate_and_collect_thiones(mol, aromatic_ring_atoms=ring_atoms)
    if len(thiones) != 1:
        raise UnsupportedStructure(
            "more than one thione group alongside a benzene-ring "
            "substituent is not supported yet"
        )
    if specified_stereocenters(mol):
        raise UnsupportedStructure(
            "a specified stereocenter alongside a benzene-ring-substituent "
            "thione chain is not supported yet"
        )
    non_ring_unsaturation = [
        b
        for b in non_single_bonds(mol)
        if b[0] not in thiones and b[1] not in thiones and b[0] not in ring_atoms and b[1] not in ring_atoms
    ]
    if non_ring_unsaturation:
        raise UnsupportedStructure(
            "chain unsaturation alongside a benzene-ring-substituent "
            "thione chain is not supported yet"
        )

    graph = adjacency(mol)
    attachment = ring_chain_attachment(graph, ring_atoms, set())
    if attachment is None:
        raise UnsupportedStructure(
            "a benzene ring with more than one exocyclic substituent "
            "alongside a chain thione is not supported yet"
        )
    ring_atom, chain_root = attachment
    (thione_sulfur,) = thiones
    (thione_carbon,) = graph[thione_sulfur]
    if thione_carbon == chain_root:
        raise UnsupportedStructure(
            "a thione carbon directly attached to the benzene ring (an "
            "aryl thione) is out of scope for this module (see the "
            "separate aromatic-ring module)"
        )

    chain, branches = longest_branched_chain_through(graph, thione_carbon, ring_atoms, thiones, halogens=halogen_substituents(mol))
    branches_by_atom = {chain[position - 1]: roots for position, roots in branches.items()}

    chain_length = len(chain)
    halogens = halogen_substituents(mol)
    best_key = None
    best_name = None
    for candidate in (chain, list(reversed(chain))):
        position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
        thione_locants = _thione_locants(position_of, thiones, graph)
        substituents = {
            position_of[atom]: [name_branch(graph, root, atom, halogens, ring_atoms, mol=mol) for root in roots]
            for atom, roots in branches_by_atom.items()
        }
        key, name = _candidate_key(chain_length, thione_locants, [], [], substituents)
        if best_key is None or key < best_key:
            best_key, best_name = key, name
    return best_name


def _name_ring_substituent_chain_thione(mol, thiones):
    """Name a thione whose C=S lies entirely on a single
    branched chain hanging off one atom of an otherwise-plain saturated
    monocyclic ring (the ring itself bears no thione) -- e.g.
    1-cyclohexylethanethione. The ring is cited as a "cyclo..."
    substituent prefix (P-29.3.3) on the chain, which is the parent
    hydride, mirroring `_name_phenyl_chain_thione` above and
    `_ketone.py`'s `_name_ring_substituent_chain_ketone`. Narrower than
    the benzene-ring case: exactly one thione."""
    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    ring_atoms = set(mol.GetRingInfo().AtomRings()[0])

    attachment = ring_chain_attachment(graph, ring_atoms, thiones)
    if attachment is None:
        raise UnsupportedStructure(
            "a ring with more than one exocyclic branch is not supported "
            "yet"
        )
    ring_atom, chain_root = attachment
    (thione_sulfur,) = thiones
    (thione_carbon,) = graph[thione_sulfur]

    chain, branches = longest_branched_chain_through(graph, thione_carbon, ring_atoms, thiones, halogens=halogen_substituents(mol))
    branches_by_atom = {
        chain[position - 1]: [r for r in roots if r != ring_atom]
        for position, roots in branches.items()
    }
    branches_by_atom = {atom: roots for atom, roots in branches_by_atom.items() if roots}

    ring_name = "cyclo" + alkyl_name(len(ring_atoms))
    chain_length = len(chain)

    best_key = None
    best_name = None
    for candidate in (chain, list(reversed(chain))):
        position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
        thione_locants = _thione_locants(position_of, thiones, graph)
        substituents = {
            position_of[atom]: [name_branch(graph, root, atom, halogens, mol=mol) for root in roots]
            for atom, roots in branches_by_atom.items()
        }
        substituents.setdefault(position_of[chain_root], []).append((ring_name, False))
        key, name = _candidate_key(chain_length, thione_locants, [], [], substituents)
        if best_key is None or key < best_key:
            best_key, best_name = key, name
    return best_name


def _name_ring_with_thione_chain_thione(mol, thiones):
    """Name a thione compound where the ring itself bears at least
    as many thiones as a single unbranched chain hanging off
    exactly one ring atom does (P-44.1.1/P-44.1.2.2). The ring is the
    parent; the chain is cited as a '(sulfanylidene...alkyl)' substituent
    prefix, mirroring `_ketone.py`'s `_name_ring_with_ketone_chain_
    ketone`."""
    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    ring_atoms = set(mol.GetRingInfo().AtomRings()[0])

    attachment = ring_chain_attachment(graph, ring_atoms, thiones)
    if attachment is None:
        raise UnsupportedStructure(
            "a ring with more than one exocyclic branch is not supported "
            "yet"
        )
    ring_atom, chain_root = attachment
    chain = ordered_chain(graph, chain_root, ring_atom, thiones)
    if chain is None:
        raise UnsupportedStructure(
            "a branched substituent chain hanging off the ring is not "
            "supported yet"
        )

    chain_set = set(chain)
    chain_thiones = {o for o in thiones if next(iter(graph[o])) in chain_set}
    ring_thiones = thiones - chain_thiones
    if len(ring_thiones) < len(chain_thiones):
        ring_name, ring_is_compound = name_branch(
            graph, ring_atom, chain_root, {**halogens, **{o: "sulfanylidene" for o in ring_thiones}}, mol=mol
        )
        chain_length = len(chain)
        best_key = None
        best_name = None
        for candidate in (chain, list(reversed(chain))):
            position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
            thione_locants = _thione_locants(position_of, chain_thiones, graph)
            substituents = {position_of[chain_root]: [(ring_name, ring_is_compound)]}
            key, name = _candidate_key(chain_length, thione_locants, [], [], substituents)
            if best_key is None or key < best_key:
                best_key, best_name = key, name
        return best_name

    chain_name, chain_is_compound = name_branch(
        graph, chain_root, ring_atom, {**halogens, **{o: "sulfanylidene" for o in chain_thiones}}, mol=mol
    )

    ring_order = ring_cycle(graph, list(ring_atoms))
    ring_size = len(ring_order)
    best_key = None
    best_name = None
    for start in range(ring_size):
        rotated = ring_order[start:] + ring_order[:start]
        for candidate in (rotated, list(reversed(rotated))):
            position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
            thione_locants = _thione_locants(position_of, ring_thiones, graph)
            substituents = {position_of[ring_atom]: [(chain_name, chain_is_compound)]}
            key = _ring_candidate_key(ring_size, thione_locants, [], [], substituents)
            if best_key is None or key < best_key:
                best_key, best_name = key, key[-1]
    return best_name


def name_thione(mol) -> str:
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() == 1:
        ring_atoms = set(ring_info.AtomRings()[0])
        if is_plain_benzene_ring(mol, ring_atoms):
            return _name_phenyl_chain_thione(mol, ring_atoms)
    thiones = _validate_and_collect_thiones(mol)
    stereo = specified_stereocenters(mol)
    graph = adjacency(mol)
    all_non_single = [b for b in non_single_bonds(mol) if b[0] not in thiones and b[1] not in thiones]
    bonds = [b for b in all_non_single if b[2] in (ENE_BOND_ORDER, YNE_BOND_ORDER)]
    if len(bonds) != len(all_non_single):
        raise UnsupportedStructure(
            "a bond order other than single, double, or triple is not "
            "supported (see P-31.1.1.1)"
        )

    ring_info = mol.GetRingInfo()
    num_rings = ring_info.NumRings()
    if num_rings == 0:
        return _name_acyclic_thione(mol, thiones, bonds, stereo)
    if num_rings == 1:
        ring_atoms = set(ring_info.AtomRings()[0])
        if any(a not in ring_atoms or b not in ring_atoms for a, b, _ in bonds):
            raise UnsupportedStructure(
                "unsaturation outside the ring alongside a cyclic thione "
                "is not supported yet (see P-31.1.3, cycloalkenes and "
                "cycloalkynes)"
            )
        if any(order == YNE_BOND_ORDER for _, _, order in bonds):
            raise UnsupportedStructure(
                "a ring triple bond (cycloalkyne) alongside a thione is "
                "not supported yet -- only a ring double bond is in scope "
                "for this first pass (see P-31.1.3)"
            )
        ring_thiones = {s for s in thiones if next(iter(graph[s])) in ring_atoms}
        if not bonds and not ring_thiones and len(thiones) == 1:
            if stereo is not None:
                raise UnsupportedStructure(
                    "a stereocenter on a substituent branch rather than "
                    "the ring itself is not supported yet (see P-92)"
                )
            return _name_ring_substituent_chain_thione(mol, thiones)
        if not bonds and ring_thiones and ring_thiones != thiones:
            if stereo is not None:
                raise UnsupportedStructure(
                    "a stereocenter alongside a ring-vs-chain thione "
                    "comparison is not supported yet (see P-92)"
                )
            return _name_ring_with_thione_chain_thione(mol, thiones)
        return _name_cyclic_thione(mol, thiones, stereo, bonds)
    raise UnsupportedStructure(
        "polycyclic and spiro thiones are not supported yet"
    )
