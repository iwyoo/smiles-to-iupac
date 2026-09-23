"""Naming of selones (the '-selone' suffix, C=Se with two carbon
substituents), the selenium analogue of ketones, on acyclic saturated or
unsaturated carbon chains and on simple monocyclic saturated rings, per
the IUPAC 2013 Recommendations ("the Blue Book"):

- P-64.6.1 (Chapter P-6, https://iupac.qmul.ac.uk/BlueBook/PDF/P6.pdf):
  "Chalcogen analogues of ketones, pseudoketones and heterones are named
  by using the following suffixes and prefixes: ... =Se '-selone' ..." --
  this module mirrors `_thione.py`'s own '-thione' (C=S) suffix machinery
  exactly, just with selenium instead of sulfur (a selone carbon, like a
  ketone/thione carbon, always has two carbon neighbors and no
  hydrogens). Confirmed via PubChem PUG REST structure match:
  `CC(=[Se])C` -> "propane-2-selone" (the same structure as the Blue
  Book's own worked example for the sulfur case, "propane-2-thione
  (PIN)"), `CC(=[Te])C`'s own sibling case in `_tellone.py` confirms the
  pattern generalizes across chalcogens. The Blue Book's own worked
  example for this exact suffix is "hexane-3-selone (PIN)"
  (`tmp/bluebook/P6.txt`), matching the structure of `CCC(=[Se])CC`.
- Two or more C=Se groups (a diselone) mirrors `_thione.py`'s own
  dithione generalization -- the locant/suffix machinery (ported
  unchanged from `_thione.py`) already generalizes over a list of selone
  locants.
- Like '-thione', '-selone' begins with a consonant, so the parent
  hydride's terminal 'e' is never elided (P-16.3.3).

- P-91.3/P-92: a
  molecule with one or more *specified* tetrahedral stereocenters -- every
  one on the principal chain/ring itself, no unspecified one alongside
  them, and no C=C/C#N double-bond E/Z element -- gets a
  "(<locant><R/S>,...)-" prefix, ascending locant order, same pattern as
  `_thione.py`/`_ketone.py` (chain and ring both). A selone's C=Se carbon
  is double-bonded to selenium exactly like a ketone's C=O carbon to
  oxygen (confirmed via RDKit's `Chem.FindPotentialStereo` to never
  itself be a potential stereocenter), so this support is unconditional.

Scope, deliberately narrow (mirrors `_thione.py`'s own scope exactly):
one or more C=Se groups on an acyclic chain or a single saturated
monocyclic ring, with no other heteroatom (in particular no ketone C=O,
thione C=S, hydroxyl -OH, or any other chalcogen) anywhere in the
molecule.

P-31.1.3: a monocyclic ring bearing a selone and exactly one C=C ring
double bond -- e.g. 'cyclohex-2-ene-1-selone', 'cyclohex-3-ene-1-selone',
both confirmed via PubChem PUG REST. Mirrors `_thione.py`'s identical
extension; deliberately narrow: a ring triple bond, and any other
substituent alongside the ring double bond, are both still explicitly
rejected pending further verification.

Explicitly out of scope (raise `UnsupportedStructure`): any oxygen or
sulfur at all, a selone carbon with fewer than two carbon neighbors (a
selenoaldehyde, a different suffix), an aromatic selone carbon; more than
one selone, or a specified stereocenter, on a von Baeyer polycyclic or
spiro skeleton (a single selone on such a skeleton is supported, P-23/
P-24 numbering integration via `_polycyclic_suffix.py` with
`elide_e=False`, mirroring `_thiol.py`); unsaturation reaching outside
the ring or a ring triple bond, and a selone on a substituent
branch off an otherwise-unsubstituted *saturated* ring. One narrow
*aromatic*-ring exception: `_name_phenyl_chain_selone` names a single
selone whose chain hangs off a plain, unsubstituted benzene ring (e.g.
'1-phenylpropane-2-selone'), mirroring `_thione.py`'s identical benzene-
ring-substituent path -- narrower than the acyclic path: exactly one
selone, no chain unsaturation, and no specified stereocenter, and a
selone carbon directly attached to the ring (an aryl selone) stays out
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
from ._bicyclic import find_bicyclic_core
from ._numerals import alkyl_name
from ._polycyclic import find_polycyclic_core
from ._polycyclic_suffix import name_monospiro_suffix, name_von_baeyer_suffix
from ._spiro import find_monospiro_atom
from ._substituents import (
    substituents_for_chain,
    branch_atom_locant,
    format_substituent_prefixes,
    name_branch,
    ring_branch_stereo_display,
    substituents_for_ring,
)

_SELENIUM = 34
_ALLOWED_ATOMIC_NUMS = {6, _SELENIUM, *HALOGEN_PREFIXES}


def has_selone_shape(mol) -> bool:
    """True iff `mol` has at least one selenium double-bonded to a carbon
    (a selone, C=Se) -- unlike a plain -SeH/-Se- selenium (`_selenol.py`'s/
    `_selenide.py`'s own loose "any selenium atom" checks), this is precise
    enough to route correctly regardless of check order."""
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != _SELENIUM:
            continue
        if atom.GetDegree() == 1:
            (bond,) = atom.GetBonds()
            if bond.GetBondTypeAsDouble() == 2.0:
                return True
    return False


def _validate_and_collect_selones(mol, aromatic_ring_atoms=frozenset()):
    """`aromatic_ring_atoms`: atom indices already independently verified
    (by the caller, before this function runs) to form a single plain
    benzene ring with exactly one exocyclic attachment -- exempted from
    the aromatic-atom rejection below so `name_selone`'s benzene-ring-
    substituent path (see `_name_phenyl_chain_selone`) can reuse this same
    validation for the rest of the molecule. Empty by default, so every
    other caller's behavior is unchanged."""
    selones = set()
    has_carbon = False
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atomic_num not in _ALLOWED_ATOMIC_NUMS:
            raise UnsupportedStructure(
                "heteroatoms other than a selone selenium (P-64.6.1) and "
                "halogen substituents (P-35.2.1) are not supported yet -- "
                "in particular, a coexisting ketone C=O or hydroxyl -OH "
                "needs Table 3.3 seniority-coexistence handling not yet "
                "implemented for selones"
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
        elif atomic_num == _SELENIUM:
            if atom.GetDegree() != 1:
                raise UnsupportedStructure(
                    "a selenium bonded to more than one heavy atom (e.g. a "
                    "selenoether/selenide) is out of scope; only an isolated "
                    "selone is supported (P-64.6.1)"
                )
            (bond,) = atom.GetBonds()
            (carbon,) = atom.GetNeighbors()
            if carbon.GetAtomicNum() != 6:
                raise UnsupportedStructure("a selone selenium must be attached to a carbon atom")
            if bond.GetBondTypeAsDouble() != 2.0:
                raise UnsupportedStructure(
                    "a selenium that isn't a selone (C=Se) is out of scope for "
                    "this module (e.g. a selenol, -SeH)"
                )
            if carbon.GetIsAromatic():
                raise UnsupportedStructure(
                    "a selone on an aromatic ring is out of scope for this module"
                )
            carbon_neighbors = [n for n in carbon.GetNeighbors() if n.GetAtomicNum() == 6]
            if len(carbon_neighbors) != 2:
                raise UnsupportedStructure(
                    "a selone carbon with fewer than two carbon neighbors "
                    "(a selenoaldehyde) is a different suffix, which "
                    "this module does not attempt to disambiguate"
                )
            selones.add(atom.GetIdx())
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
    if not selones:
        raise UnsupportedStructure("no selone (C=Se) group found; this module only handles selones")
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")
    return selones


def _name_from_substituents(chain_length, selone_locants, ene_locants, yne_locants, grouped):
    return format_substituent_prefixes(grouped) + name_from_substituents(
        chain_length, ene_locants, yne_locants, multiplied_word(len(selone_locants), "selone"), selone_locants
    )


def _candidate_key(chain_length, selone_locants, ene_locants, yne_locants, substituents):
    grouped = group_substituents(substituents)
    locant_set, total_count, citation_locants = substituent_locant_set_and_citation(grouped)
    selone_locant_set = lowest_locant_set(selone_locants)
    combined_locant_set = lowest_locant_set(ene_locants + yne_locants)
    ene_locant_set = lowest_locant_set(ene_locants)
    name = _name_from_substituents(chain_length, selone_locants, ene_locants, yne_locants, grouped)
    return (
        (
            selone_locant_set,
            combined_locant_set,
            ene_locant_set,
            -total_count,
            locant_set,
            citation_locants,
            name,
        ),
        name,
    )


def _selone_locants(position_of, selones, graph):
    locants = []
    for s in selones:
        (carbon,) = graph[s]
        if carbon not in position_of:
            return None
        locants.append(position_of[carbon])
    return locants


def _name_acyclic_selone(mol, selones, bonds, stereo=None):
    """`stereo`: None, or a list of (stereocenter_atom_idx, "R"/"S") from
    `specified_stereocenters` -- if given, only chain candidates that
    include every stereocenter are eligible (P-92: a stereocenter on a
    substituent branch rather than the principal chain is out of scope,
    mirroring `_thione.py`'s identical treatment), and the winning
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
        if _selone_locants(position_of, selones, graph) is None:
            continue
        if bonds and bond_locants(chain, bonds) is None:
            continue
        chain_set = set(chain)
        if stereo is not None and any(atom not in chain_set for atom in stereo_atoms):
            continue
        eligible.append(chain)
    if not eligible:
        if stereo is not None and any(
            _selone_locants({a: i + 1 for i, a in enumerate(c)}, selones, graph) is not None
            and (not bonds or bond_locants(c, bonds) is not None)
            for c in chains
        ):
            raise UnsupportedStructure(
                "a stereocenter on a substituent branch rather than the "
                "principal chain is not supported yet (see P-92)"
            )
        raise UnsupportedStructure(
            "not every selone-bearing carbon (and/or multiple bond) lies "
            "on a single longest carbon chain; a shorter principal chain "
            "capturing more C=Se groups is not supported yet"
        )

    best_key = None
    best_name = None
    best_position_of = None
    for chain in eligible:
        for candidate in (chain, list(reversed(chain))):
            position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
            selone_locants = _selone_locants(position_of, selones, graph)
            ene_locants, yne_locants = bond_locants(candidate, bonds) if bonds else ([], [])
            substituents = substituents_for_chain(graph, candidate, halogens, selones, mol=mol)
            key, name = _candidate_key(chain_length, selone_locants, ene_locants, yne_locants, substituents)
            if best_key is None or key < best_key:
                best_key, best_name, best_position_of = key, name, position_of

    if stereo is not None:
        labels = sorted((best_position_of[atom], code) for atom, code in stereo)
        prefix = ",".join(f"{locant}{code}" for locant, code in labels)
        return f"({prefix})-{best_name}"
    return best_name




def _ring_name_from_substituents(ring_size, selone_locants, ene_locants, yne_locants, grouped):
    total_subs = sum(len(info["locants"]) for info in grouped.values())
    prefix = format_substituent_prefixes(grouped)
    return ring_name_from_substituents(
        ring_size,
        ene_locants,
        yne_locants,
        prefix,
        total_subs,
        multiplied_word(len(selone_locants), "selone"),
        selone_locants,
    )


def _ring_candidate_key(ring_size, selone_locants, ene_locants, yne_locants, substituents):
    grouped = group_substituents(substituents)
    locant_set, _, citation_locants = substituent_locant_set_and_citation(grouped)
    selone_locant_set = lowest_locant_set(selone_locants)
    combined_locant_set = lowest_locant_set(ene_locants + yne_locants)
    ene_locant_set = lowest_locant_set(ene_locants)
    name = _ring_name_from_substituents(ring_size, selone_locants, ene_locants, yne_locants, grouped)
    return selone_locant_set, combined_locant_set, ene_locant_set, locant_set, citation_locants, name


def _ring_branch_stereo_display(graph, ring_order, selones, stereo, halogens, mol=None):
    return ring_branch_stereo_display(graph, ring_order, selones, stereo, halogens, mol=mol)


def _name_cyclic_selone(mol, selones, stereo=None, bonds=()):
    """`stereo`: None, or a list of (stereocenter_atom_idx, "R"/"S") from
    `specified_stereocenters` -- if given, every stereocenter must normally
    lie on the ring itself (P-92: a stereocenter on a substituent branch is
    out of scope, mirroring `_thione.py`'s `_name_cyclic_thione`), and the
    winning ring numbering's own locants for those atoms are used to
    format a "(<locant><R/S>,...)-" prefix onto the name, ascending
    locant order (P-91.3). The one narrow exception
    (`_ring_branch_stereo_display`, mirroring `_ketone.py`'s own case):
    exactly one stereocenter on the ring's sole substituent branch instead
    embeds a bracketed descriptor into that substituent's own name, in
    place of the usual ring-locant prefix."""
    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    ring_info = mol.GetRingInfo()
    ring_atoms = list(ring_info.AtomRings()[0])
    ring_order = ring_cycle(graph, ring_atoms)
    ring_size = len(ring_order)
    branch_stereo = None
    if stereo is not None and any(atom not in ring_order for atom, _ in stereo):
        branch_stereo = _ring_branch_stereo_display(graph, ring_order, selones, stereo, halogens, mol=mol)
        if branch_stereo is None:
            raise UnsupportedStructure(
                "a stereocenter on a substituent branch rather than the ring "
                "itself is not supported yet (see P-92)"
            )
    if bonds and any(substituents_for_ring(graph, ring_order, halogens, selones, mol=mol).values()):
        raise UnsupportedStructure(
            "a substituent alongside both a ring double/triple bond and a "
            "selone is not supported yet (see module docstring)"
        )

    best_key = None
    best_name = None
    best_position_of = None
    for start in range(ring_size):
        rotated = ring_order[start:] + ring_order[:start]
        for candidate in (rotated, list(reversed(rotated))):
            position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
            selone_locants = _selone_locants(position_of, selones, graph)
            if selone_locants is None:
                raise UnsupportedStructure(
                    "a selone not on the ring itself (e.g. on a "
                    "substituent branch) is not supported yet"
                )
            substituents = substituents_for_ring(graph, candidate, halogens, selones, mol=mol)
            if branch_stereo is not None:
                branch_ring_atom, display = branch_stereo
                substituents[position_of[branch_ring_atom]] = [(display, False)]
            ene_locants, yne_locants = ring_bond_locants(position_of, bonds, ring_size)
            key = _ring_candidate_key(ring_size, selone_locants, ene_locants, yne_locants, substituents)
            if best_key is None or key < best_key:
                best_key, best_name, best_position_of = key, key[-1], position_of

    if stereo is not None and branch_stereo is None:
        labels = sorted((best_position_of[atom], r_or_s) for atom, r_or_s in stereo)
        prefix = ",".join(f"{locant}{r_or_s}" for locant, r_or_s in labels)
        return f"({prefix})-{best_name}"
    return best_name


def _name_phenyl_chain_selone(mol, ring_atoms):
    """Name a selone whose C=Se lies entirely on a single unbranched
    chain hanging off one atom of an otherwise-plain, unsubstituted
    benzene ring -- e.g. 1-phenylpropane-2-selone. The ring is cited as a
    'phenyl' substituent prefix (via `name_branch`'s aromatic-ring
    recognition) on the chain, which is the parent hydride, mirroring
    `_thione.py`'s `_name_phenyl_chain_thione`. Narrower than the acyclic
    path above: exactly one selone, no chain unsaturation, and no
    specified stereocenter -- each is a separate follow-up rather than
    being combined with the ring case in this first slice. A selone
    carbon directly attached to the ring (an aryl selone) stays out of
    scope, same as the module's existing acyclic-carbonyl check above."""
    selones = _validate_and_collect_selones(mol, aromatic_ring_atoms=ring_atoms)
    if len(selones) != 1:
        raise UnsupportedStructure(
            "more than one selone group alongside a benzene-ring "
            "substituent is not supported yet"
        )
    if specified_stereocenters(mol):
        raise UnsupportedStructure(
            "a specified stereocenter alongside a benzene-ring-substituent "
            "selone chain is not supported yet"
        )
    non_ring_unsaturation = [
        b
        for b in non_single_bonds(mol)
        if b[0] not in selones and b[1] not in selones and b[0] not in ring_atoms and b[1] not in ring_atoms
    ]
    if non_ring_unsaturation:
        raise UnsupportedStructure(
            "chain unsaturation alongside a benzene-ring-substituent "
            "selone chain is not supported yet"
        )

    graph = adjacency(mol)
    attachment = ring_chain_attachment(graph, ring_atoms, set())
    if attachment is None:
        raise UnsupportedStructure(
            "a benzene ring with more than one exocyclic substituent "
            "alongside a chain selone is not supported yet"
        )
    ring_atom, chain_root = attachment
    (selone_selenium,) = selones
    (selone_carbon,) = graph[selone_selenium]
    if selone_carbon == chain_root:
        raise UnsupportedStructure(
            "a selone carbon directly attached to the benzene ring (an "
            "aryl selone) is out of scope for this module (see the "
            "separate aromatic-ring module)"
        )

    chain, branches = longest_branched_chain_through(graph, selone_carbon, ring_atoms, selones, halogens=halogen_substituents(mol))
    branches_by_atom = {chain[position - 1]: roots for position, roots in branches.items()}

    chain_length = len(chain)
    halogens = halogen_substituents(mol)
    best_key = None
    best_name = None
    for candidate in (chain, list(reversed(chain))):
        position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
        selone_locants = _selone_locants(position_of, selones, graph)
        substituents = {
            position_of[atom]: [name_branch(graph, root, atom, halogens, ring_atoms, mol=mol) for root in roots]
            for atom, roots in branches_by_atom.items()
        }
        key, name = _candidate_key(chain_length, selone_locants, [], [], substituents)
        if best_key is None or key < best_key:
            best_key, best_name = key, name
    return best_name


def _name_ring_substituent_chain_selone(mol, selones):
    """Name a selone whose C=S lies entirely on a single
    branched chain hanging off one atom of an otherwise-plain saturated
    monocyclic ring (the ring itself bears no selone) -- e.g.
    1-cyclohexylethaneselone. The ring is cited as a "cyclo..."
    substituent prefix (P-29.3.3) on the chain, which is the parent
    hydride, mirroring `_name_phenyl_chain_selone` above and
    `_ketone.py`'s `_name_ring_substituent_chain_ketone`. Narrower than
    the benzene-ring case: exactly one selone."""
    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    ring_atoms = set(mol.GetRingInfo().AtomRings()[0])

    attachment = ring_chain_attachment(graph, ring_atoms, selones)
    if attachment is None:
        raise UnsupportedStructure(
            "a ring with more than one exocyclic branch is not supported "
            "yet"
        )
    ring_atom, chain_root = attachment
    (selone_selenium,) = selones
    (selone_carbon,) = graph[selone_selenium]

    chain, branches = longest_branched_chain_through(graph, selone_carbon, ring_atoms, selones, halogens=halogen_substituents(mol))
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
        selone_locants = _selone_locants(position_of, selones, graph)
        substituents = {
            position_of[atom]: [name_branch(graph, root, atom, halogens, mol=mol) for root in roots]
            for atom, roots in branches_by_atom.items()
        }
        substituents.setdefault(position_of[chain_root], []).append((ring_name, False))
        key, name = _candidate_key(chain_length, selone_locants, [], [], substituents)
        if best_key is None or key < best_key:
            best_key, best_name = key, name
    return best_name


def _name_ring_with_selone_chain_selone(mol, selones):
    """Name a selone compound where the ring itself bears at least
    as many selones as a single unbranched chain hanging off
    exactly one ring atom does (P-44.1.1/P-44.1.2.2). The ring is the
    parent; the chain is cited as a '(selanylidene...alkyl)' substituent
    prefix, mirroring `_ketone.py`'s `_name_ring_with_ketone_chain_
    ketone`."""
    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    ring_atoms = set(mol.GetRingInfo().AtomRings()[0])

    attachment = ring_chain_attachment(graph, ring_atoms, selones)
    if attachment is None:
        raise UnsupportedStructure(
            "a ring with more than one exocyclic branch is not supported "
            "yet"
        )
    ring_atom, chain_root = attachment
    chain = ordered_chain(graph, chain_root, ring_atom, selones)
    if chain is None:
        raise UnsupportedStructure(
            "a branched substituent chain hanging off the ring is not "
            "supported yet"
        )

    chain_set = set(chain)
    chain_selones = {o for o in selones if next(iter(graph[o])) in chain_set}
    ring_selones = selones - chain_selones
    if len(ring_selones) < len(chain_selones):
        ring_name, ring_is_compound = name_branch(
            graph, ring_atom, chain_root, {**halogens, **{o: "selanylidene" for o in ring_selones}}, mol=mol
        )
        chain_length = len(chain)
        best_key = None
        best_name = None
        for candidate in (chain, list(reversed(chain))):
            position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
            selone_locants = _selone_locants(position_of, chain_selones, graph)
            substituents = {position_of[chain_root]: [(ring_name, ring_is_compound)]}
            key, name = _candidate_key(chain_length, selone_locants, [], [], substituents)
            if best_key is None or key < best_key:
                best_key, best_name = key, name
        return best_name

    chain_name, chain_is_compound = name_branch(
        graph, chain_root, ring_atom, {**halogens, **{o: "selanylidene" for o in chain_selones}}, mol=mol
    )

    ring_order = ring_cycle(graph, list(ring_atoms))
    ring_size = len(ring_order)
    best_key = None
    best_name = None
    for start in range(ring_size):
        rotated = ring_order[start:] + ring_order[:start]
        for candidate in (rotated, list(reversed(rotated))):
            position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
            selone_locants = _selone_locants(position_of, ring_selones, graph)
            substituents = {position_of[ring_atom]: [(chain_name, chain_is_compound)]}
            key = _ring_candidate_key(ring_size, selone_locants, [], [], substituents)
            if best_key is None or key < best_key:
                best_key, best_name = key, key[-1]
    return best_name


def name_selone(mol) -> str:
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() == 1:
        ring_atoms = set(ring_info.AtomRings()[0])
        if is_plain_benzene_ring(mol, ring_atoms):
            return _name_phenyl_chain_selone(mol, ring_atoms)
    selones = _validate_and_collect_selones(mol)
    stereo = specified_stereocenters(mol)
    graph = adjacency(mol)
    all_non_single = [b for b in non_single_bonds(mol) if b[0] not in selones and b[1] not in selones]
    bonds = [b for b in all_non_single if b[2] in (ENE_BOND_ORDER, YNE_BOND_ORDER)]
    if len(bonds) != len(all_non_single):
        raise UnsupportedStructure(
            "a bond order other than single, double, or triple is not "
            "supported (see P-31.1.1.1)"
        )

    ring_info = mol.GetRingInfo()
    num_rings = ring_info.NumRings()
    if num_rings == 0:
        return _name_acyclic_selone(mol, selones, bonds, stereo)
    if num_rings == 1:
        ring_atoms = set(ring_info.AtomRings()[0])
        if any(a not in ring_atoms or b not in ring_atoms for a, b, _ in bonds):
            raise UnsupportedStructure(
                "unsaturation outside the ring alongside a cyclic selone "
                "is not supported yet (see P-31.1.3, cycloalkenes and "
                "cycloalkynes)"
            )
        if any(order == YNE_BOND_ORDER for _, _, order in bonds):
            raise UnsupportedStructure(
                "a ring triple bond (cycloalkyne) alongside a selone is "
                "not supported yet -- only a ring double bond is in scope "
                "for this first pass (see P-31.1.3)"
            )
        ring_selones = {s for s in selones if next(iter(graph[s])) in ring_atoms}
        if not bonds and not ring_selones and len(selones) == 1:
            if stereo is not None:
                raise UnsupportedStructure(
                    "a stereocenter on a substituent branch rather than "
                    "the ring itself is not supported yet (see P-92)"
                )
            return _name_ring_substituent_chain_selone(mol, selones)
        if not bonds and ring_selones and ring_selones != selones:
            if stereo is not None:
                raise UnsupportedStructure(
                    "a stereocenter alongside a ring-vs-chain selone "
                    "comparison is not supported yet (see P-92)"
                )
            return _name_ring_with_selone_chain_selone(mol, selones)
        return _name_cyclic_selone(mol, selones, stereo, bonds)
    return _name_von_baeyer_or_spiro_selone(mol, selones, stereo, bonds)


def _name_von_baeyer_or_spiro_selone(mol, selones, stereo, bonds):
    """P-23.2.1/P-24.2.1's von Baeyer bicyclic/polycyclic/monospiro
    numbering extended with a single selone (C=Se) suffix, via
    `_polycyclic_suffix.name_von_baeyer_suffix`/`name_monospiro_suffix`'s
    `elide_e=False` ('selone' begins with a consonant, P-16.3.3, same as
    `_thiol.py`'s 'thiol' -- e.g. 'bicyclo[2.2.1]heptane-2-selone',
    PubChem CID 175998577). Mirrors `_thiol.py`'s own bicyclic/polycyclic-
    before-spiro dispatch order and restrictions: exactly one selone on
    the ring system itself, no ring unsaturation, no specified
    stereocenter."""
    if len(selones) != 1:
        raise UnsupportedStructure(
            "more than one selone on a von Baeyer bicyclic/polycyclic or "
            "monospiro ring system is not supported yet"
        )
    if bonds:
        raise UnsupportedStructure(
            "an unsaturated von Baeyer bicyclic/polycyclic or monospiro "
            "ring system is not supported yet (see P-31.1.4/P-31.1.5)"
        )
    if stereo is not None:
        raise UnsupportedStructure(
            "a specified stereocenter alongside a von Baeyer bicyclic/"
            "polycyclic or monospiro selone is not supported yet (see P-92)"
        )

    (selenium,) = selones
    graph = adjacency(mol)
    (selone_carbon,) = graph[selenium]

    bicyclic_core = find_bicyclic_core(mol)
    polycyclic_core = None
    von_baeyer_ring_count = None
    if bicyclic_core is None:
        for candidate_ring_count in (3, 4, 5, 6):
            polycyclic_core = find_polycyclic_core(mol, candidate_ring_count)
            if polycyclic_core is not None:
                von_baeyer_ring_count = candidate_ring_count
                break
    if bicyclic_core is not None or polycyclic_core is not None:
        return name_von_baeyer_suffix(
            mol,
            selone_carbon,
            selones,
            "selone",
            "selone",
            bicyclic_core,
            polycyclic_core,
            von_baeyer_ring_count,
            elide_e=False,
        )

    spiro_atom = find_monospiro_atom(mol)
    if spiro_atom is not None:
        return name_monospiro_suffix(mol, selone_carbon, selones, "selone", "selone", spiro_atom, elide_e=False)

    raise UnsupportedStructure(
        "polycyclic and fused-ring selones are not supported yet (P-23/"
        "P-25 numbering integration with a suffix group is future work)"
    )
