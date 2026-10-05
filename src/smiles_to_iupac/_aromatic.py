"""Substituted benzene (P-22.1.2, P-14.3.3): the single aromatic six-membered carbocycle with halogen, alkyl and other
branches, named with the substituent machinery shared by the other ring modules; fused systems are named by
`_fusion_name` and numbered by `_fused_numbering`."""

from ._multiplicative_text import enclose
from ._common import (
    UnsupportedStructure,
    adjacency,
    carbon_adjacency,
    group_substituents,
    halogen_substituents,
    non_single_bonds,
    ring_cycle,
    specified_stereocenters,
    substituent_locant_set_and_citation,
    validate_atoms_and_bonds,
)
from ._substituents import branch_atom_locant, format_substituent_prefixes, name_branch


def find_aromatic_fused_core(mol):
    """The benzene ring of `mol` as (atom_rings, ring_atom_sets), or None when `mol` has any other ring system."""
    atom_rings = mol.GetRingInfo().AtomRings()
    if len(atom_rings) != 1 or len(atom_rings[0]) != 6:
        return None
    if any(mol.GetAtomWithIdx(i).GetAtomicNum() != 6 or not mol.GetAtomWithIdx(i).GetIsAromatic() for i in atom_rings[0]):
        return None
    return atom_rings, [set(atom_rings[0])]


def _validate_aromatic_bonds(mol, ring_atoms=frozenset()):
    """A non-aromatic multiple bond lying entirely on an exocyclic branch
    (both atoms outside `ring_atoms`) is exempted -- P-61.2.3/P-31.1.3.4:
    a cyclic hydrocarbon substituted by an unsaturated chain is named
    with the ring as parent and the chain cited as an ordinary
    'alkenyl'/'alkynyl' substituent prefix (see
    `unbranched_unsaturated_substituent_name`), e.g. 'ethenylbenzene'
    (PIN; 'styrene' is a retained name restricted to general nomenclature
    only, P-31.1.3.4). A bond straddling the ring and a branch, or one
    entirely inside the ring itself (or a second, non-aromatic ring),
    stays rejected exactly as before."""
    for bond in mol.GetBonds():
        if bond.GetBondTypeAsDouble() == 1.0 or bond.GetIsAromatic():
            continue
        a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if a not in ring_atoms and b not in ring_atoms:
            continue
        raise UnsupportedStructure(
            "a non-aromatic multiple bond is not supported within or "
            "attached to an aromatic ring system (see P-61.2.3)"
        )


def _benzene_candidates(graph, ring_atoms):
    cycle = ring_cycle(graph, list(ring_atoms))
    n = len(cycle)
    candidates = []
    for start in range(n):
        rotated = cycle[start:] + cycle[:start]
        for seq in (rotated, list(reversed(rotated))):
            candidates.append({atom: i + 1 for i, atom in enumerate(seq)})
    return candidates


def _branch_name(graph, carbon_graph, root, ring_atom, ring_atoms, halogens, unsaturated_bonds, mol=None):
    return name_branch(graph, root, ring_atom, halogens, mol=mol)


def _ring_substituents(graph, locants, ring_atoms, halogens, mol=None, carbon_graph=None, unsaturated_bonds=()):
    """Like `_substituents.substituents_for_ring`, but `locants` (nonfusion atom
    -> integer position) covers only part of the ring skeleton -- a fusion
    carbon has no free valence and so is never itself a locant -- while
    `ring_atoms` (the full skeleton, fusion included) is what a neighboring
    atom must be excluded from to count as a substituent branch root."""
    substituents = {}
    for atom, position in locants.items():
        branch_roots = [n for n in graph[atom] if n not in ring_atoms]
        if branch_roots:
            substituents[position] = [
                _branch_name(graph, carbon_graph, root, atom, ring_atoms, halogens, unsaturated_bonds, mol=mol)
                for root in branch_roots
            ]
    return substituents


def _candidate_key(
    parent,
    locants,
    ring_atoms,
    graph,
    halogens,
    omit_single_locant,
    stereo_display=None,
    mol=None,
    carbon_graph=None,
    unsaturated_bonds=(),
):
    substituents = _ring_substituents(
        graph, locants, ring_atoms, halogens, mol=mol, carbon_graph=carbon_graph, unsaturated_bonds=unsaturated_bonds
    )
    grouped = group_substituents(substituents)
    locant_set, total_count, citation_locants = substituent_locant_set_and_citation(grouped)
    if not grouped:
        name = parent
    elif omit_single_locant and total_count == 1:
        # P-14.3.3: on benzene, every ring position is equivalent before
        # substitution, so a single substituent's locant is not essential.
        if stereo_display is not None:
            display = stereo_display
        else:
            (only_name,) = grouped
            display = enclose(only_name) if grouped[only_name]["compound"] else only_name
        name = f"{display}{parent}"
    else:
        if stereo_display is not None:
            # Not benzene, so the ring attachment locant is significant and
            # still needs citing; only the substituent's own display text is
            # replaced by the stereo-decorated one (already validated to be
            # the system's sole substituent -- see `_stereo_display`).
            (only_name,) = grouped
            grouped = {stereo_display: {"locants": grouped[only_name]["locants"], "compound": False}}
        name = format_substituent_prefixes(grouped) + parent
    return locant_set, citation_locants, name


def _stereo_display(mol, graph, n, ring_atoms, halogens):
    """If `mol` has one or more specified tetrahedral stereocenters
    (P-92), build the bracketed "[(<locant><R/S>)-<name>]" substituent
    display P-91.3 requires when the stereocenter sits on a substituent
    branch rather than the ring itself (the Blue Book's own worked
    example, since an all-carbon aromatic ring atom is never itself a
    stereocenter). Returns None if there's no specified stereocenter at
    all (the caller proceeds exactly as before). Deliberately narrow: only a plain benzene
    or naphthalene ring (n in (1, 2)) with exactly one substituent,
    carrying exactly one specified stereocenter, is supported; anything
    else raises `UnsupportedStructure`. For naphthalene the ring
    attachment locant is still significant (unlike benzene) and is cited
    by the caller's normal substituent-prefix machinery -- this function
    only builds the substituent's own decorated display text."""
    stereo = specified_stereocenters(mol)
    if stereo is None:
        return None
    if n not in (1, 2):
        raise UnsupportedStructure(
            "a specified stereocenter combined with anything other than a "
            "plain benzene or naphthalene ring parent is not supported yet "
            "(see P-91.3)"
        )
    if len(stereo) != 1:
        raise UnsupportedStructure(
            "more than one specified stereocenter on a substituent branch "
            "is not supported yet (see P-91.3)"
        )
    stereo_atom, r_or_s = stereo[0]
    if stereo_atom in ring_atoms:
        raise UnsupportedStructure(
            "a specified stereocenter on the aromatic ring itself is not "
            "supported (an all-carbon aromatic ring atom can't be a "
            "genuine stereocenter)"
        )
    branch_attachments = [
        (ring_atom, neighbor)
        for ring_atom in ring_atoms
        for neighbor in graph[ring_atom]
        if neighbor not in ring_atoms
    ]
    if len(branch_attachments) != 1:
        raise UnsupportedStructure(
            "a specified stereocenter combined with anything other than "
            "exactly one substituent on the ring is not supported yet "
            "(see P-91.3)"
        )
    ring_atom, branch_root = branch_attachments[0]
    branch_name, branch_compound = name_branch(graph, branch_root, ring_atom, halogens, mol=mol)
    site_locant = branch_atom_locant(graph, branch_root, ring_atom, stereo_atom, halogens, mol=mol)
    descriptor = f"({site_locant}{r_or_s})-{branch_name}"
    return f"[{descriptor}]" if branch_compound else f"({descriptor})"


def name_aromatic_fused(mol, core) -> str:
    validate_atoms_and_bonds(mol)

    atom_rings, ring_atom_sets = core
    n = len(atom_rings)
    ring_atoms = set(ring_atom_sets[0])
    _validate_aromatic_bonds(mol, ring_atoms)

    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    unsaturated_bonds = [b for b in non_single_bonds(mol) if b[0] not in ring_atoms and b[1] not in ring_atoms]
    carbon_graph = carbon_adjacency(mol) if unsaturated_bonds else None

    stereo_display = _stereo_display(mol, graph, n, ring_atoms, halogens)
    if stereo_display is not None and unsaturated_bonds:
        raise UnsupportedStructure(
            "a specified stereocenter combined with an unsaturated "
            "exocyclic substituent is not supported yet"
        )

    best_key = None
    best_name = None
    for locants in _benzene_candidates(graph, ring_atom_sets[0]):
        key = _candidate_key(
            "benzene",
            locants,
            ring_atoms,
            graph,
            halogens,
            True,
            stereo_display,
            mol=mol,
            carbon_graph=carbon_graph,
            unsaturated_bonds=unsaturated_bonds,
        )
        if best_key is None or key < best_key:
            best_key, best_name = key, key[-1]
    return best_name
