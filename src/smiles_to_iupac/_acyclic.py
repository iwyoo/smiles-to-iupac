"""Naming of acyclic saturated hydrocarbons (alkanes), per the IUPAC 2013
Recommendations ("the Blue Book"):

- P-44.3 (Chapter P-4, https://iupac.qmul.ac.uk/BlueBook/PDF/P4.pdf): the
  principal chain is the one with the greater number of skeletal atoms.
- P-45.2 (same chapter): remaining ties are broken, in order, by (1) the
  maximum number of substituent prefixes, (2) the lowest locant set for those
  prefixes, and (3) the lowest locants in the prefixes' order of citation.
- P-14.3.5 / P-14.4 / P-14.5 (Chapter P-1, https://iupac.qmul.ac.uk/BlueBook/PDF/P1.pdf):
  lowest-locant-set comparison, numbering, and alphanumerical order of prefixes.
- P-29.3.2.1 (Chapter P-2): unbranched substituent groups (methyl, ethyl, propyl, ...).
- P-29.4 / P-46 (Chapter P-2, P-4): branched ("compound") substituent groups,
  e.g. `(1-methylpropyl)` for a sec-butyl-like branch — see `_substituents.py`.
- P-35.2.1 (Chapter P-3): halogen substituents (fluoro, chloro, bromo, iodo)
  are never skeletal atoms (P-44.3), so the principal chain is found over
  carbon-carbon connectivity only (`carbon_adjacency`, see `_common.py`)
  while substituent detection still uses the full atom graph.

`winning_chain_from_carbon_graph` additionally exposes the winning chain
itself (not just its name) for stereodescriptor locant lookups --
`name_from_carbon_graph` delegates to it unchanged, so `_ether.py`/
`_peroxide.py`/`_nitro.py` (its other callers) see no behavior change.
`name_acyclic_alkane` uses it directly to apply the same P-91.3/P-92
stereodescriptor treatment to a halogenated chain's own stereocenters,
mirroring `_ether.py`'s pattern.
"""

from ._common import (
    UnsupportedStructure,
    adjacency,
    carbon_adjacency,
    group_substituents,
    halogen_substituents,
    longest_chains,
    non_single_bonds,
    specified_stereocenters,
    substituent_locant_set_and_citation,
    validate_atoms_and_bonds,
)
from ._numerals import alkane_name
from ._substituents import format_substituent_prefixes, name_branch


def _substituents_for_chain(graph, chain, halogens, mol=None):
    """Return {position (1-based) -> [(name, is_compound), ...]} for a
    candidate chain."""
    chain_set = set(chain)
    substituents = {}
    for position, atom in enumerate(chain, start=1):
        branch_roots = [n for n in graph[atom] if n not in chain_set]
        if not branch_roots:
            continue
        substituents[position] = [name_branch(graph, root, atom, halogens, mol=mol) for root in branch_roots]
    return substituents


def _name_from_substituents(chain_length, grouped):
    if chain_length == 1 and grouped:
        # P-14.3.4.2(a): every substituent's locant on a mononuclear parent
        # hydride is always '1' and is never cited, however many there are
        # (e.g. 'chloromethane (PIN)' for CH3Cl, 'dichlorosilane' for SiH2Cl2).
        return format_substituent_prefixes(grouped, omit_locants=True) + alkane_name(chain_length)
    total_subs = sum(len(info["locants"]) for info in grouped.values())
    if chain_length == 2 and total_subs == 1:
        (name,) = grouped
        info = grouped[name]
        if not info["compound"]:
            # P-14.3.4.2(b): a homogeneous two-carbon chain bearing exactly
            # one substituent has only one possible structure regardless of
            # numbering direction, so the locant is omittable, e.g.
            # 'chloroethane', analogous to 'ethanol (PIN)' for CH3-CH2-OH.
            return name + alkane_name(chain_length)
    prefix = format_substituent_prefixes(grouped)
    return prefix + alkane_name(chain_length)


def _candidate_key(chain_length, substituents):
    """Sort key implementing P-45.2.1-P-45.2.3, most-preferred first."""
    grouped = group_substituents(substituents)
    locant_set, total_count, citation_locants = substituent_locant_set_and_citation(grouped)
    name = _name_from_substituents(chain_length, grouped)
    # Higher substituent count and lower locants are preferred, so negate the count
    # to sort every field in ascending "most preferred first" order.
    return (-total_count, locant_set, citation_locants, name), name


def longest_chain_length(carbon_graph) -> int:
    """P-44.3's own top-level criterion (greater number of skeletal atoms in
    the chain) applied to a single candidate side, before any P-45.2
    locant-set tie-break: the length of the longest simple path in
    `carbon_graph`. Used by `_ether.py`/`_peroxide.py`/`_sulfide.py`/
    `_selenide.py`/`_telluride.py` to compare two candidate parent sides
    that tie in total skeletal-atom count but not necessarily in the length
    of the chain each would actually use as its own parent hydride (e.g. a
    neopentyl arm's 5 carbons max out at chain length 3 around its
    quaternary carbon, while an isopentyl arm's 5 carbons reach chain
    length 4 -- PubChem CID-confirmed `CC(C)(C)COCCC(C)C` ->
    '1-(2,2-dimethylpropoxy)-3-methylbutane' picks the isopentyl side as
    parent on this basis alone, before any locant-set comparison is even
    needed)."""
    return len(longest_chains(carbon_graph)[0])


def _best_candidate(full_graph, carbon_graph, terminals, mol=None):
    """Shared search behind `winning_chain_from_carbon_graph` and
    `winning_chain_with_key`: every candidate chain/direction's P-45.2 sort
    key, alongside the winning chain and name."""
    chains = longest_chains(carbon_graph)
    chain_length = len(chains[0])

    best_key = None
    best_chain = None
    best_name = None
    for chain in chains:
        for candidate in (chain, list(reversed(chain))):
            substituents = _substituents_for_chain(full_graph, candidate, terminals, mol=mol)
            key, name = _candidate_key(chain_length, substituents)
            if best_key is None or key < best_key:
                best_key, best_chain, best_name = key, candidate, name

    return best_key, best_chain, best_name


def winning_chain_from_carbon_graph(full_graph, carbon_graph, terminals, mol=None):
    """Same P-44.3/P-45.2 tie-break as `name_from_carbon_graph` below, but
    also returns the winning candidate chain itself (root-to-tip, in the
    direction that won), not just its name -- used by a caller (e.g.
    `_acetal.py`) that additionally needs to locate a specific atom's
    position on that chain, such as for a stereodescriptor's locant
    (P-91.3). Kept as the single source of truth so `name_from_carbon_graph`
    and any such caller can never disagree about which chain was chosen."""
    _, chain, name = _best_candidate(full_graph, carbon_graph, terminals, mol=mol)
    return chain, name


def winning_chain_with_key(full_graph, carbon_graph, terminals, mol=None):
    """Same as `winning_chain_from_carbon_graph`, but also exposes the
    P-45.2 sort key (lowest locant set, then alphanumerical order) used to
    pick it. A caller comparing two *different* candidate parent sides of
    equal `longest_chain_length` (e.g. an isobutyl arm vs. a tert-butyl arm,
    both built on a 3-long chain) can compare their keys directly to decide
    which side is senior, the same way this function's own inner loop
    already decides between numbering directions of one fixed chain."""
    return _best_candidate(full_graph, carbon_graph, terminals, mol=mol)


def name_from_carbon_graph(full_graph, carbon_graph, terminals, mol=None) -> str:
    """Name the acyclic saturated skeleton given by `carbon_graph` (P-44.3
    chain search), with `full_graph` used for substituent detection and
    `terminals` ({atom_idx -> prefix name}) naming any leaf substituent
    that's excluded from the chain search and never recursed into — the
    same role `halogen_substituents` plays for halogens (P-35.2.1), reused
    by `_ether.py` for an ether oxygen's precomputed 'alkoxy' prefix."""
    _, name = winning_chain_from_carbon_graph(full_graph, carbon_graph, terminals, mol=mol)
    return name


def name_acyclic_alkane(mol) -> str:
    validate_atoms_and_bonds(mol)
    if non_single_bonds(mol):
        raise UnsupportedStructure(
            "unsaturation is not supported by this module (see "
            "smiles_to_iupac._unsaturated for alkenes/alkynes)"
        )
    if mol.GetRingInfo().NumRings() > 0:
        raise UnsupportedStructure(
            "rings are not supported by this module (see smiles_to_iupac._cyclic)"
        )

    full_graph = adjacency(mol)
    chain, name = winning_chain_from_carbon_graph(full_graph, carbon_adjacency(mol), halogen_substituents(mol), mol=mol)

    stereo = specified_stereocenters(mol)
    if stereo is None:
        return name

    # P-92: every specified stereocenter must lie on the winning principal
    # chain; one on a substituent branch is out of scope, mirroring
    # `_ether.py`'s identical treatment.
    position_of = {atom: i + 1 for i, atom in enumerate(chain)}
    if any(atom not in position_of for atom, _ in stereo):
        raise UnsupportedStructure(
            "a stereocenter on a substituent branch rather than the "
            "principal chain is not supported yet (see P-92)"
        )
    labels = sorted((position_of[atom], code) for atom, code in stereo)
    prefix = ",".join(f"{locant}{code}" for locant, code in labels)
    return f"({prefix})-{name}"
