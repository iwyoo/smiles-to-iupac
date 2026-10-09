"""Naming of simple monocyclic saturated hydrocarbons (cycloalkanes), per the
IUPAC 2013 Recommendations ("the Blue Book"):

- P-22.1.1 (Chapter P-2, https://iupac.qmul.ac.uk/BlueBook/PDF/P2.pdf): the name
  is formed by attaching the nondetachable prefix 'cyclo' to the name of the
  acyclic saturated hydrocarbon with the same number of carbon atoms.
- P-14.4 / P-45.2 (Chapter P-1 / P-4): the numbering direction and starting atom
  are chosen, as for a chain, to give the lowest locant set to substituents,
  then the lowest locants in their order of citation.
- P-14.3.3 (Chapter P-1): a locant is cited only when essential to the
  structure; for a single substituent on an otherwise unsubstituted ring, every
  ring atom is equivalent before substitution, so the locant is not essential
  and is omitted (e.g. 'methylcyclohexane', not '1-methylcyclohexane').
- P-29.4 / P-46 (Chapter P-2, P-4): branched ("compound") substituent groups,
  e.g. `(1-methylpropyl)` for a sec-butyl-like ring substituent — see
  `_substituents.py`.
- P-35.2.1 (Chapter P-3): halogen substituents (fluoro, chloro, bromo, iodo)
  hang off a ring atom the same way any other substituent does; a ring's own
  atom sequence needs no carbon-only filtering here since RDKit's ring
  perception (`GetRingInfo`) never includes a monovalent atom in a ring.
- P-93.5.1.1: ring stereocenters are cited by CIP descriptors in a PIN; the general ring namer assigns them.
  A single specified stereocenter on the ring's own sole substituent branch gets `_alcohol.py`'s
  `ring_branch_stereo_display` treatment, a bracketed "[(<locant><R/S>)-<name>]" branch display (P-92).

- P-29.2 / P-31.1.4.3.4 (Chapter P-2/P-3): an exocyclic multiple bond makes the
  substituent an ylidene/enyl/ynyl prefix ('ethylidene', 'prop-1-en-1-yl',
  ...) named by `name_branch` and numbered like any detachable prefix.

Fused, bridged, and spiro ring systems are out of scope for this module and
raise NotImplementedError.
"""

from rdkit import Chem

from ._multiplicative_text import enclose
from ._common import (
    UnsupportedStructure,
    adjacency,
    group_substituents,
    halogen_substituents,
    non_single_bonds,
    ring_cycle,
    specified_stereocenters,
    substituent_locant_set_and_citation,
    validate_atoms_and_bonds,
)
from ._locant_omission import omits_all_locants
from ._numerals import alkane_name
from ._substituents import format_substituent_prefixes, ring_branch_stereo_display, substituents_for_ring




def _name_from_substituents(ring_size, grouped, omit_locants=False):
    parent = "cyclo" + alkane_name(ring_size)
    total_count = sum(len(info["locants"]) for info in grouped.values())
    if total_count == 0:
        return parent
    if total_count == 1:
        # P-14.3.3: the locant is not essential on an otherwise unsubstituted ring.
        (name,) = grouped
        display_name = enclose(name) if grouped[name]["compound"] else name
        return f"{display_name}{parent}"
    prefix = format_substituent_prefixes(grouped, omit_all=omit_locants)
    return prefix + parent


def _candidate_key(ring_size, substituents, omit=None):
    """Sort key implementing P-45.2.2/P-45.2.3, most-preferred first (the
    substituent count is fixed for a given ring, so unlike the acyclic case
    there is no P-45.2.1 dimension to break ties on)."""
    grouped = group_substituents(substituents)
    locant_set, _, citation_locants = substituent_locant_set_and_citation(grouped)
    name = _name_from_substituents(ring_size, grouped, omit is not None and omit(grouped))
    return locant_set, citation_locants, name


def name_cycloalkane(mol) -> str:
    validate_atoms_and_bonds(mol)
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() != 1:
        raise UnsupportedStructure(
            "only simple monocyclic ring systems are supported so far (see "
            "P-23/P-24/P-25 for polyalicyclic, spiro, and fused ring systems)"
        )

    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    ring_atoms = list(ring_info.AtomRings()[0])
    ring_set = set(ring_atoms)
    # `core.py` routes a ring-internal multiple bond to `_cyclic_unsaturated.py`
    # first, so one reaching here is an unsupported shape.
    if any(a in ring_set and b in ring_set for a, b, _ in non_single_bonds(mol)):
        raise UnsupportedStructure(
            "unsaturated rings are not supported yet (see P-31.1.3, "
            "cycloalkenes and cycloalkynes)"
        )
    ring_order = ring_cycle(graph, ring_atoms)
    ring_size = len(ring_order)

    stereo = specified_stereocenters(mol)
    branch_stereo = None
    if stereo is not None and any(atom not in ring_order for atom, _ in stereo):
        branch_stereo = ring_branch_stereo_display(graph, ring_order, frozenset(), stereo, halogens, mol=mol)
        if branch_stereo is None:
            raise UnsupportedStructure(
                "a stereocenter on a substituent branch rather than the "
                "ring itself is not supported yet (see P-92)"
            )

    if branch_stereo is None and any(
        e.specified == Chem.StereoSpecified.Specified for e in Chem.FindPotentialStereo(mol)
    ):
        raise UnsupportedStructure("ring stereocenters are cited by CIP descriptors, which the general ring namer assigns")

    best_key = None
    best_name = None
    for start in range(ring_size):
        rotated = ring_order[start:] + ring_order[:start]
        for candidate in (rotated, list(reversed(rotated))):
            position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
            substituents = substituents_for_ring(graph, candidate, halogens, mol=mol, unsaturated=True)
            if branch_stereo is not None:
                branch_ring_atom, display = branch_stereo
                substituents[position_of[branch_ring_atom]] = [(display, False)]
            key = _candidate_key(
                ring_size, substituents, lambda grouped, ring=candidate: omits_all_locants(mol, ring, grouped)
            )
            if best_key is None or key < best_key:
                best_key, best_name = key, key[-1]

    return best_name
