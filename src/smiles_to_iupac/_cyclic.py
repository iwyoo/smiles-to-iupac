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
- P-93.5.1.3 / P-91.2.1.2.1(b) (Chapter P-9,
  https://iupac.qmul.ac.uk/BlueBook/P9.html), as of
  `tasks/ring-cis-trans-naming.md` (2026-08-25): 'cis'/'trans' is valid
  *general* nomenclature for the relative configuration of a disubstituted
  alicyclic ring -- NOT a PIN construction (a PIN always uses full CIP R/S
  descriptors instead, e.g. '(1R,2R)-1,2-dimethylcyclohexane', per
  P-91.2.1/P-91.3), matching this project's existing precedent of
  supporting Blue Book general nomenclature alongside PINs (e.g. the
  Hantzsch-Widman/biphenyl retained names). Scoped to the simplest case
  only: exactly two ring carbons, adjacent (1,2-), each bearing exactly
  one identical substituent -- symmetric enough that there's no reference-
  substituent selection question (P-93.5.1.3's full r/c/t rules, needed
  for 3+ substituents or a 1,2-pair of *different* substituents, are out
  of scope). The relation is determined by direct 3D geometry (embed a
  conformer, check whether the two substituents sit on the same side of
  the ring's mean plane) rather than by CIP R/S: which R/S combination
  corresponds to cis vs. trans is known to flip between different
  substitution patterns (e.g. 1,2- vs 1,3-), so reusing R/S labels here
  would need separate, unperformed verification for each pattern -- pure
  geometry avoids that entirely. Confirmed against PubChem's own isomeric
  SMILES for cis-/trans-1,2-dimethylcyclohexane (CID 16628/23313).
  A ring stereocenter combination other than exactly this shape (a single
  specified stereocenter, more than two, an asymmetric 1,2-pair, or any
  left unspecified alongside a specified one) raises `UnsupportedStructure`
  explicitly rather than silently dropping the stereochemistry, mirroring
  `_alcohol.py`'s R/S and `_unsaturated.py`'s E/Z handling.

Fused, bridged, and spiro ring systems are out of scope for this module and
raise NotImplementedError.
"""

from rdkit import Chem
from rdkit.Chem import AllChem

from ._common import (
    UnsupportedStructure,
    adjacency,
    halogen_substituents,
    lowest_locant_set,
    non_single_bonds,
    validate_atoms_and_bonds,
)
from ._numerals import alkane_name
from ._substituents import alpha_sort_key, format_substituent_prefixes, name_branch


def _ring_cycle(graph, ring_atoms):
    """Order a simple ring's atoms into a cyclic sequence by walking its bonds."""
    ring_set = set(ring_atoms)
    order = [ring_atoms[0]]
    previous = None
    while len(order) < len(ring_atoms):
        current = order[-1]
        next_atom = next(n for n in graph[current] if n in ring_set and n != previous)
        order.append(next_atom)
        previous = current
    return order


def _substituents_for_ring(graph, ring_order, halogens):
    ring_set = set(ring_order)
    substituents = {}
    for position, atom in enumerate(ring_order, start=1):
        branch_roots = [n for n in graph[atom] if n not in ring_set]
        if not branch_roots:
            continue
        substituents[position] = [name_branch(graph, root, atom, halogens) for root in branch_roots]
    return substituents


def _group(substituents):
    grouped = {}
    for position, entries in substituents.items():
        for name, is_compound in entries:
            info = grouped.setdefault(name, {"locants": [], "compound": is_compound})
            info["locants"].append(position)
    return grouped


def _name_from_substituents(ring_size, grouped):
    parent = "cyclo" + alkane_name(ring_size)
    total_count = sum(len(info["locants"]) for info in grouped.values())
    if total_count == 0:
        return parent
    if total_count == 1:
        # P-14.3.3: the locant is not essential on an otherwise unsubstituted ring.
        (name,) = grouped
        display_name = f"({name})" if grouped[name]["compound"] else name
        return f"{display_name}{parent}"
    prefix = format_substituent_prefixes(grouped)
    return prefix + parent


def _candidate_key(ring_size, substituents):
    """Sort key implementing P-45.2.2/P-45.2.3, most-preferred first (the
    substituent count is fixed for a given ring, so unlike the acyclic case
    there is no P-45.2.1 dimension to break ties on)."""
    grouped = _group(substituents)
    locant_set = lowest_locant_set(loc for info in grouped.values() for loc in info["locants"])
    citation_locants = tuple(
        loc
        for name in sorted(grouped, key=alpha_sort_key)
        for loc in sorted(grouped[name]["locants"])
    )
    name = _name_from_substituents(ring_size, grouped)
    return locant_set, citation_locants, name


def _ring_stereocenters(mol, ring_set):
    """If `mol` has stereo elements that are exactly two specified
    tetrahedral atom stereocenters (and nothing else at all -- no third
    element, no unspecified one alongside them), both on ring atoms,
    return the two atom indices. If there's no specified stereo element at
    all, return None (the caller proceeds exactly as before, no cis-/
    trans- prefix). Any other combination raises `UnsupportedStructure`
    explicitly (see module docstring)."""
    elements = Chem.FindPotentialStereo(mol)
    specified = [e for e in elements if e.specified == Chem.StereoSpecified.Specified]
    if not specified:
        return None
    if (
        len(specified) != 2
        or len(elements) != 2
        or any(e.type != Chem.StereoType.Atom_Tetrahedral for e in specified)
    ):
        raise UnsupportedStructure(
            "ring stereochemistry beyond exactly two specified tetrahedral "
            "ring stereocenters is not supported yet (see P-93.5.1.3)"
        )
    atoms = [e.centeredOn for e in specified]
    if not all(a in ring_set for a in atoms):
        raise UnsupportedStructure(
            "a specified stereocenter not on the ring itself is not "
            "supported yet"
        )
    return atoms


def _ring_normal(coords):
    """Newell's method: unit polygon normal for a (possibly non-planar)
    ring, given `coords` as a list of (x, y, z) tuples in ring order."""
    nx = ny = nz = 0.0
    n = len(coords)
    for i in range(n):
        x1, y1, z1 = coords[i]
        x2, y2, z2 = coords[(i + 1) % n]
        nx += (y1 - y2) * (z1 + z2)
        ny += (z1 - z2) * (x1 + x2)
        nz += (x1 - x2) * (y1 + y2)
    length = (nx * nx + ny * ny + nz * nz) ** 0.5
    return nx / length, ny / length, nz / length


def _substituents_same_face(mol, ring_atoms, ring_atom_a, sub_atom_a, ring_atom_b, sub_atom_b):
    """Embed a 3D conformer and check whether the two substituent atoms
    sit on the same side of the ring's mean plane (see module docstring
    for why this is done geometrically rather than via CIP R/S)."""
    mol_h = Chem.AddHs(mol)
    cids = AllChem.EmbedMultipleConfs(mol_h, numConfs=5, randomSeed=0xC0FFEE)
    if not cids:
        raise UnsupportedStructure(
            "could not generate a 3D conformer to determine cis/trans (see "
            "P-93.5.1.3)"
        )
    props = AllChem.MMFFGetMoleculeProperties(mol_h)
    best_cid, best_energy = None, None
    for cid in cids:
        ff = AllChem.MMFFGetMoleculeForceField(mol_h, props, confId=cid)
        ff.Minimize()
        energy = ff.CalcEnergy()
        if best_energy is None or energy < best_energy:
            best_energy, best_cid = energy, cid
    conf = mol_h.GetConformer(best_cid)

    coords = [tuple(conf.GetAtomPosition(a)) for a in ring_atoms]
    nx, ny, nz = _ring_normal(coords)

    pos_a = conf.GetAtomPosition(ring_atom_a)
    pos_sub_a = conf.GetAtomPosition(sub_atom_a)
    pos_b = conf.GetAtomPosition(ring_atom_b)
    pos_sub_b = conf.GetAtomPosition(sub_atom_b)
    side_a = (pos_sub_a.x - pos_a.x) * nx + (pos_sub_a.y - pos_a.y) * ny + (pos_sub_a.z - pos_a.z) * nz
    side_b = (pos_sub_b.x - pos_b.x) * nx + (pos_sub_b.y - pos_b.y) * ny + (pos_sub_b.z - pos_b.z) * nz
    return (side_a > 0) == (side_b > 0)


def _cis_trans_prefix(mol, graph, ring_atoms, ring_set, halogens, stereo_atoms):
    """(prefix, is_compound) if `stereo_atoms` is exactly the minimal
    symmetric 1,2-disubstituted shape this module supports (see module
    docstring), else raises `UnsupportedStructure`."""
    atom_a, atom_b = stereo_atoms
    if atom_b not in graph[atom_a]:
        raise UnsupportedStructure(
            "specified ring stereocenters that aren't adjacent (1,2-) to "
            "each other are not supported yet (see P-93.5.1.3)"
        )
    subs_a = [n for n in graph[atom_a] if n not in ring_set]
    subs_b = [n for n in graph[atom_b] if n not in ring_set]
    if len(subs_a) != 1 or len(subs_b) != 1:
        raise UnsupportedStructure(
            "a ring stereocenter without exactly one substituent is not "
            "supported yet (see P-93.5.1.3)"
        )
    (sub_a,), (sub_b,) = subs_a, subs_b
    name_a, compound_a = name_branch(graph, sub_a, atom_a, halogens)
    name_b, compound_b = name_branch(graph, sub_b, atom_b, halogens)
    if (name_a, compound_a) != (name_b, compound_b):
        raise UnsupportedStructure(
            "an asymmetric 1,2-disubstituted pair (two different "
            "substituents) is not supported yet -- P-93.5.1.3's reference-"
            "substituent selection rules are needed (see P-93.5.1.3)"
        )
    same_face = _substituents_same_face(mol, ring_atoms, atom_a, sub_a, atom_b, sub_b)
    return ("cis-" if same_face else "trans-"), False


def name_cycloalkane(mol) -> str:
    validate_atoms_and_bonds(mol)
    if non_single_bonds(mol):
        raise UnsupportedStructure(
            "unsaturated rings are not supported yet (see P-31.1.3, "
            "cycloalkenes and cycloalkynes)"
        )
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
    ring_order = _ring_cycle(graph, ring_atoms)
    ring_size = len(ring_order)

    stereo_atoms = _ring_stereocenters(mol, ring_set)
    cis_trans = None
    if stereo_atoms is not None:
        cis_trans = _cis_trans_prefix(mol, graph, ring_atoms, ring_set, halogens, stereo_atoms)

    best_key = None
    best_name = None
    for start in range(ring_size):
        rotated = ring_order[start:] + ring_order[:start]
        for candidate in (rotated, list(reversed(rotated))):
            substituents = _substituents_for_ring(graph, candidate, halogens)
            key = _candidate_key(ring_size, substituents)
            if best_key is None or key < best_key:
                best_key, best_name = key, key[-1]

    if cis_trans is not None:
        prefix, _ = cis_trans
        return f"{prefix}{best_name}"
    return best_name
