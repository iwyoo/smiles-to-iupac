"""Naming of thioate anions (R-CO-S(-) <-> R-CS-O(-)), the chalcogen
analogue of a carboxylate anion with one oxygen replaced by sulfur,
restricted to a single such group on an acyclic carbon chain, per the
IUPAC 2013 Recommendations ("the Blue Book"):

- P-72.2.2.2.1.1 (Chapter P-7, https://iupac.qmul.ac.uk/BlueBook/P7.html):
  the anion formed by removing a hydron from either chalcogen of a thioic
  acid (`_thioic_acid.py`, P-65.1.5) is named by replacing 'thioic acid'
  with 'thioate'. Unlike the neutral acid, whose two tautomers are
  distinct, fully-determined structures each needing its own O-/S- letter
  locant, the anion's negative charge is delocalized across both
  chalcogens -- the primary source draws this explicitly with a
  resonance arrow and gives both drawn forms the *same* name, with no
  letter at all: 'CH3-CH2-CO-S(-) <-> CH3-CH2-CS-O(-)' -> 'propanethioate
  (PIN)', 'CH3-CO-S(-) <-> CH3-CS-O(-)' -> 'ethanethioate (PIN)'. This
  project's own PubChem cross-check confirms the same collapse at the
  structure level, not just the name: 'CCC(=O)[S-]' and 'CCC(=S)[O-]'
  both resolve to the identical CID (22717040) and name "propanethioate"
  ('CC(=O)[S-]'/'CC(=S)[O-]' likewise both CID 3815167, "ethanethioate")
  -- so, unlike `_thioic_acid.py`, this module never distinguishes which
  chalcogen is drawn double- vs. single-bonded.
- Otherwise mirrors `_carboxylate.py` exactly: acyclic only, the
  thioate carbon is always chain terminus C1 (P-14.3.3, its own locant
  never cited), 'ene'/'yne' unsaturation and halogen substituents are
  supported the same way, and only exactly one thioate group is in
  scope (no 'bis(thioate)').
- Formic acid's analogue (chain length 1, R = H) is named with the
  'methane' stem, mirroring `_thioic_acid.py`'s identical treatment.

Explicitly out of scope (raise `UnsupportedStructure`):
- A ring anywhere in the molecule (mirrors `_carboxylate.py`'s acyclic-
  only scope).
- More than one thioate group, or any oxygen/sulfur that isn't part of
  the single thioate's carbonyl/anion pair.
- Any charged or radical atom other than the single anionic chalcogen.
- The selenium analogue is handled separately by `_selenoate.py`;
  tellurium (telluroate) is a separate follow-up module, mirroring how
  `_selenoic_acid.py`/`_telluroic_acid.py` followed `_thioic_acid.py`.
- Any other heteroatom.

P-91.3/P-92: a molecule
with one or more *specified* tetrahedral stereocenters -- every one on the
principal chain itself, no unspecified one alongside them, and no
C=C/C#N double-bond E/Z element -- gets a "(<locant><R/S>,...)-" prefix,
ascending locant order, same pattern as `_carboxylate.py` (the thioate
carbon is always C1, so it never affects numbering). The thioate carbon
itself is never a potential stereocenter, confirmed via RDKit
`FindPotentialStereo` on both drawn tautomers
(`CC[C@@H](C)C(=O)[S-]`/`CC[C@@H](C)C(=S)[O-]`).
"""

from rdkit import Chem

from ._common import (
    elides_before,
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
    longest_branched_chain,
    longest_chains,
    lowest_locant_set,
    multiplied_word,
    non_single_bonds,
    ring_chain_attachment,
    ring_chain_attachment_with_halogens,
    specified_stereocenters,
    substituent_locant_set_and_citation,
)
from ._numerals import alkane_name
from ._substituents import (
    format_substituent_prefixes,
    name_branch,
    plain_alkyl_ring_substituents,
)

_ALLOWED_ATOMIC_NUMS = {6, 8, 16, *HALOGEN_PREFIXES}
_CHALCOGENS = (8, 16)


def _thioate_matches(mol):
    """Carbons shaped like a thioate group: a neutral, monovalent,
    double-bonded O or S plus a formal-charge -1, monovalent, single-
    bonded chalcogen of the *other* element. Returns a list of
    (carbon, carbonyl_atom, anion_atom) triples."""
    matches = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6:
            continue
        chalcogens = [n for n in atom.GetNeighbors() if n.GetAtomicNum() in _CHALCOGENS]
        if len(chalcogens) != 2:
            continue
        carbonyls = [
            o
            for o in chalcogens
            if o.GetDegree() == 1
            and o.GetFormalCharge() == 0
            and mol.GetBondBetweenAtoms(atom.GetIdx(), o.GetIdx()).GetBondTypeAsDouble() == 2.0
        ]
        anions = [
            o
            for o in chalcogens
            if o.GetDegree() == 1
            and o.GetFormalCharge() == -1
            and mol.GetBondBetweenAtoms(atom.GetIdx(), o.GetIdx()).GetBondTypeAsDouble() == 1.0
        ]
        if len(carbonyls) != 1 or len(anions) != 1:
            continue
        if {carbonyls[0].GetAtomicNum(), anions[0].GetAtomicNum()} != {8, 16}:
            continue
        matches.append((atom, carbonyls[0], anions[0]))
    return matches


def has_thioate_shape(mol) -> bool:
    return bool(_thioate_matches(mol))


def _find_thioate_group(mol):
    matches = _thioate_matches(mol)
    if len(matches) != 1:
        raise UnsupportedStructure(
            "exactly one thioate (-COS-/-CSO-) group is required; zero or "
            "multiple such groups are not supported yet (P-72.2.2.2.1.1)"
        )
    thioate_carbon, carbonyl_atom, anion_atom = matches[0]
    total_chalcogens = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() in _CHALCOGENS)
    if total_chalcogens != 2:
        raise UnsupportedStructure(
            "a chalcogen outside the single thioate group's carbonyl/anion "
            "pair is out of scope for this module"
        )
    carbon_neighbors = [n for n in thioate_carbon.GetNeighbors() if n.GetAtomicNum() == 6]
    if len(carbon_neighbors) > 1:
        raise UnsupportedStructure(
            "a thioate carbon with more than one carbon neighbor is not a "
            "valid thioate group"
        )
    return thioate_carbon, carbonyl_atom, anion_atom


def _suffix_body(ene_locants, yne_locants):
    segments = []
    if ene_locants:
        segments.append((sorted(ene_locants), multiplied_word(len(ene_locants), "ene")))
    if yne_locants:
        segments.append((sorted(yne_locants), multiplied_word(len(yne_locants), "yne")))

    words = [word for _, word in segments] + ["thioate"]
    for i in range(len(words) - 1):
        if words[i].endswith("e") and elides_before(words[i + 1]):
            words[i] = words[i][:-1]

    if segments:
        locant_parts = [
            f"{','.join(str(loc) for loc in locants)}-{word}"
            for (locants, _), word in zip(segments, words[:-1])
        ]
        body = "-".join(locant_parts) + words[-1]
    else:
        body = words[-1]
    elide_stem = words[0][0] in "aeiouy"
    return body, elide_stem


def _name_from_substituents(chain_length, ene_locants, yne_locants, grouped):
    has_unsaturation = bool(ene_locants or yne_locants)
    prefix = format_substituent_prefixes(grouped)
    if has_unsaturation:
        stem = alkane_name(chain_length)[:-3]
        needs_stem_a = (len(ene_locants) >= 2) if ene_locants else (len(yne_locants) >= 2)
    else:
        stem = alkane_name(chain_length)
        needs_stem_a = False

    body, elide_stem = _suffix_body(ene_locants, yne_locants)
    if not has_unsaturation and elide_stem:
        stem = stem[:-1]
    separator = "-" if (ene_locants or yne_locants) else ""
    return prefix + stem + ("a" if needs_stem_a else "") + separator + body


def _candidate_key(chain_length, ene_locants, yne_locants, substituents):
    grouped = group_substituents(substituents)
    locant_set, total_count, citation_locants = substituent_locant_set_and_citation(grouped)
    combined_locant_set = lowest_locant_set(ene_locants + yne_locants)
    ene_locant_set = lowest_locant_set(ene_locants)
    name = _name_from_substituents(chain_length, ene_locants, yne_locants, grouped)
    return (
        (
            combined_locant_set,
            ene_locant_set,
            -total_count,
            locant_set,
            citation_locants,
            name,
        ),
        name,
    )


def _substituents_for_chain(graph, chain, halogens, excluded_atoms, mol=None):
    chain_set = set(chain)
    substituents = {}
    for position, atom in enumerate(chain, start=1):
        branch_roots = [n for n in graph[atom] if n not in chain_set and n not in excluded_atoms]
        if branch_roots:
            substituents[position] = [name_branch(graph, root, atom, halogens, mol=mol) for root in branch_roots]
    return substituents


def _name_acyclic_thioate(mol, thioate_carbon_idx, excluded_atoms, bonds, stereo=None):
    """`stereo`: None, or a list of (stereocenter_atom_idx, "R"/"S") from
    `specified_stereocenters` -- if given, only chain candidates that
    include every stereocenter are eligible (P-92: a stereocenter on a
    substituent branch rather than the principal chain is out of scope,
    mirroring `_carboxylate.py`), and the winning candidate's own locants
    are used to format a "(<locant><R/S>,...)-" prefix onto the final
    name."""
    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    chains = longest_chains(carbon_adjacency(mol))
    chain_length = len(chains[0])
    stereo_atoms = [atom for atom, _ in stereo] if stereo is not None else []

    eligible = []
    for chain in chains:
        if thioate_carbon_idx not in chain:
            continue
        if bonds and bond_locants(chain, bonds) is None:
            continue
        chain_set = set(chain)
        if stereo is not None and any(atom not in chain_set for atom in stereo_atoms):
            continue
        eligible.append(chain)
    if not eligible:
        if stereo is not None and any(
            thioate_carbon_idx in c and (not bonds or bond_locants(c, bonds) is not None) for c in chains
        ):
            raise UnsupportedStructure(
                "a stereocenter on a substituent branch rather than the "
                "principal chain is not supported yet (see P-92)"
            )
        raise UnsupportedStructure(
            "the thioate carbon (and/or multiple bonds) does not lie on a "
            "single longest carbon chain"
        )

    best_key = None
    best_name = None
    best_position_of = None
    for chain in eligible:
        for candidate in (chain, list(reversed(chain))):
            if candidate[0] != thioate_carbon_idx:
                continue
            position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
            ene_locants, yne_locants = bond_locants(candidate, bonds) if bonds else ([], [])
            substituents = _substituents_for_chain(graph, candidate, halogens, excluded_atoms, mol=mol)
            key, name = _candidate_key(chain_length, ene_locants, yne_locants, substituents)
            if best_key is None or key < best_key:
                best_key, best_name, best_position_of = key, name, position_of

    if stereo is not None:
        labels = sorted((best_position_of[atom], code) for atom, code in stereo)
        prefix = ",".join(f"{locant}{code}" for locant, code in labels)
        return f"({prefix})-{best_name}"
    return best_name


def _validate_and_prepare_thioate(mol, aromatic_ring_atoms=frozenset()):
    """(thioate_carbon, excluded_atoms, bonds, stereo) after checking the
    molecule fits this module's scope (see module docstring).

    `aromatic_ring_atoms`: atom indices already independently verified (by
    the caller, before this function runs) to form a single plain benzene
    ring with exactly one exocyclic attachment -- exempted from the
    aromatic-atom rejection below so `name_thioate`'s benzene-ring-
    substituent path (see `_name_phenyl_chain_thioate`) can reuse this
    same validation for the rest of the molecule. Empty by default, so
    every other caller's behavior is unchanged. Mirrors
    `_carboxylate.py`'s equivalent aromatic-exemption pattern."""
    thioate_carbon, carbonyl_atom, anion_atom = _find_thioate_group(mol)
    excluded_atoms = {carbonyl_atom.GetIdx(), anion_atom.GetIdx()}
    stereo = specified_stereocenters(mol)

    has_carbon = False
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atomic_num not in _ALLOWED_ATOMIC_NUMS:
            raise UnsupportedStructure(
                "heteroatoms other than the thioate's own chalcogens "
                "(P-72.2.2.2.1.1) and halogen substituents (P-35.2.1) are "
                "not supported yet"
            )
        if atom.GetIdx() in excluded_atoms:
            continue
        if atom.GetFormalCharge() != 0 or atom.GetIsotope() != 0:
            raise UnsupportedStructure(
                "a charged or isotopically modified atom other than the "
                "single thioate anion chalcogen is not supported yet"
            )
        if atomic_num == 6:
            has_carbon = True
            if atom.GetIsAromatic() and atom.GetIdx() not in aromatic_ring_atoms:
                raise UnsupportedStructure(
                    "aromatic rings are out of scope for this module (see "
                    "the separate aromatic-ring module)"
                )
        elif atomic_num not in _CHALCOGENS and atom.GetDegree() != 1:
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

    all_non_single = [
        b
        for b in non_single_bonds(mol)
        if b[0] not in excluded_atoms
        and b[1] not in excluded_atoms
        and (b[0] not in aromatic_ring_atoms or b[1] not in aromatic_ring_atoms)
    ]
    bonds = [b for b in all_non_single if b[2] in (ENE_BOND_ORDER, YNE_BOND_ORDER)]
    if len(bonds) != len(all_non_single):
        raise UnsupportedStructure(
            "a bond order other than single, double, or triple is not "
            "supported (see P-31.1.1.1)"
        )
    return thioate_carbon, excluded_atoms, bonds, stereo


def _name_phenyl_chain_thioate(mol, ring_atoms):
    """Name a thioate whose -CO-S(-) lies entirely on a single unbranched
    chain hanging off one atom of an otherwise-plain, unsubstituted
    benzene ring -- e.g. 3-phenylpropanethioate. The ring is cited as a
    'phenyl' substituent prefix (via `name_branch`'s aromatic-ring
    recognition) on the chain, which is the parent hydride, mirroring
    `_carboxylate.py`'s `_name_phenyl_chain_carboxylate`. Narrower than
    the acyclic path above: no chain unsaturation and no specified
    stereocenter."""
    thioate_carbon, excluded_atoms, bonds, stereo = _validate_and_prepare_thioate(
        mol, aromatic_ring_atoms=ring_atoms
    )
    if stereo:
        raise UnsupportedStructure(
            "a specified stereocenter alongside a benzene-ring-substituent "
            "thioate chain is not supported yet"
        )
    non_ring_unsaturation = [b for b in bonds if b[0] not in ring_atoms or b[1] not in ring_atoms]
    if non_ring_unsaturation:
        raise UnsupportedStructure(
            "chain unsaturation alongside a benzene-ring-substituent "
            "thioate chain is not supported yet"
        )

    graph = adjacency(mol)
    halogens = {**halogen_substituents(mol), **plain_alkyl_ring_substituents(mol, graph, ring_atoms)}
    attachment = ring_chain_attachment_with_halogens(graph, ring_atoms, set(), halogens)
    if attachment is None:
        raise UnsupportedStructure(
            "a benzene ring with more than one non-halogen, non-alkyl "
            "exocyclic substituent alongside a chain thioate is not "
            "supported yet"
        )
    chain, branches = longest_branched_chain(graph, thioate_carbon.GetIdx(), ring_atoms, excluded_atoms, halogens=halogen_substituents(mol))
    if len(chain) < 2:
        raise UnsupportedStructure(
            "a thioate group directly attached to the benzene ring uses a "
            "separate construction, out of scope for this "
            "acyclic-chain-parent module"
        )

    chain_length = len(chain)
    substituents = {
        position: [name_branch(graph, root, chain[position - 1], halogens, ring_atoms, mol=mol) for root in roots]
        for position, roots in branches.items()
    }
    grouped = group_substituents(substituents)
    return _name_from_substituents(chain_length, [], [], grouped)


def name_thioate(mol) -> str:
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() == 1:
        ring_atoms = set(ring_info.AtomRings()[0])
        if is_plain_benzene_ring(mol, ring_atoms):
            return _name_phenyl_chain_thioate(mol, ring_atoms)
    if mol.GetRingInfo().NumRings() > 0:
        raise UnsupportedStructure(
            "a thioate group on/in a ring is out of scope for this "
            "acyclic-only module"
        )

    thioate_carbon, excluded_atoms, bonds, stereo = _validate_and_prepare_thioate(mol)
    return _name_acyclic_thioate(mol, thioate_carbon.GetIdx(), excluded_atoms, bonds, stereo)
