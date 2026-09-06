"""Naming of selenoate anions (R-CO-Se(-) <-> R-CSe-O(-)), the selenium
analogue of a thioate anion (`_thioate.py`), restricted to a single such
group on an acyclic carbon chain, per the IUPAC 2013 Recommendations
("the Blue Book"):

- P-72.2.2.2.1.1 (Chapter P-7, https://iupac.qmul.ac.uk/BlueBook/P7.html):
  the same 'ic acid' -> 'ate' replacement rule that produced 'thioate'
  from a thioic acid applies uniformly to O, S, Se, and Te (the rule
  text names all four chalcogens together). Like `_thioate.py`, and
  unlike the neutral selenoic acid's two distinct tautomers (each needing
  its own O-/Se- letter locant, `_selenoic_acid.py`), the anion's charge
  is delocalized across both chalcogens, so both drawn forms get the
  *same* name with no letter at all. Confirmed directly against PubChem:
  'CC(=O)[Se-]' and 'CC(=[Se])[O-]' both resolve to the identical CID
  (136810730) and name "ethaneselenoate" -- the longer-chain
  'propaneselenoate' analogue isn't registered in PubChem (CID 0), so
  this project's own two-tautomer 'ethaneselenoate' collapse is the
  direct verification; everything else here mirrors the
  already-double-verified `_thioate.py` mechanically (same pattern PubChem
  itself uses for the naming rule, per P-72.2.2.2.1.1's own O/S/Se/Te
  wording).
- Otherwise mirrors `_thioate.py`/`_carboxylate.py` exactly: acyclic
  only, the selenoate carbon is always chain terminus C1 (P-14.3.3, its
  own locant never cited), 'ene'/'yne' unsaturation and halogen
  substituents are supported the same way, and only exactly one
  selenoate group is in scope (no 'bis(selenoate)').
- Formic acid's analogue (chain length 1, R = H) is named with the
  'methane' stem, mirroring `_thioate.py`'s identical treatment.

Explicitly out of scope (raise `UnsupportedStructure`):
- A ring anywhere in the molecule (mirrors `_carboxylate.py`'s acyclic-
  only scope).
- More than one selenoate group, or any oxygen/selenium that isn't part
  of the single selenoate's carbonyl/anion pair.
- Any charged or radical atom other than the single anionic chalcogen.
- The tellurium analogue (telluroate) -- a separate follow-up module,
  mirroring how `_telluroic_acid.py` followed `_selenoic_acid.py`.
- Any other heteroatom.

P-91.3/P-92: a molecule
with one or more *specified* tetrahedral stereocenters -- every one on the
principal chain itself, no unspecified one alongside them, and no
C=C/C#N double-bond E/Z element -- gets a "(<locant><R/S>,...)-" prefix,
ascending locant order, same pattern as `_carboxylate.py`/`_thioate.py`
(the selenoate carbon is always C1, so it never affects numbering). The
selenoate carbon itself is never a potential stereocenter, confirmed via
RDKit `FindPotentialStereo` on both drawn tautomers
(`CC[C@@H](C)C(=O)[Se-]`/`CC[C@@H](C)C(=[Se])[O-]`).
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
    longest_branched_chain,
    longest_chains,
    lowest_locant_set,
    multiplied_word,
    non_single_bonds,
    ring_chain_attachment,
    specified_stereocenters,
)
from ._numerals import alkane_name
from ._substituents import alpha_sort_key, format_substituent_prefixes, name_branch

_ALLOWED_ATOMIC_NUMS = {6, 8, 34, *HALOGEN_PREFIXES}
_CHALCOGENS = (8, 34)


def _selenoate_matches(mol):
    """Carbons shaped like a selenoate group: a neutral, monovalent,
    double-bonded O or Se plus a formal-charge -1, monovalent, single-
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
        if {carbonyls[0].GetAtomicNum(), anions[0].GetAtomicNum()} != {8, 34}:
            continue
        matches.append((atom, carbonyls[0], anions[0]))
    return matches


def has_selenoate_shape(mol) -> bool:
    return bool(_selenoate_matches(mol))


def _find_selenoate_group(mol):
    matches = _selenoate_matches(mol)
    if len(matches) != 1:
        raise UnsupportedStructure(
            "exactly one selenoate (-COSe-/-CSeO-) group is required; zero or "
            "multiple such groups are not supported yet (P-72.2.2.2.1.1)"
        )
    selenoate_carbon, carbonyl_atom, anion_atom = matches[0]
    total_chalcogens = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() in _CHALCOGENS)
    if total_chalcogens != 2:
        raise UnsupportedStructure(
            "a chalcogen outside the single selenoate group's carbonyl/anion "
            "pair is out of scope for this module"
        )
    carbon_neighbors = [n for n in selenoate_carbon.GetNeighbors() if n.GetAtomicNum() == 6]
    if len(carbon_neighbors) > 1:
        raise UnsupportedStructure(
            "a selenoate carbon with more than one carbon neighbor is not a "
            "valid selenoate group"
        )
    return selenoate_carbon, carbonyl_atom, anion_atom


def _suffix_body(ene_locants, yne_locants):
    segments = []
    if ene_locants:
        segments.append((sorted(ene_locants), multiplied_word(len(ene_locants), "ene")))
    if yne_locants:
        segments.append((sorted(yne_locants), multiplied_word(len(yne_locants), "yne")))

    words = [word for _, word in segments] + ["selenoate"]
    for i in range(len(words) - 1):
        if words[i].endswith("e") and words[i + 1][0] in "aeiouy":
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
    total_count = sum(len(info["locants"]) for info in grouped.values())
    locant_set = lowest_locant_set(loc for info in grouped.values() for loc in info["locants"])
    citation_locants = tuple(
        loc
        for name in sorted(grouped, key=alpha_sort_key)
        for loc in sorted(grouped[name]["locants"])
    )
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


def _substituents_for_chain(graph, chain, halogens, excluded_atoms):
    chain_set = set(chain)
    substituents = {}
    for position, atom in enumerate(chain, start=1):
        branch_roots = [n for n in graph[atom] if n not in chain_set and n not in excluded_atoms]
        if branch_roots:
            substituents[position] = [name_branch(graph, root, atom, halogens) for root in branch_roots]
    return substituents


def _name_acyclic_selenoate(mol, selenoate_carbon_idx, excluded_atoms, bonds, stereo=None):
    """`stereo`: None, or a list of (stereocenter_atom_idx, "R"/"S") from
    `specified_stereocenters` -- if given, only chain candidates that
    include every stereocenter are eligible (P-92: a stereocenter on a
    substituent branch rather than the principal chain is out of scope,
    mirroring `_thioate.py`/`_carboxylate.py`), and the winning candidate's
    own locants are used to format a "(<locant><R/S>,...)-" prefix onto
    the final name."""
    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    chains = longest_chains(carbon_adjacency(mol))
    chain_length = len(chains[0])
    stereo_atoms = [atom for atom, _ in stereo] if stereo is not None else []

    eligible = []
    for chain in chains:
        if selenoate_carbon_idx not in chain:
            continue
        if bonds and bond_locants(chain, bonds) is None:
            continue
        chain_set = set(chain)
        if stereo is not None and any(atom not in chain_set for atom in stereo_atoms):
            continue
        eligible.append(chain)
    if not eligible:
        if stereo is not None and any(
            selenoate_carbon_idx in c and (not bonds or bond_locants(c, bonds) is not None) for c in chains
        ):
            raise UnsupportedStructure(
                "a stereocenter on a substituent branch rather than the "
                "principal chain is not supported yet (see P-92)"
            )
        raise UnsupportedStructure(
            "the selenoate carbon (and/or multiple bonds) does not lie on a "
            "single longest carbon chain"
        )

    best_key = None
    best_name = None
    best_position_of = None
    for chain in eligible:
        for candidate in (chain, list(reversed(chain))):
            if candidate[0] != selenoate_carbon_idx:
                continue
            position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
            ene_locants, yne_locants = bond_locants(candidate, bonds) if bonds else ([], [])
            substituents = _substituents_for_chain(graph, candidate, halogens, excluded_atoms)
            key, name = _candidate_key(chain_length, ene_locants, yne_locants, substituents)
            if best_key is None or key < best_key:
                best_key, best_name, best_position_of = key, name, position_of

    if stereo is not None:
        labels = sorted((best_position_of[atom], code) for atom, code in stereo)
        prefix = ",".join(f"{locant}{code}" for locant, code in labels)
        return f"({prefix})-{best_name}"
    return best_name


def _validate_and_prepare_selenoate(mol, aromatic_ring_atoms=frozenset()):
    """(selenoate_carbon, excluded_atoms, bonds, stereo) after checking the
    molecule fits this module's scope (see module docstring).

    `aromatic_ring_atoms`: atom indices already independently verified (by
    the caller, before this function runs) to form a single plain benzene
    ring with exactly one exocyclic attachment -- exempted from the
    aromatic-atom rejection below so `name_selenoate`'s benzene-ring-
    substituent path (see `_name_phenyl_chain_selenoate`) can reuse this
    same validation for the rest of the molecule. Empty by default, so
    every other caller's behavior is unchanged. Mirrors
    `_carboxylate.py`'s equivalent aromatic-exemption pattern."""
    selenoate_carbon, carbonyl_atom, anion_atom = _find_selenoate_group(mol)
    excluded_atoms = {carbonyl_atom.GetIdx(), anion_atom.GetIdx()}
    stereo = specified_stereocenters(mol)

    has_carbon = False
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atomic_num not in _ALLOWED_ATOMIC_NUMS:
            raise UnsupportedStructure(
                "heteroatoms other than the selenoate's own chalcogens "
                "(P-72.2.2.2.1.1) and halogen substituents (P-35.2.1) are "
                "not supported yet"
            )
        if atom.GetIdx() in excluded_atoms:
            continue
        if atom.GetFormalCharge() != 0 or atom.GetIsotope() != 0:
            raise UnsupportedStructure(
                "a charged or isotopically modified atom other than the "
                "single selenoate anion chalcogen is not supported yet"
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
    return selenoate_carbon, excluded_atoms, bonds, stereo


def _name_phenyl_chain_selenoate(mol, ring_atoms):
    """Name a selenoate whose -CO-Se(-) lies entirely on a single
    unbranched chain hanging off one atom of an otherwise-plain,
    unsubstituted benzene ring -- e.g. 3-phenylpropaneselenoate. The ring
    is cited as a 'phenyl' substituent prefix (via `name_branch`'s
    aromatic-ring recognition) on the chain, which is the parent hydride,
    mirroring `_carboxylate.py`'s `_name_phenyl_chain_carboxylate`.
    Narrower than the acyclic path above: no chain unsaturation and no
    specified stereocenter."""
    selenoate_carbon, excluded_atoms, bonds, stereo = _validate_and_prepare_selenoate(
        mol, aromatic_ring_atoms=ring_atoms
    )
    if stereo:
        raise UnsupportedStructure(
            "a specified stereocenter alongside a benzene-ring-substituent "
            "selenoate chain is not supported yet"
        )
    non_ring_unsaturation = [b for b in bonds if b[0] not in ring_atoms or b[1] not in ring_atoms]
    if non_ring_unsaturation:
        raise UnsupportedStructure(
            "chain unsaturation alongside a benzene-ring-substituent "
            "selenoate chain is not supported yet"
        )

    graph = adjacency(mol)
    attachment = ring_chain_attachment(graph, ring_atoms, set())
    if attachment is None:
        raise UnsupportedStructure(
            "a benzene ring with more than one exocyclic substituent "
            "alongside a chain selenoate is not supported yet"
        )
    chain, branches = longest_branched_chain(graph, selenoate_carbon.GetIdx(), ring_atoms, excluded_atoms)
    if len(chain) < 2:
        raise UnsupportedStructure(
            "a selenoate group directly attached to the benzene ring uses "
            "a separate construction, out of scope for this "
            "acyclic-chain-parent module"
        )

    chain_length = len(chain)
    halogens = halogen_substituents(mol)
    substituents = {
        position: [name_branch(graph, root, chain[position - 1], halogens, ring_atoms) for root in roots]
        for position, roots in branches.items()
    }
    grouped = group_substituents(substituents)
    return _name_from_substituents(chain_length, [], [], grouped)


def name_selenoate(mol) -> str:
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() == 1:
        ring_atoms = set(ring_info.AtomRings()[0])
        if is_plain_benzene_ring(mol, ring_atoms):
            return _name_phenyl_chain_selenoate(mol, ring_atoms)
    if mol.GetRingInfo().NumRings() > 0:
        raise UnsupportedStructure(
            "a selenoate group on/in a ring is out of scope for this "
            "acyclic-only module"
        )

    selenoate_carbon, excluded_atoms, bonds, stereo = _validate_and_prepare_selenoate(mol)
    return _name_acyclic_selenoate(mol, selenoate_carbon.GetIdx(), excluded_atoms, bonds, stereo)
