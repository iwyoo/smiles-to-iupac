"""Naming of sulfonate anions (the 'sulfonate' suffix, -S(=O)(=O)-O(-)) on
acyclic saturated or unsaturated carbon chains, per the IUPAC 2013
Recommendations ("the Blue Book"):

- P-72.2.2.2.1.1 (Chapter P-7, https://iupac.qmul.ac.uk/BlueBook/P7.html):
  the preferred IUPAC name of an anion formed by removing a hydron from the
  chalcogen atom of an acid is formed by replacing the acid name's 'ic
  acid' ending with 'ate' -- 'sulfonic acid' -> 'sulfonate'. Otherwise
  mirrors `_sulfonic_acid.py`'s acyclic path: the same chain-search/
  numbering, the same 'ene'/'yne' mechanics, and the same halogen-
  substituent/stereocenter handling -- PubChem PUG REST confirms the
  parallel construction directly: 'CS(=O)(=O)[O-]' -> 'methanesulfonate'
  (CID 85257), 'CCS(=O)(=O)[O-]' -> 'ethanesulfonate' (CID 3717105),
  'CCCS(=O)(=O)[O-]' -> 'propane-1-sulfonate' (CID 4431756),
  'C=CCS(=O)(=O)[O-]' -> 'prop-2-ene-1-sulfonate' (CID 4447605),
  'CC[C@@H](C)S(=O)(=O)[O-]' -> '(2R)-butane-2-sulfonate' (CID 49867167).
- P-14.3.4.2(a)/(b) locant omission: a mononuclear parent never cites the
  suffix locant (e.g. 'methanesulfonate'), and on a *two-carbon* chain the
  sulfonate carbon is always forced to be C1 (the principal characteristic
  group gets the lowest locant), so its own locant is never cited there
  either -- regardless of how many other substituents are present or
  where they sit, confirmed via PubChem: 'CC(Cl)S(=O)(=O)[O-]' ->
  '1-chloroethanesulfonate' (CID 19003700), 'ClCS(=O)(=O)[O-]' ->
  'chloromethanesulfonate', 'ClCC(Cl)S(=O)(=O)[O-]' ->
  '1,2-dichloroethanesulfonate'. This is broader than
  `_sulfonic_acid.py`'s equivalent check (which only omits the locant when
  no other substituent is present) -- that module has the same latent
  gap (confirmed live: 'CC(Cl)S(=O)(=O)O' -> the module's current
  '1-chloroethane-1-sulfonic acid' vs. PubChem's
  '1-chloroethanesulfonic acid'), left as a follow-up there rather than
  fixed as part of this module's own scope.
- Unlike `_sulfonic_acid.py`, this first pass only supports an acyclic
  chain (no ring/benzene-ring path yet -- follow-up work, mirroring how
  `_carboxylate.py`/`_thioate.py` started acyclic-first relative to their
  neutral-acid counterparts) and only exactly one -SO3(-) group.

Explicitly out of scope (raise `UnsupportedStructure`):
- Any ring anywhere in the molecule.
- More than one sulfonate group, or any oxygen that isn't part of the
  single sulfonate's own three oxygens.
- Any charged or radical atom other than the single anionic oxygen's
  formal charge -1.
- A sulfonate carbon that is also part of a C=C/C#C bond (mirrors
  `_sulfonic_acid.py`'s identical `_reject_enesulfonic_carbon` rejection).
- Any other heteroatom (N, ...). Combination with a cation into a full
  salt name (e.g. 'sodium methanesulfonate') is a separate follow-up for
  `_salt.py`, tracked in the coverage index.

P-91.3/P-92: a molecule with one or more *specified* tetrahedral
stereocenters -- every one on the principal chain itself, no unspecified
one alongside them, and no C=C/C#N double-bond E/Z element -- gets a
"(<locant><R/S>,...)-" prefix, ascending locant order, same pattern as
`_sulfonic_acid.py`.
"""

from rdkit import Chem

from ._common import (
    HALOGEN_PREFIXES,
    UnsupportedStructure,
    adjacency,
    bond_locant,
    bond_locants,
    carbon_adjacency,
    elides_before,
    group_substituents,
    halogen_substituents,
    longest_chains,
    lowest_locant_set,
    multiplied_word,
    non_single_bonds,
    specified_stereocenters,
    substituent_locant_set_and_citation,
)
from ._numerals import alkane_name
from ._substituents import format_substituent_prefixes, name_branch

_ENE_ORDER = 2.0
_YNE_ORDER = 3.0


def _sulfonate_sulfur_atoms(mol):
    """Sulfur atoms shaped like a sulfonate group: bonded to exactly one
    carbon, two double-bonded (terminal) oxygens, and one single-bonded,
    monovalent, formal-charge -1 anionic oxygen (mirrors
    `_sulfonic_acid.py`'s `_sulfonic_sulfur_atoms`, with the hydroxyl
    oxygen replaced by an anionic one, no H)."""
    matches = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 16 or atom.GetDegree() != 4:
            continue
        neighbors = atom.GetNeighbors()
        carbons = [n for n in neighbors if n.GetAtomicNum() == 6]
        oxygens = [n for n in neighbors if n.GetAtomicNum() == 8]
        if len(carbons) != 1 or len(oxygens) != 3:
            continue
        double_os = [o for o in oxygens if mol.GetBondBetweenAtoms(atom.GetIdx(), o.GetIdx()).GetBondTypeAsDouble() == 2.0]
        anion_os = [o for o in oxygens if mol.GetBondBetweenAtoms(atom.GetIdx(), o.GetIdx()).GetBondTypeAsDouble() == 1.0]
        if len(double_os) != 2 or len(anion_os) != 1:
            continue
        if any(o.GetDegree() != 1 for o in double_os):
            continue
        (anion_o,) = anion_os
        if anion_o.GetDegree() != 1 or anion_o.GetFormalCharge() != -1:
            continue
        matches.append(atom)
    return matches


def has_sulfonate_shape(mol) -> bool:
    return bool(_sulfonate_sulfur_atoms(mol))


def _validate_and_collect_sulfonate(mol):
    sulfur_atoms = _sulfonate_sulfur_atoms(mol)
    if not sulfur_atoms:
        raise UnsupportedStructure(
            "no sulfonate (-SO3-) group found; this module only handles "
            "sulfonate anions"
        )
    if len(sulfur_atoms) > 1:
        raise UnsupportedStructure(
            "more than one sulfonate group is out of scope for this module"
        )
    (sulfur,) = sulfur_atoms
    sulfonate_atom_idxs = {sulfur.GetIdx()}
    sulfonate_atom_idxs.update(n.GetIdx() for n in sulfur.GetNeighbors() if n.GetAtomicNum() == 8)

    has_carbon = False
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atomic_num == 6:
            has_carbon = True
            if atom.GetIsAromatic():
                raise UnsupportedStructure(
                    "aromatic rings are out of scope for this module"
                )
        elif atomic_num in HALOGEN_PREFIXES:
            if atom.GetDegree() != 1:
                raise UnsupportedStructure(
                    "a halogen atom must be a monovalent substituent (P-35.2.1)"
                )
        elif atom.GetIdx() not in sulfonate_atom_idxs:
            raise UnsupportedStructure(
                "heteroatoms other than a sulfonate group (P-72.2.2.2.1.1) "
                "and halogen substituents (P-35.2.1) are not supported yet"
            )
        if atom.GetIdx() not in sulfonate_atom_idxs and (
            atom.GetFormalCharge() != 0 or atom.GetIsotope() != 0
        ):
            raise UnsupportedStructure("charged or isotopically modified atoms are not supported yet")
    if not has_carbon:
        raise UnsupportedStructure(
            "a structure with no carbon atom has no hydrocarbon parent "
            "hydride to substitute"
        )
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")

    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() > 0:
        raise UnsupportedStructure(
            "a ring anywhere in the molecule is out of scope for this "
            "acyclic-only module"
        )

    (carbon,) = (n for n in sulfur.GetNeighbors() if n.GetAtomicNum() == 6)
    return sulfur.GetIdx(), carbon.GetIdx()


def _reject_enesulfonate_carbon(so3_carbon, bonds):
    unsaturated_atoms = {a for a, b, _ in bonds} | {b for a, b, _ in bonds}
    if so3_carbon in unsaturated_atoms:
        raise UnsupportedStructure(
            "a sulfonate on a carbon that is also part of a C=C/C#C bond "
            "is out of scope for this module"
        )


def _suffix_body(ene_locants, yne_locants, so3_locant):
    segments = []
    if ene_locants:
        segments.append((sorted(ene_locants), multiplied_word(len(ene_locants), "ene")))
    if yne_locants:
        segments.append((sorted(yne_locants), multiplied_word(len(yne_locants), "yne")))
    segments.append(([so3_locant], "sulfonate"))

    words = [word for _, word in segments]
    for i in range(len(words) - 1):
        if words[i].endswith("e") and elides_before(words[i + 1]):
            words[i] = words[i][:-1]

    parts = [
        f"{','.join(str(loc) for loc in locants)}-{word}"
        for (locants, _), word in zip(segments, words)
    ]
    return "-".join(parts)


def _name_from_substituents(chain_length, so3_locant, ene_locants, yne_locants, grouped):
    has_unsaturation = bool(ene_locants or yne_locants)

    if chain_length == 1:
        # P-14.3.4.2(a): a mononuclear parent's locant is always '1' and
        # never cited.
        return format_substituent_prefixes(grouped, omit_locants=True) + alkane_name(1) + "sulfonate"

    if chain_length == 2 and not has_unsaturation:
        # P-14.3.4.2(b): on a two-carbon chain the sulfonate carbon is
        # always forced to be C1 (the principal characteristic group gets
        # the lowest locant), so citing '1' is never informative -- the
        # locant is omitted regardless of how many other substituents are
        # present or where they sit, e.g. 'ethanesulfonate',
        # '1-chloroethanesulfonate', '2-chloroethanesulfonate',
        # '1,2-dichloroethanesulfonate' (all confirmed via PubChem PUG
        # REST). This is broader than `_sulfonic_acid.py`'s equivalent
        # check (`total_subs == 0` only) -- that module has the same
        # latent gap, confirmed live against
        # '1-chloroethanesulfonic acid'/'2-chloroethanesulfonic acid'/
        # '1,2-dichloroethanesulfonic acid', tracked as a follow-up there
        # rather than fixed in this module's own scope.
        return format_substituent_prefixes(grouped) + alkane_name(2) + "sulfonate"

    prefix = format_substituent_prefixes(grouped)
    if has_unsaturation:
        stem = alkane_name(chain_length)[:-3]
        needs_stem_a = (len(ene_locants) >= 2) if ene_locants else (len(yne_locants) >= 2)
    else:
        stem = alkane_name(chain_length)
        needs_stem_a = False

    body = _suffix_body(ene_locants, yne_locants, so3_locant)
    return prefix + stem + ("a" if needs_stem_a else "") + "-" + body


def _candidate_key(chain_length, so3_locant, ene_locants, yne_locants, substituents):
    grouped = group_substituents(substituents)
    locant_set, total_count, citation_locants = substituent_locant_set_and_citation(grouped)
    combined_locant_set = lowest_locant_set(ene_locants + yne_locants)
    ene_locant_set = lowest_locant_set(ene_locants)
    name = _name_from_substituents(chain_length, so3_locant, ene_locants, yne_locants, grouped)
    return (
        (
            so3_locant,
            combined_locant_set,
            ene_locant_set,
            -total_count,
            locant_set,
            citation_locants,
            name,
        ),
        name,
    )


def _substituents_for_chain(graph, chain, halogens, excluded, mol=None):
    chain_set = set(chain)
    substituents = {}
    for position, atom in enumerate(chain, start=1):
        branch_roots = [n for n in graph[atom] if n not in chain_set and n not in excluded]
        if branch_roots:
            substituents[position] = [name_branch(graph, root, atom, halogens, mol=mol) for root in branch_roots]
    return substituents


def _name_acyclic_sulfonate(mol, sulfur_idx, so3_carbon, bonds, stereo=None):
    """`stereo`: None, or a list of (stereocenter_atom_idx, "R"/"S") from
    `specified_stereocenters` -- if given, only chain candidates that
    include every stereocenter are eligible (P-92), and the winning
    candidate's own locants are used to format a "(<locant><R/S>,...)-"
    prefix onto the final name, mirroring `_sulfonic_acid.py`'s
    `_name_acyclic_sulfonic_acid`."""
    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    excluded = {sulfur_idx}
    chains = longest_chains(carbon_adjacency(mol))
    chain_length = len(chains[0])
    stereo_atoms = [atom for atom, _ in stereo] if stereo is not None else []

    eligible = []
    for chain in chains:
        chain_set = set(chain)
        if so3_carbon not in chain_set:
            continue
        if bonds and bond_locants(chain, bonds) is None:
            continue
        if stereo is not None and any(atom not in chain_set for atom in stereo_atoms):
            continue
        eligible.append(chain)
    if not eligible:
        if stereo is not None and any(
            so3_carbon in c and (not bonds or bond_locants(c, bonds) is not None) for c in chains
        ):
            raise UnsupportedStructure(
                "a stereocenter on a substituent branch rather than the "
                "principal chain is not supported yet (see P-92)"
            )
        raise UnsupportedStructure(
            "the sulfonate-bearing carbon (and/or a multiple bond) does "
            "not lie on a single longest carbon chain; a shorter "
            "principal chain is not supported yet"
        )

    best_key = None
    best_name = None
    best_position_of = None
    for chain in eligible:
        for candidate in (chain, list(reversed(chain))):
            position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
            so3_locant = position_of[so3_carbon]
            ene_locants, yne_locants = bond_locants(candidate, bonds) if bonds else ([], [])
            substituents = _substituents_for_chain(graph, candidate, halogens, excluded, mol=mol)
            key, name = _candidate_key(chain_length, so3_locant, ene_locants, yne_locants, substituents)
            if best_key is None or key < best_key:
                best_key, best_name, best_position_of = key, name, position_of

    if stereo is not None:
        labels = sorted((best_position_of[atom], code) for atom, code in stereo)
        prefix = ",".join(f"{locant}{code}" for locant, code in labels)
        return f"({prefix})-{best_name}"
    return best_name


def name_sulfonate(mol) -> str:
    sulfur_idx, so3_carbon = _validate_and_collect_sulfonate(mol)
    stereo = specified_stereocenters(mol)
    all_non_single = non_single_bonds(mol)
    bonds = [b for b in all_non_single if b[2] in (_ENE_ORDER, _YNE_ORDER) and sulfur_idx not in (b[0], b[1])]
    if len(bonds) != len(all_non_single) - 2:
        # The two S=O double bonds are always present and excluded above;
        # anything else non-single must be a chain ene/yne bond.
        raise UnsupportedStructure(
            "a bond order other than single, double, or triple is not "
            "supported (see P-31.1.1.1)"
        )
    _reject_enesulfonate_carbon(so3_carbon, bonds)
    return _name_acyclic_sulfonate(mol, sulfur_idx, so3_carbon, bonds, stereo)
