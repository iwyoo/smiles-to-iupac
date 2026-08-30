"""Naming of hydrazides (the '-hydrazide' suffix, terminal
-C(=O)NH-NH2) on acyclic saturated or unsaturated carbon chains, per the
IUPAC 2013 Recommendations ("the Blue Book"):

- P-66.3.1.1 (Chapter P-6a, https://iupac.qmul.ac.uk/BlueBook/PDF/P6a.pdf):
  'hydrazide' is the suffix for a chain-terminal -CO-NH-NH2 group.
  Structurally and mechanically this mirrors `_amide.py`'s -CONH2 exactly
  (P-66.1): the hydrazide carbon is always a chain terminus (after its
  carbonyl and its first hydrazide nitrogen, at most one more substituent
  -- another chain carbon -- fits), so its own suffix locant is never
  cited (P-14.3.3) and the ene/yne/halogen mechanics are unchanged.
  Confirmed via PubChem PUG REST: CID 77079 (`CCCC(=O)NN`) ->
  "butanehydrazide" (matching the Blue Book's own worked example
  'butanehydrazide (PIN)' directly), CID 93202 (`CCCCC(=O)NN`) ->
  "pentanehydrazide" (also a direct worked-example match), CID 79890
  (`CCC(=O)NN`) -> "propanehydrazide", CID 429824 (`ClCCC(=O)NN`) ->
  "3-chloropropanehydrazide", CID 219399 (`CC=CC(=O)NN`) ->
  "but-2-enehydrazide", CID 269313 (`OCCC(=O)NN`) ->
  "3-hydroxypropanehydrazide" (a coexisting standalone hydroxyl is cited
  as the 'hydroxy' prefix, same as `_amide.py`).
- P-66.3.1.2.1: **unlike amide** (where the systematic name, e.g.
  'ethanamide', is the actual PIN and the common name 'acetamide' is not),
  the Blue Book explicitly carves out five retained names as the
  preferred IUPAC names for hydrazide: 'cyanohydrazide', 'formohydrazide',
  'acetohydrazide', 'benzohydrazide', 'oxalohydrazide' -- of these, the
  two that would otherwise be plain acyclic chain cases are the
  mononuclear (methane, formohydrazide) and dinuclear (ethane,
  acetohydrazide) chain lengths, both implemented directly as their
  retained names rather than the systematic '-hydrazide' suffix this
  module otherwise produces -- confirmed directly from the primary text,
  not inferred from PubChem's own name choice (which happens to agree
  with the Blue Book here, unlike some other retained-name cases
  elsewhere in this project). Confirmed via PubChem PUG REST: CID 12229
  (`NNC=O`) -> "formohydrazide" (mononuclear, no room for a substituent),
  CID 14039 (`CC(=O)NN`) -> "acetohydrazide", CID 101883 (`ClCC(=O)NN`)
  -> "2-chloroacetohydrazide" (a substituent on the dinuclear case's
  terminal carbon is cited as an ordinary prefix, same mechanism as the
  systematic chain-length cases below). Chain lengths of three carbons
  and up are unaffected: P-66.3.1.2.3 confirms 'butanehydrazide
  (PIN)... not butyrohydrazide'
  directly, i.e. the systematic name is the real PIN there, same as this
  module already assumes.
- P-91.3/P-92 (`tasks/hydrazide-stereocenter-naming.md`): a molecule with
  one or more *specified* tetrahedral stereocenters -- every one on the
  principal chain itself (chain length >= 3; the 1-/2-carbon retained-name
  cases structurally can't have a genuine stereocenter), no unspecified
  one alongside them, and no C=C/C#N double-bond E/Z element -- gets a
  "(<locant><R/S>,...)-" prefix, ascending locant order, e.g.
  '(2R)-2-methylbutanehydrazide' (PubChem CID 30066157), same pattern as
  `_amide.py`. Both hydrazide nitrogens are never stereocenters (verified
  via RDKit's `Chem.FindPotentialStereo`), so this support is
  unconditional.
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
    longest_chains,
    lowest_locant_set,
    multiplied_word,
    non_single_bonds,
    specified_stereocenters,
)
from ._numerals import alkane_name
from ._substituents import alpha_sort_key, format_substituent_prefixes, name_branch

_ALLOWED_ATOMIC_NUMS = {6, 7, 8, *HALOGEN_PREFIXES}


def has_hydrazide_shape(mol) -> bool:
    """True if some carbon carries a doubly-bonded, monovalent carbonyl
    oxygen and a singly-bonded nitrogen (degree 2, one H) that is itself
    singly bonded to a terminal, two-H nitrogen -- a -CO-NH-NH2 pattern --
    regardless of whether the rest of the molecule is in scope. Used by
    `core.py` to route ahead of the amide/aldehyde/ketone dispatch, since
    a hydrazide carbon would otherwise look amide-shaped (P-66.1's own
    -CONH2 check doesn't look past the first nitrogen)."""
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6:
            continue
        has_carbonyl = any(
            n.GetAtomicNum() == 8
            and n.GetDegree() == 1
            and mol.GetBondBetweenAtoms(atom.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 2.0
            for n in atom.GetNeighbors()
        )
        if not has_carbonyl:
            continue
        for n1 in atom.GetNeighbors():
            if (
                n1.GetAtomicNum() == 7
                and n1.GetDegree() == 2
                and n1.GetTotalNumHs() == 1
                and mol.GetBondBetweenAtoms(atom.GetIdx(), n1.GetIdx()).GetBondTypeAsDouble() == 1.0
            ):
                (n2,) = [n for n in n1.GetNeighbors() if n.GetIdx() != atom.GetIdx()]
                if (
                    n2.GetAtomicNum() == 7
                    and n2.GetDegree() == 1
                    and n2.GetTotalNumHs() == 2
                    and mol.GetBondBetweenAtoms(n1.GetIdx(), n2.GetIdx()).GetBondTypeAsDouble() == 1.0
                ):
                    return True
    return False


def _validate_and_collect_hydrazide(mol):
    """Check the molecule fits this module's scope (see module docstring)
    and return (hydrazide_carbon, hydrazide_oxygen, n1, n2, hydroxyls):
    the single -CO-NH-NH2 carbon/oxygen/n1(-NH-)/n2(-NH2) atom indices,
    and the set of any coexisting standalone hydroxyl-oxygen indices."""
    has_carbon = False
    hydrazide_carbons = set()
    oxygen_by_carbon = {}
    n1_by_carbon = {}
    hydroxyls = set()
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atomic_num not in _ALLOWED_ATOMIC_NUMS:
            raise UnsupportedStructure(
                "heteroatoms other than a hydrazide's oxygen/nitrogens "
                "(P-66.3.1.1) and halogen substituents (P-35.2.1) are not "
                "supported yet"
            )
        if atom.GetFormalCharge() != 0 or atom.GetIsotope() != 0:
            raise UnsupportedStructure("charged or isotopically modified atoms are not supported yet")
        if atomic_num == 6:
            has_carbon = True
            if atom.GetIsAromatic():
                raise UnsupportedStructure(
                    "aromatic rings are out of scope for this module (see "
                    "the separate aromatic-ring module)"
                )
        elif atomic_num == 8:
            if atom.GetDegree() != 1:
                raise UnsupportedStructure(
                    "an oxygen bonded to more than one heavy atom (e.g. an "
                    "ether or ester) is out of scope; only an isolated "
                    "hydrazide carbonyl or a standalone hydroxyl is "
                    "supported (P-66.3.1.1)"
                )
            (bond,) = atom.GetBonds()
            (carbon,) = atom.GetNeighbors()
            if carbon.GetAtomicNum() != 6:
                raise UnsupportedStructure("a hydrazide/hydroxyl oxygen must be attached to a carbon atom")
            bond_order = bond.GetBondTypeAsDouble()
            if bond_order == 1.0 and atom.GetTotalNumHs() == 1:
                hydroxyls.add(atom.GetIdx())
                continue
            if bond_order != 2.0:
                raise UnsupportedStructure(
                    "an oxygen that isn't a carbonyl (=O) or hydroxyl (-OH) "
                    "is out of scope for this module"
                )
            if carbon.GetIsAromatic():
                raise UnsupportedStructure(
                    "a carbonyl on an aromatic ring is out of scope for "
                    "this module"
                )
            oxygen_by_carbon.setdefault(carbon.GetIdx(), []).append(atom.GetIdx())
        elif atomic_num == 7:
            neighbors = atom.GetNeighbors()
            carbon_neighbors = [n for n in neighbors if n.GetAtomicNum() == 6]
            nitrogen_neighbors = [n for n in neighbors if n.GetAtomicNum() == 7]
            is_n1_shape = (
                len(carbon_neighbors) == 1
                and len(nitrogen_neighbors) == 1
                and atom.GetDegree() == 2
                and atom.GetTotalNumHs() == 1
                and mol.GetBondBetweenAtoms(atom.GetIdx(), carbon_neighbors[0].GetIdx()).GetBondTypeAsDouble()
                == 1.0
                and mol.GetBondBetweenAtoms(atom.GetIdx(), nitrogen_neighbors[0].GetIdx()).GetBondTypeAsDouble()
                == 1.0
            )
            is_n2_shape = (
                not carbon_neighbors
                and len(nitrogen_neighbors) == 1
                and atom.GetDegree() == 1
                and atom.GetTotalNumHs() == 2
            )
            if is_n1_shape:
                (n2,) = nitrogen_neighbors
                if not (
                    n2.GetDegree() == 1
                    and n2.GetTotalNumHs() == 2
                    and not [x for x in n2.GetNeighbors() if x.GetAtomicNum() == 6]
                ):
                    raise UnsupportedStructure(
                        "a hydrazide nitrogen chain with any substituent "
                        "other than its hydrogens (N-substituted "
                        "hydrazide) is out of scope for this module"
                    )
                n1_by_carbon.setdefault(carbon_neighbors[0].GetIdx(), []).append((atom.GetIdx(), n2.GetIdx()))
                continue
            if is_n2_shape:
                # Already validated as part of its N1 partner's own check
                # above; nothing further to do for this atom in isolation.
                continue
            raise UnsupportedStructure(
                "a nitrogen shaped like neither a hydrazide's -NH- nor its "
                "terminal -NH2 (P-66.3.1.1) is out of scope for this "
                "module (e.g. a plain primary amide/amine nitrogen bonded "
                "directly to a carbon; see _amide.py for a plain primary "
                "amide)"
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

    for carbon_idx in set(oxygen_by_carbon) | set(n1_by_carbon):
        oxygens = oxygen_by_carbon.get(carbon_idx, [])
        n1_pairs = n1_by_carbon.get(carbon_idx, [])
        if len(oxygens) != 1 or len(n1_pairs) != 1:
            raise UnsupportedStructure(
                "an oxygen/nitrogen pattern that isn't exactly one "
                "carbonyl oxygen and one -NH-NH2 chain on the same carbon "
                "is a more/less senior characteristic group than a plain "
                "hydrazide (Table 3.3/4.4), which this module does not "
                "attempt to disambiguate"
            )
        carbon = mol.GetAtomWithIdx(carbon_idx)
        carbon_neighbors = [n for n in carbon.GetNeighbors() if n.GetAtomicNum() == 6]
        if len(carbon_neighbors) > 1:
            raise UnsupportedStructure(
                "a hydrazide carbon with more than one carbon neighbor is "
                "not a valid terminal hydrazide carbon (a ketone-shaped "
                "carbon is out of scope for this module)"
            )
        hydrazide_carbons.add(carbon_idx)

    if not hydrazide_carbons:
        raise UnsupportedStructure(
            "no hydrazide (-CO-NH-NH2) group found; this module only "
            "handles hydrazides"
        )
    if len(hydrazide_carbons) > 1:
        raise UnsupportedStructure(
            "more than one hydrazide group is out of scope for this module"
        )
    (hydrazide_carbon,) = hydrazide_carbons
    (hydrazide_oxygen,) = oxygen_by_carbon[hydrazide_carbon]
    ((n1, n2),) = n1_by_carbon[hydrazide_carbon]
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")
    return hydrazide_carbon, hydrazide_oxygen, n1, n2, hydroxyls


def _suffix_body(ene_locants, yne_locants):
    """Locant-and-suffix string for the combined 'ene'/'yne'/'hydrazide'
    ending (e.g. '2-enehydrazide'); the hydrazide group's own locant is
    never cited (P-14.3.3, see module docstring)."""
    segments = []
    if ene_locants:
        segments.append((sorted(ene_locants), multiplied_word(len(ene_locants), "ene")))
    if yne_locants:
        segments.append((sorted(yne_locants), multiplied_word(len(yne_locants), "yne")))

    words = [word for _, word in segments] + ["hydrazide"]
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
    separator = "-" if has_unsaturation else ""
    return prefix + stem + ("a" if needs_stem_a else "") + separator + body


def _candidate_key(chain_length, ene_locants, yne_locants, substituents):
    """Sort key implementing P-44.4.1.10 (ene/yne locants) ahead of P-45.2
    (substituent-prefix locants), most-preferred first. The hydrazide
    group's own locant isn't part of this key: candidates are
    pre-filtered so the hydrazide carbon always sits at C1 (see
    `_name_acyclic_hydrazide`)."""
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


def _substituents_for_chain(graph, chain, halogens, excluded):
    chain_set = set(chain)
    substituents = {}
    for position, atom in enumerate(chain, start=1):
        branch_roots = [n for n in graph[atom] if n not in chain_set and n not in excluded]
        if branch_roots:
            substituents[position] = [name_branch(graph, root, atom, halogens) for root in branch_roots]
    return substituents


def _name_acyclic_hydrazide(mol, hydrazide_carbon, excluded, hydroxyls, bonds, stereo=None):
    """`stereo`: None, or a list of (stereocenter_atom_idx, "R"/"S") from
    `specified_stereocenters` -- if given, only chain candidates that
    include every stereocenter are eligible (P-92: a stereocenter on a
    substituent branch rather than the principal chain is out of scope,
    mirroring `_amide.py`'s identical treatment), and the winning
    candidate's own locants are used to format a "(<locant><R/S>,...)-"
    prefix onto the final name -- applied outermost, since a
    stereodescriptor always sits at the very front of the complete name
    (P-91.3). A 1- or 2-carbon chain (the retained-name special cases)
    structurally can't have a genuine stereocenter, so `stereo` is only
    ever non-None here for chain_length >= 3."""
    graph = adjacency(mol)
    halogens = {**halogen_substituents(mol), **{o: "hydroxy" for o in hydroxyls}}
    chains = longest_chains(carbon_adjacency(mol))
    chain_length = len(chains[0])

    if chain_length < 3:
        if bonds:
            raise UnsupportedStructure(
                "unsaturation alongside a 1- or 2-carbon hydrazide chain "
                "is not supported"
            )
        (chain,) = [c for c in chains if hydrazide_carbon in c]
        if chain[0] != hydrazide_carbon:
            chain = list(reversed(chain))
        if chain_length == 1:
            # P-66.3.1.2.1: the mononuclear case's retained name
            # ('formohydrazide') is itself the PIN; there's no room on a
            # single carbon for any substituent.
            return "formohydrazide"
        # P-66.3.1.2.1: the dinuclear case's retained name
        # ('acetohydrazide') is the PIN, with any substituent on its
        # terminal carbon cited as an ordinary prefix (e.g.
        # '2-chloroacetohydrazide', PubChem CID 101883).
        substituents = _substituents_for_chain(graph, chain, halogens, excluded)
        grouped = group_substituents(substituents)
        prefix = format_substituent_prefixes(grouped)
        return prefix + "acetohydrazide"

    stereo_atoms = [atom for atom, _ in stereo] if stereo is not None else []

    eligible = []
    for chain in chains:
        if hydrazide_carbon not in chain:
            continue
        if bonds and bond_locants(chain, bonds) is None:
            continue
        chain_set = set(chain)
        if stereo is not None and any(atom not in chain_set for atom in stereo_atoms):
            continue
        eligible.append(chain)
    if not eligible:
        if stereo is not None and any(
            hydrazide_carbon in c and (not bonds or bond_locants(c, bonds) is not None) for c in chains
        ):
            raise UnsupportedStructure(
                "a stereocenter on a substituent branch rather than the "
                "principal chain is not supported yet (see P-92)"
            )
        raise UnsupportedStructure(
            "the hydrazide-bearing carbon (and/or a multiple bond) does "
            "not lie on a single longest carbon chain; a shorter "
            "principal chain is not supported yet"
        )

    best_key = None
    best_name = None
    best_position_of = None
    for chain in eligible:
        for candidate in (chain, list(reversed(chain))):
            if candidate[0] != hydrazide_carbon:
                # The hydrazide carbon must sit at C1 (see module
                # docstring); a direction that doesn't start there is
                # never valid.
                continue
            ene_locants, yne_locants = bond_locants(candidate, bonds) if bonds else ([], [])
            substituents = _substituents_for_chain(graph, candidate, halogens, excluded)
            key, name = _candidate_key(chain_length, ene_locants, yne_locants, substituents)
            if best_key is None or key < best_key:
                position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
                best_key, best_name, best_position_of = key, name, position_of

    if stereo is not None:
        labels = sorted((best_position_of[atom], code) for atom, code in stereo)
        prefix = ",".join(f"{locant}{code}" for locant, code in labels)
        return f"({prefix})-{best_name}"
    return best_name


def name_hydrazide(mol) -> str:
    hydrazide_carbon, hydrazide_oxygen, n1, n2, hydroxyls = _validate_and_collect_hydrazide(mol)
    stereo = specified_stereocenters(mol)
    if mol.GetRingInfo().NumRings() > 0:
        raise UnsupportedStructure(
            "a hydrazide on/in a ring is out of scope for this "
            "acyclic-only module"
        )

    excluded = {hydrazide_oxygen, n1, n2}
    all_non_single = [b for b in non_single_bonds(mol) if b[0] not in excluded and b[1] not in excluded]
    bonds = [b for b in all_non_single if b[2] in (ENE_BOND_ORDER, YNE_BOND_ORDER)]
    if len(bonds) != len(all_non_single):
        raise UnsupportedStructure(
            "a bond order other than single, double, or triple is not "
            "supported (see P-31.1.1.1)"
        )
    graph = adjacency(mol)
    ene_yne_carbons = {a for a, b, _ in bonds} | {b for a, b, _ in bonds}
    for o in hydroxyls:
        (carbon,) = graph[o]
        if carbon in ene_yne_carbons:
            raise UnsupportedStructure(
                "a hydroxyl on a carbon that is also part of a C=C/C#C bond "
                "(an enol) is a tautomer of a more senior carbonyl form and "
                "is out of scope for this module (P-31.1.4.2.4)"
            )

    return _name_acyclic_hydrazide(mol, hydrazide_carbon, excluded, hydroxyls, bonds, stereo)
