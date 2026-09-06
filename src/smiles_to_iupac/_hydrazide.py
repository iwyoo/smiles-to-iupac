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
- P-91.3/P-92: a molecule with
  one or more *specified* tetrahedral stereocenters -- every one on the
  principal chain itself (chain length >= 3; the 1-/2-carbon retained-name
  cases structurally can't have a genuine stereocenter), no unspecified
  one alongside them, and no C=C/C#N double-bond E/Z element -- gets a
  "(<locant><R/S>,...)-" prefix, ascending locant order, e.g.
  '(2R)-2-methylbutanehydrazide' (PubChem CID 30066157), same pattern as
  `_amide.py`. Both hydrazide nitrogens are never stereocenters (verified
  via RDKit's `Chem.FindPotentialStereo`), so this support is
  unconditional.
- A plain, unbranched, unsubstituted alkyl substituent on either
  hydrazide nitrogen is cited as
  an "N-"/"N'-" prefix (the carbonyl-adjacent nitrogen is "N", the
  terminal one "N'", mirroring `_amide.py`'s "N-"/`_urea.py`'s
  "N-"/"N'-" convention), placed directly ahead of the acyl stem, in
  alphabetical order by substituent name; identically-named substituents
  (whether on the same nitrogen or split across both) share one
  multiplying prefix with their locants listed together, e.g. "N,N'-di"/
  "N,N',N'-tri". Confirmed via PubChem PUG REST: 'N-methylacetohydrazide'
  (`CC(=O)N(C)N`, CID 19051), "N'-methylacetohydrazide" (`CC(=O)NNC`, CID
  122488), "N,N'-dimethylacetohydrazide" (`CC(=O)N(C)NC`, CID 12596294),
  "N',N'-dimethylacetohydrazide" (`CC(=O)NN(C)C`, CID 80385),
  "N,N',N'-trimethylacetohydrazide" (`CC(=O)N(C)N(C)C`, CID 12362621),
  "N-ethyl-N'-methylacetohydrazide" (`CC(=O)N(CC)NC`, CID 89037658),
  "N'-ethyl-N-methylacetohydrazide" (`CC(=O)N(C)NCC`, CID 119096326), and
  the systematic (non-retained-name) chain lengths behave the same way:
  "N-methylbutanehydrazide" (`CCCC(=O)N(C)N`, CID 53649946),
  "N'-methylbutanehydrazide" (`CCCC(=O)NNC`, CID 88556637). The
  carbonyl-adjacent nitrogen can carry at most one alkyl substituent
  (it's already bonded to the carbonyl carbon and the other nitrogen);
  the terminal nitrogen can carry up to two.
- P-66.3.3.3: a *diacylhydrazide* (R-CO-NH-NH-CO-R', both hydrazide
  nitrogens bearing no other substituent) is named by keeping the senior
  acyl group's own hydrazide name as the parent and citing the other as
  an "N'-<acyl>" prefix, e.g. 'N'-benzoylbenzohydrazide' (not
  '1,2-dibenzoylhydrazine'). Determining which acyl group is senior
  needs the general acid-seniority rules (Table 3.3/4.4); this module
  only handles the *symmetric* case (both acyl groups name identically,
  so which one plays which role is arbitrary), mirroring `_imide.py`'s
  identical narrowing for the same reason. Confirmed via PubChem PUG
  REST: 'N'-acetylacetohydrazide' (`CC(=O)NNC(=O)C`, CID 72884, dinuclear
  retained-name case on both sides) and 'N'-propanoylpropanehydrazide'
  (`CCC(=O)NNC(=O)CC`, CID 73715, systematic three-carbon case on both
  sides). An unsymmetric diacylhydrazide, an acyl group on the
  carbonyl-adjacent nitrogen instead of the terminal one, and any other
  substituent alongside either acyl chain (halogen, branching, an
  aromatic acyl group) are out of scope for this first pass.
- One narrow exception to the "acyclic-only" scope below:
  `_name_phenyl_chain_hydrazide` names a hydrazide's chain hanging off a
  single plain, unsubstituted benzene ring (e.g. '3-phenylpropanehydrazide',
  or '2-phenylacetohydrazide' for the dinuclear retained-name case,
  mirroring the module's own '2-chloroacetohydrazide'
  substituent-on-retained-name pattern), mirroring `_amide.py`'s
  identical benzene-ring-substituent path. Narrower than the acyclic path: no N-/N'-alkyl substitution, no
  coexisting standalone hydroxyl, no chain unsaturation, and no specified
  stereocenter, and a hydrazide directly on the ring (benzohydrazide-
  style) stays out of scope for this chain-parent module.
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
    linear_branch,
    longest_chains,
    lowest_locant_set,
    multiplied_word,
    non_single_bonds,
    ordered_chain,
    ring_chain_attachment,
    ring_chain_attachment_with_halogens,
    specified_stereocenters,
)
from ._numerals import alkane_name, alkyl_name, multiplying_prefix
from ._substituents import (
    alpha_sort_key,
    format_substituent_prefixes,
    name_branch,
    plain_alkyl_ring_substituents,
)

_ALLOWED_ATOMIC_NUMS = {6, 7, 8, *HALOGEN_PREFIXES}


def _is_carbonyl_carbon(mol, carbon_atom):
    return any(
        o.GetAtomicNum() == 8
        and o.GetDegree() == 1
        and mol.GetBondBetweenAtoms(carbon_atom.GetIdx(), o.GetIdx()).GetBondTypeAsDouble() == 2.0
        for o in carbon_atom.GetNeighbors()
    )


def has_hydrazide_shape(mol) -> bool:
    """True if some carbon carries a doubly-bonded, monovalent carbonyl
    oxygen and a singly-bonded nitrogen that is itself singly bonded to
    exactly one other nitrogen -- a -CO-N(-)-N(-) pattern, regardless of
    how many (if any) plain alkyl substituents sit on either nitrogen --
    regardless of whether the rest of the molecule is in scope. Used by
    `core.py` to route ahead of the amide/aldehyde/ketone dispatch, since
    a hydrazide carbon would otherwise look amide-shaped (P-66.1's own
    -CONH2 check doesn't look past the first nitrogen). Doesn't check
    substituent count/H totals (unlike the old, unsubstituted-only
    version of this check) so N-/N'-substituted hydrazides route here
    too; `_validate_and_collect_hydrazide` below does the real scope
    enforcement."""
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6:
            continue
        if not _is_carbonyl_carbon(mol, atom):
            continue
        for n1 in atom.GetNeighbors():
            if n1.GetAtomicNum() != 7:
                continue
            if mol.GetBondBetweenAtoms(atom.GetIdx(), n1.GetIdx()).GetBondTypeAsDouble() != 1.0:
                continue
            nitrogen_neighbors = [n for n in n1.GetNeighbors() if n.GetAtomicNum() == 7]
            if len(nitrogen_neighbors) != 1:
                continue
            (n2,) = nitrogen_neighbors
            if mol.GetBondBetweenAtoms(n1.GetIdx(), n2.GetIdx()).GetBondTypeAsDouble() == 1.0:
                return True
    return False


def _diacyl_hydrazide_core(mol):
    """(n1, n2, c1, o1, c2, o2) if the molecule is a plain
    R-CO-NH-NH-CO-R' diacylhydrazide -- each bridging nitrogen bears
    exactly one H and no substituent besides its own acyl carbon and the
    N-N bond (P-66.3.3.3). None otherwise."""
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 7 or atom.GetDegree() != 2 or atom.GetTotalNumHs() != 1:
            continue
        neighbors = list(atom.GetNeighbors())
        carbons = [n for n in neighbors if n.GetAtomicNum() == 6]
        nitrogens = [n for n in neighbors if n.GetAtomicNum() == 7]
        if len(carbons) != 1 or len(nitrogens) != 1:
            continue
        c1, n2 = carbons[0], nitrogens[0]
        if not _is_carbonyl_carbon(mol, c1):
            continue
        if mol.GetBondBetweenAtoms(atom.GetIdx(), c1.GetIdx()).GetBondTypeAsDouble() != 1.0:
            continue
        if mol.GetBondBetweenAtoms(atom.GetIdx(), n2.GetIdx()).GetBondTypeAsDouble() != 1.0:
            continue
        if n2.GetDegree() != 2 or n2.GetTotalNumHs() != 1:
            continue
        n2_other = [n for n in n2.GetNeighbors() if n.GetIdx() != atom.GetIdx()]
        if len(n2_other) != 1 or n2_other[0].GetAtomicNum() != 6:
            continue
        c2 = n2_other[0]
        if not _is_carbonyl_carbon(mol, c2):
            continue
        if mol.GetBondBetweenAtoms(n2.GetIdx(), c2.GetIdx()).GetBondTypeAsDouble() != 1.0:
            continue

        def _carbonyl_oxygen(carbon):
            oxygens = [
                o
                for o in carbon.GetNeighbors()
                if o.GetAtomicNum() == 8
                and o.GetDegree() == 1
                and mol.GetBondBetweenAtoms(carbon.GetIdx(), o.GetIdx()).GetBondTypeAsDouble() == 2.0
            ]
            return oxygens[0].GetIdx() if len(oxygens) == 1 else None

        o1, o2 = _carbonyl_oxygen(c1), _carbonyl_oxygen(c2)
        if o1 is None or o2 is None:
            continue
        return atom.GetIdx(), n2.GetIdx(), c1.GetIdx(), o1, c2.GetIdx(), o2
    return None


def has_diacyl_hydrazide_shape(mol) -> bool:
    return _diacyl_hydrazide_core(mol) is not None


def _validate_and_collect_diacyl_hydrazide(mol):
    n1, n2, c1, o1, c2, o2 = _diacyl_hydrazide_core(mol)
    core_atoms = {n1, n2, o1, o2}

    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atomic_num not in (6, 7, 8):
            raise UnsupportedStructure(
                "heteroatoms other than a diacylhydrazide's own two "
                "bridging nitrogens and carbonyl oxygens (P-66.3.3.3) are "
                "out of scope for this module"
            )
        if atom.GetFormalCharge() != 0 or atom.GetIsotope() != 0:
            raise UnsupportedStructure("charged or isotopically modified atoms are not supported yet")
        if atomic_num == 6:
            if atom.GetIsAromatic():
                raise UnsupportedStructure(
                    "an aromatic acyl group (e.g. benzoyl) is out of scope "
                    "for this module"
                )
        elif atomic_num == 7:
            if atom.GetIdx() not in (n1, n2):
                raise UnsupportedStructure(
                    "a nitrogen other than the diacylhydrazide's own two "
                    "bridging nitrogens needs seniority handling not yet "
                    "implemented here"
                )
        elif atomic_num == 8:
            if atom.GetIdx() not in core_atoms:
                raise UnsupportedStructure(
                    "an oxygen other than the diacylhydrazide's own two "
                    "carbonyl oxygens (e.g. a coexisting hydroxyl or ester) "
                    "is out of scope for this module"
                )
    if mol.GetRingInfo().NumRings() > 0:
        raise UnsupportedStructure(
            "a cyclic diacylhydrazide is a structurally different shape "
            "and out of scope for this module"
        )
    if len(non_single_bonds(mol)) != 2:
        # Exactly the two acyl C=O bonds are always present; anything else
        # is chain unsaturation, out of scope for this narrow first pass.
        raise UnsupportedStructure(
            "chain unsaturation is out of scope for this module"
        )
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")
    for carbon_idx in (c1, c2):
        carbon_neighbors = [n for n in mol.GetAtomWithIdx(carbon_idx).GetNeighbors() if n.GetAtomicNum() == 6]
        if len(carbon_neighbors) > 1:
            raise UnsupportedStructure(
                "an acyl carbon with more than one carbon neighbor is not "
                "a valid diacylhydrazide acyl carbon"
            )
    return c1, c2


def _acyl_chain_length(mol, acyl_carbon):
    length = linear_branch(carbon_adjacency(mol), acyl_carbon, None)
    if length is None:
        raise UnsupportedStructure("a branched acyl chain is out of scope for this module")
    return length


def _diacyl_prefix_name(chain_length):
    if chain_length == 1:
        return "formyl"
    if chain_length == 2:
        return "acetyl"
    return alkane_name(chain_length)[:-1] + "oyl"


def _diacyl_parent_name(chain_length):
    if chain_length == 1:
        return "formohydrazide"
    if chain_length == 2:
        return "acetohydrazide"
    return alkane_name(chain_length) + "hydrazide"


def name_diacyl_hydrazide(mol) -> str:
    c1, c2 = _validate_and_collect_diacyl_hydrazide(mol)
    length1, length2 = _acyl_chain_length(mol, c1), _acyl_chain_length(mol, c2)
    if length1 != length2:
        raise UnsupportedStructure(
            "an unsymmetric diacylhydrazide (the two acyl groups name "
            "differently) needs Table 3.3 acid-seniority handling not yet "
            "implemented here"
        )
    return f"N'-{_diacyl_prefix_name(length1)}{_diacyl_parent_name(length2)}"


def _validate_and_collect_hydrazide(mol, aromatic_ring_atoms=frozenset()):
    """Check the molecule fits this module's scope (see module docstring)
    and return (hydrazide_carbon, hydrazide_oxygen, n1, n2, n1_alkyl,
    n2_alkyl, hydroxyls): the single -CO-N(R)-N(R')(R'') carbon/oxygen/
    n1/n2 atom indices, n1's 0-1 and n2's 0-2 N-alkyl substituent carbon
    indices, and the set of any coexisting standalone hydroxyl-oxygen
    indices.

    `aromatic_ring_atoms`: atom indices already independently verified (by
    the caller, before this function runs) to form a single plain benzene
    ring with exactly one exocyclic attachment -- exempted from the
    aromatic-atom rejection below so `name_hydrazide`'s benzene-ring-
    substituent path (see `_name_phenyl_chain_hydrazide`) can reuse this
    same validation for the rest of the molecule. Empty by default, so
    every other caller's behavior is unchanged."""
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
            if atom.GetIsAromatic() and atom.GetIdx() not in aromatic_ring_atoms:
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
            neighbors = list(atom.GetNeighbors())
            if len(neighbors) > 3:
                raise UnsupportedStructure(
                    "a hydrazide nitrogen with more than two substituents "
                    "besides its ring/chain neighbor is not a valid "
                    "hydrazide nitrogen"
                )
            if any(n.GetAtomicNum() not in (6, 7) for n in neighbors):
                raise UnsupportedStructure(
                    "a hydrazide nitrogen bonded to anything other than "
                    "carbon or its partner nitrogen is out of scope for "
                    "this module"
                )
            for bond in atom.GetBonds():
                if bond.GetBondTypeAsDouble() != 1.0:
                    raise UnsupportedStructure(
                        "a hydrazide nitrogen must be singly bonded to all its neighbors"
                    )
            carbon_neighbors = [n for n in neighbors if n.GetAtomicNum() == 6]
            nitrogen_neighbors = [n for n in neighbors if n.GetAtomicNum() == 7]
            carbonyl_neighbors = [n for n in carbon_neighbors if _is_carbonyl_carbon(mol, n)]
            alkyl_neighbors = [n for n in carbon_neighbors if n not in carbonyl_neighbors]
            is_n1_shape = len(carbonyl_neighbors) == 1 and len(nitrogen_neighbors) == 1 and len(alkyl_neighbors) <= 1
            is_n2_shape = not carbonyl_neighbors and len(nitrogen_neighbors) == 1
            if is_n1_shape:
                (n2,) = nitrogen_neighbors
                n2_neighbors = list(n2.GetNeighbors())
                if len(n2_neighbors) > 3:
                    raise UnsupportedStructure(
                        "a hydrazide's terminal nitrogen with more than "
                        "two alkyl substituents is not a valid hydrazide "
                        "nitrogen"
                    )
                if any(n.GetIdx() != atom.GetIdx() and n.GetAtomicNum() != 6 for n in n2_neighbors):
                    raise UnsupportedStructure(
                        "a hydrazide's terminal nitrogen bonded to "
                        "anything other than carbon (besides its own N-N "
                        "bond) is out of scope for this module"
                    )
                for bond in n2.GetBonds():
                    if bond.GetBondTypeAsDouble() != 1.0:
                        raise UnsupportedStructure(
                            "a hydrazide nitrogen must be singly bonded to all its neighbors"
                        )
                n1_alkyl = tuple(n.GetIdx() for n in alkyl_neighbors)
                n2_alkyl = tuple(n.GetIdx() for n in n2_neighbors if n.GetIdx() != atom.GetIdx())
                n1_by_carbon.setdefault(carbonyl_neighbors[0].GetIdx(), []).append(
                    (atom.GetIdx(), n2.GetIdx(), n1_alkyl, n2_alkyl)
                )
                continue
            if is_n2_shape:
                # Already validated as part of its N1 partner's own check
                # above; nothing further to do for this atom in isolation.
                continue
            raise UnsupportedStructure(
                "a nitrogen shaped like neither a hydrazide's -N(R)- nor "
                "its terminal -N(R')(R'') (P-66.3.1.1) is out of scope for "
                "this module (e.g. a plain primary amide/amine nitrogen "
                "bonded directly to a carbon; see _amide.py for a plain "
                "primary amide)"
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
    ((n1, n2, n1_alkyl, n2_alkyl),) = n1_by_carbon[hydrazide_carbon]
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")
    return hydrazide_carbon, hydrazide_oxygen, n1, n2, n1_alkyl, n2_alkyl, hydroxyls


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


def _collect_n_alkyl(full_carbon_graph, mol, hydroxyls, graph, n1_alkyl, n2_alkyl):
    """Validate n1's (0-1) and n2's (0-2) alkyl substituents -- each must
    be a plain, unbranched, unsubstituted, saturated alkyl chain, the same
    restriction `_amide.py` places on its own N-substituents -- and return
    (entries, n_alkyl_atoms): `entries` is a list of ("N"/"N'", name)
    pairs (module docstring's N/N' convention) ready for `_format_n_prefix`,
    and `n_alkyl_atoms` is the full set of atom indices spanned by every
    N-substituent, to exclude from the principal-chain search below."""
    entries = []
    n_alkyl_atoms = set()
    for locant_label, alkyl_roots in (("N", n1_alkyl), ("N'", n2_alkyl)):
        for root in alkyl_roots:
            length = linear_branch(full_carbon_graph, root, None)
            if length is None:
                raise UnsupportedStructure("a branched N-substituent is not supported yet")
            atoms = set()
            previous, current = None, root
            while current is not None:
                atoms.add(current)
                neighbors = [n for n in full_carbon_graph[current] if n != previous]
                previous, current = current, (neighbors[0] if neighbors else None)
            if any(b[0] in atoms or b[1] in atoms for b in non_single_bonds(mol)):
                raise UnsupportedStructure("an unsaturated N-substituent is not supported yet")
            if any(next(iter(graph[o])) in atoms for o in hydroxyls):
                raise UnsupportedStructure(
                    "a substituted N-substituent (e.g. bearing a hydroxyl) "
                    "is not supported yet; only a plain, unsubstituted "
                    "alkyl N-substituent is in scope"
                )
            entries.append((locant_label, alkyl_name(length)))
            n_alkyl_atoms |= atoms
    return entries, n_alkyl_atoms


def _format_n_prefix(entries):
    """entries: [("N"/"N'", name), ...] -> the assembled "N-"/"N'-" prefix
    string (module docstring), grouping identically-named substituents
    (whether on the same nitrogen or split across both) under one shared
    multiplying prefix, e.g. [("N", "methyl"), ("N'", "methyl")] ->
    "N,N'-dimethyl". '' if entries is empty."""
    if not entries:
        return ""
    grouped = {}
    for locant, name in entries:
        grouped.setdefault(name, []).append(locant)
    parts = []
    for name in sorted(grouped, key=alpha_sort_key):
        locants = sorted(grouped[name])
        multiplier = multiplying_prefix(len(locants)) if len(locants) > 1 else ""
        parts.append(f"{','.join(locants)}-{multiplier}{name}")
    return "-".join(parts)


def _name_acyclic_hydrazide(mol, hydrazide_carbon, excluded, hydroxyls, bonds, n1_alkyl, n2_alkyl, stereo=None):
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
    full_carbon_graph = carbon_adjacency(mol)
    n_entries, n_alkyl_atoms = _collect_n_alkyl(full_carbon_graph, mol, hydroxyls, graph, n1_alkyl, n2_alkyl)
    n_prefix = _format_n_prefix(n_entries)

    # N-alkyl substituent carbons hang off the (excluded) hydrazide
    # nitrogens, not off any acyl-chain carbon, so they form their own
    # isolated component(s) in the carbon-only graph; the whole subtree
    # must be removed before picking the longest chain, or a longer
    # N-substituent would be mistaken for the acyl chain itself (same
    # issue `_amide.py` guards against).
    carbon_graph = {k: v for k, v in full_carbon_graph.items() if k not in n_alkyl_atoms}
    chains = longest_chains(carbon_graph)
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
            # single carbon for any substituent. Unlike the dinuclear case
            # below, a substituent on either nitrogen also breaks the
            # retained name itself (PubChem switches to 'formamide' +
            # amino-substituent naming for `O=CN(C)N`/`O=CNNC`, a
            # different naming paradigm this module doesn't attempt), so
            # N-substitution is out of scope here specifically.
            if n_entries:
                raise UnsupportedStructure(
                    "an N-/N'-substituted formohydrazide is not supported "
                    "yet (the retained name itself doesn't survive "
                    "substitution here, unlike the dinuclear case)"
                )
            return "formohydrazide"
        # P-66.3.1.2.1: the dinuclear case's retained name
        # ('acetohydrazide') is the PIN, with any substituent on its
        # terminal carbon cited as an ordinary prefix (e.g.
        # '2-chloroacetohydrazide', PubChem CID 101883).
        substituents = _substituents_for_chain(graph, chain, halogens, excluded)
        grouped = group_substituents(substituents)
        prefix = format_substituent_prefixes(grouped)
        name = prefix + "acetohydrazide"
        if n_prefix:
            separator = "-" if name[0].isdigit() else ""
            name = f"{n_prefix}{separator}{name}"
        return name

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

    if n_prefix:
        separator = "-" if best_name[0].isdigit() else ""
        best_name = f"{n_prefix}{separator}{best_name}"

    if stereo is not None:
        labels = sorted((best_position_of[atom], code) for atom, code in stereo)
        prefix = ",".join(f"{locant}{code}" for locant, code in labels)
        return f"({prefix})-{best_name}"
    return best_name


def _name_phenyl_chain_hydrazide(mol, ring_atoms):
    """Name a hydrazide whose -CO-NH-NH2 lies entirely on a single
    unbranched chain hanging off one atom of an otherwise-plain,
    unsubstituted benzene ring -- e.g. 3-phenylpropanehydrazide, or
    '2-phenylacetohydrazide' for the dinuclear retained-name case
    (mirroring the module's own '2-chloroacetohydrazide' pattern). The
    ring is cited as a 'phenyl' substituent prefix on the chain, which is
    the parent hydride, mirroring `_amide.py`'s
    `_name_phenyl_chain_amide`. Narrower than the acyclic path above: no
    N-/N'-alkyl substitution, no coexisting standalone hydroxyl, no chain
    unsaturation, and no specified stereocenter -- each is a separate
    follow-up (see
    tasks/phenyl-substituent-on-hydrazide-chain.md's scope note)."""
    hydrazide_carbon, hydrazide_oxygen, n1, n2, n1_alkyl, n2_alkyl, hydroxyls = _validate_and_collect_hydrazide(
        mol, aromatic_ring_atoms=ring_atoms
    )
    if n1_alkyl or n2_alkyl:
        raise UnsupportedStructure(
            "an N-/N'-alkyl-substituted hydrazide alongside a benzene-ring "
            "substituent is not supported yet"
        )
    if hydroxyls:
        raise UnsupportedStructure(
            "a standalone hydroxyl alongside a benzene-ring-substituent "
            "hydrazide chain is not supported yet"
        )
    if specified_stereocenters(mol):
        raise UnsupportedStructure(
            "a specified stereocenter alongside a benzene-ring-substituent "
            "hydrazide chain is not supported yet"
        )
    excluded = {hydrazide_oxygen, n1, n2}
    non_ring_unsaturation = [
        b
        for b in non_single_bonds(mol)
        if b[0] not in excluded and b[1] not in excluded and b[0] not in ring_atoms and b[1] not in ring_atoms
    ]
    if non_ring_unsaturation:
        raise UnsupportedStructure(
            "chain unsaturation alongside a benzene-ring-substituent "
            "hydrazide chain is not supported yet"
        )

    graph = adjacency(mol)
    halogens = {**halogen_substituents(mol), **plain_alkyl_ring_substituents(mol, graph, ring_atoms)}
    attachment = ring_chain_attachment_with_halogens(graph, ring_atoms, set(), halogens)
    if attachment is None:
        raise UnsupportedStructure(
            "a benzene ring with more than one non-halogen, non-alkyl "
            "exocyclic substituent alongside a chain hydrazide is not "
            "supported yet"
        )
    ring_atom, chain_root = attachment
    chain = ordered_chain(graph, chain_root, ring_atom, excluded)
    if chain is None:
        raise UnsupportedStructure(
            "a branched chain hanging off the benzene ring alongside a "
            "hydrazide is not supported yet"
        )
    if chain[-1] != hydrazide_carbon:
        raise UnsupportedStructure(
            "the hydrazide carbon must be the chain's far terminus from "
            "the benzene ring for this benzene-substituent path"
        )
    if len(chain) < 2:
        raise UnsupportedStructure(
            "a hydrazide directly attached to the benzene ring "
            "(benzohydrazide-style naming) uses a separate construction, "
            "out of scope for this acyclic-chain-parent module"
        )

    ordered = list(reversed(chain))
    chain_length = len(ordered)
    position_of = {atom: i + 1 for i, atom in enumerate(ordered)}
    substituents = {
        position_of[chain_root]: [name_branch(graph, ring_atom, chain_root, halogens, ring_atoms)]
    }
    grouped = group_substituents(substituents)
    if chain_length == 2:
        # P-66.3.1.2.1: the dinuclear case's retained name ('acetohydrazide')
        # is the PIN, with the substituent cited as an ordinary prefix
        # (see module docstring's '2-chloroacetohydrazide').
        prefix = format_substituent_prefixes(grouped)
        return prefix + "acetohydrazide"
    return _name_from_substituents(chain_length, [], [], grouped)


def name_hydrazide(mol) -> str:
    if has_diacyl_hydrazide_shape(mol):
        return name_diacyl_hydrazide(mol)
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() == 1:
        ring_atoms = set(ring_info.AtomRings()[0])
        if is_plain_benzene_ring(mol, ring_atoms):
            return _name_phenyl_chain_hydrazide(mol, ring_atoms)
    hydrazide_carbon, hydrazide_oxygen, n1, n2, n1_alkyl, n2_alkyl, hydroxyls = _validate_and_collect_hydrazide(mol)
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

    return _name_acyclic_hydrazide(mol, hydrazide_carbon, excluded, hydroxyls, bonds, n1_alkyl, n2_alkyl, stereo)
