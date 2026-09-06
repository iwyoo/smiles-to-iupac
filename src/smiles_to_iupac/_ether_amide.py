"""Naming of a molecule combining a plain ether (-O-R') with exactly one
separate unsubstituted primary amide (-CONH2) on an acyclic saturated
carbon skeleton, per the IUPAC 2013 Recommendations ("the Blue Book"):

- P-41 Table 4.1: ethers have no principal-characteristic-group suffix at
  all (class 41, junior even to plain carbon compounds) -- exactly the
  same reasoning `_ether_amine.py`/`_ether_thiol.py`/`_ether_ketone.py`/
  `_ether_aldehyde.py` document in full: an ether oxygen is *always* the
  'R-oxy' substituent prefix (P-63.2.2.1.1), never the parent, so there is
  no seniority competition to resolve. This module mirrors
  `_ether_aldehyde.py`'s structure, swapping in `_amide.py`'s own
  chain-naming machinery.
- `COCC(N)=O` -> PubChem's own '2-methoxyacetamide' confirms the amide is
  always the suffix parent, the ether always the 'alkoxy' prefix (this
  project keeps its existing 'ethanamide' stem convention rather than
  PubChem's retained '-acetamide', matching `_aldehyde.py`'s analogous
  'ethanal' vs. 'acetaldehyde' choice).
- `_amide.py`'s `_name_acyclic_amide` already had a mechanism to remove a
  substituent's own carbon component before the principal-chain search
  (its N-alkyl substituents' `n_substituent_atoms`); it gained one new
  optional `extra_excluded_carbons` parameter (merged into that same
  exclusion set) so this module can remove the ether's alkoxy-branch
  component the identical way, rather than the whole-graph-replacement
  `carbon_graph` parameter `_thiol.py`/`_ketone.py`/`_aldehyde.py` each
  gained (PR #428-430) -- `_amide.py` already builds its own graph, so
  extending its existing exclusion set was the smaller change.

Scope, deliberately narrow (mirrors `_ether_aldehyde.py`): exactly one
plain ether oxygen (both sides acyclic saturated carbon) plus exactly one
unsubstituted primary amide (-CONH2, no N-alkyl), on one acyclic
*saturated* skeleton, halogens allowed. Explicitly out of scope (raise
`UnsupportedStructure`): more than one ether oxygen or amide, any
N-alkyl-substituted amide, any other heteroatom (including a coexisting
hydroxyl), any ring, any chain unsaturation (ene/yne) besides the amide's
own C=O, and any specified stereocenter.
"""

from rdkit import Chem

from ._amide import _is_carbonyl_carbon, _name_acyclic_amide
from ._common import (
    HALOGEN_PREFIXES,
    UnsupportedStructure,
    adjacency,
    bfs,
    carbon_adjacency,
    non_single_bonds,
    specified_stereocenters,
)
from ._ether import _oxy_prefix
from ._substituents import name_branch

_ALLOWED_ATOMIC_NUMS = {6, 7, 8, *HALOGEN_PREFIXES}


def _find_ether_oxygens(mol):
    """The molecule's sole ether-shaped oxygen among possibly several --
    see `_ether_ketone.py`'s identical helper (an amide molecule has a
    second oxygen, its own carbonyl)."""
    return [
        atom
        for atom in mol.GetAtoms()
        if atom.GetAtomicNum() == 8
        and atom.GetDegree() == 2
        and all(n.GetAtomicNum() == 6 for n in atom.GetNeighbors())
    ]


def _find_amide(mol):
    """Return (amide_carbon, amide_oxygen, amide_nitrogen) for a single
    unsubstituted primary amide (-CONH2), or None -- mirrors
    `_amide_amine.py`'s identical helper, plus one extra check: any
    carbonyl-carbon neighbor other than the carbonyl oxygen and this
    amide nitrogen must be carbon (or absent, the formamide shape).
    Without that check, a carbamate's R-O-C(=O)-NH2 carbon would also
    match (it has a carbonyl oxygen and a qualifying -NH2 too), since
    `_is_carbonyl_carbon` doesn't look at the carbon's other substituents
    -- the R-O- ether oxygen there is exactly the coexisting-ether shape
    this module looks for, so without this check every carbamate would be
    misidentified as an ether/amide pair."""
    found = None
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6 or not _is_carbonyl_carbon(mol, atom):
            continue
        for n in atom.GetNeighbors():
            if (
                n.GetAtomicNum() == 7
                and n.GetDegree() == 1
                and n.GetTotalNumHs() == 2
                and mol.GetBondBetweenAtoms(atom.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 1.0
            ):
                (oxygen,) = [
                    o
                    for o in atom.GetNeighbors()
                    if o.GetAtomicNum() == 8
                    and o.GetDegree() == 1
                    and mol.GetBondBetweenAtoms(atom.GetIdx(), o.GetIdx()).GetBondTypeAsDouble() == 2.0
                ]
                others = [
                    nb for nb in atom.GetNeighbors() if nb.GetIdx() not in (oxygen.GetIdx(), n.GetIdx())
                ]
                if any(nb.GetAtomicNum() != 6 for nb in others):
                    continue
                if found is not None:
                    return None
                found = (atom.GetIdx(), oxygen.GetIdx(), n.GetIdx())
    return found


def has_ether_amide_shape(mol) -> bool:
    if mol.GetRingInfo().NumRings() > 0:
        return False
    ethers = _find_ether_oxygens(mol)
    if len(ethers) != 1:
        return False
    amide = _find_amide(mol)
    if amide is None:
        return False
    total_oxygens = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() == 8)
    return total_oxygens == 2


def _validate(mol, ether_oxygen, amide_oxygen, amide_nitrogen):
    has_carbon = False
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atomic_num not in _ALLOWED_ATOMIC_NUMS:
            raise UnsupportedStructure(
                "heteroatoms other than a plain ether oxygen (P-63.2.1) "
                "and an unsubstituted primary amide's oxygen/nitrogen "
                "(P-66.1) are not supported yet"
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
            if atom.GetIdx() not in (ether_oxygen, amide_oxygen):
                raise UnsupportedStructure(
                    "an oxygen that isn't the single plain ether oxygen or "
                    "the single amide carbonyl is out of scope for this "
                    "module"
                )
        elif atomic_num == 7:
            if atom.GetIdx() != amide_nitrogen:
                raise UnsupportedStructure(
                    "a nitrogen that isn't the amide's own unsubstituted "
                    "-CONH2 nitrogen is out of scope for this module"
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
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")


def name_ether_amide(mol) -> str:
    if mol.GetRingInfo().NumRings() > 0:
        raise UnsupportedStructure(
            "an ether/amide combination on/in a ring uses a different "
            "naming construction, out of scope for this acyclic-only module"
        )
    ethers = _find_ether_oxygens(mol)
    if len(ethers) != 1:
        raise UnsupportedStructure(
            "exactly one plain ether oxygen coexisting with an amide is "
            "supported here (see _ether.py for a plain ether)"
        )
    (ether_oxygen,) = ethers
    ether_oxygen_idx = ether_oxygen.GetIdx()

    amide = _find_amide(mol)
    if amide is None:
        raise UnsupportedStructure(
            "exactly one unsubstituted primary amide (-CONH2) coexisting "
            "with the single ether is supported here (see _amide.py for a "
            "plain amide)"
        )
    amide_carbon, amide_oxygen, amide_nitrogen = amide
    _validate(mol, ether_oxygen_idx, amide_oxygen, amide_nitrogen)

    excluded = {amide_oxygen, amide_nitrogen}
    other_unsaturation = [b for b in non_single_bonds(mol) if b[0] not in excluded and b[1] not in excluded]
    if other_unsaturation:
        raise UnsupportedStructure(
            "chain unsaturation (ene/yne) combined with an ether/amide "
            "coexistence is out of scope for this module (1st-pass scope: "
            "saturated only)"
        )
    if specified_stereocenters(mol):
        raise UnsupportedStructure(
            "a specified stereocenter alongside an ether/amide "
            "coexistence is not supported yet"
        )

    full_graph = adjacency(mol)
    full_carbon_graph = carbon_adjacency(mol)
    ether_carbons = [n.GetIdx() for n in ether_oxygen.GetNeighbors()]
    amide_component, _ = bfs(full_carbon_graph, amide_carbon)
    main_side = [c for c in ether_carbons if c in amide_component]
    other_side = [c for c in ether_carbons if c not in amide_component]
    if len(main_side) != 1 or len(other_side) != 1:
        raise UnsupportedStructure(
            "the ether oxygen must sit between the amide's own carbon "
            "skeleton and a separate alkoxy branch for this module"
        )
    (r_prime_carbon,) = other_side

    r_prime_component, _ = bfs(full_carbon_graph, r_prime_carbon)

    sub_name, sub_compound = name_branch(full_graph, r_prime_carbon, ether_oxygen_idx, {}, mol=mol)
    if sub_compound:
        sub_name = f"({sub_name})"
    oxy_term = _oxy_prefix(sub_name)
    extra_names = {ether_oxygen_idx: oxy_term}

    return _name_acyclic_amide(
        mol,
        amide_carbon,
        amide_nitrogen,
        excluded,
        (),
        set(),
        (),
        extra_names=extra_names,
        extra_excluded_carbons=frozenset(r_prime_component),
    )
