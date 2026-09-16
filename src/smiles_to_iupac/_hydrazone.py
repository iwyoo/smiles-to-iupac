"""Naming of simple hydrazones (R2C=N-NH2, terminal nitrogen unsubstituted)
on acyclic saturated carbon chains, per the IUPAC 2013 Recommendations
("the Blue Book"):

- P-68.3.1.2.2 (Chapter P-6a, https://iupac.qmul.ac.uk/BlueBook/PDF/P6a.pdf):
  a hydrazone is not given its own suffix here -- the PIN treats it as
  `_hydrazine.py`'s hydrazine parent (H2N-NH2) bearing a "-ylidene"
  (double-bond) substituent prefix on one nitrogen instead of an ordinary
  "-yl" alkyl prefix, reusing the same carbon-side chain search
  `_imine.py` already has for its own C=N group (longest chain containing
  the double-bond carbon, tried in both directions, lowest achievable
  locant wins). Confirmed via PubChem PUG REST: CID 81125 (`C=NN`) ->
  "methylidenehydrazine", CID 53651639 (`CC=NN`) -> "ethylidenehydrazine",
  CID 53726800 (`CCC=NN`) -> "propylidenehydrazine", CID 78937
  (`CC(C)=NN`) -> "propan-2-ylidenehydrazine", CID 71356052
  (`CCC(C)=NN` / `CC(=NN)CC`) -> "butan-2-ylidenehydrazine".
- Unlike `_imine.py`'s own "-imine" suffix (which always cites a locant
  for a 3+-carbon chain, even at position 1, e.g. "hexan-1-imine"), the
  "-ylidene" prefix here follows the ordinary substituent-prefix rule
  instead (P-29.2: the free valence is always locant 1 and never cited):
  when the double-bond carbon lands at position 1 of its chain (an
  "aldehyde-shaped" attachment, 0 or 1 carbon neighbors), this module
  delegates straight to `_substituents.py`'s `name_branch` -- the exact
  same old-style (non-systematic, e.g. "propyl" not "propan-1-yl")
  numbering and branch-naming `name_branch` already produces for an
  ordinary substituent rooted at its own attachment point -- and appends
  "idene" (turning the "-yl" ending into "-ylidene"). Confirmed via
  PubChem PUG REST: CID 71356051 (`CC(C)C=NN`, a branched aldehyde
  hydrazone) -> "2-methylpropylidenehydrazine" (the same "2-methylpropyl"
  old-style compound name `name_branch` already produces elsewhere, with
  "idene" appended). When the double-bond carbon instead lands at
  position 2+ (a "ketone-shaped" attachment, 2 carbon neighbors, so it
  can never be a chain endpoint), this module falls back to the
  systematic "alkan-N-ylidene" form directly, mirroring `_imine.py`'s own
  chain-substituent/locant assembly with the suffix swapped.
- P-35.2.1: a halogen substituent on the carbon chain coexists freely,
  reusing `halogen_substituents`/`name_branch` unchanged. Confirmed via
  PubChem PUG REST: CID 163551743 (`ClCC=NN`) ->
  "2-chloroethylidenehydrazine", CID 174934626 (`ClCCC=NN`) ->
  "3-chloropropylidenehydrazine".
- P-91.2(e)/P-93.1 (E/Z stereo): when this module's own C=N bond has its
  geometry specified in the input, a bare `"(E)-"`/`"(Z)-"` prefix (no
  locant, confirmed against PubChem's own auto-generated name even for a
  ketone-shaped attachment whose own name already carries one, e.g.
  `"(E)-butan-2-ylidenehydrazine"`) is added via
  `_common.specified_double_bond_stereo`, the same helper `_imine.py`
  uses for its own C=N. A prior revision of this docstring claimed
  RDKit's `Chem.FindPotentialStereo` *always* flags this bond as an
  unspecified potential stereo element regardless of input -- that claim
  no longer holds against the currently pinned RDKit (confirmed:
  `C/C=N/N` now reports it as specified) and was corrected here, mirroring
  the identical correction already made in `_imine.py`. A specified
  tetrahedral chain stereocenter coexisting with the C=N bond (specified
  or not) is still rejected, since combining the two kinds of descriptor
  is out of scope for this module (P-92/P-93). `_azine.py`/`_amidine.py`
  share the identical (now-stale) rejection pattern independently -- not
  updated by this change.

Explicitly out of scope (raise `UnsupportedStructure`):
- An N-substituted hydrazone (=N-NH-R, the terminal nitrogen bearing an
  alkyl group instead of being a plain -NH2): a completely different
  naming scheme, confirmed via PubChem PUG REST to *not* follow the
  ylidenehydrazine pattern -- CID 86520 (`CC=NNC`) ->
  "N-(ethylideneamino)methanamine", an amine-parent name. This needs
  `_amine.py` support and is a separate, unverified follow-up.
- An azine (R2C=N-N=CR2, both hydrazine nitrogens double-bonded to their
  own carbon): also a different scheme -- CID 79085 (`CC(C)=NN=C(C)C`)
  -> "N-(propan-2-ylideneamino)propan-2-imine", an imine-parent name.
  Separate, unverified follow-up.
- Semicarbazide (H2N-NH-C(=O)-NH2) and other hydrazone-adjacent
  P-68.3.1.2.x shapes involving a carbonyl: this project has no urea
  parent-hydride module yet (PubChem's own PIN for the parent,
  `NC(=O)NN`, is "aminourea"), so this is a separate, much larger
  follow-up.
- More than one C=N double bond, any ring, any aromatic atom, any
  unsaturation besides the hydrazone's own C=N, any heteroatom other than
  the hydrazone's own two nitrogens and halogen substituents, and charged
  or isotopically modified atoms -- all unverified in this first pass,
  mirroring `_imine.py`'s own equivalent exclusions.
"""

from rdkit import Chem

from ._common import (
    HALOGEN_PREFIXES,
    UnsupportedStructure,
    adjacency,
    carbon_adjacency,
    group_substituents,
    halogen_substituents,
    longest_chains,
    lowest_locant_set,
    specified_double_bond_stereo,
)
from ._numerals import alkane_name
from ._substituents import alpha_sort_key, format_substituent_prefixes, name_branch

_ALLOWED_ATOMIC_NUMS = {6, 7, *HALOGEN_PREFIXES}


def _find_hydrazone_bond(mol):
    return [
        bond
        for bond in mol.GetBonds()
        if bond.GetBondTypeAsDouble() == 2.0
        and {bond.GetBeginAtom().GetAtomicNum(), bond.GetEndAtom().GetAtomicNum()} == {6, 7}
    ]


def has_hydrazone_shape(mol) -> bool:
    """A C=N double bond whose nitrogen has exactly one other, single-bonded
    neighbor which is itself a plain terminal nitrogen (degree 1) -- the
    R2C=N-NH2 shape this module handles. A more complete validation (atom
    types, rings, extra unsaturation, ...) happens in `name_hydrazone`."""
    for bond in _find_hydrazone_bond(mol):
        a, b = bond.GetBeginAtom(), bond.GetEndAtom()
        nitrogen = a if a.GetAtomicNum() == 7 else b
        if nitrogen.GetDegree() != 2:
            continue
        carbon_idx = (a if a.GetAtomicNum() == 6 else b).GetIdx()
        (other,) = [n for n in nitrogen.GetNeighbors() if n.GetIdx() != carbon_idx]
        if other.GetAtomicNum() != 7:
            continue
        nn_bond = mol.GetBondBetweenAtoms(nitrogen.GetIdx(), other.GetIdx())
        if nn_bond.GetBondTypeAsDouble() == 1.0 and other.GetDegree() == 1:
            return True
    return False


def _validate_and_find_hydrazone(mol):
    hydrazone_bonds = _find_hydrazone_bond(mol)
    if not hydrazone_bonds:
        raise UnsupportedStructure(
            "no hydrazone (C=N-NH2) group found; this module only handles hydrazones"
        )
    if len(hydrazone_bonds) > 1:
        raise UnsupportedStructure("more than one C=N double bond is not supported yet")
    (bond,) = hydrazone_bonds
    a, b = bond.GetBeginAtom(), bond.GetEndAtom()
    carbon = a if a.GetAtomicNum() == 6 else b
    imine_n = a if a.GetAtomicNum() == 7 else b

    if imine_n.GetDegree() != 2:
        raise UnsupportedStructure(
            "the hydrazone's C=N nitrogen must have exactly one other substituent"
        )
    (other,) = [n for n in imine_n.GetNeighbors() if n.GetIdx() != carbon.GetIdx()]
    if other.GetAtomicNum() != 7:
        raise UnsupportedStructure(
            "no hydrazone (C=N-NH2) group found; this module only handles hydrazones"
        )
    nn_bond = mol.GetBondBetweenAtoms(imine_n.GetIdx(), other.GetIdx())
    if nn_bond.GetBondTypeAsDouble() != 1.0:
        raise UnsupportedStructure("the hydrazone's N-N bond must be a single bond")
    if other.GetDegree() != 1:
        raise UnsupportedStructure(
            "an N-substituted hydrazone (=N-NH-R, P-68.3.1.2.2's more general "
            "case) is out of scope for this module"
        )
    if imine_n.GetFormalCharge() != 0 or other.GetFormalCharge() != 0:
        raise UnsupportedStructure("charged atoms are not supported yet")

    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in _ALLOWED_ATOMIC_NUMS:
            raise UnsupportedStructure(
                "heteroatoms other than the hydrazone's own two nitrogens "
                "and halogen substituents (P-35.2.1) are not supported yet"
            )
        if atom.GetFormalCharge() != 0 or atom.GetIsotope() != 0:
            raise UnsupportedStructure("charged or isotopically modified atoms are not supported yet")
        if atom.GetAtomicNum() == 6 and atom.GetIsAromatic():
            raise UnsupportedStructure("an aromatic carbon skeleton is out of scope for this module")
    if mol.GetRingInfo().NumRings() != 0:
        raise UnsupportedStructure("a ring anywhere in the molecule is out of scope for this module")
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")

    other_non_single = [
        b for b in mol.GetBonds() if b.GetBondTypeAsDouble() != 1.0 and b.GetIdx() != bond.GetIdx()
    ]
    if other_non_single:
        raise UnsupportedStructure(
            "a bond order other than single (besides the hydrazone C=N itself) "
            "is not supported yet"
        )

    carbon_neighbors = [n for n in carbon.GetNeighbors() if n.GetAtomicNum() == 6]
    if len(carbon_neighbors) not in (0, 1, 2):
        raise UnsupportedStructure(
            "a hydrazone carbon must have zero, one, or two carbon neighbors"
        )

    return carbon.GetIdx(), imine_n.GetIdx()


def _substituents_for_chain(graph, chain, halogens, exclude, mol=None):
    chain_set = set(chain)
    substituents = {}
    for position, atom in enumerate(chain, start=1):
        branch_roots = [n for n in graph[atom] if n not in chain_set and n not in exclude]
        if branch_roots:
            substituents[position] = [name_branch(graph, root, atom, halogens, mol=mol) for root in branch_roots]
    return substituents


def _ylidene_name(chain_length, locant, substituents):
    grouped = group_substituents(substituents)
    prefix = format_substituent_prefixes(grouped)
    stem = alkane_name(chain_length)[:-1]
    locant_set = lowest_locant_set(loc for info in grouped.values() for loc in info["locants"])
    citation_locants = tuple(
        loc
        for name in sorted(grouped, key=alpha_sort_key)
        for loc in sorted(grouped[name]["locants"])
    )
    name = f"{prefix}{stem}-{locant}-ylidene"
    return (locant, locant_set, citation_locants, name), name


def _name_hydrazone_carbon(mol, carbon_idx, imine_n_idx):
    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    chains = longest_chains(carbon_adjacency(mol))

    eligible = [chain for chain in chains if carbon_idx in chain]
    if not eligible:
        raise UnsupportedStructure(
            "the hydrazone carbon does not lie on any longest carbon chain; "
            "a shorter principal chain (P-44.1.1) is not supported yet"
        )

    best_locant = min(
        position_of[carbon_idx]
        for chain in eligible
        for position_of in (
            {atom: i + 1 for i, atom in enumerate(chain)},
            {atom: i + 1 for i, atom in enumerate(reversed(chain))},
        )
    )

    if best_locant == 1:
        # P-29.2: the free valence is always locant 1 of a substituent
        # prefix and never cited -- the exact numbering/branch-naming
        # `name_branch` already applies for an ordinary '-yl' substituent
        # rooted at its own attachment point (see module docstring).
        branch_name, _ = name_branch(graph, carbon_idx, imine_n_idx, halogens, mol=mol)
        return branch_name + "idene"

    chain_length = len(chains[0])
    best_key = None
    best_name = None
    for chain in eligible:
        for candidate in (chain, list(reversed(chain))):
            position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
            locant = position_of[carbon_idx]
            substituents = _substituents_for_chain(graph, candidate, halogens, {imine_n_idx}, mol=mol)
            key, name = _ylidene_name(chain_length, locant, substituents)
            if best_key is None or key < best_key:
                best_key, best_name = key, name
    return best_name


def name_hydrazone(mol) -> str:
    carbon_idx, imine_n_idx = _validate_and_find_hydrazone(mol)
    stereo = specified_double_bond_stereo(mol)
    name = _name_hydrazone_carbon(mol, carbon_idx, imine_n_idx) + "hydrazine"
    if stereo is not None:
        # PubChem registers this stereo separately and its own
        # auto-generated name already cites it as a bare "(E)-"/"(Z)-"
        # prefix with no locant, even for a "ketone-shaped" attachment
        # whose own name already carries one (e.g. "(E)-butan-2-
        # ylidenehydrazine") -- confirmed directly (module docstring).
        ((_, code),) = stereo
        name = f"({code})-{name}"
    return name
