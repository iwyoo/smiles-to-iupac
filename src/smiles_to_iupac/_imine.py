"""Naming of imines (the '-imine' suffix, C=N-H or C=N-R with a carbon on
the other end of the double bond) on acyclic saturated carbon chains, per
the IUPAC 2013 Recommendations ("the Blue Book"):

- P-62.3 (Chapter P-6a, https://iupac.qmul.ac.uk/BlueBook/PDF/P6a.pdf):
  'imine' is the suffix for a C=N-H or C=N-R group, generically called
  'aldimine' (R-CH=N-R') when the imine carbon has one carbon neighbor, or
  'ketimine' (R(R')C=N-R'') when it has two. Confirmed via PubChem PUG
  REST: CID 123139 (`C=N`) -> "methanimine", CID 140746 (`CC=N`) ->
  "ethanimine", CID 142304 (`CC(C)=N`) -> "propan-2-imine", CID 12355452
  (`CCC(C)=N`) -> "butan-2-imine".
- Unlike a ketone carbonyl (`_ketone.py`), an aldimine carbon (one carbon
  neighbor) is a perfectly ordinary case for this suffix -- there is no
  senior aldehyde-shaped group to defer to, since '-imine' has no such
  competing aldehyde-like sibling. P-62.3.1.1's own worked example cites
  the locant even at chain position 1 on a 3+-carbon chain: CID 5289472
  (`CCC=N`) -> "propan-1-imine", CID 14112901 (`CCCCCC=N`) ->
  "hexan-1-imine" (matching the Blue Book's own "hexan-1-imine (PIN)"
  worked example directly). This is a real difference from
  '-al'/'-onitrile'/etc., whose sole-possible-position locant is always
  omitted regardless of chain length.
- P-14.3.4.2(a)/(b): a mononuclear (one-carbon) chain never cites a
  locant (CID 123139, "methanimine", not "methan-1-imine"), and a
  homogeneous two-carbon chain omits the locant too, regardless of other
  substituents present: CID 54110962 (`ClCC=N`) -> "2-chloroethanimine",
  not "2-chloroethan-1-imine". Same rule the shared `_common.py`
  `name_from_substituents` hook now implements project-wide.
- P-62.3.1.1: an N-substituent (imine nitrogen bonded to one additional
  carbon group instead of H) is cited as an "N-" prefix with no locant
  (nitrogen is never part of the numbered chain), same style already used
  by `_imide.py`. Confirmed via PubChem PUG REST: CID 144069 (`CC=NC`) ->
  "N-methylethanimine", CID 138743 (`CC(C)=NC`) ->
  "N-methylpropan-2-imine", CID 137199 (`C=NC`) -> "N-methylmethanimine".
  A nitrogen's remaining valence after the C=N double bond is only one
  more bond, so at most one N-substituent is ever possible (never
  N,N-disubstituted, unlike an amine nitrogen).
- P-35.2.1: halogen substituents on the carbon chain are prefix-only and
  coexist freely with the imine suffix, reusing
  `halogen_substituents`/`format_substituent_prefixes` unchanged (same as
  `_ketone.py`).
- P-68.3.1.1.2 (oximes): the PIN for R2C=N-OH is itself defined as the
  N-hydroxy derivative of the imine named by this module -- confirmed by
  the worked example 'N-hydroxypentan-2-imine (PIN)' for pentan-2-one
  oxime. This module therefore also accepts the imine nitrogen's one
  remaining substituent being an oxygen instead of a carbon, formatted the
  same "N-" way as an N-alkyl substituent (P-62.3.1.1): plain -OH gives
  'N-hydroxy...', and an O-alkyl ether oxime (=N-O-R, an unbranched alkyl
  R) gives 'N-<alkoxy>...' (confirmed via PubChem PUG REST: CID 54150571,
  `CCC=NOCC` -> "N-ethoxypropan-1-imine", matching the Blue Book's own
  'N-ethoxypropan-1-imine (PIN)' worked example directly). The plain -OH
  case has a PubChem-autoname-vs-PIN mismatch worth flagging: PubChem's
  own auto-generated name for e.g. `CC(=NO)CCC` (CID 136433) is
  "N-pentan-2-ylidenehydroxylamine" (a hydroxylamine-parent 'ylidene'
  prefix pattern) rather than "N-hydroxypentan-2-imine" -- implemented per
  the Blue Book's own direct PIN citation instead, consistent with this
  project's established practice elsewhere (see e.g. `_hydroxylamine.py`,
  `_sulfoxide.py`, `_nitro.py`).
- P-91.2(e)/P-93.1 (E/Z stereo): when this module's own C=N bond has its
  geometry specified in the input (`/`/`\\`), a locanted `"(nE)-"`/`"(nZ)-"`
  prefix is added to the whole name (including any `N-` substituent
  prefix), the same way `_unsaturated.py` handles an acyclic C=C bond, via
  `_common.specified_double_bond_stereo`. A prior revision of this
  docstring claimed RDKit's `Chem.FindPotentialStereo` *always* flags this
  bond as an unspecified potential stereo element regardless of input --
  that claim no longer holds against the currently pinned RDKit
  (confirmed: `C/C=N/C` now reports it as specified) and was corrected
  here; a specified tetrahedral chain stereocenter coexisting with the
  imine bond (specified or not) is still rejected, since combining the two
  kinds of descriptor is out of scope for this module (P-92/P-93).
  `_hydrazone.py`/`_azine.py`/`_amidine.py` share the identical
  (now-stale) rejection pattern independently -- not updated by this
  change.

Explicitly out of scope (raise `UnsupportedStructure`):
- Any ring anywhere in the molecule (cyclic/aromatic imines, e.g.
  'thiolan-2-imine', are a separate, unverified case here).
- More than one C=N imine bond (polyimines, P-62.3.1.1's multiplying-
  prefix case) -- unverified in this first pass.
- Any heteroatom other than the imine nitrogen, halogen substituents, and
  (for an oxime) its own N-O(-R) oxygen (no coexisting -OH elsewhere,
  other C=N/C=O, etc. -- a suffix-seniority competition this module
  doesn't attempt).
- A branched imine nitrogen substituent (carbon or O-alkyl), an aromatic
  N-substituent, or any bond order other than single/double, or a second
  double/triple bond elsewhere on the chain (P-31 unsaturation combined
  with imine is unverified in this first pass).
- Nitrolic/nitrosolic acids (P-68.3.1.1.3, an oxime combined with an
  adjacent nitro/nitroso group) -- a separate, unverified follow-up.
- Charged or isotopically modified atoms.
"""

from rdkit import Chem

from ._common import (
    HALOGEN_PREFIXES,
    UnsupportedStructure,
    adjacency,
    carbon_adjacency,
    group_substituents,
    halogen_substituents,
    is_plain_benzene_ring,
    longest_branched_chain_through,
    longest_chains,
    name_from_substituents,
    non_single_bonds,
    ring_chain_attachment,
    ring_branch_attachments,
    separate_aromatic_monocycles,
    specified_double_bond_stereo,
    substituent_locant_set_and_citation,
)
from ._substituents import (
    substituents_for_chain,
    format_substituent_prefixes,
    name_branch,
)

_ALLOWED_ATOMIC_NUMS = {6, 7, *HALOGEN_PREFIXES}
_OXIME_ALLOWED_ATOMIC_NUMS = {6, 7, 8, *HALOGEN_PREFIXES}
_OXY_PREFIX = {"methyl": "methoxy", "ethyl": "ethoxy", "propyl": "propoxy", "butyl": "butoxy"}


def has_simple_imine_shape(mol) -> bool:
    return any(
        bond.GetBondTypeAsDouble() == 2.0
        and {bond.GetBeginAtom().GetAtomicNum(), bond.GetEndAtom().GetAtomicNum()} == {6, 7}
        for bond in mol.GetBonds()
    )


def _validate_and_find_imine(mol, aromatic_ring_atoms=frozenset()):
    """Check the molecule fits this module's scope (see module docstring)
    and return (imine_carbon_idx, imine_nitrogen_idx, n_substituent_root,
    oxime_oxygen_idx, oxime_alkyl_root). The last two are only set for an
    oxime/O-alkyl oxime ether (oxime_alkyl_root stays None for a plain
    -OH); n_substituent_root and the oxime pair are mutually exclusive,
    since the imine nitrogen has only one substituent position free.

    `aromatic_ring_atoms`: atom indices already independently verified (by
    the caller, before this function runs) to form a single plain benzene
    ring with exactly one exocyclic attachment -- exempted from the
    aromatic-atom rejection below so `name_imine`'s benzene-ring-
    substituent path (see `_name_phenyl_chain_imine`) can reuse this same
    validation for the rest of the molecule. Empty by default, so every
    other caller's behavior is unchanged. Mirrors `_sulfonic_acid.py`'s
    equivalent aromatic-exemption pattern."""
    imine_bonds = [
        bond
        for bond in mol.GetBonds()
        if bond.GetBondTypeAsDouble() == 2.0
        and {bond.GetBeginAtom().GetAtomicNum(), bond.GetEndAtom().GetAtomicNum()} == {6, 7}
    ]
    if not imine_bonds:
        raise UnsupportedStructure("no imine (C=N) group found; this module only handles imines")
    if len(imine_bonds) > 1:
        raise UnsupportedStructure("more than one C=N imine bond (polyimines) is not supported yet")
    (imine_bond,) = imine_bonds
    carbon = imine_bond.GetBeginAtom() if imine_bond.GetBeginAtom().GetAtomicNum() == 6 else imine_bond.GetEndAtom()
    nitrogen = imine_bond.GetBeginAtom() if imine_bond.GetBeginAtom().GetAtomicNum() == 7 else imine_bond.GetEndAtom()

    n_substituent_root = None
    oxime_oxygen_idx = None
    if nitrogen.GetDegree() == 2:
        (other,) = [n for n in nitrogen.GetNeighbors() if n.GetIdx() != carbon.GetIdx()]
        if other.GetAtomicNum() == 8:
            oxime_oxygen_idx = other.GetIdx()
        elif other.GetAtomicNum() == 6 and not other.GetIsAromatic():
            n_substituent_root = other.GetIdx()
        else:
            raise UnsupportedStructure(
                "an imine nitrogen substituent other than a plain carbon "
                "group or an oxime oxygen (P-68.3.1.1.2) is not supported "
                "yet"
            )
    elif nitrogen.GetDegree() != 1:
        raise UnsupportedStructure(
            "an imine nitrogen must have exactly one substituent (or none, "
            "i.e. =N-H)"
        )

    allowed = _OXIME_ALLOWED_ATOMIC_NUMS if oxime_oxygen_idx is not None else _ALLOWED_ATOMIC_NUMS
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atomic_num not in allowed:
            raise UnsupportedStructure(
                "heteroatoms other than the imine nitrogen (P-62.3), a "
                "single oxime oxygen (P-68.3.1.1.2), and halogen "
                "substituents (P-35.2.1) are not supported yet"
            )
        if atomic_num == 8 and atom.GetIdx() != oxime_oxygen_idx:
            raise UnsupportedStructure(
                "an oxygen atom not shaped like a plain oxime N-OH/N-O-R "
                "is out of scope for this module"
            )
        if atom.GetFormalCharge() != 0 or atom.GetIsotope() != 0:
            raise UnsupportedStructure("charged or isotopically modified atoms are not supported yet")
        if atomic_num == 6 and atom.GetIsAromatic() and atom.GetIdx() not in aromatic_ring_atoms:
            raise UnsupportedStructure("aromatic rings are out of scope for this module")
        if atomic_num == 7 and atom.GetIsAromatic():
            raise UnsupportedStructure("an aromatic nitrogen is out of scope for this module")
    if not aromatic_ring_atoms and mol.GetRingInfo().NumRings() != 0:
        raise UnsupportedStructure("a ring anywhere in the molecule is out of scope for this module")
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")

    other_non_single = [
        bond
        for bond in mol.GetBonds()
        if bond.GetBondTypeAsDouble() != 1.0
        and bond.GetIdx() != imine_bond.GetIdx()
        and not (
            bond.GetBeginAtomIdx() in aromatic_ring_atoms and bond.GetEndAtomIdx() in aromatic_ring_atoms
        )
    ]
    if other_non_single:
        raise UnsupportedStructure(
            "a bond order other than single (besides the imine C=N itself) "
            "is not supported yet"
        )

    carbon_neighbors = [n for n in carbon.GetNeighbors() if n.GetAtomicNum() == 6]
    if len(carbon_neighbors) not in (0, 1, 2):
        raise UnsupportedStructure(
            "an imine carbon must have zero (mononuclear, e.g. "
            "methanimine), one (aldimine), or two (ketimine) carbon "
            "neighbors"
        )

    oxime_alkyl_root = None
    if oxime_oxygen_idx is not None:
        oxygen = mol.GetAtomWithIdx(oxime_oxygen_idx)
        if oxygen.GetDegree() == 2:
            (alkyl,) = [n for n in oxygen.GetNeighbors() if n.GetIdx() != nitrogen.GetIdx()]
            if alkyl.GetAtomicNum() != 6 or alkyl.GetIsAromatic():
                raise UnsupportedStructure(
                    "an oxime O-substituent other than a plain carbon "
                    "group is not supported yet"
                )
            oxime_alkyl_root = alkyl.GetIdx()
        elif oxygen.GetDegree() != 1:
            raise UnsupportedStructure("an oxime oxygen must be -OH or -O-R (degree 1 or 2)")

    return carbon.GetIdx(), nitrogen.GetIdx(), n_substituent_root, oxime_oxygen_idx, oxime_alkyl_root


def _imine_locant(position_of, imine_carbon):
    return position_of.get(imine_carbon)


def _n_substituent_prefix(graph, imine_nitrogen, n_substituent_root, oxime_oxygen_idx, oxime_alkyl_root, mol=None):
    """Build the 'N-...' prefix (without the leading 'N-' itself, e.g.
    'methyl'/'hydroxy'/'ethoxy') for an N-substituted imine or an oxime,
    or None for a plain =N-H imine. Shared by the acyclic and
    benzene-ring-substituent-chain naming paths so the two stay in sync."""
    if n_substituent_root is not None:
        name, _ = name_branch(graph, n_substituent_root, imine_nitrogen, mol=mol)
        return name
    if oxime_oxygen_idx is not None:
        if oxime_alkyl_root is None:
            return "hydroxy"
        alkyl_name_, is_compound = name_branch(graph, oxime_alkyl_root, oxime_oxygen_idx, mol=mol)
        if is_compound:
            raise UnsupportedStructure("a branched oxime O-substituent is not supported yet")
        return _OXY_PREFIX.get(alkyl_name_, alkyl_name_ + "oxy")
    return None


def _name_from_substituents(chain_length, imine_locant, grouped):
    # P-14.3.4.2(a): a mononuclear parent's substituent locants (not just
    # the imine's own) are always '1' and never cited either -- unreachable
    # before the benzene-ring-substituent chain-length-1 case (a bare
    # methanimine carbon has no room for any other substituent), so this
    # branch was never previously exercised.
    prefix = format_substituent_prefixes(grouped, omit_locants=(chain_length == 1))
    return prefix + name_from_substituents(chain_length, [], [], "imine", [imine_locant])


def _candidate_key(chain_length, imine_locant, substituents):
    grouped = group_substituents(substituents)
    locant_set, _, citation_locants = substituent_locant_set_and_citation(grouped)
    name = _name_from_substituents(chain_length, imine_locant, grouped)
    return (imine_locant, locant_set, citation_locants, name), name


def _name_acyclic_imine_with_locant(mol, imine_carbon, exclude):
    """Like `_name_acyclic_imine`, but also returns the winning candidate's
    imine locant (needed to place a `(nE)-`/`(nZ)-` stereodescriptor)."""
    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    chains = longest_chains(carbon_adjacency(mol))
    chain_length = len(chains[0])

    eligible = [chain for chain in chains if imine_carbon in chain]
    if not eligible:
        raise UnsupportedStructure(
            "the imine carbon does not lie on any longest carbon chain; a "
            "shorter principal chain (P-44.1.1) is not supported yet"
        )

    best_key = None
    best_name = None
    for chain in eligible:
        for candidate in (chain, list(reversed(chain))):
            position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
            imine_locant = _imine_locant(position_of, imine_carbon)
            substituents = substituents_for_chain(graph, candidate, halogens, exclude, mol=mol)
            key, name = _candidate_key(chain_length, imine_locant, substituents)
            if best_key is None or key < best_key:
                best_key, best_name = key, name
    return best_name, best_key[0]


def _name_acyclic_imine(mol, imine_carbon, exclude):
    name, _ = _name_acyclic_imine_with_locant(mol, imine_carbon, exclude)
    return name


def _name_phenyl_chain_imine(mol, ring_atoms):
    """Name an imine (plain =N-H, N-alkyl-substituted, or a plain -OH
    oxime) whose C=N lies entirely on a single unbranched chain hanging
    off one atom of an otherwise-plain, unsubstituted benzene ring -- e.g.
    3-phenylpropan-1-imine, N-methyl-3-phenylpropan-1-imine,
    N-hydroxy-3-phenylpropan-1-imine. The ring is cited as a 'phenyl'
    substituent prefix (via `name_branch`'s aromatic-ring recognition) on
    the chain, which is the parent hydride, mirroring `_sulfonic_acid.py`'s
    `_name_phenyl_chain_sulfonic_acid`. Unlike a carboxylic-acid-style
    chain terminus, the imine carbon's own locant is a genuine choice
    (P-62.3.1.1), same as the sulfonic acid group. An O-alkyl oxime ether
    (=N-O-R) is out of scope here -- no PubChem-registered example exists
    to verify the combination against (module docstring's own oxime
    section notes this project already diverges from PubChem's own
    auto-generated oxime name in the non-benzene case, so an unverifiable
    third variant isn't a safe mechanical extension)."""
    imine_carbon, imine_nitrogen, n_substituent_root, oxime_oxygen_idx, oxime_alkyl_root = (
        _validate_and_find_imine(mol, aromatic_ring_atoms=ring_atoms)
    )
    if oxime_alkyl_root is not None:
        raise UnsupportedStructure(
            "an oxime O-alkyl ether alongside a benzene-ring substituent "
            "chain is not supported yet"
        )
    stereo = specified_double_bond_stereo(mol)

    exclude = {imine_nitrogen}
    non_ring_unsaturation = [
        b
        for b in non_single_bonds(mol)
        if b[0] not in exclude and b[1] not in exclude and b[0] not in ring_atoms and b[1] not in ring_atoms
    ]
    if non_ring_unsaturation:
        raise UnsupportedStructure(
            "chain unsaturation alongside a benzene-ring-substituent "
            "imine chain is not supported yet"
        )

    graph = adjacency(mol)
    halogens = halogen_substituents(mol)
    rings = separate_aromatic_monocycles(mol, graph) or [set(ring_atoms)]
    attachment = ring_branch_attachments(mol, graph, rings)
    if not attachment:
        raise UnsupportedStructure(
            "a benzene ring with more than one non-halogen, non-alkyl "
            "exocyclic substituent alongside a chain imine is not "
            "supported yet"
        )
    chain, branches = longest_branched_chain_through(graph, imine_carbon, ring_atoms, exclude, halogens=halogen_substituents(mol))
    branches_by_atom = {chain[position - 1]: roots for position, roots in branches.items()}

    chain_length = len(chain)
    best_key = None
    best_name = None
    for candidate in (chain, list(reversed(chain))):
        position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
        imine_locant = _imine_locant(position_of, imine_carbon)
        substituents = {
            position_of[atom]: [name_branch(graph, root, atom, halogens, ring_atoms, mol=mol) for root in roots]
            for atom, roots in branches_by_atom.items()
        }
        key, name = _candidate_key(chain_length, imine_locant, substituents)
        if best_key is None or key < best_key:
            best_key, best_name = key, name

    n_name = _n_substituent_prefix(graph, imine_nitrogen, n_substituent_root, oxime_oxygen_idx, oxime_alkyl_root, mol=mol)
    if n_name is not None:
        separator = "-" if best_name[0].isdigit() else ""
        best_name = f"N-{n_name}{separator}{best_name}"
    if stereo is not None:
        ((_, code),) = stereo
        best_name = f"({best_key[0]}{code})-{best_name}"
    return best_name


def name_imine(mol) -> str:
    aromatic_rings = separate_aromatic_monocycles(mol, adjacency(mol))
    if aromatic_rings is not None:
        return _name_phenyl_chain_imine(mol, set().union(*aromatic_rings))
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() == 1:
        ring_atoms = set(ring_info.AtomRings()[0])
        if is_plain_benzene_ring(mol, ring_atoms):
            return _name_phenyl_chain_imine(mol, ring_atoms)

    imine_carbon, imine_nitrogen, n_substituent_root, oxime_oxygen_idx, oxime_alkyl_root = _validate_and_find_imine(
        mol
    )
    stereo = specified_double_bond_stereo(mol)
    graph = adjacency(mol)
    exclude = {imine_nitrogen}
    name, imine_locant = _name_acyclic_imine_with_locant(mol, imine_carbon, exclude)

    n_name = _n_substituent_prefix(graph, imine_nitrogen, n_substituent_root, oxime_oxygen_idx, oxime_alkyl_root, mol=mol)

    if n_name is not None:
        separator = "-" if name[0].isdigit() else ""
        name = f"N-{n_name}{separator}{name}"
    if stereo is not None:
        # `_validate_and_find_imine` already rejects a second C=N/other
        # non-single bond anywhere else in the molecule, so this module's
        # own imine bond is the only stereo element `stereo` can ever
        # contain -- no locant-sorting/grouping needed (P-91.3).
        ((_, code),) = stereo
        name = f"({imine_locant}{code})-{name}"
    return name
