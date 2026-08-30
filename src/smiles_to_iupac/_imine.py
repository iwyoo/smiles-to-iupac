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
  homogeneous two-carbon chain omits the locant too -- but, unlike this
  project's own `_alcohol.py`/`_ketone.py`, that omission is *not* gated
  on having zero other substituents: CID 54110962 (`ClCC=N`) ->
  "2-chloroethanimine", not "2-chloroethan-1-imine". This module
  therefore omits the locant for any two-carbon aldimine chain outright
  (this looks like a more complete implementation of P-14.3.4.2(b) than
  `_alcohol.py`/`_ketone.py`'s own "total substituents == 0" gate, which
  is a known, already-documented limitation there -- see
  `_alcohol.py`'s "2-cyclohexylethan-1-ol" test comment -- but fixing
  those other modules is out of scope here).
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
- P-92 stereocenters (`tasks/imine-azine-stereocenter-naming.md`): this
  module's own C=N bond is *always* flagged by RDKit's
  `Chem.FindPotentialStereo` as an unspecified potential Bond_Double
  stereo element, regardless of substituents, N-substitution, or oxime
  form -- same conclusion as `_amidine.py`/`_hydrazone.py`. Any specified
  chain tetrahedral stereocenter therefore always coexists with this
  unspecified C=N bond, and `_common.specified_stereocenters` correctly
  rejects the combination as partially specified (P-92/P-93) rather than
  silently dropping either one. This module only ever explicitly rejects
  a specified stereocenter rather than attempting to cite one; `_azine.py`
  (which reuses this module's internal chain-assembly helper directly)
  carries the identical restriction independently.

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
    longest_chains,
    lowest_locant_set,
    specified_stereocenters,
)
from ._numerals import alkane_name
from ._substituents import alpha_sort_key, format_substituent_prefixes, name_branch

_ALLOWED_ATOMIC_NUMS = {6, 7, *HALOGEN_PREFIXES}
_OXIME_ALLOWED_ATOMIC_NUMS = {6, 7, 8, *HALOGEN_PREFIXES}
_OXY_PREFIX = {"methyl": "methoxy", "ethyl": "ethoxy", "propyl": "propoxy", "butyl": "butoxy"}


def has_simple_imine_shape(mol) -> bool:
    return any(
        bond.GetBondTypeAsDouble() == 2.0
        and {bond.GetBeginAtom().GetAtomicNum(), bond.GetEndAtom().GetAtomicNum()} == {6, 7}
        for bond in mol.GetBonds()
    )


def _validate_and_find_imine(mol):
    """Check the molecule fits this module's scope (see module docstring)
    and return (imine_carbon_idx, imine_nitrogen_idx, n_substituent_root,
    oxime_oxygen_idx, oxime_alkyl_root). The last two are only set for an
    oxime/O-alkyl oxime ether (oxime_alkyl_root stays None for a plain
    -OH); n_substituent_root and the oxime pair are mutually exclusive,
    since the imine nitrogen has only one substituent position free."""
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
        if atomic_num == 6 and atom.GetIsAromatic():
            raise UnsupportedStructure("aromatic rings are out of scope for this module")
        if atomic_num == 7 and atom.GetIsAromatic():
            raise UnsupportedStructure("an aromatic nitrogen is out of scope for this module")
    if mol.GetRingInfo().NumRings() != 0:
        raise UnsupportedStructure("a ring anywhere in the molecule is out of scope for this module")
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")

    other_non_single = [
        bond
        for bond in mol.GetBonds()
        if bond.GetBondTypeAsDouble() != 1.0 and bond.GetIdx() != imine_bond.GetIdx()
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


def _substituents_for_chain(graph, chain, halogens, exclude):
    chain_set = set(chain)
    substituents = {}
    for position, atom in enumerate(chain, start=1):
        branch_roots = [n for n in graph[atom] if n not in chain_set and n not in exclude]
        if branch_roots:
            substituents[position] = [name_branch(graph, root, atom, halogens) for root in branch_roots]
    return substituents


def _name_from_substituents(chain_length, imine_locant, grouped):
    prefix = format_substituent_prefixes(grouped)
    # 'imine' always starts with a vowel, so the alkane stem's trailing 'e'
    # is always elided (P-16.3.3 / P-16.6), whether or not a locant lands
    # between them -- 'methanimine', 'propan-2-imine', not 'methaneimine'/
    # 'propane-2-imine'.
    stem = alkane_name(chain_length)[:-1]

    if chain_length in (1, 2):
        # P-14.3.4.2(a): a mononuclear chain never cites a locant. (b): a
        # two-carbon aldimine chain has only one possible imine position
        # (there is no room for a second carbon neighbor, so it's always
        # an aldimine) -- this project's own `_alcohol.py`/`_ketone.py`
        # additionally require zero other substituents for this omission,
        # but the real rule doesn't gate on that (see module docstring's
        # "2-chloroethanimine" citation), so it's omitted here regardless
        # of any halogen prefix already present.
        return prefix + stem + "imine"

    return f"{prefix}{stem}-{imine_locant}-imine"


def _candidate_key(chain_length, imine_locant, substituents):
    grouped = group_substituents(substituents)
    locant_set = lowest_locant_set(loc for info in grouped.values() for loc in info["locants"])
    citation_locants = tuple(
        loc
        for name in sorted(grouped, key=alpha_sort_key)
        for loc in sorted(grouped[name]["locants"])
    )
    name = _name_from_substituents(chain_length, imine_locant, grouped)
    return (imine_locant, locant_set, citation_locants, name), name


def _name_acyclic_imine(mol, imine_carbon, exclude):
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
            substituents = _substituents_for_chain(graph, candidate, halogens, exclude)
            key, name = _candidate_key(chain_length, imine_locant, substituents)
            if best_key is None or key < best_key:
                best_key, best_name = key, name
    return best_name


def name_imine(mol) -> str:
    imine_carbon, imine_nitrogen, n_substituent_root, oxime_oxygen_idx, oxime_alkyl_root = _validate_and_find_imine(
        mol
    )
    if specified_stereocenters(mol) is not None:
        # This module's own C=N bond is always flagged by RDKit's
        # `Chem.FindPotentialStereo` as an unspecified potential
        # Bond_Double stereo element, regardless of substituents,
        # N-substitution, or oxime form (module docstring) -- so any
        # specified chain stereocenter always coexists with it, and
        # `specified_stereocenters` correctly rejects the combination
        # (P-92/P-93) instead of the silent drop this project's
        # stereodescriptor safety net exists to fix, same conclusion as
        # `_amidine.py` (PR #209) / `_hydrazone.py` (PR #210).
        raise UnsupportedStructure(
            "a specified stereocenter alongside this module's own "
            "always-unspecified C=N bond is not supported yet (see "
            "P-92/P-93, module docstring)"
        )
    graph = adjacency(mol)
    exclude = {imine_nitrogen}
    name = _name_acyclic_imine(mol, imine_carbon, exclude)

    n_name = None
    if n_substituent_root is not None:
        n_name, _ = name_branch(graph, n_substituent_root, imine_nitrogen)
    elif oxime_oxygen_idx is not None:
        if oxime_alkyl_root is None:
            n_name = "hydroxy"
        else:
            alkyl_name_, is_compound = name_branch(graph, oxime_alkyl_root, oxime_oxygen_idx)
            if is_compound:
                raise UnsupportedStructure(
                    "a branched oxime O-substituent is not supported yet"
                )
            n_name = _OXY_PREFIX.get(alkyl_name_, alkyl_name_ + "oxy")

    if n_name is not None:
        separator = "-" if name[0].isdigit() else ""
        name = f"N-{n_name}{separator}{name}"
    return name
