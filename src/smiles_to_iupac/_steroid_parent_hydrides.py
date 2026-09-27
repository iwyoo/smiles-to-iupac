"""Naming of the seven 1989 IUPAC steroid parent ring hydrides (Rule 2.1
and Rule 3S-2.2/2.3/2.4 Table 1, https://iupac.qmul.ac.uk/steroid/3S02a.html,
the document Blue Book P-6 defers to for steroid parent names), each by
hardcoded retained-name recognition of its exact unsubstituted structure.

Each skeleton in this lineage is the previous one plus one defined
extension, built on the plain perhydrocyclopenta[a]phenanthrene ring system
(three fused six-membered rings plus one five-membered ring, angularly
ortho-fused):

- gonane (Rule 2.1): the bare ring system, no angular methyls, no C17 side
  chain. PubChem CID 6857523.
- androstane (Rule 3S-2.3): gonane + angular methyls at C10/C13. CID 6857536.
- estrane (Rule 3S-2.2): gonane + only the C13 angular methyl (no C10
  methyl). CID 5460658.
- pregnane (Rule 3S-2.4, Table 1): androstane + a plain two-carbon ethyl
  side chain at C17 (C20-C21). CID 439513.
- cholane: pregnane's side chain extended to five carbons (C20-C24) with a
  methyl branch at C20. CID 6857459.
- cholestane: cholane's side chain extended by a symmetric gem-dimethyl
  terminus (C25-C27). CID 6857534.
- ergostane: cholestane's side chain extended by one more methyl branch at
  C24 (C28). CID 6857535.
- campestane: ergostane's C24 epimer -- identical constitution (a
  5,6-dimethylheptan-2-yl side chain), opposite configuration at that
  carbon. CID 6857532.
- poriferastane: a 5-ethyl-6-methylheptan-2-yl side chain (ergostane's
  extra C24 methyl replaced by an ethyl group). CID 6857528.
- stigmastane: poriferastane's C24 epimer -- identical constitution,
  opposite configuration at that carbon. CID 6857438.

Each entry below was independently confirmed by removing the next
skeleton's side-chain/methyl carbons from its canonical SMILES and
re-canonicalizing, reproducing the previous skeleton's own canonical key
exactly (e.g. removing androstane's two angular methyls reproduces gonane's
key; removing ergostane's nine-carbon C17 side chain reproduces
androstane's key).

Since each module's only job is recognizing one exact unsubstituted parent,
a whole-molecule canonical-SMILES match is sufficient: any substituent,
extra methyl, longer side chain, or unsaturation changes the canonical
SMILES and so is correctly left unmatched, falling through to
`_polycyclic.py`'s general von Baeyer engine (which already names the bare
gonane skeleton as "tetracyclo[8.7.0.0^2,7.0^11,15]heptadecane" if nothing
here intercepts it first -- confirmed by direct testing).

`Chem.MolToSmiles` includes stereo markers (@/@@) when present on the input
mol, so a stereo-specified ring-fusion or side-chain input only matches if
it's the exact natural configuration listed below (Rule 3S-2.2/2.3/2.4's
own worked structures, PubChem's own isomeric SMILES for each CID above)
-- any other stereoisomer (a ring-fusion epimer, a partially-specified
input, or a side-chain epimer at C20/C24 not matching the natural series)
still falls through to the general engine unmatched, since it produces a
different canonical SMILES string. None of Rule 2.1/3S-2.2/3S-2.3 assign
ring-fusion stereochemistry beyond this one natural configuration per
skeleton -- a differently configured stereoisomer (e.g. 5-beta) needs its
own descriptor and lookup row, out of scope here (tracked separately).

Each skeleton's own PubChem CID above sometimes leaves the C5 ring-fusion
stereocenter unspecified in its own isomeric SMILES record (a PubChem
data-depiction quirk, not a chemistry fact -- C5 is a genuine stereocenter
in every skeleton here) -- confirmed for gonane/androstane/estrane/
cholestane/ergostane/campestane/poriferastane/stigmastane (only pregnane
and cholane's own CID records happen to already be fully specified). Each
of those 8 skeletons therefore also has a second natural-configuration
entry below, from that skeleton's own dedicated "5alpha-<name>" PubChem
name/CID, which does fully specify C5 -- additive alongside the
already-registered partially-specified one, so a real input matching
either still resolves correctly.
"""

from rdkit import Chem
from rdkit.Chem import BondType, RWMol, rdCIPLabeler

from ._parent_hydride_stripping import strip_substituents
from ._polycyclic_suffix import _suffixed_parent

_PARENT_HYDRIDES = {
    "C1CCCC2CCC3C(C12)CCC4C3CCC4": "gonane",
    "CC12CCCC1C3CCC4CCCCC4(C3CC2)C": "androstane",
    "CC12CCCC1C1CCC3CCCCC3C1CC2": "estrane",
    "CCC1CCC2C1(CCC3C2CCC4C3(CCCC4)C)C": "pregnane",
    "CCCC(C)C1CCC2C1(CCC3C2CCC4C3(CCCC4)C)C": "cholane",
    "CC(C)CCCC(C)C1CCC2C1(CCC3C2CCC4C3(CCCC4)C)C": "cholestane",
    "CC(C)C(C)CCC(C)C1CCC2C1(CCC3C2CCC4C3(CCCC4)C)C": "ergostane",
    # Natural-configuration isomeric SMILES, one per skeleton above (same
    # PubChem CIDs), additive alongside the stereo-free entries -- these
    # are the shape essentially every real-world instance of these
    # compounds actually has.
    "C1CCC2CC[C@H]3[C@@H]4CCC[C@H]4CC[C@@H]3[C@H]2C1": "gonane",
    "C[C@@]12CCC[C@H]1[C@@H]3CCC4CCCC[C@@]4([C@H]3CC2)C": "androstane",
    "C[C@@]12CCC[C@H]1[C@@H]3CCC4CCCC[C@@H]4[C@H]3CC2": "estrane",
    "CC[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC[C@H]4[C@@]3(CCCC4)C)C": "pregnane",
    "CCC[C@@H](C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC[C@H]4[C@@]3(CCCC4)C)C": "cholane",
    "C[C@H](CCCC(C)C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CCC4[C@@]3(CCCC4)C)C": "cholestane",
    "C[C@H](CC[C@H](C)C(C)C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CCC4[C@@]3(CCCC4)C)C": "ergostane",
    "C[C@H](CC[C@@H](C)C(C)C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CCC4[C@@]3(CCCC4)C)C": "campestane",
    "CC[C@@H](CC[C@@H](C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CCC4[C@@]3(CCCC4)C)C)C(C)C": "poriferastane",
    "CC[C@H](CC[C@@H](C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CCC4[C@@]3(CCCC4)C)C)C(C)C": "stigmastane",
    # Fully C5-specified natural-configuration entries (see module
    # docstring), from each skeleton's own "5alpha-<name>" PubChem CID.
    "C1CC[C@H]2[C@H](C1)CC[C@@H]3[C@@H]2CC[C@H]4[C@H]3CCC4": "gonane",
    "C[C@@]12CCC[C@H]1[C@@H]3CC[C@H]4CCCC[C@@]4([C@H]3CC2)C": "androstane",
    "C[C@@]12CCC[C@H]1[C@@H]3CC[C@H]4CCCC[C@@H]4[C@H]3CC2": "estrane",
    "C[C@H](CCCC(C)C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC[C@@H]4[C@@]3(CCCC4)C)C": "cholestane",
    "C[C@H](CC[C@H](C)C(C)C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC[C@@H]4[C@@]3(CCCC4)C)C": "ergostane",
    "C[C@H](CC[C@@H](C)C(C)C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC[C@@H]4[C@@]3(CCCC4)C)C": "campestane",
    "CC[C@@H](CC[C@@H](C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC[C@@H]4[C@@]3(CCCC4)C)C)C(C)C": "poriferastane",
    "CC[C@H](CC[C@@H](C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC[C@@H]4[C@@]3(CCCC4)C)C)C(C)C": "stigmastane",
}
_CANONICAL_TO_NAME = {Chem.CanonSmiles(smiles): name for smiles, name in _PARENT_HYDRIDES.items()}
assert len(_CANONICAL_TO_NAME) == len(_PARENT_HYDRIDES), (
    "two entries above canonicalized to the same key -- a real name "
    "collision, not just a duplicate row (see module docstring)"
)

# The plain (no-stereo) entries above, keyed by name instead of SMILES --
# used by both `_nor_steroid.py`'s skeletal-modification prefixes and
# `_ketone.py`'s steroid-suffix integration (#1027) to recognize a
# constitution match irrespective of stereochemistry: a suffix-bearing
# steroid's own stereo, if any, is already rejected upstream (P-92
# descriptor citation on a retained steroid name is separate, unimplemented
# future work), so only the ring skeleton's constitution -- not its
# configuration -- is ever checked here.
_RAW_PARENTS = {name: smiles for smiles, name in _PARENT_HYDRIDES.items() if "@" not in smiles}
_PLAIN_CANONICAL_TO_NAME = {Chem.CanonSmiles(smiles): name for name, smiles in _RAW_PARENTS.items()}


def _build_locant_query(extra_locants):
    """A hand-built query graph encoding the standard steroid numbering
    (Table 10.1: ring positions 1-17, angular methyls 18/19), used to map
    a recognized parent's own atom indices onto their steroid locants via
    substructure match -- not the raw SMILES entries' own atom order,
    which isn't guaranteed to follow locant order. `extra_locants`: which
    of the two angular methyls (18 at C13, 19 at C10) the target parent
    also has (gonane has neither, estrane has only 18, every androstane-
    family parent has both)."""
    rw = RWMol()
    atom_idx = {}
    for locant in range(1, 18):
        atom_idx[locant] = rw.AddAtom(Chem.Atom(6))
    for locant in extra_locants:
        atom_idx[locant] = rw.AddAtom(Chem.Atom(6))

    def bond(a, b):
        rw.AddBond(atom_idx[a], atom_idx[b], BondType.SINGLE)

    bond(1, 2), bond(2, 3), bond(3, 4), bond(4, 5), bond(5, 6), bond(6, 7)
    bond(7, 8), bond(8, 9), bond(9, 10), bond(10, 1), bond(5, 10)
    bond(9, 11), bond(11, 12), bond(12, 13), bond(13, 14), bond(14, 8)
    bond(13, 17), bond(17, 16), bond(16, 15), bond(15, 14)
    if 18 in extra_locants:
        bond(13, 18)
    if 19 in extra_locants:
        bond(10, 19)

    mol = rw.GetMol()
    Chem.SanitizeMol(mol)
    return mol, atom_idx


_RING_QUERY, _RING_QUERY_IDX = _build_locant_query(())
_ESTRANE_QUERY, _ESTRANE_QUERY_IDX = _build_locant_query((18,))
_ANDROSTANE_QUERY, _ANDROSTANE_QUERY_IDX = _build_locant_query((18, 19))

_QUERY_BY_PARENT = {
    "gonane": (_RING_QUERY, _RING_QUERY_IDX),
    "estrane": (_ESTRANE_QUERY, _ESTRANE_QUERY_IDX),
}


def _locant_map(name, mol):
    """`{locant: atom_idx}` for `mol`'s own occurrence of the named
    parent's ring skeleton (`mol` must actually contain it, e.g. already
    confirmed via `_CANONICAL_TO_NAME`/`_PLAIN_CANONICAL_TO_NAME`)."""
    query, query_idx = _QUERY_BY_PARENT.get(name, (_ANDROSTANE_QUERY, _ANDROSTANE_QUERY_IDX))
    matches = mol.GetSubstructMatches(query, uniquify=False)
    assert len(matches) == 1, f"{name}'s steroid-numbering query matched {len(matches)} times, expected 1"
    match = matches[0]
    return {locant: match[idx] for locant, idx in query_idx.items()}


def steroid_suffix_name(mol, suffix_atom, attachment_carbon, suffix_word, elide_e=True):
    """Suffix-agnostic: if the ring system, with `suffix_atom` (the suffix
    group's own heteroatom -- a ketone's =O or an alcohol's -OH oxygen)
    removed, exactly matches one of the seven bare steroid parent skeletons
    above (constitution only, stereochemistry ignored -- see
    `_PLAIN_CANONICAL_TO_NAME`'s own docstring note), return the retained
    steroid name with `attachment_carbon`'s fixed steroid locant and
    `suffix_word` suffixed on (e.g. 'androstan-3-one', 'androstan-3-ol').
    Otherwise return None -- a non-steroid polycyclic parent, or a steroid
    skeleton not among the seven recognized here."""
    stripped, old_to_new = strip_substituents(mol, {suffix_atom})
    if stripped is None:
        return None
    new_attachment_carbon = old_to_new[attachment_carbon]
    Chem.RemoveStereochemistry(stripped)
    name = _PLAIN_CANONICAL_TO_NAME.get(Chem.MolToSmiles(stripped))
    if name is None:
        return None

    locant_of_atom = {atom: locant for locant, atom in _locant_map(name, stripped).items()}
    locant = locant_of_atom.get(new_attachment_carbon)
    if locant is None:
        return None
    return _suffixed_parent(name, locant, suffix_word, elide_e=elide_e)


# P-101.2.6.1.1's alpha/beta symbolism (`tmp/bluebook/P10.txt` lines
# 222-246): "an atom or group attached to the ring is called 'alpha' if it
# lies below or 'beta' if it lies above the plane of the paper" in the
# standard steroid projection. The same primary text's own worked example
# ("5beta,9beta,10alpha-pregnane") states which configuration each of
# C8/C9/C10/C13/C14 is in *when unspecified/matching the fundamental
# parent structure* (8/10/13 = beta, 9/14 = alpha) -- those five are
# therefore only cited when *inverted* from that implied configuration.
# C5 is different: the same text says C5's configuration "when relevant,
# is indicated by alpha, beta, or xi" with no default implied by the bare
# parent name at all, so it's always cited once specified -- consistent
# with this module's own pre-existing convention of registering a second,
# fully C5-specified natural-configuration entry per skeleton under a
# "5alpha-<name>" PubChem name (see the module docstring), confirming
# 5-alpha is the label used for the natural/commonly-drawn configuration.
_ALPHA_BETA_NATURAL_LABEL = {5: "alpha", 8: "beta", 9: "alpha", 10: "beta", 13: "beta", 14: "alpha"}
_ALPHA_BETA_OPPOSITE = {"alpha": "beta", "beta": "alpha"}
_ALPHA_BETA_ALWAYS_CITE = frozenset({5})


def _natural_cip_by_locant(name):
    """The androstane-numbered ring-fusion locants (a subset of
    5/8/9/10/13/14, whichever exist as stereocenters for this particular
    skeleton -- gonane/estrane lack one or both angular methyls) mapped to
    their real CIP code in `name`'s own natural-configuration reference
    entry above, computed fresh per family rather than reused across
    families: the same spatial (alpha/beta) configuration gets a
    *different* CIP letter in different families, since CIP priority
    depends on the whole local substituent environment (e.g. the C17 side
    chain), not just spatial arrangement -- confirmed by direct
    computation (androstane's natural C13 is CIP 'S', pregnane's is 'R',
    same spatial 'beta' configuration both times). Picks whichever
    stereo-specified entry for `name` has the most fusion locants with a
    specified CIP code (the fully C5-specified second-block entry, when
    one exists per the module docstring's own note about which skeletons
    have it)."""
    best = {}
    for smiles, entry_name in _PARENT_HYDRIDES.items():
        if entry_name != name or "@" not in smiles:
            continue
        mol = Chem.MolFromSmiles(smiles)
        rdCIPLabeler.AssignCIPLabels(mol)
        locant_map = _locant_map(name, mol)
        codes = {}
        for locant in _ALPHA_BETA_NATURAL_LABEL:
            atom_idx = locant_map.get(locant)
            if atom_idx is None:
                continue
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.HasProp("_CIPCode"):
                codes[locant] = atom.GetProp("_CIPCode")
        if len(codes) > len(best):
            best = codes
    return best


_NATURAL_CIP_BY_LOCANT_BY_NAME = {name: _natural_cip_by_locant(name) for name in _PLAIN_CANONICAL_TO_NAME.values()}


def _alpha_beta_stereo_prefix(name, mol):
    """`"<locants>-"` alpha/beta stereodescriptor prefix (P-101.2.6.1.1)
    for `mol`, an occurrence of the named steroid parent skeleton `name`
    (already confirmed via `_PLAIN_CANONICAL_TO_NAME`), or `""` if `mol`
    has no specified stereochemistry at all, or `None` if it has some
    specified ring-fusion stereocenters but not every locant this family
    has a natural-configuration reference for -- a partially-specified
    input is out of scope, same all-or-nothing policy the module's prior
    single-locant mechanism already had. Locant 5, when specified, is
    always cited (the bare parent name never implies its configuration);
    each of 8/9/10/13/14 is cited only when *inverted* from the implied
    natural configuration, omitted when it matches."""
    natural = _NATURAL_CIP_BY_LOCANT_BY_NAME.get(name, {})
    if not natural:
        return ""
    if not any(atom.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED for atom in mol.GetAtoms()):
        return ""
    rdCIPLabeler.AssignCIPLabels(mol)
    locant_map = _locant_map(name, mol)
    parts = []
    for locant, natural_cip in sorted(natural.items()):
        atom_idx = locant_map.get(locant)
        atom = mol.GetAtomWithIdx(atom_idx) if atom_idx is not None else None
        if atom is None or not atom.HasProp("_CIPCode"):
            return None
        actual_cip = atom.GetProp("_CIPCode")
        natural_label = _ALPHA_BETA_NATURAL_LABEL[locant]
        if actual_cip == natural_cip:
            if locant in _ALPHA_BETA_ALWAYS_CITE:
                parts.append((locant, natural_label))
        else:
            parts.append((locant, _ALPHA_BETA_OPPOSITE[natural_label]))
    if not parts:
        return ""
    return ",".join(f"{locant}{label}" for locant, label in parts) + "-"


def _stereo_specified_parent_hydride_name(mol):
    """If `mol`'s constitution (ignoring stereo) matches one of the seven
    bare steroid parent skeletons, and its ring-fusion stereocenters are
    either all unspecified or all specified, return the alpha/beta-
    prefixed name (e.g. "5beta,9beta,10alpha-pregnane", or plain
    "androstane" if every fusion locant matches the natural configuration
    and none needs citing). Returns None for a constitution mismatch or a
    partially-specified stereocenter (out of scope, same as before)."""
    stripped = Chem.Mol(mol)
    Chem.RemoveStereochemistry(stripped)
    name = _PLAIN_CANONICAL_TO_NAME.get(Chem.MolToSmiles(stripped))
    if name is None:
        return None
    prefix = _alpha_beta_stereo_prefix(name, mol)
    if prefix is None:
        return None
    return prefix + name


def _steroid_unsaturated_name(mol):
    """If `mol` is exactly one of the seven bare steroid parent skeletons
    (constitution only, per `_PLAIN_CANONICAL_TO_NAME`) plus exactly one
    ring C=C double bond at a standard, sequential steroid locant pair
    (n, n+1) -- not a ring-fusion pair like C5-C10, whose own compound-
    locant citation (e.g. 'estra-5(10)-ene') is separate, unimplemented
    machinery -- return the retained name with an '-ene' suffix at the
    double bond's lower locant (e.g. 'androst-5-ene'). Otherwise None.

    Mirrors `steroid_suffix_name`'s own strip-then-match-then-locate
    approach, except no atom is removed here (only a bond order change),
    so the original atom indices stay valid for `_locant_map` directly."""
    double_bonds = [
        (bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
        for bond in mol.GetBonds()
        if bond.GetBondTypeAsDouble() == 2.0
    ]
    if len(double_bonds) != 1:
        return None
    (a, b) = double_bonds[0]
    if mol.GetAtomWithIdx(a).GetAtomicNum() != 6 or mol.GetAtomWithIdx(b).GetAtomicNum() != 6:
        return None

    rw = Chem.RWMol(mol)
    rw.GetBondBetweenAtoms(a, b).SetBondType(BondType.SINGLE)
    for idx in (a, b):
        atom = rw.GetAtomWithIdx(idx)
        atom.SetNoImplicit(False)
        atom.SetNumExplicitHs(0)
    stripped = rw.GetMol()
    try:
        Chem.SanitizeMol(stripped)
    except Chem.rdchem.KekulizeException:
        return None
    Chem.RemoveStereochemistry(stripped)
    name = _PLAIN_CANONICAL_TO_NAME.get(Chem.MolToSmiles(stripped))
    if name is None:
        return None

    locant_of_atom = {atom: locant for locant, atom in _locant_map(name, stripped).items()}
    locant_a, locant_b = locant_of_atom.get(a), locant_of_atom.get(b)
    if locant_a is None or locant_b is None:
        return None
    lower, upper = sorted((locant_a, locant_b))
    if upper != lower + 1:
        return None
    return name[:-3] + f"-{lower}-ene"


def _steroid_aromatic_a_ring_name(mol):
    """If `mol` has exactly one aromatic ring, a 6-membered benzo ring at
    the standard steroid A-ring locants (1,2,3,4,5,10), and de-aromatizing
    it (to plain single bonds, recomputing each ring atom's H count) yields
    exactly one of the seven bare steroid parent skeletons, return the
    retained name with the standard mancude-aromatic-A-ring suffix, e.g.
    'estra-1,3,5(10)-triene'. The A-ring's locants are fixed by P-31
    steroid numbering for every recognized parent, so the triene locant
    set is always this one literal string, never computed per-molecule.
    Otherwise None -- any other aromatic ring shape, position, or count is
    out of scope here."""
    aromatic_atoms = {a.GetIdx() for a in mol.GetAtoms() if a.GetIsAromatic()}
    if len(aromatic_atoms) != 6:
        return None
    ring_info = mol.GetRingInfo()
    aromatic_rings = [
        ring
        for ring in ring_info.AtomRings()
        if len(ring) == 6 and all(idx in aromatic_atoms for idx in ring)
    ]
    if len(aromatic_rings) != 1:
        return None
    (ring,) = aromatic_rings
    if any(mol.GetAtomWithIdx(idx).GetAtomicNum() != 6 for idx in ring):
        return None

    rw = Chem.RWMol(mol)
    ring_set = set(ring)
    for bond in list(rw.GetBonds()):
        if bond.GetBeginAtomIdx() in ring_set and bond.GetEndAtomIdx() in ring_set:
            bond.SetBondType(BondType.SINGLE)
            bond.SetIsAromatic(False)
    for idx in ring:
        atom = rw.GetAtomWithIdx(idx)
        atom.SetIsAromatic(False)
        atom.SetNoImplicit(False)
        atom.SetNumExplicitHs(0)
    stripped = rw.GetMol()
    try:
        Chem.SanitizeMol(stripped)
    except Chem.rdchem.KekulizeException:
        return None
    Chem.RemoveStereochemistry(stripped)
    name = _PLAIN_CANONICAL_TO_NAME.get(Chem.MolToSmiles(stripped))
    if name is None:
        return None

    locant_of_atom = {atom: locant for locant, atom in _locant_map(name, stripped).items()}
    if {locant_of_atom.get(idx) for idx in ring} != {1, 2, 3, 4, 5, 10}:
        return None
    return name[:-3] + "a-1,3,5(10)-triene"


def has_steroid_parent_hydride_name(mol) -> bool:
    if Chem.MolToSmiles(mol) in _CANONICAL_TO_NAME:
        return True
    return _stereo_specified_parent_hydride_name(mol) is not None


def name_steroid_parent_hydride(mol) -> str:
    key = Chem.MolToSmiles(mol)
    if key in _CANONICAL_TO_NAME:
        return _CANONICAL_TO_NAME[key]
    return _stereo_specified_parent_hydride_name(mol)


def has_steroid_unsaturated_name(mol) -> bool:
    return _steroid_unsaturated_name(mol) is not None


def name_steroid_unsaturated(mol) -> str:
    return _steroid_unsaturated_name(mol)


def has_steroid_aromatic_a_ring_name(mol) -> bool:
    return _steroid_aromatic_a_ring_name(mol) is not None


def name_steroid_aromatic_a_ring(mol) -> str:
    return _steroid_aromatic_a_ring_name(mol)
