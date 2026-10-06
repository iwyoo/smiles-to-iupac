"""Naming of a molecule combining exactly one primary amide (-C(=O)NH2) with
one or more internal ketone (=O) carbonyls on the same acyclic saturated
carbon chain, per the IUPAC 2013 Recommendations ("the Blue Book"):

- P-41, Table 3.3: 'amide' (`_amide.py`) outranks 'one' (`_ketone.py`), so a
  coexisting ketone is demoted to the 'oxo' substituent prefix instead of
  its own '-one' suffix, e.g. 'CC(=O)CC(N)=O' (acetoacetamide) ->
  '3-oxobutanamide' (a well-known worked example). This mirrors
  `_aldehyde_ketone.py`'s aldehyde+ketone demotion, reusing the same
  {atom_idx -> 'oxo'} injection into the chain's substituent-prefix
  machinery.
- Otherwise mirrors `_amide.py` exactly: the amide carbon is always a chain
  terminus and always becomes C1 with its own locant never cited
  (P-14.3.3); a ketone carbon's locant, by contrast, is always cited.

Scope, deliberately narrow (first extension of the general suffix-vs-suffix
demotion problem beyond aldehyde+ketone): a single primary amide plus one
or more ketones, all on one acyclic saturated chain, with halogen substituents
allowed. Explicitly out of scope (raise `UnsupportedStructure`): any chain
unsaturation (ene/yne), any ring, an N-substituted amide, more than one
amide, a coexisting hydroxyl/ether/other heteroatom, and any ketone not
captured by a single longest chain.
"""

from ._amide import _name_acyclic_amide
from ._coexisting_groups import name_via_senior_acyclic
from ._common import (
    UnsupportedStructure,
    adjacency,
    group_substituents,
    halogen_substituents,
    is_plain_benzene_ring,
    longest_branched_chain,
    name_from_substituents,
    non_single_bonds,
    ring_chain_attachment,
    validate_allowed_atoms,
)
from ._substituents import format_substituent_prefixes, substituents_for_chain


def _find_amide_carbon(mol):
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6:
            continue
        carbonyls = [
            n
            for n in atom.GetNeighbors()
            if n.GetAtomicNum() == 8
            and n.GetDegree() == 1
            and mol.GetBondBetweenAtoms(atom.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 2.0
        ]
        amide_nitrogens = [
            n
            for n in atom.GetNeighbors()
            if n.GetAtomicNum() == 7
            and n.GetDegree() == 1
            and n.GetTotalNumHs() == 2
            and mol.GetBondBetweenAtoms(atom.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 1.0
        ]
        carbon_neighbors = [n for n in atom.GetNeighbors() if n.GetAtomicNum() == 6]
        if len(carbonyls) == 1 and len(amide_nitrogens) == 1 and len(carbon_neighbors) <= 1:
            return atom, carbonyls[0], amide_nitrogens[0]
    return None


def _extra_ketones(mol, excluded_oxygens):
    ketones = set()
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 8 or atom.GetIdx() in excluded_oxygens:
            continue
        if atom.GetDegree() != 1:
            continue
        (bond,) = atom.GetBonds()
        if bond.GetBondTypeAsDouble() != 2.0:
            continue
        (carbon,) = atom.GetNeighbors()
        if carbon.GetAtomicNum() != 6:
            continue
        carbon_neighbors = [n for n in carbon.GetNeighbors() if n.GetAtomicNum() == 6]
        if len(carbon_neighbors) == 2:
            ketones.add(atom.GetIdx())
    return ketones


def has_ketone_amide_shape(mol) -> bool:
    found = _find_amide_carbon(mol)
    if found is None:
        return False
    _, amide_oxygen, _ = found
    return bool(_extra_ketones(mol, {amide_oxygen.GetIdx()}))


def _validate(mol, excluded_oxygens, amide_nitrogen, aromatic_ring_atoms=frozenset()):
    validate_allowed_atoms(
        mol,
        "heteroatoms other than the amide's own oxygen/nitrogen and "
        "ketone carbonyl oxygens (P-41, Table 3.3) and halogen "
        "substituents (P-35.2.1) are not supported yet",
        [
            (
                8,
                excluded_oxygens,
                "an oxygen that isn't the amide's own carbonyl oxygen or a "
                "ketone-shaped carbonyl is out of scope for this module "
                "(e.g. a coexisting hydroxyl or ether)",
            ),
            (
                7,
                {amide_nitrogen.GetIdx()},
                "more than one nitrogen, or a nitrogen that isn't the "
                "amide's own, is out of scope for this module",
            ),
        ],
        aromatic_ring_atoms=aromatic_ring_atoms,
    )


def _name_from_substituents(chain_length, grouped):
    return format_substituent_prefixes(grouped) + name_from_substituents(chain_length, [], [], "amide")


def _name_phenyl_chain_ketone_amide(mol, ring_atoms):
    """Name a ketone+amide combination whose -CONH2 lies entirely on a
    single unbranched chain hanging off one atom of an otherwise-plain,
    unsubstituted benzene ring -- e.g. 3-oxo-4-phenylbutanamide. The ring
    is cited as a 'phenyl' substituent prefix (via `name_branch`'s
    aromatic-ring recognition) on the chain, which is the parent hydride,
    mirroring `_aldehyde_ketone.py`'s `_name_phenyl_chain_aldehyde_ketone`.
    Narrower than the acyclic path above: no chain unsaturation."""
    found = _find_amide_carbon(mol)
    if found is None:
        raise UnsupportedStructure(
            "no primary amide (-CONH2) group found; this module only "
            "handles primary amides"
        )
    amide_carbon, amide_oxygen, amide_nitrogen = found
    ketones = _extra_ketones(mol, {amide_oxygen.GetIdx()})
    if not ketones:
        raise UnsupportedStructure(
            "no coexisting ketone found; this module only handles an amide "
            "combined with at least one ketone (see _amide.py for a plain "
            "amide)"
        )
    own_excluded = {amide_oxygen.GetIdx(), amide_nitrogen.GetIdx()}
    excluded = own_excluded | ketones
    _validate(mol, excluded, amide_nitrogen, aromatic_ring_atoms=ring_atoms)

    all_non_single = non_single_bonds(mol)
    carbonyl_bonds = [b for b in all_non_single if b[0] in excluded or b[1] in excluded]
    ring_bonds = [b for b in all_non_single if b[0] in ring_atoms and b[1] in ring_atoms]
    if len(carbonyl_bonds) + len(ring_bonds) != len(all_non_single):
        raise UnsupportedStructure(
            "chain unsaturation (ene/yne) alongside a benzene-ring-"
            "substituent ketone/amide combination is out of scope for "
            "this module"
        )

    graph = adjacency(mol)
    names = {**halogen_substituents(mol), **{o: "oxo" for o in ketones}}
    amide_carbon_idx = amide_carbon.GetIdx()
    attachment = ring_chain_attachment(graph, ring_atoms, set())
    if attachment is None:
        raise UnsupportedStructure(
            "a benzene ring with more than one exocyclic substituent "
            "alongside a chain ketone/amide combination is not supported "
            "yet"
        )
    chain, _ = longest_branched_chain(graph, amide_carbon_idx, ring_atoms, own_excluded | set(names), halogens=halogen_substituents(mol))
    if any(graph[o][0] not in chain for o in ketones):
        raise UnsupportedStructure(
            "not every ketone-bearing carbon lies on the chain hanging "
            "off the benzene ring for this benzene-substituent path"
        )
    if len(chain) < 2:
        raise UnsupportedStructure(
            "an amide group directly attached to the benzene ring uses a "
            "separate naming construction, out of scope for this "
            "acyclic-chain-parent module"
        )

    chain_length = len(chain)
    substituents = substituents_for_chain(graph, chain, names, own_excluded, mol=mol, aromatic_atoms=ring_atoms)
    grouped = group_substituents(substituents)
    return _name_from_substituents(chain_length, grouped)


def name_ketone_amide(mol) -> str:
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() == 1:
        ring_atoms = set(ring_info.AtomRings()[0])
        if is_plain_benzene_ring(mol, ring_atoms):
            return _name_phenyl_chain_ketone_amide(mol, ring_atoms)
    if mol.GetRingInfo().NumRings() > 0:
        raise UnsupportedStructure(
            "an amide on/in a ring (a lactam) is out of scope for this "
            "acyclic-only module"
        )
    found = _find_amide_carbon(mol)
    if found is None:
        raise UnsupportedStructure(
            "no primary amide (-CONH2) group found; this module only "
            "handles primary amides"
        )
    amide_carbon, amide_oxygen, amide_nitrogen = found
    ketones = _extra_ketones(mol, {amide_oxygen.GetIdx()})
    if not ketones:
        raise UnsupportedStructure(
            "no coexisting ketone found; this module only handles an amide "
            "combined with at least one ketone (see _amide.py for a plain "
            "amide)"
        )
    own_excluded = {amide_oxygen.GetIdx(), amide_nitrogen.GetIdx()}
    excluded = own_excluded | ketones
    _validate(mol, excluded, amide_nitrogen)

    all_non_single = non_single_bonds(mol)
    carbonyl_bonds = [b for b in all_non_single if b[0] in excluded or b[1] in excluded]
    if len(carbonyl_bonds) != len(all_non_single):
        raise UnsupportedStructure(
            "chain unsaturation (ene/yne) combined with a ketone-on-amide "
            "demotion is out of scope for this module"
        )

    graph = adjacency(mol)
    ketone_carbons = {graph[o][0] for o in ketones}
    return name_via_senior_acyclic(
        _name_acyclic_amide,
        "amide",
        "ketone",
        (mol, amide_carbon.GetIdx(), amide_nitrogen.GetIdx(), own_excluded, (), set(), []),
        {o: "oxo" for o in ketones},
        required_atoms=ketone_carbons,
    )
