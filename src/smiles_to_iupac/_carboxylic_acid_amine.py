"""Naming of a molecule combining exactly one carboxylic acid (-COOH) with
exactly one primary amine (-NH2) on the same acyclic saturated carbon chain,
per the IUPAC 2013 Recommendations ("the Blue Book"):

- P-41, Table 3.3: 'oic acid' (`_carboxylic_acid.py`) far outranks 'amine'
  (`_amine.py`), so a coexisting primary amine is demoted to the 'amino'
  substituent prefix instead of its own '-amine' suffix, e.g.
  'NCCCC(=O)O' (4-aminobutanoic acid, GABA) -> '4-aminobutanoic acid'
  (PubChem CID 119's own IUPACName). This mirrors
  `_aldehyde_carboxylic_acid.py`'s carbonyl demotion, except the demoted
  group here is a substituent atom (the amine nitrogen itself), not a
  chain carbon,
  reusing `_carboxylic_acid.py`'s own `{atom_idx -> 'hydroxy'}` injection
  pattern for a coexisting standalone hydroxyl, just with 'amino' instead.
- Otherwise mirrors `_carboxylic_acid.py` exactly: the -COOH carbon is
  always a chain terminus and always becomes C1 with its own locant never
  cited (P-14.3.3); the amine nitrogen's locant, by contrast, is always
  cited via the 'amino' prefix.
- Glycine (H2N-CH2-COOH) is a real-world worked example, but its actual PIN
  ('2-aminoacetic acid', PubChem CID 750) uses the retained acid name
  'acetic acid' rather than the systematic 'ethanoic acid' this project's
  own `_carboxylic_acid.py` always produces (see that module's own tests,
  e.g. 'CC(=O)O' -> 'ethanoic acid', not 'acetic acid') -- this module
  follows that same already-established project convention, so it produces
  '2-aminoethanoic acid' for glycine's structure instead, a pre-existing
  divergence from PubChem inherited from `_carboxylic_acid.py`, not a new
  one introduced here.

Scope, deliberately narrow: a single carboxylic acid plus a single primary
amine, both on one acyclic *saturated* chain, with halogen substituents
allowed.
Explicitly out of scope (raise `UnsupportedStructure`): any chain
unsaturation (ene/yne), any ring, more than one carboxylic acid or amine, a
secondary/tertiary amine, a coexisting hydroxyl/ether/other heteroatom, and
any acid/amine not captured by a single longest chain.

Any *specified* tetrahedral stereocenter (e.g. the alpha carbon of an amino
acid such as alanine/valine) is labeled via `_common.specified_stereocenters`
using the same P-91.3/P-92 mechanism as `_carboxylic_acid.py` (P-92: a
stereocenter on a substituent branch rather than the principal chain is out
of scope).
"""

from rdkit import Chem

from ._coexisting_groups import name_via_senior_acyclic
from ._common import (
    UnsupportedStructure,
    adjacency,
    find_primary_amines,
    group_substituents,
    halogen_substituents,
    is_plain_benzene_ring,
    longest_branched_chain,
    name_from_substituents,
    non_single_bonds,
    ring_chain_attachment,
    specified_stereo_elements,
    specified_stereocenters,
    substituent_locant_set_and_citation,
    validate_allowed_atoms,
)
from ._carboxylic_acid import _name_acyclic_carboxylic_acid
from ._substituents import format_substituent_prefixes, name_branch, substituents_for_chain



def _find_carboxylic_acid_carbon(mol):
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6:
            continue
        oxygens = [n for n in atom.GetNeighbors() if n.GetAtomicNum() == 8]
        if len(oxygens) < 2:
            continue
        carbonyls = [
            o
            for o in oxygens
            if o.GetDegree() == 1 and mol.GetBondBetweenAtoms(atom.GetIdx(), o.GetIdx()).GetBondTypeAsDouble() == 2.0
        ]
        hydroxyls = [
            o
            for o in oxygens
            if o.GetDegree() == 1
            and mol.GetBondBetweenAtoms(atom.GetIdx(), o.GetIdx()).GetBondTypeAsDouble() == 1.0
            and o.GetTotalNumHs() == 1
        ]
        carbon_neighbors = [n for n in atom.GetNeighbors() if n.GetAtomicNum() == 6]
        if len(carbonyls) == 1 and len(hydroxyls) == 1 and len(carbon_neighbors) <= 1:
            return atom, carbonyls[0], hydroxyls[0]
    return None

def _amine_on_a_different_carbon(mol, acid_carbon_idx, amines):
    """False if some amine nitrogen is bonded directly to the acid carbon
    itself (H2N-COOH, carbamic acid -- a distinct retained functional class,
    `_carbamate.py`'s territory, not this module's "amine elsewhere on the
    chain" shape)."""
    return all(next(iter(mol.GetAtomWithIdx(n).GetNeighbors())).GetIdx() != acid_carbon_idx for n in amines)


def has_carboxylic_acid_amine_shape(mol) -> bool:
    found = _find_carboxylic_acid_carbon(mol)
    if found is None:
        return False
    amines = find_primary_amines(mol)
    return bool(amines) and _amine_on_a_different_carbon(mol, found[0].GetIdx(), amines)


def _validate(mol, excluded_oxygens, amines, aromatic_ring_atoms=frozenset()):
    validate_allowed_atoms(
        mol,
        "heteroatoms other than the acid's own oxygens, a primary amine "
        "nitrogen (P-41, Table 3.3), and halogen substituents (P-35.2.1) "
        "are not supported yet",
        [
            (
                8,
                excluded_oxygens,
                "an oxygen that isn't part of the single carboxylic acid "
                "group is out of scope for this module (e.g. a coexisting "
                "hydroxyl, ether, or carbonyl)",
            ),
            (
                7,
                amines,
                "a nitrogen that isn't a plain primary amine (-NH2) is out "
                "of scope for this module (secondary/tertiary amines, "
                "imines, and nitriles are not supported)",
            ),
        ],
        aromatic_ring_atoms=aromatic_ring_atoms,
    )


def _name_from_substituents(chain_length, grouped):
    return format_substituent_prefixes(grouped) + name_from_substituents(chain_length, [], [], "oic acid")


def _candidate_key(chain_length, grouped):
    locant_set, _, citation_locants = substituent_locant_set_and_citation(grouped)
    name = _name_from_substituents(chain_length, grouped)
    return (locant_set, citation_locants, name), name


def _name_phenyl_chain_carboxylic_acid_amine(mol, ring_atoms):
    """Name a carboxylic-acid+amine combination (e.g. phenylalanine's
    2-amino-3-phenylpropanoic acid shape) whose -COOH lies entirely on a
    single unbranched chain hanging off one atom of an otherwise-plain,
    unsubstituted benzene ring. The ring is cited as a 'phenyl'
    substituent prefix (via `name_branch`'s aromatic-ring recognition) on
    the chain, which is the parent hydride, mirroring
    `_carboxylic_acid.py`'s `_name_phenyl_chain_carboxylic_acid`. Unlike
    `_aldehyde_carboxylic_acid.py`, the demoted amine is a substituent
    atom (like a halogen), not a chain carbon, so it can sit anywhere
    along the chain -- no aldehyde-style terminus conflict with the ring.
    Narrower than the acyclic path above: no chain unsaturation."""
    found = _find_carboxylic_acid_carbon(mol)
    if found is None:
        raise UnsupportedStructure(
            "no carboxylic acid (-COOH) group found; this module only "
            "handles carboxylic acids"
        )
    acid_carbon, carbonyl_oxygen, hydroxyl_oxygen = found
    excluded_acid_oxygens = {carbonyl_oxygen.GetIdx(), hydroxyl_oxygen.GetIdx()}
    amines = find_primary_amines(mol)
    if len(amines) != 1:
        raise UnsupportedStructure(
            "exactly one primary amine (-NH2) coexisting with the single "
            "carboxylic acid is supported here (see _carboxylic_acid.py "
            "for a plain carboxylic acid, _amine.py for a plain amine)"
        )
    if not _amine_on_a_different_carbon(mol, acid_carbon.GetIdx(), amines):
        raise UnsupportedStructure(
            "a primary amine bonded directly to the acid carbon itself "
            "(H2N-COOH, carbamic acid) is a distinct retained functional "
            "class, out of scope here (see _carbamate.py)"
        )
    _validate(mol, excluded_acid_oxygens, amines, aromatic_ring_atoms=ring_atoms)

    all_non_single = non_single_bonds(mol)
    carbonyl_bonds = [b for b in all_non_single if b[0] in excluded_acid_oxygens or b[1] in excluded_acid_oxygens]
    ring_bonds = [b for b in all_non_single if b[0] in ring_atoms and b[1] in ring_atoms]
    if len(carbonyl_bonds) + len(ring_bonds) != len(all_non_single):
        raise UnsupportedStructure(
            "chain unsaturation (ene/yne) alongside a benzene-ring-"
            "substituent carboxylic-acid/amine combination is out of "
            "scope for this module"
        )

    graph = adjacency(mol)
    names = {**halogen_substituents(mol), **{n: "amino" for n in amines}}
    acid_carbon_idx = acid_carbon.GetIdx()
    stereo = specified_stereocenters(mol)
    if stereo:
        raise UnsupportedStructure(
            "a specified stereocenter alongside a benzene-ring-substituent "
            "carboxylic-acid/amine chain is not supported yet"
        )

    attachment = ring_chain_attachment(graph, ring_atoms, set())
    if attachment is None:
        raise UnsupportedStructure(
            "a benzene ring with more than one exocyclic substituent "
            "alongside a chain carboxylic-acid/amine combination is not "
            "supported yet"
        )
    chain, _ = longest_branched_chain(graph, acid_carbon_idx, ring_atoms, excluded_acid_oxygens | set(names), halogens=halogen_substituents(mol))
    if any(graph[n][0] not in chain for n in amines):
        raise UnsupportedStructure(
            "not every amine-bearing carbon lies on the chain hanging "
            "off the benzene ring for this benzene-substituent path"
        )
    if len(chain) < 2:
        raise UnsupportedStructure(
            "a -COOH group directly attached to the benzene ring uses a "
            "different naming construction, out of scope for this "
            "acyclic-chain-parent module"
        )

    chain_length = len(chain)
    substituents = substituents_for_chain(graph, chain, names, excluded_acid_oxygens, mol=mol, aromatic_atoms=ring_atoms)
    grouped = group_substituents(substituents)
    return _name_from_substituents(chain_length, grouped)


def name_carboxylic_acid_amine(mol) -> str:
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() == 1:
        ring_atoms = set(ring_info.AtomRings()[0])
        if is_plain_benzene_ring(mol, ring_atoms):
            return _name_phenyl_chain_carboxylic_acid_amine(mol, ring_atoms)
    if mol.GetRingInfo().NumRings() > 0:
        raise UnsupportedStructure(
            "a -COOH and/or -NH2 group on/in a ring uses a different "
            "naming construction, out of scope for this acyclic-only module"
        )
    found = _find_carboxylic_acid_carbon(mol)
    if found is None:
        raise UnsupportedStructure(
            "no carboxylic acid (-COOH) group found; this module only "
            "handles carboxylic acids"
        )
    acid_carbon, carbonyl_oxygen, hydroxyl_oxygen = found
    excluded_acid_oxygens = {carbonyl_oxygen.GetIdx(), hydroxyl_oxygen.GetIdx()}
    amines = find_primary_amines(mol)
    if len(amines) != 1:
        raise UnsupportedStructure(
            "exactly one primary amine (-NH2) coexisting with the single "
            "carboxylic acid is supported here (see _carboxylic_acid.py "
            "for a plain carboxylic acid, _amine.py for a plain amine)"
        )
    if not _amine_on_a_different_carbon(mol, acid_carbon.GetIdx(), amines):
        raise UnsupportedStructure(
            "a primary amine bonded directly to the acid carbon itself "
            "(H2N-COOH, carbamic acid) is a distinct retained functional "
            "class, out of scope here (see _carbamate.py)"
        )
    _validate(mol, excluded_acid_oxygens, amines)

    all_non_single = non_single_bonds(mol)
    carbonyl_bonds = [b for b in all_non_single if b[0] in excluded_acid_oxygens or b[1] in excluded_acid_oxygens]
    if len(carbonyl_bonds) != len(all_non_single):
        raise UnsupportedStructure(
            "chain unsaturation (ene/yne) combined with a carboxylic-acid/"
            "amine seniority demotion is out of scope for this module "
            "(1st-pass scope: saturated only)"
        )

    graph = adjacency(mol)
    acid_carbon_idx = acid_carbon.GetIdx()
    stereo = specified_stereo_elements(mol)
    amine_carbons = {graph[n][0] for n in amines}
    return name_via_senior_acyclic(
        _name_acyclic_carboxylic_acid,
        "carboxylic_acid",
        "amine",
        (mol, {acid_carbon_idx}, excluded_acid_oxygens, set(), []),
        {n: "amino" for n in amines},
        stereo=stereo,
        required_atoms=amine_carbons,
    )
