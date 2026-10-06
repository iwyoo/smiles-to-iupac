"""Naming of acetals and ketals (RR'C(O-R'')(O-R'''), two alkoxy
substituents on the same carbon) on acyclic saturated hydrocarbon chains,
per the IUPAC 2013 Recommendations ("the Blue Book"):

- P-66.6.5.1 (Chapter P-6a, https://iupac.qmul.ac.uk/BlueBook/PDF/P6a.pdf):
  an acetal/ketal has no principal-characteristic-group suffix of its
  own -- its preferred IUPAC name is a plain alkane parent hydride
  substituted by two 'alkoxy' prefixes (P-63.2.2.1.1), exactly the same
  substitutive mechanism `_ether.py` already uses for a single ether
  oxygen, just with two such oxygens landing on the same carbon.
  Confirmed via PubChem PUG REST: CID 20858 (`CCC(OCC)OCC`) ->
  "1,1-diethoxypropane" (matching the Blue Book's own worked example
  '1,1-diethoxypropane (PIN)' directly), CID 10795 (`COC(OC)C`) ->
  "1,1-dimethoxyethane", CID 6495 (`COC(C)(C)OC`, a ketal -- neither R nor
  R' is hydrogen) -> "2,2-dimethoxypropane".
- This module reuses `_acyclic.py`'s `name_from_carbon_graph` exactly the
  way `_ether.py`/`_nitro.py`/`_peroxide.py` do: each acetal oxygen is
  passed in as a precomputed 'terminals' entry (a leaf substituent
  excluded from the carbon-only chain search), so two identical alkoxy
  groups combine through the same multiplying-prefix machinery already
  used for 'dichloro'/'dinitro' (P-14.5.2), giving 'di-' + one shared
  locant set (e.g. '1,1-diethoxypropane'), and two different alkoxy
  groups cite each individually in alphabetical order (matching the
  ether-style module's own reuse of that mechanism).
- P-63.2.2.1.1: each alkoxy prefix reuses `_ether.py`'s own
  `_oxy_prefix`/contracted-name table ('methoxy', 'ethoxy', 'propoxy',
  'butoxy' for the four shortest unbranched chains).
- P-91.3/P-92: a molecule with one
  or more *specified* tetrahedral stereocenters -- every one on the
  principal chain itself (the acetal carbon included, when its two
  alkoxy groups differ enough to make it a genuine stereocenter -- when
  they're identical it's symmetric and never flagged by RDKit's
  `Chem.FindPotentialStereo`), no unspecified one alongside them, and no
  C=C/C#N double-bond E/Z element -- gets a "(<locant><R/S>,...)-"
  prefix, ascending locant order, e.g. '(2R)-1,1-dimethoxy-2-methylbutane'
  (PubChem CID 90324169), '(1R)-1-ethoxy-1-methoxypropane' (PubChem CID
  97550416). Uses `_acyclic.py`'s `winning_chain_from_carbon_graph`
  (added alongside `name_from_carbon_graph`, which every other caller of
  that module keeps using unchanged) to locate the winning chain's own
  locant for each stereocenter.

Explicitly out of scope (raise `UnsupportedStructure`), mirroring
`_ether.py`'s own first-pass scope:
- A branched O-substituent (R'' or R'''), for the same enclosure-
  interaction reason `_ether.py` excludes a branched R'.
- More than one acetal/ketal-shaped carbon in the same molecule.
- A cyclic acetal/ketal (P-66.6.5.1.2, a separate heterocyclic/spiro
  naming problem) or any ring anywhere in the molecule.
- Any heteroatom other than the two acetal oxygens (in particular, no
  verified halogen coexistence yet, unlike `_ether.py`).
- Any unsaturation, or an aromatic ring.
- A hemiacetal/hemiketal (RR'C(OH)(O-R''), P-66.6.5.2) -- a different,
  unverified shape with one hydroxyl instead of a second alkoxy group.
"""

from rdkit import Chem

from ._acyclic import winning_chain_from_carbon_graph
from ._common import (
    UnsupportedStructure,
    adjacency,
    carbon_adjacency,
    component_subgraph,
    non_single_bonds,
    specified_stereocenters,
)
from ._ether import _oxy_prefix
from ._substituents import name_branch


def has_acetal_shape(mol) -> bool:
    """True iff some carbon carries exactly two singly-bonded, degree-2
    oxygens each also bonded to a carbon (a plain acetal/ketal -O-C
    pattern on both sides, P-66.6.5.1) -- regardless of whether the rest
    of the molecule is in scope."""
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6:
            continue
        oxy_neighbors = [
            n
            for n in atom.GetNeighbors()
            if n.GetAtomicNum() == 8
            and n.GetDegree() == 2
            and mol.GetBondBetweenAtoms(atom.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 1.0
        ]
        if len(oxy_neighbors) == 2:
            return True
    return False


def _validate_and_collect_acetal(mol):
    """Check the molecule fits this module's scope (see module docstring)
    and return (acetal_carbon, oxygen_1, oxygen_2): the single acetal
    carbon's atom index and its two alkoxy-oxygen atom indices."""
    has_carbon = False
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atomic_num not in (6, 8):
            raise UnsupportedStructure(
                "heteroatoms other than the two acetal/ketal alkoxy "
                "oxygens (P-66.6.5.1) are not supported yet"
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
        else:
            if atom.GetDegree() != 2 or any(n.GetAtomicNum() != 6 for n in atom.GetNeighbors()):
                raise UnsupportedStructure(
                    "an oxygen that isn't a plain ether-type -O- bonded to "
                    "two carbons is out of scope for this module (e.g. a "
                    "hydroxyl or hemiacetal -OH, P-66.6.5.2)"
                )
            if any(
                mol.GetBondBetweenAtoms(atom.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() != 1.0
                for n in atom.GetNeighbors()
            ):
                raise UnsupportedStructure("an acetal/ketal oxygen must be singly bonded to both carbons")
    if not has_carbon:
        raise UnsupportedStructure(
            "a structure with no carbon atom has no hydrocarbon parent "
            "hydride to substitute"
        )
    if mol.GetRingInfo().NumRings() > 0:
        raise UnsupportedStructure(
            "a cyclic acetal/ketal is out of scope for this acyclic-only "
            "module (P-66.6.5.1.2 names it as a separate heterocyclic/"
            "spiro compound)"
        )
    if non_single_bonds(mol):
        raise UnsupportedStructure("unsaturation is not supported by this module yet")

    acetal_carbons = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6:
            continue
        oxy_neighbors = [n for n in atom.GetNeighbors() if n.GetAtomicNum() == 8]
        if len(oxy_neighbors) == 2:
            acetal_carbons.append((atom.GetIdx(), [o.GetIdx() for o in oxy_neighbors]))
    if len(acetal_carbons) != 1:
        raise UnsupportedStructure(
            "more than one acetal/ketal group is out of scope for this "
            "module"
        )
    (acetal_carbon, oxygens) = acetal_carbons[0]
    oxygen_1, oxygen_2 = oxygens
    total_oxygens = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if total_oxygens != 2:
        raise UnsupportedStructure(
            "an oxygen elsewhere in the molecule besides the single "
            "acetal/ketal group's own two alkoxy oxygens is out of scope "
            "for this module"
        )
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")
    return acetal_carbon, oxygen_1, oxygen_2


def name_acetal(mol) -> str:
    acetal_carbon, oxygen_1, oxygen_2 = _validate_and_collect_acetal(mol)
    stereo = specified_stereocenters(mol)
    full_graph = adjacency(mol)
    carbon_graph = component_subgraph(carbon_adjacency(mol), acetal_carbon)

    terminals = {}
    for oxygen_idx in (oxygen_1, oxygen_2):
        (sub_root,) = (idx for idx in full_graph[oxygen_idx] if idx != acetal_carbon)
        sub_name, sub_compound = name_branch(full_graph, sub_root, oxygen_idx, {}, mol=mol)
        if sub_compound:
            raise UnsupportedStructure(
                "a branched alkoxy substituent's enclosing marks are not "
                "supported yet (see P-63.2.2.1.1, mirroring _ether.py's "
                "own restriction)"
            )
        terminals[oxygen_idx] = _oxy_prefix(sub_name)

    chain, name = winning_chain_from_carbon_graph(full_graph, carbon_graph, terminals, mol=mol)
    if stereo is None:
        return name

    # P-92: every specified stereocenter -- the acetal carbon itself (only
    # a genuine stereocenter when the two alkoxy groups differ) and/or any
    # chain carbon -- must lie on the winning parent chain; one on a
    # substituent branch is out of scope, mirroring
    # `_carboxylic_acid.py`'s identical treatment.
    position_of = {atom: i + 1 for i, atom in enumerate(chain)}
    if any(atom not in position_of for atom, _ in stereo):
        raise UnsupportedStructure(
            "a stereocenter on a substituent branch rather than the "
            "principal chain is not supported yet (see P-92)"
        )
    labels = sorted((position_of[atom], code) for atom, code in stereo)
    prefix = ",".join(f"{locant}{code}" for locant, code in labels)
    return f"({prefix})-{name}"
