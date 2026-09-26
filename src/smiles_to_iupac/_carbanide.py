"""Naming of a plain acyclic hydrocarbon carbanion (a formal-charge -1
carbon on an otherwise unfunctionalized carbon chain, the '-ide' suffix),
per the IUPAC 2013 Recommendations ("the Blue Book"):

- P-72.2.2.1 (Chapter P-7, https://iupac.qmul.ac.uk/BlueBook/P7.html): "An
  anion derived formally by the removal of one or more hydrons from any
  position of a neutral parent hydride is preferably named by using the
  suffix '-ide', with elision of the final letter 'e' of the parent
  hydride" -- e.g. 'methanide (PIN)' for H3C- (line 889), and P-71.1's own
  text confirms 'ethanide' (not 'ethyl anion') for CH3-CH2- (line 131). The
  suffix locant is cited the same way `_alcohol.py`'s '-ol' is: omitted at
  a mononuclear (one-carbon) parent and at a saturated two-carbon parent
  (P-14.3.4.2(a)/(b), `_common.should_omit_mononuclear_locants`, reused via
  `name_from_substituents`), otherwise always cited -- confirmed by
  P-72.2.2.1's own worked '2-(...)propan-2-ide' examples (e.g. line 2872),
  which always cite the '2-' locant on a 3+-carbon chain.
- A branched skeleton is named via the same longest-chain-through-the-
  anion-carbon search and substituent-branch machinery `_alkoxide.py` uses
  for its own oxygen substituent, e.g. `C[CH-]C` (isopropyl anion) ->
  'propan-2-ide'.

Explicitly out of scope (raise `UnsupportedStructure`), mirroring
`_alkoxide.py`'s own acyclic-only pilot scope but narrower still:
- Any ring anywhere in the molecule.
- Any halogen substituent or other heteroatom.
- Any unsaturation (a C=C/C#C bond) anywhere in the molecule.
- More than one carbanion center, or a charged/isotopically modified atom
  other than the single carbanion carbon's own formal charge -1.
- A "functionalized" carbanion where the anion carbon is also part of a
  characteristic group (e.g. the acetyl anion, '1-oxoethan-1-ide') --
  that combination needs its own, separate module.
"""

from rdkit import Chem

from ._common import (
    UnsupportedStructure,
    adjacency,
    carbon_adjacency,
    group_substituents,
    longest_chains,
    name_from_substituents,
    substituent_locant_set_and_citation,
)
from ._substituents import format_substituent_prefixes, substituents_for_chain


def _is_isocyanide_carbon(mol, atom):
    # An isocyanide carbon (-N+#C-, P-61.9) also carries a formal charge
    # -1, but it's a triple-bonded terminus paired with a +1 nitrogen, not
    # a plain hydrocarbon carbanion -- `_isocyanide.py` already handles
    # this shape (dispatched separately, later in `core.py`), so it must
    # be excluded here rather than misrouted into this module.
    if atom.GetDegree() != 1:
        return False
    (bond,) = atom.GetBonds()
    if bond.GetBondTypeAsDouble() != 3.0:
        return False
    (neighbor,) = atom.GetNeighbors()
    return neighbor.GetAtomicNum() == 7 and neighbor.GetFormalCharge() == 1


def _find_carbanide_carbons(mol):
    return [
        atom
        for atom in mol.GetAtoms()
        if atom.GetAtomicNum() == 6 and atom.GetFormalCharge() == -1 and not _is_isocyanide_carbon(mol, atom)
    ]


def has_carbanide_shape(mol) -> bool:
    return bool(_find_carbanide_carbons(mol))


def _validate_and_find_carbanide(mol):
    matches = _find_carbanide_carbons(mol)
    if len(matches) != 1:
        raise UnsupportedStructure(
            "exactly one carbanion (-1 charge) carbon is required; zero or "
            "multiple such centers are not supported yet (P-72.2.2.1)"
        )
    (anion,) = matches
    if anion.GetIsotope() != 0:
        raise UnsupportedStructure(
            "an isotopically modified carbanion center is not supported yet"
        )
    if mol.GetRingInfo().NumRings() > 0:
        raise UnsupportedStructure("a ring is out of scope for this acyclic-only module")

    for atom in mol.GetAtoms():
        if atom.GetIdx() == anion.GetIdx():
            continue
        if atom.GetAtomicNum() != 6:
            raise UnsupportedStructure(
                "heteroatoms and halogens are out of scope for this pilot "
                "carbanide module (P-72.2.2.1)"
            )
        if atom.GetFormalCharge() != 0 or atom.GetIsotope() != 0 or atom.GetIsAromatic():
            raise UnsupportedStructure(
                "a charged, isotopically modified, or aromatic carbon "
                "alongside the carbanion center is not supported yet"
            )
    if len(Chem.GetMolFrags(mol)) > 1:
        raise UnsupportedStructure("multi-fragment structures are not supported yet")
    for bond in mol.GetBonds():
        if bond.GetBondTypeAsDouble() != 1.0:
            raise UnsupportedStructure(
                "unsaturation is out of scope for this pilot carbanide module"
            )
    return anion


def _candidate_key(chain_length, locant, substituents):
    grouped = group_substituents(substituents)
    locant_set, total_count, citation_locants = substituent_locant_set_and_citation(grouped)
    prefix = format_substituent_prefixes(grouped, omit_locants=chain_length == 1)
    name = prefix + name_from_substituents(chain_length, [], [], "ide", [locant])
    return (locant, -total_count, locant_set, citation_locants, name), name


def name_carbanide(mol) -> str:
    anion = _validate_and_find_carbanide(mol)
    anion_idx = anion.GetIdx()
    graph = adjacency(mol)
    chains = longest_chains(carbon_adjacency(mol))
    chain_length = len(chains[0])

    eligible = [chain for chain in chains if anion_idx in chain]
    if not eligible:
        raise UnsupportedStructure(
            "the carbanion center does not lie on a single longest carbon "
            "chain"
        )

    best_key = None
    best_name = None
    for chain in eligible:
        for candidate in (chain, list(reversed(chain))):
            position_of = {atom: i + 1 for i, atom in enumerate(candidate)}
            locant = position_of[anion_idx]
            substituents = substituents_for_chain(graph, candidate, {}, mol=mol)
            key, name = _candidate_key(chain_length, locant, substituents)
            if best_key is None or key < best_key:
                best_key, best_name = key, name
    return best_name
