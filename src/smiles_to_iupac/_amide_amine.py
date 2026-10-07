"""Naming of a molecule combining a primary amide (-CONH2, unsubstituted)
with exactly one separate primary amine (-NH2) on the same acyclic
saturated carbon chain, per the IUPAC 2013 Recommendations ("the Blue
Book"):

- P-41 (`_seniority.py`): 'amide' (`_amide.py`, Table 4.4 class 16) far
  outranks 'amine' (class 51), so a coexisting primary amine is demoted to
  the 'amino' substituent prefix instead of its own '-amine' suffix, e.g.
  'NCCC(N)=O' -> '3-aminopropanamide' (PubChem's own IUPACName, an exact
  match -- no retained-name divergence for this class).
  Reuses `_amide.py`'s own `_name_acyclic_amide` chain-search/numbering
  function directly (via `_coexisting_groups.py`), the same pattern
  `_ketone_amine.py`/`_aldehyde_amine.py` used for their own senior
  modules.
- Otherwise mirrors `_amide.py`'s acyclic path exactly: the amide carbon is
  always chain C1 (P-14.3.3, no locant), and the amine nitrogen's own
  locant is always cited via the 'amino' prefix.

Scope, deliberately narrow (mirrors `_ketone_amine.py`/`_aldehyde_amine.py`):
a single unsubstituted primary amide (-CONH2, no N-alkyl) plus a single
primary amine, both on one acyclic *saturated* chain, with halogen
substituents allowed.
Explicitly out of scope (raise `UnsupportedStructure`): any N-alkyl amide
substitution, any chain unsaturation (ene/yne), any specified
stereocenter, any ring, more than one amide or amine, a secondary/tertiary
amine, a coexisting standalone hydroxyl/other heteroatom, and any
amide/amine not captured by a single longest chain. A benzene-ring-
substituent chain variant (mirroring `_amide.py`'s
`_name_phenyl_chain_amide`) is a separate follow-up.
"""


from ._amide import _is_carbonyl_carbon, _name_acyclic_amide
from ._coexisting_groups import name_via_senior_acyclic
from ._common import (
    UnsupportedStructure,
    adjacency,
    find_primary_amines,
    non_single_bonds,
    specified_stereocenters,
    validate_allowed_atoms,
)



def _find_amide(mol):
    """Return (amide_carbon, amide_oxygen, amide_nitrogen) for a single
    unsubstituted primary amide (-CONH2), or None if there isn't exactly
    one such group. N-alkyl-substituted amides are deliberately excluded
    here (out of scope, see module docstring) rather than partially
    handled."""
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
                oxygens = [o for o in atom.GetNeighbors() if o.GetAtomicNum() == 8]
                if len(oxygens) != 1:
                    continue
                if found is not None:
                    return None
                found = (atom.GetIdx(), oxygens[0].GetIdx(), n.GetIdx())
    return found


def has_amide_amine_shape(mol) -> bool:
    amide = _find_amide(mol)
    if amide is None:
        return False
    return len(find_primary_amines(mol, {amide[2]})) == 1


def _validate(mol, amide_atoms, amines):
    amide_carbon, amide_oxygen, amide_nitrogen = amide_atoms
    validate_allowed_atoms(
        mol,
        "heteroatoms other than a primary amide's oxygen/nitrogen "
        "(P-66.1), a primary amine nitrogen (P-41, Table 3.3), and "
        "halogen substituents (P-35.2.1) are not supported yet",
        [
            (
                8,
                {amide_oxygen},
                "an oxygen that isn't the single amide carbonyl is out of "
                "scope for this module (e.g. a coexisting hydroxyl, ether, "
                "or a second carbonyl)",
            ),
            (
                7,
                {amide_nitrogen} | amines,
                "a nitrogen that isn't the amide's own -CONH2 nitrogen or "
                "a plain primary amine (-NH2) is out of scope for this "
                "module (N-alkyl amides, secondary/tertiary amines, "
                "imines, and nitriles are not supported)",
            ),
        ],
    )


def name_amide_amine(mol) -> str:
    if mol.GetRingInfo().NumRings() > 0:
        raise UnsupportedStructure(
            "an amide/amine combination on/in a ring uses a different "
            "naming construction, out of scope for this acyclic-only module"
        )
    amide = _find_amide(mol)
    if amide is None:
        raise UnsupportedStructure(
            "exactly one unsubstituted primary amide (-CONH2) coexisting "
            "with a primary amine is supported here (see _amide.py for a "
            "plain amide)"
        )
    amide_carbon, amide_oxygen, amide_nitrogen = amide
    amines = find_primary_amines(mol, {amide_nitrogen})
    if len(amines) != 1:
        raise UnsupportedStructure(
            "exactly one primary amine (-NH2) coexisting with the single "
            "amide is supported here (see _amine.py for a plain amine)"
        )
    _validate(mol, amide, amines)

    excluded = {amide_oxygen, amide_nitrogen}
    all_non_single = [b for b in non_single_bonds(mol) if b[0] not in excluded and b[1] not in excluded]
    if all_non_single:
        raise UnsupportedStructure(
            "chain unsaturation (ene/yne) combined with an amide/amine "
            "seniority demotion is out of scope for this module (1st-pass "
            "scope: saturated only)"
        )
    if specified_stereocenters(mol):
        raise UnsupportedStructure(
            "a specified stereocenter alongside an amide/amine seniority "
            "demotion is not supported yet"
        )

    graph = adjacency(mol)
    (amine_nitrogen,) = amines
    (amine_carbon,) = graph[amine_nitrogen]
    return name_via_senior_acyclic(
        _name_acyclic_amide,
        "amide",
        "amine",
        (mol, amide_carbon, amide_nitrogen, excluded, (), set(), []),
        {amine_nitrogen: "amino"},
        required_atoms={amine_carbon},
    )
