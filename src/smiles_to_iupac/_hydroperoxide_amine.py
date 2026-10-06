"""Naming of a molecule combining a hydroperoxide (-OOH) with exactly one
separate primary amine (-NH2) on the same acyclic saturated carbon chain,
per the IUPAC 2013 Recommendations ("the Blue Book"):

- P-41 Table 4.1 (`_seniority.py`, the Blue Book): hydroperoxides
  (class 18) outrank amines (class 19), so a coexisting primary amine is
  demoted to the 'amino' substituent prefix instead of its own '-amine'
  suffix, e.g. 'NCCOO' -> '2-aminoethane-1-peroxol'.
- **Note on PubChem divergence**: PubChem's own auto-generated IUPACName
  for 'NCCOO' is "2-hydroperoxyethanamine" (amine as the suffix parent,
  hydroperoxide demoted) -- the *opposite* seniority conclusion from
  Table 4.1's explicit text above. PubChem's automated namer is known
  elsewhere in this project to deviate from strict IUPAC rules (see the
  natural-products scoping notes); this module follows the Blue Book
  primary source directly rather than PubChem's tool output for this one
  pair, per a re-investigation of the primary text.
- Otherwise mirrors `_hydroperoxide.py`'s acyclic path exactly: the
  -OOH-bearing carbon gets the lowest available locant ahead of the
  amine's own 'amino' prefix locant (P-44.1.1), and 'peroxol' never elides
  the parent hydride's final 'e' (P-16.3.3, consonant suffix) the way
  '-ol' does.

Scope, deliberately narrow (mirrors `_alcohol_amine.py`/`_thiol_amine.py`):
a single -OOH group plus a single primary amine, both on one acyclic
*saturated* chain, with halogen substituents allowed.
Explicitly out of scope (raise `UnsupportedStructure`): any chain
unsaturation (ene/yne), any specified stereocenter, any ring, more than
one hydroperoxide or amine, a secondary/tertiary amine, and a coexisting
standalone hydroxyl/ether/other heteroatom.
"""


from ._coexisting_groups import name_via_senior_acyclic
from ._common import (
    UnsupportedStructure,
    adjacency,
    find_primary_amines,
    non_single_bonds,
    specified_stereocenters,
    validate_allowed_atoms,
)
from ._hydroperoxide import _hydroperoxide_oxygens, _name_acyclic_hydroperoxide


def has_hydroperoxide_amine_shape(mol) -> bool:
    if mol.GetRingInfo().NumRings() > 0:
        return False
    oxygens = _hydroperoxide_oxygens(mol)
    if oxygens is None:
        return False
    attach, terminal = oxygens
    amines = find_primary_amines(mol, {attach.GetIdx(), terminal.GetIdx()})
    return len(amines) == 1


def _validate(mol, exclude, amine_nitrogen):
    validate_allowed_atoms(
        mol,
        "heteroatoms other than a hydroperoxide's own two oxygens "
        "(P-56.1), a primary amine nitrogen (P-41, Table 4.1), and "
        "halogen substituents (P-35.2.1) are not supported yet",
        [
            (
                8,
                exclude,
                "an oxygen that isn't part of the single hydroperoxide's "
                "-O-O-H pair is out of scope for this module (e.g. a "
                "coexisting hydroxyl or ether)",
            ),
            (
                7,
                {amine_nitrogen},
                "a nitrogen that isn't a single plain primary amine (-NH2) "
                "is out of scope for this module (secondary/tertiary "
                "amines, imines, and nitriles are not supported)",
            ),
        ],
    )


def name_hydroperoxide_amine(mol) -> str:
    if mol.GetRingInfo().NumRings() > 0:
        raise UnsupportedStructure(
            "a hydroperoxide/amine combination on/in a ring uses a "
            "different naming construction, out of scope for this "
            "acyclic-only module"
        )
    oxygens = _hydroperoxide_oxygens(mol)
    if oxygens is None:
        raise UnsupportedStructure(
            "not a plain hydroperoxide (-OOH, P-56.1) shape coexisting "
            "with a primary amine"
        )
    attach, terminal = oxygens
    exclude = {attach.GetIdx(), terminal.GetIdx()}
    (site,) = (n for n in attach.GetNeighbors() if n.GetIdx() != terminal.GetIdx())
    site_idx = site.GetIdx()

    amines = find_primary_amines(mol, exclude)
    if len(amines) != 1:
        raise UnsupportedStructure(
            "exactly one primary amine (-NH2) coexisting with the single "
            "hydroperoxide is supported here (see _amine.py for a plain "
            "amine)"
        )
    (amine_nitrogen,) = amines
    _validate(mol, exclude, amine_nitrogen)

    all_non_single = [b for b in non_single_bonds(mol) if b[0] not in exclude and b[1] not in exclude]
    if all_non_single:
        raise UnsupportedStructure(
            "chain unsaturation (ene/yne) combined with a hydroperoxide/"
            "amine seniority demotion is out of scope for this module "
            "(1st-pass scope: saturated only)"
        )
    if specified_stereocenters(mol):
        raise UnsupportedStructure(
            "a specified stereocenter alongside a hydroperoxide/amine "
            "seniority demotion is not supported yet"
        )

    graph = adjacency(mol)
    (amine_carbon,) = graph[amine_nitrogen]
    return name_via_senior_acyclic(
        _name_acyclic_hydroperoxide,
        "hydroperoxide",
        "amine",
        (mol, site_idx, exclude),
        {amine_nitrogen: "amino"},
        required_atoms={amine_carbon},
    )
