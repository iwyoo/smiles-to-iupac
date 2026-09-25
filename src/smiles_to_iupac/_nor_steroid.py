"""Naming of the 'nor' skeletal-modification prefix (P-101.3.1,
https://iupac.qmul.ac.uk/BlueBook/PDF/P10.pdf) on the bare steroid parent
hydrides already recognized by `_steroid_parent_hydrides.py`: a single
skeletal atom -- a ring atom at a non-fusion ("atomic connector") position,
or one of the two angular methyls (C18/C19) -- removed, with its attached
hydrogens, from an otherwise-unmodified parent skeleton. Named
`<locant>-nor-<parent>` where `<locant>` is that atom's position in the
standard steroid numbering (Table 10.1: ring positions 1-17, angular
methyls 18/19).

Only non-fusion ring positions are handled here, not the six true
ring-fusion carbons (5, 8, 9, 10, 13, 14) or a side-chain-bearing C17:
removing a fusion atom merges two rings into one rather than simply
shrinking the ring system by one atom, a structurally different case the
Blue Book's own worked examples (all at plain "atomic connector"
positions, e.g. C4 in `4-nor-5beta-pregnane`) don't cover -- out of scope
here. This is detected automatically rather than via a hardcoded locant
list: a steroid parent's own ring/methyl atom is eligible for `nor` only
if its degree in that parent's full molecular graph (side chain included)
is at most 2, which is exactly the non-fusion, non-branched condition.

Reference locant numbering is established via substructure matching
against a hand-built query graph (not the existing dict's own SMILES atom
order, which isn't guaranteed to follow locant order): one query encodes
the bare 17-atom ring skeleton (locants 1-17), a second adds the C13
methyl (18), a third adds both angular methyls (18, 19) -- matched against
whichever of gonane/estrane/the androstane-family parents applies. Each
query was confirmed (see module's own test) to match its target with
exactly one substructure match (no other skeletal automorphism), so the
locant assignment is unambiguous.

Some single-atom removals coincide with an already-listed retained-name
parent rather than a genuine 'nor' case -- confirmed empirically:
removing androstane's C19 methyl reproduces estrane exactly, removing
estrane's remaining C18 methyl reproduces gonane, and removing ergostane's
extra side-chain methyl reproduces cholestane. These are excluded from
the lookup below (not asserted as errors): `core.py` dispatches to
`_steroid_parent_hydrides.py`'s own exact-parent match before this
module, so a real occurrence of one of these structures is already named
correctly as the retained name, and this module never needs to (and must
not) re-claim it under a `nor` name.
"""

from rdkit import Chem
from rdkit.Chem import RWMol

from ._steroid_parent_hydrides import _CANONICAL_TO_NAME, _PARENT_HYDRIDES, _RAW_PARENTS, _locant_map


def _nor_smiles(name, smiles, locant, atom_idx):
    mol = Chem.MolFromSmiles(smiles)
    if mol.GetAtomWithIdx(atom_idx).GetDegree() > 2:
        return None
    rw = RWMol(mol)
    rw.RemoveAtom(atom_idx)
    try:
        Chem.SanitizeMol(rw)
    except (Chem.rdchem.KekulizeException, Chem.rdchem.AtomValenceException):
        return None
    return Chem.MolToSmiles(rw)


def _build_nor_lookup():
    lookup = {}
    for name, smiles in _RAW_PARENTS.items():
        mol = Chem.MolFromSmiles(smiles)
        locants = _locant_map(name, mol)
        for locant, atom_idx in locants.items():
            nor_smiles = _nor_smiles(name, smiles, locant, atom_idx)
            if nor_smiles is None or nor_smiles in _CANONICAL_TO_NAME:
                continue
            if nor_smiles in lookup:
                assert lookup[nor_smiles] == (locant, name), (
                    f"nor-{name} at {locant} collides with an unrelated nor entry"
                )
                continue
            lookup[nor_smiles] = (locant, name)
    return lookup


_NOR_LOOKUP = _build_nor_lookup()


def has_nor_steroid_shape(mol) -> bool:
    return Chem.MolToSmiles(mol) in _NOR_LOOKUP


def name_nor_steroid(mol) -> str:
    locant, parent = _NOR_LOOKUP[Chem.MolToSmiles(mol)]
    return f"{locant}-nor-{parent}"
