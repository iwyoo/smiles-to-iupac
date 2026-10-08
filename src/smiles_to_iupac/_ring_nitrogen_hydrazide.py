"""An unsubstituted hydrazide group on a ring nitrogen (P-66.3.1.1): 'piperidine-1-carbohydrazide'. The group is named as
the carboxamide it ranks beside, then the ending is exchanged."""

from rdkit import Chem

from ._common import UnsupportedStructure

_RING_N_HYDRAZIDE = Chem.MolFromSmarts("[#7;R;X3;+0]-[CX3](=O)-[NX3;H1]-[NX3;H2]")


def has_ring_nitrogen_hydrazide_shape(mol) -> bool:
    return len(mol.GetSubstructMatches(_RING_N_HYDRAZIDE)) == 1


def name_ring_nitrogen_hydrazide(mol) -> str:
    if not has_ring_nitrogen_hydrazide_shape(mol):
        raise UnsupportedStructure("not a hydrazide group on a ring nitrogen")
    terminal = mol.GetSubstructMatch(_RING_N_HYDRAZIDE)[-1]
    healed = Chem.RWMol(mol)
    healed.RemoveAtom(terminal)
    healed = healed.GetMol()
    Chem.SanitizeMol(healed)
    from .core import _smiles_to_iupac_unabridged

    name = _smiles_to_iupac_unabridged(Chem.MolToSmiles(healed))
    if not name.endswith("-carboxamide") or name.count("carboxamide") != 1:
        raise UnsupportedStructure("the amide analogue is not named by a carboxamide suffix")
    return name[: -len("carboxamide")] + "carbohydrazide"
