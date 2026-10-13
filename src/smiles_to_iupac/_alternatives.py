"""General-nomenclature names that the Blue Book allows next to the preferred IUPAC name (P-15.2.2, P-15.6,
P-72.4, P-73.4). `smiles_to_iupac` returns the PIN; this module returns the alternatives, each labelled by kind."""

import warnings
from typing import NamedTuple

from rdkit import Chem, rdBase

from ._functional_modifier import functional_modifier_names
from ._ionic_prefix import ionic_prefix_names
from .core import _parse_smiles, smiles_to_iupac


class AlternativeName(NamedTuple):
    kind: str
    name: str
    rule: str


_GENERATORS = (("functional-modifier", functional_modifier_names), ("ionic-prefix", ionic_prefix_names))


def smiles_to_alternative_names(smiles: str) -> list[AlternativeName]:
    if not isinstance(smiles, str):
        raise TypeError(f"smiles must be a str, not {type(smiles).__name__}")
    legacy = Chem.GetUseLegacyStereoPerception()
    Chem.SetUseLegacyStereoPerception(False)
    try:
        with rdBase.BlockLogs(), warnings.catch_warnings():
            warnings.simplefilter("ignore")
            mol = _parse_smiles(smiles) if smiles.strip() else None
            if mol is None:
                raise ValueError(f"invalid SMILES: {smiles!r}")
            found = []
            for kind, generator in _GENERATORS:
                try:
                    found += [AlternativeName(kind, name, rule) for name, rule in generator(mol, smiles_to_iupac)]
                except NotImplementedError:
                    continue
            return found
    finally:
        Chem.SetUseLegacyStereoPerception(legacy)
