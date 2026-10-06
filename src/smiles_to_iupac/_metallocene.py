"""Naming of the 7 retained metallocene ("-ocene") names (P-69.2.7, Chapter
P-6a, the Blue Book): "bis(eta5-cyclopenta-2,4-
dien-1-yl)metal" sandwiches of Fe/Ru/Os/Ni/Cr/Co/V are retained PINs
(ferrocene/ruthenocene/osmocene/nickelocene/chromocene/cobaltocene/
vanadocene), exempted from the general P-69.2 coordination-nomenclature
machinery.

Real structures use either of two equally valid charge conventions for
the same conceptual sandwich -- a neutral metal with two cyclopentadienyl
biradical rings, or a +2 metal cation with two cyclopentadienide anion
rings -- so both are matched, mirroring `_steroid_parent_hydrides.py`'s/
`_fullerene.py`'s/`_nucleoside.py`'s identical retained-name-by-exact-
structure pattern. No substituent support, no hapticity citation, no
general coordination mechanism: any substituted ring, wrong metal, or
wrong ring count is correctly left unmatched.
"""

from rdkit import Chem

_METAL_SYMBOLS = {
    "ferrocene": "Fe",
    "ruthenocene": "Ru",
    "osmocene": "Os",
    "nickelocene": "Ni",
    "chromocene": "Cr",
    "cobaltocene": "Co",
    "vanadocene": "V",
}

_CANONICAL_TO_NAME = {}
for _name, _symbol in _METAL_SYMBOLS.items():
    _neutral = f"C1=C[CH]C=C1.C1=C[CH]C=C1.[{_symbol}]"
    _anionic = f"[CH-]1C=CC=C1.[CH-]1C=CC=C1.[{_symbol}+2]"
    _CANONICAL_TO_NAME[Chem.CanonSmiles(_neutral)] = _name
    _CANONICAL_TO_NAME[Chem.CanonSmiles(_anionic)] = _name
assert len(_CANONICAL_TO_NAME) == 2 * len(_METAL_SYMBOLS), (
    "two entries above canonicalized to the same key -- a real name "
    "collision, not just a duplicate row"
)


def has_metallocene_name(mol) -> bool:
    return Chem.MolToSmiles(mol) in _CANONICAL_TO_NAME


def name_metallocene(mol) -> str:
    return _CANONICAL_TO_NAME[Chem.MolToSmiles(mol)]
