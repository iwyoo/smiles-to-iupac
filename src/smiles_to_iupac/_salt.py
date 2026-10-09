"""Fragment-level ion naming shared by the salt writer in `_acid_salts.py` (P-72, P-73, P-77): fixed cation names,
the retained polyatomic anions and the organic anion namers."""

from rdkit import Chem

from ._alkoxide import has_alkoxide_shape, name_alkoxide
from ._ammonium import has_ammonium_shape, name_ammonium
from ._anion import has_general_anion_shape, name_anion
from ._carbanide import name_carbanide
from ._carboxylate import has_carboxylate_shape, name_carboxylate
from ._common import HALOGEN_PREFIXES, UnsupportedStructure
from ._selenoate import has_selenoate_shape, name_selenoate
from ._thioate import has_thioate_shape, name_thioate


def _has_halide_anion_shape(frag):
    if frag.GetNumAtoms() != 1:
        return False
    atom = frag.GetAtomWithIdx(0)
    return atom.GetAtomicNum() in HALOGEN_PREFIXES and atom.GetFormalCharge() == -1


def _name_halide_anion(frag) -> str:
    prefix = HALOGEN_PREFIXES[frag.GetAtomWithIdx(0).GetAtomicNum()]
    return prefix[:-1] + "ide"


def _has_hydroxide_anion_shape(frag):
    return frag.GetNumAtoms() == 1 and frag.GetAtomWithIdx(0).GetAtomicNum() == 8 and frag.GetAtomWithIdx(0).GetFormalCharge() == -1


def _name_hydroxide_anion(frag) -> str:
    return "hydroxide"


def _has_carbanide_anion_shape(frag):
    if sum(a.GetFormalCharge() for a in frag.GetAtoms()) != -1:
        return False
    try:
        name_carbanide(frag)
    except UnsupportedStructure:
        return False
    return True


_ANION_KINDS = [
    (has_carboxylate_shape, name_carboxylate),
    (has_alkoxide_shape, name_alkoxide),
    (has_thioate_shape, name_thioate),
    (has_selenoate_shape, name_selenoate),
    (_has_halide_anion_shape, _name_halide_anion),
    (_has_hydroxide_anion_shape, _name_hydroxide_anion),
    (_has_carbanide_anion_shape, name_carbanide),
    (has_general_anion_shape, name_anion),
]


def _is_bare_polyatomic_ion(frag, center_atomic_num, num_oxygens, total_charge):
    """True iff `frag` is a single central atom (`center_atomic_num`, formal
    charge 0) bonded to exactly `num_oxygens` monovalent oxygens and
    nothing else, with the fragment's combined formal charge exactly
    `total_charge` -- accepts any resonance depiction of the charge split
    across the oxygens (e.g. sulfate's two anionic + two neutral
    doubly-bonded oxygens), since only the stoichiometry and net charge
    matter for naming these fixed-formula ions (P-65.6.2)."""
    if frag.GetNumAtoms() != 1 + num_oxygens:
        return False
    centers = [a for a in frag.GetAtoms() if a.GetAtomicNum() == center_atomic_num]
    if len(centers) != 1:
        return False
    (center,) = centers
    # The central atom's own formal charge varies by resonance depiction
    # (e.g. nitrate's N is routinely written charge +1, '[O-][N+](=O)[O-]',
    # to satisfy valence -- unlike carbonate's neutral-C depiction) --
    # only the fragment's *combined* charge below is checked, not any
    # individual atom's.
    if center.GetIsotope() != 0 or center.GetDegree() != num_oxygens:
        return False
    others = [a for a in frag.GetAtoms() if a.GetIdx() != center.GetIdx()]
    if any(a.GetAtomicNum() != 8 or a.GetDegree() != 1 or a.GetIsotope() != 0 for a in others):
        return False
    return sum(a.GetFormalCharge() for a in frag.GetAtoms()) == total_charge


def has_sulfate_shape(frag) -> bool:
    return _is_bare_polyatomic_ion(frag, 16, 4, -2)


def name_sulfate(frag) -> str:
    return "sulfate"


def has_carbonate_shape(frag) -> bool:
    return _is_bare_polyatomic_ion(frag, 6, 3, -2)


def name_carbonate(frag) -> str:
    return "carbonate"


def has_nitrate_shape(frag) -> bool:
    return _is_bare_polyatomic_ion(frag, 7, 3, -1)


def name_nitrate(frag) -> str:
    return "nitrate"


def has_phosphate_shape(frag) -> bool:
    return _is_bare_polyatomic_ion(frag, 15, 4, -3)


def name_phosphate(frag) -> str:
    return "phosphate"


_POLYATOMIC_ANION_KINDS = [
    (has_sulfate_shape, name_sulfate),
    (has_carbonate_shape, name_carbonate),
    (has_nitrate_shape, name_nitrate),
    (has_phosphate_shape, name_phosphate),
    (has_general_anion_shape, name_anion),
]


def _polyatomic_anion(frag):
    """(namer, charge_magnitude) for `frag` if it matches one of the fixed-
    formula polyatomic anions above, else None."""
    for has_shape, namer in _POLYATOMIC_ANION_KINDS:
        if has_shape(frag):
            magnitude = -sum(atom.GetFormalCharge() for atom in frag.GetAtoms())
            return namer, magnitude
    return None

_SINGLY_CHARGED_CATION_ONLY_ANIONS = {name_alkoxide, name_thioate, name_selenoate, name_carbanide, name_anion}

_MONOATOMIC_CATION_NAMES = {
    ("Li", 1): "lithium",
    ("Na", 1): "sodium",
    ("K", 1): "potassium",
    ("Rb", 1): "rubidium",
    ("Cs", 1): "cesium",
    ("Be", 2): "beryllium",
    ("Mg", 2): "magnesium",
    ("Ca", 2): "calcium",
    ("Sr", 2): "strontium",
    ("Ba", 2): "barium",
    ("Al", 3): "aluminium",
}


def _monoatomic_cation(frag):
    if frag.GetNumAtoms() != 1:
        return None
    atom = frag.GetAtomWithIdx(0)
    charge = atom.GetFormalCharge()
    name = _MONOATOMIC_CATION_NAMES.get((atom.GetSymbol(), charge))
    if name is None:
        return None
    return name, charge


def _cation(frag):
    cation = _monoatomic_cation(frag)
    if cation is not None:
        return cation
    if has_ammonium_shape(frag):
        # `has_ammonium_shape` only looks at the charged nitrogen itself,
        # ignoring the rest of the fragment -- when called (as below) on a
        # would-be *single*-fragment "cation" that is actually a whole
        # zwitterion, `name_ammonium` then tries to name the entire
        # fragment as if it were just the ammonium compound and fails
        # deeper in (e.g. on a coexisting carboxylate oxygen). That failure
        # means this fragment was never a real bare cation, not that
        # `_split_cation_anions` itself should crash -- `has_zwitterion_shape`
        # (see `_zwitterion.py`) is routed ahead of this module for the
        # zwitterion case specifically, but this catch stays regardless as
        # a defensive fallback for any other shape that superficially
        # matches on the nitrogen alone.
        try:
            return name_ammonium(frag), 1
        except UnsupportedStructure:
            pass
    return _organic_cation(frag)


def _organic_cation(frag):
    """Name and charge of an organic cation whose positive centres all carry one charge, named as a whole by the
    substitutive engine (an '-ium' name, P-73.1, P-73.5)."""
    from .core import smiles_to_iupac

    charges = [a.GetFormalCharge() for a in frag.GetAtoms() if a.GetFormalCharge()]
    if not charges or set(charges) != {1}:
        return None
    try:
        name = smiles_to_iupac(Chem.MolToSmiles(frag))
    except UnsupportedStructure:
        return None
    return (name, len(charges)) if name.endswith(("ium", "ium)", "ium(1+)", "ium(2+)")) else None
