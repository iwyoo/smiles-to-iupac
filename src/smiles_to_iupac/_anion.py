"""Anions on characteristic groups (P-72.2.2.2.1-.3, P-72.7): the anionic
group is always the parent suffix and every neutral group is cited as a prefix.

Each anionic chalcogen/nitrogen is protonated and marked, `_polyfunctional`
names the neutral molecule with only the marked groups eligible as the
principal group, and the neutral suffix is then swapped for its anion form.
"""

import contextlib
import re

from rdkit import Chem

from ._common import UnsupportedStructure, specified_stereo_elements

ANION_PROP = "_anion"

_CHALCOGENS = {8, 16, 34, 52}
_MULTIPLIED = {"di": "bis", "tri": "tris", "tetra": "tetrakis", "penta": "pentakis", "hexa": "hexakis"}
_COUNT = "(di|tri|tetra|penta|hexa)?"

_SWAPS = (
    (re.compile(r"phenol$"), lambda m: "phenoxide"),
    (re.compile(r"imine$"), lambda m: "iminide"),
    (re.compile(r"imidic acid$"), lambda m: "imidate"),
    (re.compile(r"aniline$"), lambda m: "benzenaminide"),
    (re.compile(r"(di|tri|tetra|penta|hexa)(ol|thiol|amine)$"), lambda m: _bis(m.group(1), m.group(2))),
    (re.compile(r"thiol$"), lambda m: "thiolate"),
    (re.compile(r"ol$"), lambda m: "olate"),
    (re.compile(r"amine$"), lambda m: "aminide"),
    (re.compile(r"oic acid$"), lambda m: "oate"),
    (re.compile(r"carboxylic acid$"), lambda m: "carboxylate"),
    (re.compile(r"sulfonic acid$"), lambda m: "sulfonate"),
)


def _bis(count, suffix):
    word = {"ol": "olate", "thiol": "thiolate", "amine": "aminide"}[suffix]
    return f"{_MULTIPLIED[count]}({word})"


def anion_weight(atom):
    return int(atom.GetProp(ANION_PROP))


def swap_suffix(name):
    for pattern, replace in _SWAPS:
        match = pattern.search(name)
        if match:
            return name[: match.start()] + replace(match)
    from ._acid_derivatives import anion_name

    try:
        return anion_name(name)
    except UnsupportedStructure:
        raise UnsupportedStructure("this anionic group has no anion suffix form yet") from None


def _is_nitro_oxygen(atom):
    return atom.GetAtomicNum() != 6 and any(n.GetFormalCharge() > 0 for n in atom.GetNeighbors())


def anion_atoms(mol):
    return [a for a in mol.GetAtoms() if a.GetFormalCharge() < 0 and not _is_nitro_oxygen(a)]


def _is_carbanion(atom):
    return atom.GetAtomicNum() == 6 and atom.GetFormalCharge() in (-1, -2) and not atom.GetNumRadicalElectrons()


def _is_acyl_like(carbon):
    return any(
        b.GetBondTypeAsDouble() >= 2.0 and b.GetOtherAtom(carbon).GetAtomicNum() in (7, 8, 16, 34, 52)
        for b in carbon.GetBonds()
    )


def _is_peroxy_host(host, atom):
    return (
        host.GetAtomicNum() in _CHALCOGENS
        and host.GetDegree() == 2
        and any(n.GetAtomicNum() == 6 for n in host.GetNeighbors() if n.GetIdx() != atom.GetIdx())
    )


def _demoted_peroxy(atom, centers):
    """A peroxolate-type center beside a more senior anionic group is cited as a prefix (P-72.7(e))."""
    if not (_is_group_anion(atom) and atom.GetAtomicNum() in _CHALCOGENS and _is_peroxy_host(atom.GetNeighbors()[0], atom)):
        return False
    return any(
        c is not atom
        and (_is_group_anion(c) or _is_carbanion(c))
        and not _is_junior_to_peroxol(c)
        and not _is_peroxy_anion(c)
        for c in centers
    )


def _is_peroxy_anion(atom):
    return atom.GetAtomicNum() in _CHALCOGENS and atom.GetDegree() == 1 and _is_peroxy_host(atom.GetNeighbors()[0], atom)


def _is_junior_to_peroxol(atom):
    if atom.GetAtomicNum() == 7:
        return True
    if atom.GetAtomicNum() in (16, 34, 52) and atom.GetDegree() == 1:
        host = atom.GetNeighbors()[0]
        return host.GetAtomicNum() == 6 and not any(
            b.GetBondTypeAsDouble() == 2.0 and b.GetOtherAtom(host).GetAtomicNum() in (7, 8, 16, 34, 52) for b in host.GetBonds()
        )
    return False


def _is_secondary_aminide(atom):
    return (
        atom.GetAtomicNum() == 7
        and atom.GetFormalCharge() == -1
        and atom.GetDegree() == 2
        and atom.GetTotalNumHs() == 0
        and not atom.IsInRing()
        and not atom.GetNumRadicalElectrons()
        and all(
            n.GetAtomicNum() == 6 and not _is_acyl_like(n) and atom.GetOwningMol().GetBondBetweenAtoms(atom.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 1.0
            for n in atom.GetNeighbors()
        )
    )


def _is_group_anion(atom):
    if _is_secondary_aminide(atom):
        return True
    charge = -atom.GetFormalCharge()
    if charge not in (1, 2) or atom.GetDegree() != 1 or atom.GetNumRadicalElectrons():
        return False
    z, hydrogens = atom.GetAtomicNum(), atom.GetTotalNumHs()
    host = atom.GetNeighbors()[0]
    order = atom.GetBonds()[0].GetBondTypeAsDouble()
    if z in _CHALCOGENS:
        return charge == 1 and order == 1.0 and hydrogens == 0 and (host.GetAtomicNum() == 6 or (host.GetAtomicNum() == 16 and host.GetDegree() >= 3)) or _is_peroxy_host(host, atom)
    if z == 7 and host.GetAtomicNum() == 6:
        if order == 2.0:
            return charge == 1 and hydrogens == 0
        return order == 1.0 and not _is_acyl_like(host) and hydrogens == 2 - charge
    return False


def _is_ammonium_prefix(atom):
    """An acyclic ammonium nitrogen that an alcoholate or thiolate parent cites as an 'azaniumyl' prefix (P-74.1.3)."""
    return (
        atom.GetAtomicNum() == 7
        and atom.GetFormalCharge() == 1
        and not atom.IsInRing()
        and all(n.GetAtomicNum() == 6 for n in atom.GetNeighbors())
    )


def marked_neutral(mol):
    from ._anion_center import center_kind, is_center_atom

    if len(Chem.GetMolFrags(mol)) != 1:
        raise UnsupportedStructure("a multi-fragment anionic structure is not supported here")
    centers = anion_atoms(mol)
    cations = [a for a in mol.GetAtoms() if a.GetFormalCharge() > 0 and not _is_nitro_nitrogen(a)]
    prefix_cations = [a for a in cations if _is_ammonium_prefix(a)] if centers and all(_is_group_anion(a) for a in centers) else []
    cations = [a for a in cations if a not in prefix_cations]
    if any(not (a.IsInRing() and a.GetFormalCharge() == 1) for a in cations) or len(cations) > 1:
        raise UnsupportedStructure("cationic centers beside an anionic group are not supported here")
    others = [a for a in centers if not (_is_group_anion(a) or _is_carbanion(a)) or _demoted_peroxy(a, centers)]
    if not centers or not all(is_center_atom(a) for a in others):
        raise UnsupportedStructure("this anionic center is not on a supported characteristic group")
    editable = Chem.RWMol(mol)
    for center in centers:
        atom = editable.GetAtomWithIdx(center.GetIdx())
        charge = -atom.GetFormalCharge()
        hydrogens = atom.GetTotalNumHs()
        atom.SetFormalCharge(0)
        if center in others:
            word, lam = center_kind(center, preferred="ide" if center.IsInRing() else "uide")
            atom.SetNoImplicit(True)
            atom.SetNumExplicitHs(hydrogens)
            atom.SetProp("_anion_word", word)
            atom.SetProp("_anion_charge", str(charge))
            if lam:
                atom.SetProp("_anion_lambda", str(lam))
        else:
            atom.SetNumExplicitHs(hydrogens + charge)
            atom.SetNoImplicit(True)
            atom.SetProp(ANION_PROP, str(charge))
    neutral = editable.GetMol()
    neutral.UpdatePropertyCache(strict=False)
    Chem.SanitizeMol(neutral, Chem.SANITIZE_ALL ^ Chem.SANITIZE_PROPERTIES)
    for atom in neutral.GetAtoms():
        if atom.HasProp("_anion_word"):
            atom.SetNumRadicalElectrons(0)
    return neutral


def _is_nitro_nitrogen(atom):
    return atom.GetAtomicNum() == 7 and any(n.GetFormalCharge() < 0 for n in atom.GetNeighbors())


_CARBANIDE_PARENT = re.compile(r"-\d+(?:,\d+)*-(?:di|tri|tetra)?ide$")
_THIOIC_TAIL = re.compile(r"((?:di|tri|tetra)?(?:carbo)?)thioic acid$")
_ANION_TOKEN = re.compile(r"(?:ide|ate|ite)(?![a-z])|ato|ido|idyl|uid-|id-|ide-")
_MULTIPLE_ANION = re.compile(r"(?:di|tri|tetra|bis|tris|tetrakis)\(?[A-Za-z-]*(?:ide|uide|ate|ite)")


def name_anion(mol):
    name = _name_anion_unchecked(mol)
    centers = len(anion_atoms(mol))
    tokens = len(_ANION_TOKEN.findall(name))
    from ._anion_oxoacid import oxoacid_center

    if oxoacid_center(mol) is not None:
        return name
    if tokens == 0 or (centers > 1 and tokens < 2 and not _MULTIPLE_ANION.search(name)):
        raise UnsupportedStructure("the anionic centers of this structure are not all cited by a supported name")
    return name


def _name_anion_unchecked(mol):
    centers = anion_atoms(mol)
    if len(centers) < 2 or any(_is_group_anion(a) for a in centers):
        return _name_substitutive(mol)
    error = None
    try:
        name = _name_substitutive(mol)
    except UnsupportedStructure as caught:
        name, error = None, caught
    if name is None or "idyl" in name or not _ANION_TOKEN.search(name):
        multiplicative = _multiplicative_anion(mol)
        if multiplicative is not None:
            return multiplicative
    if name is None:
        raise error
    return name


@contextlib.contextmanager
def _anion_stereo(mol, neutral):
    """The stereo elements of the anion carried over to its protonated copy, where two rings that differ only in
    the charged atom would otherwise leave a centre unrecognised (P-92.1.4.4)."""
    from ._polyfunctional import STEREO_OF_ISOTOPOLOGUE

    elements = specified_stereo_elements(mol) or []
    if len(elements) == len(specified_stereo_elements(neutral) or []):
        yield
        return
    located = []
    for kind, idx, code in elements:
        if kind == "bond":
            bond = mol.GetBondWithIdx(idx)
            located.append(("bond", (bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()), code))
        else:
            located.append((kind, idx, code))
    token = STEREO_OF_ISOTOPOLOGUE.set(located)
    try:
        yield
    finally:
        STEREO_OF_ISOTOPOLOGUE.reset(token)


def _name_substitutive(mol):
    from ._anion_center import has_center_anion_shape, name_center_anion
    from ._anion_oxoacid import name_oxoacid_anion, oxoacid_center
    from ._polyfunctional import (
        chalcogen_acid_variant,
        name_polyfunctional,
        peroxy_carbonyl_thio,
        peroxy_variant,
    )

    oxoacid = oxoacid_center(mol)
    if oxoacid is not None:
        return name_oxoacid_anion(mol, oxoacid)
    if has_center_anion_shape(mol) and (not _has_group_or_carbon(mol) or _ring_mixed_centers(mol)):
        return name_center_anion(mol)
    neutral = marked_neutral(mol)
    ring_cations = [a for a in neutral.GetAtoms() if a.GetFormalCharge() > 0 and a.GetAtomicNum() != 7 and not _is_nitro_nitrogen(a)]
    if ring_cations:
        from ._polycation import _name_ring_polycation

        return _name_ring_polycation(neutral, ring_cations)
    with _anion_stereo(mol, neutral):
        try:
            name = name_polyfunctional(neutral)
        except UnsupportedStructure:
            try:
                name = _multiplicative_neutral_name(neutral)
            except UnsupportedStructure:
                name = _pipeline_neutral_name(neutral)
    variant = peroxy_variant(neutral)
    if variant is not None:
        return _peroxy_swap(name, variant, peroxy_carbonyl_thio(neutral))
    chalco = chalcogen_acid_variant(neutral)
    if chalco is not None and _THIOIC_TAIL.search(name):
        return _THIOIC_TAIL.sub(lambda m: f"{m.group(1)}{chalco}ate", name)
    if not any(a.HasProp(ANION_PROP) and a.GetAtomicNum() != 6 for a in neutral.GetAtoms()):
        return _added_hydrogen(name)
    if any(a.HasProp(ANION_PROP) and anion_weight(a) == 2 for a in neutral.GetAtoms()):
        return _swap_aminediide(name)
    if _CARBANIDE_PARENT.search(name) and any(a.HasProp(ANION_PROP) and a.GetAtomicNum() == 6 for a in neutral.GetAtoms()):
        return name
    return swap_suffix(name)


def has_general_anion_shape(mol):
    if not any(a.GetFormalCharge() < 0 for a in mol.GetAtoms()):
        return False
    try:
        name_anion(mol)
    except UnsupportedStructure:
        return False
    return True


def _swap_aminediide(name):
    if name.endswith("aniline"):
        return name[: -len("aniline")] + "benzenaminediide"
    if name.endswith("amine") and not name.endswith("diamine"):
        return name[: -len("amine")] + "aminediide"
    raise UnsupportedStructure("this dianionic nitrogen has no aminediide name yet")


def _is_peroxy_center(atom):
    return atom.GetDegree() == 1 and atom.GetNeighbors()[0].GetDegree() == 2 and atom.GetNeighbors()[0].GetAtomicNum() in _CHALCOGENS


_PEROXY_TAIL = re.compile(r"(di|tri|tetra)?peroxo(ic acid|l)$")


def _peroxy_swap(name, variant, carbonyl_thio=""):
    match = _PEROXY_TAIL.search(name)
    if match is None:
        raise UnsupportedStructure("this peroxy anion has no anion suffix form yet")
    mixed = variant.startswith("(")
    core = (variant[1:] if mixed else variant) + "peroxo"
    word = core + (carbonyl_thio + "ate" if match.group(2) == "ic acid" else "late")
    stem = name[: match.start()]
    if match.group(1):
        return f"{stem}{_MULTIPLIED[match.group(1)]}({word})"
    return f"{stem}({word})" if mixed else stem + word


def _has_group_or_carbon(mol):
    return any(_is_group_anion(a) or _is_carbanion(a) for a in anion_atoms(mol))


def _multiplicative_anion(mol):
    from ._chain_multiplicative import anion_multiplicative_name

    centers = [a.GetIdx() for a in anion_atoms(mol)]
    return anion_multiplicative_name(mol, centers)


_HYDRO_IDE = re.compile(
    r"^(?P<pre>.*?)(?P<a>\d+),(?P<b>\d+)-dihydro(?P<parent>[a-z]+?)e?-(?P<c>\d+)-(?P<word>u?ide)$"
)


def _added_hydrogen(name):
    """'1,2-dihydropyridin-2-ide' -> 'pyridin-2(1H)-ide' when the other hydro position is not a center (P-72.2.2.1.1)."""
    match = _HYDRO_IDE.match(name)
    if match is None:
        return name
    a, b, c = (int(match.group(k)) for k in ("a", "b", "c"))
    if c not in (a, b):
        return name
    other = b if c == a else a
    pre = match.group("pre")
    if pre.endswith("-"):
        pre = pre[:-1]
    return f"{pre}{match.group('parent')}-{c}({other}H)-{match.group('word')}"


def _multiplicative_neutral_name(neutral):
    """Multiplicative name of the protonated structure when every principal group is an anion (P-72.5.1)."""
    from ._multiplicative import name_if_multiplicative
    from ._multiplicative_groups import classify

    groups = classify(neutral)
    if not groups:
        raise UnsupportedStructure("this anionic structure has no multiplicative name yet")
    top = min(g.rank for g in groups)
    marked = {a.GetIdx() for a in neutral.GetAtoms() if a.HasProp(ANION_PROP)}
    if any(not (g.atoms & marked) for g in groups if g.rank == top):
        raise UnsupportedStructure("a neutral group of the principal class remains beside the anions")
    name = name_if_multiplicative(neutral)
    if name is None:
        raise UnsupportedStructure("this anionic structure has no multiplicative name yet")
    return name


def _ring_mixed_centers(mol):
    centers = anion_atoms(mol)
    return (
        all(a.IsInRing() and not _is_group_anion(a) for a in centers)
        and any(a.GetAtomicNum() != 6 for a in centers)
        and any(a.GetAtomicNum() == 6 for a in centers)
    )


def _pipeline_neutral_name(neutral):
    """Name of the protonated structure by the full pipeline when every principal group is an anion."""
    from ._multiplicative_groups import classify
    from .core import smiles_to_iupac

    groups = classify(neutral)
    if not groups:
        raise UnsupportedStructure("this anionic structure has no neutral parent name yet")
    top = min(g.rank for g in groups)
    marked = {a.GetIdx() for a in neutral.GetAtoms() if a.HasProp(ANION_PROP)}
    if any(not (g.atoms & marked) for g in groups if g.rank == top):
        raise UnsupportedStructure("a neutral group of the principal class remains beside the anions")
    plain = Chem.Mol(neutral)
    for atom in plain.GetAtoms():
        for prop in (ANION_PROP, "_anion_word", "_anion_charge", "_anion_lambda"):
            if atom.HasProp(prop):
                atom.ClearProp(prop)
    return smiles_to_iupac(Chem.MolToSmiles(plain))
