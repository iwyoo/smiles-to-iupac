"""Radical ions whose radical and ionic centres are carbon atoms of one skeleton (P-75.2): the anionic or cationic
suffixes follow the parent hydride, then the radical suffixes; radical centres take the lowest locants."""

import re

from rdkit import Chem

from ._common import UnsupportedStructure
from ._numerals import alkane_name, multiplying_prefix

_MAX_ATOMS = 60


def _centres(mol):
    radicals = [a for a in mol.GetAtoms() if a.GetNumRadicalElectrons()]
    ions = [a for a in mol.GetAtoms() if a.GetFormalCharge()]
    if (
        not radicals
        or not ions
        or any(a.GetNumRadicalElectrons() != 1 or a.GetAtomicNum() != 6 for a in radicals)
        or any(a.GetAtomicNum() != 6 or abs(a.GetFormalCharge()) != 1 for a in ions)
        or len({a.GetFormalCharge() for a in ions}) != 1
        or any(a.GetIsotope() for a in mol.GetAtoms())
        or len(Chem.GetMolFrags(mol)) != 1
        or mol.GetNumAtoms() > _MAX_ATOMS
    ):
        return None
    return radicals, ions


def has_skeleton_radical_ion_shape(mol) -> bool:
    return _centres(mol) is not None


def _suffix(count, singular):
    return singular if count == 1 else multiplying_prefix(count) + singular


def _ending(charge):
    return "id" if charge < 0 else "ylium"


def _chain_or_ring(mol):
    if any(a.GetAtomicNum() != 6 or a.GetDegree() > 2 or a.GetIsAromatic() for a in mol.GetAtoms()) or any(
        b.GetBondTypeAsDouble() != 1.0 for b in mol.GetBonds()
    ):
        return None
    rings = mol.GetRingInfo().NumRings()
    if rings > 1 or (rings == 1 and mol.GetNumBonds() != mol.GetNumAtoms()):
        return None
    return rings


def _orders(mol, ring):
    graph = {a.GetIdx(): [n.GetIdx() for n in a.GetNeighbors()] for a in mol.GetAtoms()}
    size = mol.GetNumAtoms()
    walk = [next(i for i, n in graph.items() if len(n) == 1)] if not ring else [next(iter(graph))]
    previous = None
    while len(walk) < size:
        nxt = [n for n in graph[walk[-1]] if n != previous]
        previous = walk[-1]
        walk.append(nxt[0])
    if not ring:
        return [walk, walk[::-1]]
    return [[walk[(s + step * k) % size] for k in range(size)] for s in range(size) for step in (1, -1)]


def _simple_name(mol, ring, radicals, ions):
    """Unbranched chain or monocycle; a radical and an ionic centre may share an atom."""
    best = None
    for order in _orders(mol, bool(ring)):
        locant = {atom: i + 1 for i, atom in enumerate(order)}
        key = (sorted(locant[a.GetIdx()] for a in radicals), sorted(locant[a.GetIdx()] for a in ions))
        if best is None or key < best:
            best = key
    return _assemble(("cyclo" if ring else "") + alkane_name(mol.GetNumAtoms()), best[0], best[1], ions[0].GetFormalCharge())


def _assemble(base, radical_locants, ion_locants, charge, prefixes=""):
    ion_text = _suffix(len(ion_locants), _ending(charge))
    radical_text = _suffix(len(radical_locants), "yl")
    if base.endswith("ane") and ion_text[0] in "aeiouy":
        base = base[:-1]
    elif base.endswith("an") and ion_text[0] not in "aeiouy":
        base += "e"
    return (
        f"{prefixes}{base}-{','.join(map(str, ion_locants))}-{ion_text}-{','.join(map(str, radical_locants))}-{radical_text}"
    )


def _automorphisms(mol, centres):
    centre_set = {a.GetIdx() for a in centres}
    maps = mol.GetSubstructMatches(Chem.Mol(mol), uniquify=False, useChirality=False, maxMatches=5000)
    return [m for m in maps if {m[i] for i in centre_set} == centre_set]


def _general_name(mol, radicals, ions):
    from ._anion import name_anion
    from ._polyfunctional import LAST_POSITIONS

    charge = ions[0].GetFormalCharge()
    analogue = Chem.RWMol(mol)
    for atom in radicals:
        target = analogue.GetAtomWithIdx(atom.GetIdx())
        target.SetNumRadicalElectrons(0)
        target.SetFormalCharge(-1)
    for atom in ions:
        analogue.GetAtomWithIdx(atom.GetIdx()).SetFormalCharge(-1)
    analogue = analogue.GetMol()
    Chem.SanitizeMol(analogue)
    LAST_POSITIONS.set(None)
    name = name_anion(analogue)
    found = LAST_POSITIONS.get()
    tail = re.search(r"-?(\d+(?:,\d+)*)-(di|tri|tetra)?ide$", name)
    centres = {a.GetIdx() for a in radicals + ions}
    if tail is None or found is None or any(c not in found[1] for c in centres):
        raise UnsupportedStructure("the skeleton of this radical ion is not named as a polyanion parent")
    positions = found[1]
    parent = set(positions)
    flip = mol.GetRingInfo().NumRings() == 0
    size = len(parent)
    attachments = {
        a for atom in parent for a in (n.GetIdx() for n in mol.GetAtomWithIdx(atom).GetNeighbors()) if a not in parent
    }
    best = None
    for image in _automorphisms(mol, radicals + ions):
        if any(image[i] not in positions for i in positions):
            continue
        numbering = {i: positions[image[i]] for i in positions}
        for numbers in ([numbering, {i: size + 1 - n for i, n in numbering.items()}] if flip else [numbering]):
            rad = sorted(numbers[a.GetIdx()] for a in radicals)
            ion = sorted(numbers[a.GetIdx()] for a in ions)
            substituted = sorted(
                numbers[n.GetIdx()] for a in attachments for n in mol.GetAtomWithIdx(a).GetNeighbors() if n.GetIdx() in parent
            )
            key = (sorted(set(rad + ion)), rad, ion, substituted)
            if best is None or key < best[0]:
                best = (key, numbers)
    if best is None:
        raise UnsupportedStructure("no numbering of this radical ion is available")
    numbers = best[1]
    head = name[: tail.start()]
    prefixes, base = _split_parent(head)
    remap = {positions[i]: numbers[i] for i in positions}
    return _assemble(base, best[0][1], best[0][2], charge, _renumber(prefixes, remap))


def _renumber(prefixes, remap):
    depth, out, buffer = 0, [], ""
    for ch in prefixes + "\0":
        if depth == 0 and (ch.isdigit() or (ch == "," and buffer)):
            buffer += ch
            continue
        if buffer:
            if ch == "-":
                numbers = sorted(remap.get(int(x), int(x)) for x in buffer.split(","))
                out.append(",".join(map(str, numbers)))
            else:
                out.append(buffer)
            buffer = ""
        if ch == "\0":
            break
        out.append(ch)
        depth += ch in "([{"
        depth -= ch in ")]}"
    return "".join(out)


def _split_parent(head):
    match = re.search(r"([a-z]*(?:cyclo)?[a-z]*an[e]?)$", head)
    if match is None:
        raise UnsupportedStructure("the parent hydride of this radical ion is not delimited")
    return head[: match.start()], match.group(1)


def name_skeleton_radical_ion(mol) -> str:
    centres = _centres(mol)
    if centres is None:
        raise UnsupportedStructure("this radical ion is not a set of carbon centres on one skeleton")
    radicals, ions = centres
    ring = _chain_or_ring(mol)
    if ring is not None and mol.GetNumAtoms() > 1:
        return _simple_name(mol, ring, radicals, ions)
    if mol.GetNumAtoms() == 1:
        return "methanidyl" if ions[0].GetFormalCharge() < 0 else "methyliumyl"
    return _general_name(mol, radicals, ions)
