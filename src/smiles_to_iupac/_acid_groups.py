"""Recognition of acid groups and their functional-replacement analogues on a
carbon, sulfur, selenium or tellurium centre (P-65.1.3-P-65.1.5, P-65.3.1)."""

from dataclasses import dataclass

from ._acid_lexicon import AcidSpec, make_spec

_SYMBOL = {8: "O", 16: "S", 34: "Se", 52: "Te"}
_CENTER_SYMBOL = {6: "C", 16: "S", 34: "Se", 52: "Te"}


@dataclass
class AcidGroup:
    center: int
    spec: AcidSpec
    owned: set
    n_substituents: list
    n_atom: int = None


def _bond(mol, a, b):
    return mol.GetBondBetweenAtoms(a, b).GetBondTypeAsDouble()


def _terminal(atom, hydrogens):
    return atom.GetDegree() == 1 and atom.GetTotalNumHs() == hydrogens and not atom.GetFormalCharge()


def _oxo_slot(mol, center, atom):
    """('O'|'S'|'Se'|'Te'|'NH'|'NNH2', owned atoms, N-substituent atoms) for a double-bonded neighbor, else None."""
    z = atom.GetAtomicNum()
    if _bond(mol, center, atom.GetIdx()) != 2.0 or atom.GetFormalCharge():
        return None
    if z in _SYMBOL:
        return (_SYMBOL[z], {atom.GetIdx()}, []) if atom.GetDegree() == 1 else None
    if z != 7:
        return None
    others = [n for n in atom.GetNeighbors() if n.GetIdx() != center]
    if not others:
        return ("NH", {atom.GetIdx()}, []) if atom.GetTotalNumHs() == 1 else None
    if len(others) == 1 and others[0].GetAtomicNum() == 7 and _terminal(others[0], 2) and atom.GetTotalNumHs() == 0:
        return ("NNH2", {atom.GetIdx(), others[0].GetIdx()}, [])
    if atom.GetTotalNumHs() == 0 and all(_bond(mol, atom.GetIdx(), n.GetIdx()) == 1.0 for n in others):
        return ("NH", {atom.GetIdx()}, [n.GetIdx() for n in others])
    return None


def _y_chain(mol, center, atom):
    """(atoms from the centre outward, owned atoms) for a -Y-H or -Y-Y-H chain on the centre, else None."""
    z = atom.GetAtomicNum()
    if z not in _SYMBOL or _bond(mol, center, atom.GetIdx()) != 1.0 or atom.GetFormalCharge():
        return None
    if _terminal(atom, 1):
        return [_SYMBOL[z]], {atom.GetIdx()}
    onward = [n for n in atom.GetNeighbors() if n.GetIdx() != center]
    if len(onward) == 1 and onward[0].GetAtomicNum() in _SYMBOL and _terminal(onward[0], 1) and _bond(mol, atom.GetIdx(), onward[0].GetIdx()) == 1.0 and atom.GetTotalNumHs() == 0:
        return [_SYMBOL[z], _SYMBOL[onward[0].GetAtomicNum()]], {atom.GetIdx(), onward[0].GetIdx()}
    return None


def _y_anion(mol, center, atom):
    """(atoms from the centre outward, owned atoms) for a -Y(-) or -Y-Y(-) chain on the centre, else None."""
    z = atom.GetAtomicNum()
    if z not in _SYMBOL or _bond(mol, center, atom.GetIdx()) != 1.0 or atom.GetFormalCharge():
        if z in _SYMBOL and atom.GetFormalCharge() == -1 and atom.GetDegree() == 1 and _bond(mol, center, atom.GetIdx()) == 1.0:
            return [_SYMBOL[z]], {atom.GetIdx()}
        return None
    onward = [n for n in atom.GetNeighbors() if n.GetIdx() != center]
    if len(onward) == 1 and onward[0].GetAtomicNum() in _SYMBOL and onward[0].GetFormalCharge() == -1 and onward[0].GetDegree() == 1:
        return [_SYMBOL[z], _SYMBOL[onward[0].GetAtomicNum()]], {atom.GetIdx(), onward[0].GetIdx()}
    return None


def acid_group_at(mol, idx, hetero_attach=False):
    """The acid group centred on atom `idx` with exactly one carbon attachment, or None."""
    atom = mol.GetAtomWithIdx(idx)
    center = _CENTER_SYMBOL.get(atom.GetAtomicNum())
    if center is None or atom.GetFormalCharge() or atom.GetIsotope():
        return None
    oxo, y, owned, n_sub, n_atom, anion = [], None, set(), [], None, False
    attachments = []
    for n in atom.GetNeighbors():
        slot = _oxo_slot(mol, idx, n)
        if slot is not None:
            oxo.append(slot[0])
            owned |= slot[1]
            if slot[2]:
                n_sub, n_atom = slot[2], n.GetIdx()
            continue
        chain = _y_chain(mol, idx, n)
        if chain is None:
            chain = _y_anion(mol, idx, n)
            anion = anion or chain is not None
        if chain is not None and y is None:
            y = chain[0]
            owned |= chain[1]
            continue
        attachments.append(n)
    formic = center == "C" and not attachments and atom.GetTotalNumHs() == 1
    if y is None or not formic and (
        len(attachments) != 1
        or not (attachments[0].GetAtomicNum() == 6 or attachments[0].IsInRing() or hetero_attach)
        or _bond(mol, idx, attachments[0].GetIdx()) != 1.0
    ):
        return None
    if center == "C":
        if len(oxo) != 1:
            return None
    elif len(oxo) not in (1, 2) or atom.GetDegree() != len(oxo) + 2:
        return None
    return AcidGroup(idx, make_spec(center, oxo, y, anion), owned, n_sub, n_atom)
