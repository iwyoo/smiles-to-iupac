"""Prefix names of acid groups and derivatives attached through carbon (P-65.1.3-P-65.1.5, P-65.2.1.4-P-65.2.1.7, P-65.5.4):
the group C(=X)Z is the divalent acyl group of =X ('carbonyl', 'carbonothioyl', 'carbonimidoyl') with the prefix of Z in
front, except the retained 'carboxy', 'carbamoyl' ... and the halide and pseudohalide infix names ('carbonochloridoyl').
"""

from ._multiplicative_text import enclose

_X_INFIX = {"O": "", "S": "thio", "Se": "seleno", "Te": "telluro", "NH": "imido", "NNH2": "hydrazono"}
_SYMBOL = {8: "O", 16: "S", 34: "Se", 52: "Te"}
_HALIDE_INFIX = {9: "fluorid", 17: "chlorid", 35: "bromid", 53: "iodid"}
_CARBAM = {"O": "carbamoyl", "S": "carbamothioyl", "Se": "carbamoselenoyl", "Te": "carbamotelluroyl", "NH": "carbamimidoyl", "NNH2": "carbamohydrazonoyl"}
_FORM = {"O": "formyl", "S": "methanethioyl", "Se": "methaneselenoyl", "Te": "methanetelluroyl", "NH": "methanimidoyl", "NNH2": "methanehydrazonoyl"}


def _bond(mol, a, b):
    return mol.GetBondBetweenAtoms(a, b).GetBondTypeAsDouble()


def _vowel(text):
    return text[0] in "aeiouy"


def acyl_stem(x):
    """'carbonyl', 'carbonothioyl', 'carbonimidoyl' ... for the divalent acyl group of the =X atom."""
    infix = _X_INFIX[x]
    if not infix:
        return "carbonyl"
    word = infix + "yl"
    return "carbon" + ("" if _vowel(word) else "o") + word


def halide_acyl_name(infix, x):
    """'carbonochloridoyl', 'carbonochloridothioyl', 'carbonochloridimidoyl' (P-65.2.1.5)."""
    join = "o" if not _vowel(infix) else ""
    word = _X_INFIX[x]
    if not word:
        return "carbon" + join + infix + "oyl"
    tail = word + "yl"
    return "carbon" + join + infix + ("" if _vowel(tail) else "o") + tail


def _oxo(mol, graph, root, coming_from):
    """((atom, 'O'|'S'|'Se'|'Te'|'NH'|'NNH2', N-substituent atoms), other neighbors) or (None, others)."""
    others = [n for n in graph[root] if n != coming_from]
    slots = []
    for n in others:
        atom = mol.GetAtomWithIdx(n)
        if _bond(mol, root, n) != 2.0 or atom.GetFormalCharge():
            continue
        z = atom.GetAtomicNum()
        far = [m.GetIdx() for m in atom.GetNeighbors() if m.GetIdx() != root]
        if z in _SYMBOL and atom.GetDegree() == 1:
            slots.append((n, _SYMBOL[z], []))
        elif z == 7 and atom.GetTotalNumHs() == 0 and len(far) == 1 and _is_hydrazine_tail(mol, far[0], n):
            slots.append((n, "NNH2", []))
        elif z == 7 and all(_bond(mol, n, m) == 1.0 for m in far) and atom.GetTotalNumHs() + len(far) == 1:
            slots.append((n, "NH", far))
    if len(slots) != 1:
        return None, others
    return slots[0], [n for n in others if n != slots[0][0]]


def _is_hydrazine_tail(mol, idx, parent):
    atom = mol.GetAtomWithIdx(idx)
    return atom.GetAtomicNum() == 7 and atom.GetDegree() == 1 and atom.GetTotalNumHs() == 2 and not atom.GetFormalCharge()


def _cited(entries):
    """Prefix text for the Z and N substituents of an imidic-type group: locants C and N cited (P-65.1.3.1.2)."""
    from ._substituents import format_substituent_prefixes

    grouped = {}
    for locant, name, compound in entries:
        grouped.setdefault(name, {"locants": [], "compound": compound})["locants"].append(locant)
    return format_substituent_prefixes(grouped)


def acid_group_prefix(mol, graph, root, coming_from, halogens, aromatic_atoms, name_branch, enclose_mark=enclose):
    """(name, is_compound) of the acid-derived group on carbon `root`, or None."""
    atom = mol.GetAtomWithIdx(root)
    if atom.IsInRing() or atom.GetFormalCharge() or atom.GetIsotope() or _bond(mol, root, coming_from) != 1.0:
        return None
    slot, rest = _oxo(mol, graph, root, coming_from)
    if slot is None:
        return None
    x, n_subs = slot[1], slot[2]
    n_entries = [("N", *name_branch(graph, m, slot[0], halogens, aromatic_atoms, mol=mol)) for m in n_subs]
    if not rest:
        if n_entries:
            return None
        return _FORM[x], x != "O"
    if len(rest) != 1 or _bond(mol, root, rest[0]) != 1.0:
        return None
    z_idx = rest[0]
    zn = mol.GetAtomWithIdx(z_idx).GetAtomicNum()
    z_atom = mol.GetAtomWithIdx(z_idx)
    if z_atom.GetFormalCharge() or z_atom.HasProp("_anion") or z_atom.HasProp("_anion_word"):
        return None
    if zn in _HALIDE_INFIX:
        if n_entries:
            return None
        return halide_acyl_name(_HALIDE_INFIX[zn], x), False
    if zn == 6:
        if x == "O" or n_entries:
            return None
        from ._hetero_prefixes import _acyl_from_acid_name
        from ._retained_acids import is_compound_acyl

        name = _acyl_from_acid_name(mol, graph, root, coming_from)
        return name, is_compound_acyl(name)
    if zn == 7:
        subs = [n for n in graph[z_idx] if n != root]
        if mol.GetAtomWithIdx(z_idx).GetFormalCharge() or any(_bond(mol, z_idx, n) != 1.0 for n in subs):
            return None
        from ._hetero_prefixes import _amino, _amino_stem, _group_names, _is_amino_nitrogen

        if not subs and not n_entries:
            return _CARBAM[x], False
        if x == "O" and len(subs) == 1 and _is_amino_nitrogen(mol, subs[0], z_idx):
            return "hydrazinecarbonyl", True
        stem = _CARBAM[x]
        if subs and mol.GetAtomWithIdx(z_idx).IsInRing() and x in ("O", "S"):
            from ._hetero_prefixes import _ring_nitrogen_acyl

            return _ring_nitrogen_acyl(graph, z_idx, root, halogens, aromatic_atoms, mol, {"O": "carbonyl", "S": "carbothioyl"}[x])
        amino = _amino(_group_names(graph, mol, subs, z_idx, halogens, aromatic_atoms)) if subs else "amino"
        if n_entries:
            entries = [("N", n, c) for _, n, c in n_entries]
            entries += [("N'", *name_branch(graph, m, z_idx, halogens, aromatic_atoms, mol=mol)) for m in subs]
            return _cited(entries) + stem, True
        return _amino_stem(amino) + stem, True
    if zn not in _SYMBOL:
        return None
    z_name, z_compound = name_branch(graph, z_idx, root, halogens, aromatic_atoms, mol=mol)
    if x in ("NH", "NNH2"):
        entries = [("C", z_name, z_compound), *n_entries]
        return _cited(entries) + acyl_stem(x), True
    if z_name == "hydroxy":
        if x == "O":
            return "carboxy", False
        return "hydroxy" + acyl_stem(x), True
    if z_name == "hydroperoxy" and x == "O":
        return "carbonoperoxoyl", False
    if z_name == "ylooxidanyl" and x == "O":
        return "oxylcarbonyl", True
    return (enclose_mark(z_name) if z_compound else z_name) + acyl_stem(x), True
