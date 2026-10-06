"""A principal characteristic group on an acyclic substituent decides the parent hydride (P-101.7.1.1.4, P-44.1.1)."""

import re

from rdkit import Chem

from ._common import UnsupportedStructure

_PLACEHOLDERS = ("Cl", "Br", "I", "F")
_SUFFIXES = ("oic acid", "amide", "amine", "ol", "al", "one")


def chain_template(cand, view, branches):
    """(location, root, surrogate name with a placeholder prefix) for the one acyclic chain that holds the
    principal group when the skeleton holds none, else None."""
    from .core import smiles_to_iupac

    mapped = set(cand.mapping.values())
    placeholder = next((p for p in _PLACEHOLDERS if p not in set(view.elem.values())), None)
    if placeholder is None:
        return None
    found = []
    for loc, root in branches:
        if view.elem[root] != "C":
            continue
        atoms = _branch_atoms(view, root, mapped)
        if any(a in view.rings for a in atoms):
            continue
        editable = Chem.RWMol(view.mol)
        keep = sorted(atoms)
        for idx in sorted((i for i in view.adj if i not in atoms), reverse=True):
            editable.RemoveAtom(idx)
        new = editable.AddAtom(Chem.Atom(placeholder))
        editable.AddBond(keep.index(root), new, Chem.BondType.SINGLE)
        surrogate = editable.GetMol()
        for atom in surrogate.GetAtoms():
            atom.SetNoImplicit(False)
            atom.SetNumExplicitHs(0)
        try:
            Chem.SanitizeMol(surrogate)
            name = smiles_to_iupac(Chem.MolToSmiles(surrogate))
        except Exception:
            continue
        stem = placeholder_prefix(placeholder)
        match = re.fullmatch(rf"(?:(\d+)-)?{stem}([a-z]+?)(?:-(\d+)-)?({'|'.join(_SUFFIXES)})", name)
        if match is None and not re.fullmatch(rf"{stem}[a-z]+", name):
            continue
        if any(name.endswith(suffix) for suffix in _SUFFIXES):
            found.append((loc, root, name, stem))
    if len(found) > 1:
        raise UnsupportedStructure("several acyclic substituents carry principal groups")
    return found[0] if found else None


def placeholder_prefix(placeholder):
    return {"Cl": "chloro", "Br": "bromo", "I": "iodo", "F": "fluoro"}[placeholder]


def _branch_atoms(view, root, mapped):
    seen, stack = {root}, [root]
    while stack:
        for n in view.adj[stack.pop()]:
            if n not in seen and n not in mapped:
                seen.add(n)
                stack.append(n)
    return seen
