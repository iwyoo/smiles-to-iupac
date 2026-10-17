"""Oxygen nuclides of a polycarbonic acid chain (P-82.2.1, P-82.6.4): the chain is numbered with the carbons at the odd
positions and the bridging oxygens at the even ones; a bridging oxygen is cited with its chain locant, any other oxygen
with the letter O and the locant of its carbon."""

import re

from ._common import UnsupportedStructure
from ._isotope_labels import descriptor
from ._substituents import ISOTOPE_LABELS


def _chain_oxygen_nuclides(mol, order, bridges):
    """{oxygen atom: locant text} of the labelled oxygens of the chain numbered along `order`: a bridging oxygen takes its
    chain locant, any other oxygen of a carbon the letter O and the locant of that carbon (P-82.6.4)."""
    context = ISOTOPE_LABELS.get()
    if not context or context["mol"] is not mol:
        return {}
    open_labels = {a: e for a, e in context["labels"].items() if a not in context["consumed"]}
    located = {}
    for i, center in enumerate(order):
        for n in mol.GetAtomWithIdx(center).GetNeighbors():
            if n.GetIdx() in open_labels and n.GetAtomicNum() == 8 and n.GetDegree() == 1:
                located[n.GetIdx()] = f"O{2 * i + 1}"
    for i in range(len(order) - 1):
        bridge = bridges[frozenset((order[i], order[i + 1]))]
        if isinstance(bridge, int) and bridge in open_labels and mol.GetAtomWithIdx(bridge).GetAtomicNum() == 8:
            located[bridge] = str(2 * i + 2)
    return located


def chain_nuclide_key(mol, order, bridges):
    """The locant numbers of the labelled oxygens along `order`, sorted, so that the numbering giving them the lowest
    locants is taken (P-31.1.4.3.4)."""
    return sorted(int(text.lstrip("O")) for text in _chain_oxygen_nuclides(mol, order, bridges).values())


def chain_nuclides_text(mol, order, bridges):
    """'(2,O1-18O)': the descriptor of the oxygen nuclides of the chain, '' without any (P-82.2.1, P-82.6.4)."""
    located = _chain_oxygen_nuclides(mol, order, bridges)
    if not located:
        return ""
    context = ISOTOPE_LABELS.get()
    entries = [context["labels"][a] for a in located]
    if any(e["H"] or not e["skeleton"] for e in entries):
        raise UnsupportedStructure("a hydrogen nuclide on a polycarbonic acid is not supported yet")
    context["consumed"].update(located)
    bridging = {a: e for a, e in zip(located, entries) if located[a][0] != "O"}
    return descriptor(
        bridging,
        {a: int(located[a]) for a in bridging},
        False,
        extra=[(e["skeleton"], located[a], 1, False) for a, e in zip(located, entries) if a not in bridging],
    )


def insert_chain_nuclides(name, text):
    """`name` with the descriptor `text` in front of the stem of the polycarbonic acid."""
    if not text:
        return name
    found = re.search(r"(?:di|tri|tetra)carbon", name)
    if found is None:
        raise UnsupportedStructure("the stem of the polycarbonic acid is not delimited")
    return name[: found.start()] + text + name[found.start():]

