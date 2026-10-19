"""Acyl and ester words shared by the retained steroid parents and the stereoparent builder (P-65.6.3.3, P-101.7.2):
anions of O-acyl groups, alkyl groups of several esters and the expression of a junior acyl group as prefixes."""

from collections import Counter

from ._common import UnsupportedStructure, alpha_sort_key
from ._diester_anions import acid_anions
from ._multiplicative_text import enclose
from ._numerals import multiplying_prefix


def _carbonyl_oxygen(carbon):
    return next(
        (
            n
            for n in carbon.GetNeighbors()
            if n.GetAtomicNum() == 8 and n.GetDegree() == 1 and carbon.GetOwningMol().GetBondBetweenAtoms(carbon.GetIdx(), n.GetIdx()).GetBondTypeAsDouble() == 2.0
        ),
        None,
    )


def o_acyl_anions(mol, groups):
    """[(anion name, substituted)] for O-acyl groups given as (ester oxygen, acyl carbon) pairs."""
    matches = []
    for oxygen, acyl in groups:
        carbon = mol.GetAtomWithIdx(acyl)
        owner = next(n for n in mol.GetAtomWithIdx(oxygen).GetNeighbors() if n.GetIdx() != acyl)
        matches.append((carbon, _carbonyl_oxygen(carbon), mol.GetAtomWithIdx(oxygen), owner))
    return acid_anions(mol, matches)


def cite_alkyls(entries):
    """Alkyl words in front of the name of an ester of one acid parent (P-65.6.3.3.2.1): identical groups are multiplied,
    different groups follow alphanumerically, each with the locants of the acid groups it esterifies. `entries` are
    (locant text, name, compound)."""
    groups = {}
    for locant, name, compound in entries:
        entry = groups.setdefault(name, {"compound": compound, "locants": []})
        entry["locants"].append(locant)
    words = []
    for name in sorted(groups, key=alpha_sort_key):
        entry = groups[name]
        count = len(entry["locants"])
        word = name
        if count > 1:
            word = multiplying_prefix(count, compound=entry["compound"]) + (enclose(name) if entry["compound"] else name)
        elif entry["compound"] and len(groups) > 1:
            word = enclose(name)
        words.append(f"{','.join(str(loc) for loc in entry['locants'])}-{word}" if len(groups) > 1 else word)
    return " ".join(words)


def demote_acyl(cls, members, classes, branches):
    """Express acyl groups that rank below the principal group as prefixes: a substituent group is cited whole, a
    terminal one of the parent as its 'oxo' with 'hydroxy', 'amino' or 'alkoxy' on the same atom (P-65.6.3.2.3)."""
    for position, atoms, extra in members:
        if "root" in extra:
            branches.append((position, extra["root"]))
            continue
        classes.setdefault("ketone", []).append((position, atoms, {}))
        if cls == "acid":
            classes.setdefault("alcohol", []).append((position, atoms, {}))
        elif cls == "amide":
            classes.setdefault("amine", []).append((position, atoms, {}))
        elif cls == "ester":
            branches.append((position, extra["ester"][0]))
        elif cls != "aldehyde":
            raise UnsupportedStructure("a terminal acyl group of this class below the principal group is not supported")


def majority_kind(members, tie):
    """The kind ('o' terminal of the parent, 'c' carbon outside it) held by more of the members; `tie` wins a draw."""
    counts = Counter(extra.get("kind", "o") for _, _, extra in members)
    return max(counts, key=lambda kind: (counts[kind], kind == tie))
