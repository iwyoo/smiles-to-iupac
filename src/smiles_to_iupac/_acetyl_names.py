"""Retained 'acetate', 'acetyl' and 'formyl' in names of anions (P-65.1.1.1, P-65.1.7.2.1, P-72.2.2.2.1.1):
the systematic 'ethanoate'/'ethanoyl' of the acid modules become the preferred retained forms, and
substituent locants on acetate (always 2) are omitted: '2-chloroethanoate' -> 'chloroacetate'.
"""

import re

_OPEN = "([{"
_CLOSE = ")]}"
_TWO_LOCANT = re.compile(r"2(?:,2)*-")


def acetyl_names(name):
    name = name.replace("ethanoyl", "acetyl").replace("methanoyl", "formyl")
    position = name.find("ethanoate")
    while position != -1:
        start = _prefix_start(name, position)
        prefix = _drop_two_locants(name[start:position])
        name = name[:start] + prefix + "acetate" + name[position + len("ethanoate"):]
        position = name.find("ethanoate", start + len(prefix) + len("acetate"))
    return name


def _prefix_start(name, end):
    depth = 0
    i = end
    while i > 0:
        ch = name[i - 1]
        if ch in _CLOSE:
            depth += 1
        elif ch in _OPEN:
            if depth == 0:
                break
            depth -= 1
        elif ch == " " and depth == 0:
            break
        i -= 1
    return i


def _drop_two_locants(prefix):
    if not prefix:
        return prefix
    units, depth, i, current = [], 0, 0, None
    while i < len(prefix):
        ch = prefix[i]
        if depth == 0 and (i == 0 or prefix[i - 1] == "-"):
            match = _TWO_LOCANT.match(prefix[i:])
            if match:
                if current is not None:
                    units.append(current.rstrip("-"))
                current = ""
                i += match.end()
                continue
        if ch in _OPEN:
            depth += 1
        elif ch in _CLOSE:
            depth -= 1
        if current is None:
            return prefix
        current += ch
        i += 1
    if current is None:
        return prefix
    units.append(current.rstrip("-"))
    return units[0] + "".join(unit if unit.startswith(_OPEN) else f"({unit})" for unit in units[1:])
