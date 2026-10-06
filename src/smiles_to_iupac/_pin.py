"""Names for which the Blue Book defines no preferred IUPAC name (P-69.0:
Group 1-12 organometallics other than ocenes) are returned with a warning."""

import threading
import warnings

_state = threading.local()


class NonPreferredNameWarning(UserWarning):
    pass


def enter():
    depth = getattr(_state, "depth", 0)
    if depth == 0:
        _state.reasons = []
    _state.depth = depth + 1


def mark(name, reason):
    if getattr(_state, "depth", 0):
        _state.reasons.append(reason)
    return name


def nested():
    return getattr(_state, "depth", 0) > 0


def reason_count():
    return len(_state.reasons)


def reasons_since(start):
    return _state.reasons[start:]


def replay(reasons):
    _state.reasons.extend(reasons)


def leave(name):
    _state.depth -= 1
    if _state.depth:
        return
    reasons, _state.reasons = _state.reasons, []
    if name is not None and reasons:
        warnings.warn(
            f"{name!r} is a valid name but not a preferred IUPAC name: {reasons[0]}",
            NonPreferredNameWarning,
            stacklevel=3,
        )
