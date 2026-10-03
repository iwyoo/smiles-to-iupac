"""Names for which the Blue Book defines no preferred IUPAC name (P-69.0:
Group 1-12 organometallics other than ocenes) are returned with a warning."""

import threading
import warnings

_state = threading.local()


class NonPreferredNameWarning(UserWarning):
    pass


def begin():
    _state.reasons = []


def mark(name, reason):
    if hasattr(_state, "reasons"):
        _state.reasons.append(reason)
    return name


def finish(name):
    reasons = getattr(_state, "reasons", [])
    _state.reasons = []
    if reasons:
        warnings.warn(
            f"{name!r} is a valid name but not a preferred IUPAC name: {reasons[0]}",
            NonPreferredNameWarning,
            stacklevel=3,
        )
    return name
