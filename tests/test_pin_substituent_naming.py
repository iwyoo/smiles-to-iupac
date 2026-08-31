"""Guards against a recurring regression class: PR #237 fixed `name_branch`
always fixing a compound substituent's free valence at locant 1 (CAS-style
'1-methylethyl', '1,1-dimethylethyl', '1-methylpropyl', ...) instead of
numbering through a branch point per P-29.3.2.2 ('propan-2-yl',
'tert-butyl', 'butan-2-yl', ...). Several test files had accepted those
CAS-style names as "expected" for months before this was noticed (see the
PR #237 description). This scans every string literal in tests/ for that
exact discredited shape, so a future regression (in `name_branch` itself,
or in a new hand-written test) fails immediately instead of quietly
reintroducing pre-PIN names as "expected" output.
"""

import io
import pathlib
import re
import tokenize

import pytest

_STEM = r"(?:meth|eth|prop|but|pent|hex|hept|oct|non|dec)"
_BAD_COMPOUND_SUBSTITUENT_RE = re.compile(
    rf"\b1-{_STEM}yl{_STEM}yl\b|\b1,1-di{_STEM}yl{_STEM}yl\b"
)

_TESTS_DIR = pathlib.Path(__file__).parent


def _string_literals(path):
    with open(path, "rb") as f:
        tokens = tokenize.tokenize(f.readline)
        for tok in tokens:
            if tok.type == tokenize.STRING:
                yield tok.start[0], tok.string


def _violations():
    for path in sorted(_TESTS_DIR.glob("*.py")):
        if path == pathlib.Path(__file__):
            continue
        for lineno, literal in _string_literals(path):
            if _BAD_COMPOUND_SUBSTITUENT_RE.search(literal):
                yield f"{path.name}:{lineno}: {literal}"


def test_no_pre_pin_branch_point_substituent_names():
    # P-29.3.2.2: a compound substituent whose own free valence is a
    # branch point must be numbered through that point, never cited as
    # locant 1 of a chain fused directly onto another alkyl word (that
    # shape -- e.g. '1-methylethyl', '1,1-dimethylethyl', '1-methylpropyl'
    # -- is always the discredited pre-PIN CAS style; the PIN form always
    # reads '<alkane stem>-<N>-yl' or, for tert-butyl, the P-29.6.1
    # retained name).
    violations = list(_violations())
    assert not violations, (
        "found pre-PIN CAS-style compound substituent name(s) in test "
        "expected values -- renumber through the branch point per "
        "P-29.3.2.2 instead:\n" + "\n".join(violations)
    )
