"""R'-O- substituent prefixes (P-63.2.2.1.1, P-63.2.2.2). The retained methoxy, ethoxy, propoxy, butoxy, tert-butoxy and
phenoxy are simple prefixes. Every other prefix is a compound prefix: a substituted retained one ('2,2-dimethylpropoxy',
'cyclohexylmethoxy', '4-chlorophenoxy'), a glycosyl one ('beta-D-glucopyranosyloxy'), or R' followed by 'oxy' with a
compound R' enclosed ('(propan-2-yl)oxy', 'pentyloxy')."""

import re

from ._multiplicative_text import enclose

_RETAINED = {
    "methyl": "methoxy",
    "ethyl": "ethoxy",
    "propyl": "propoxy",
    "butyl": "butoxy",
    "tert-butyl": "tert-butoxy",
    "phenyl": "phenoxy",
}
_ISOTOPE_ONLY = re.compile(r"(\(\d[^()-]*\))([a-z-]+)")
_SUBSTITUTED_END = re.compile(r"(?:methyl|ethyl|propyl|butyl)$")
_NOT_SUBSTITUTED_END = ("cyclopropyl", "cyclobutyl", "tert-butyl")


def alkoxy_prefix(rname, compound=False):
    """(prefix text, is_compound) for the group R'-O- whose R' substituent prefix is `rname`."""
    if rname in _RETAINED:
        return _RETAINED[rname], False
    modified = _ISOTOPE_ONLY.fullmatch(rname)
    if modified and modified.group(2) in _RETAINED:
        return modified.group(1) + _RETAINED[modified.group(2)], False
    if rname.endswith("phenyl"):
        return rname[: -len("phenyl")] + "phenoxy", True
    if _SUBSTITUTED_END.search(rname) and not rname.endswith(_NOT_SUBSTITUTED_END):
        return rname[:-2] + "oxy", True
    if rname.endswith("osyl"):
        return rname + "oxy", True
    if compound or rname[0].isdigit() or rname[0] == "(":
        return enclose(rname) + "oxy", True
    return rname + "oxy", True
