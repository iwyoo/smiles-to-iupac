"""Unbranched phosphoric-acid chains hanging off an ester oxygen (P-106.2, P-106.3.2, P-106.3.5; P-67.2.2).

Each phosphorus carries =O/=S and OH/SH; a bridging -O-, -S-, -NH- or -CH2- joins the next one. Replacement
prefixes (thio, imido, carba) are numbered over the chain atoms P1, X2, P3, ... from the esterified phosphorus.
"""

from dataclasses import dataclass, field

from ._common import UnsupportedStructure
from ._numerals import numerical_term

_BRIDGE_PREFIX = {16: "thio", 7: "imido", 6: "carba"}


@dataclass
class PhosphorusChain:
    phosphorus: list = field(default_factory=list)
    bridges: list = field(default_factory=list)
    sulfur_per_phosphorus: list = field(default_factory=list)
    hydrogens: int = 0
    atoms: set = field(default_factory=set)
    extra: tuple = None

    @property
    def length(self):
        return len(self.phosphorus)

    @property
    def plain_monophosphate(self):
        return self.length == 1 and not self.sulfur_per_phosphorus[0] and self.extra is None


def _bond_order(mol, a, b):
    return mol.GetBondBetweenAtoms(a, b).GetBondTypeAsDouble()


def _terminal_chalcogen(mol, idx):
    atom = mol.GetAtomWithIdx(idx)
    return atom.GetAtomicNum() in (8, 16) and atom.GetDegree() == 1 and not atom.GetFormalCharge()


def _bridge_kind(mol, graph, idx, came_from):
    """('O'|'S'|'N'|'C', next atom) for a bridging atom between two phosphorus atoms, else None."""
    atom = mol.GetAtomWithIdx(idx)
    if atom.GetFormalCharge() or atom.GetIsotope() or atom.GetDegree() != 2:
        return None
    onward = [n for n in graph[idx] if n != came_from]
    if len(onward) != 1 or mol.GetAtomWithIdx(onward[0]).GetAtomicNum() != 15:
        return None
    hydrogens = atom.GetTotalNumHs()
    expected = {8: 0, 16: 0, 7: 1, 6: 2}.get(atom.GetAtomicNum())
    if expected is None or hydrogens != expected:
        return None
    return atom.GetAtomicNum(), onward[0]


def parse_chain(mol, graph, root, parent, allow_extra=False):
    """The chain starting at phosphorus `root` bonded to `parent`; with `allow_extra`, a final
    phosphorus may carry one more -O-R group, returned as `extra` = (oxygen, R atom)."""
    chain = PhosphorusChain()
    p, came_from = root, parent
    while True:
        atom = mol.GetAtomWithIdx(p)
        if (
            atom.GetAtomicNum() != 15
            or atom.GetDegree() != 4
            or atom.GetFormalCharge()
            or atom.GetIsotope()
            or atom.GetTotalNumHs()
        ):
            raise UnsupportedStructure("unsupported phosphorus atom in a phosphate chain")
        chain.phosphorus.append(p)
        chain.atoms.add(p)
        oxo, terminals, onward, extra = [], [], [], []
        for n in graph[p]:
            if n == came_from:
                continue
            na = mol.GetAtomWithIdx(n)
            order = _bond_order(mol, p, n)
            if _terminal_chalcogen(mol, n):
                if order == 2.0:
                    oxo.append(n)
                elif na.GetTotalNumHs() == 1:
                    terminals.append(n)
                else:
                    raise UnsupportedStructure("charged or substituted terminal group on a phosphate chain")
            elif order == 1.0 and _bridge_kind(mol, graph, n, p) is not None:
                onward.append(n)
            elif order == 1.0 and allow_extra and na.GetAtomicNum() == 8 and na.GetDegree() == 2:
                extra.append(n)
            else:
                raise UnsupportedStructure("unsupported group on a phosphate chain")
        if len(oxo) != 1 or len(extra) > 1 or len(onward) > 1:
            raise UnsupportedStructure("unsupported phosphorus substitution pattern")
        if extra and (onward or chain.extra):
            raise UnsupportedStructure("a phosphate chain with an extra ester group inside it")
        chain.atoms.update(oxo + terminals)
        chain.sulfur_per_phosphorus.append(
            sum(1 for n in oxo + terminals if mol.GetAtomWithIdx(n).GetAtomicNum() == 16)
        )
        chain.hydrogens += len(terminals)
        if extra:
            oxygen = extra[0]
            (r_atom,) = [n for n in graph[oxygen] if n != p]
            chain.extra = (oxygen, r_atom)
            chain.atoms.add(oxygen)
            if len(terminals) != 1:
                raise UnsupportedStructure("unsupported terminal pattern on a phosphate chain")
            return chain
        if onward:
            if len(terminals) != 1:
                raise UnsupportedStructure("unsupported chain phosphorus")
            bridge = onward[0]
            kind, p_next = _bridge_kind(mol, graph, bridge, p)
            chain.bridges.append(kind)
            chain.atoms.add(bridge)
            came_from, p = bridge, p_next
            if p_next in chain.atoms:
                raise UnsupportedStructure("cyclic phosphate chain")
            continue
        if len(terminals) != 2:
            raise UnsupportedStructure("unsupported terminal pattern on a phosphate chain")
        return chain


def _replacements(chain):
    found = []
    for i, count in enumerate(chain.sulfur_per_phosphorus):
        found += [(2 * i + 1, "thio")] * count
    for j, kind in enumerate(chain.bridges):
        if kind != 8:
            found.append((2 * j + 2, _BRIDGE_PREFIX[kind]))
    return sorted(found)


def _replacement_text(found):
    by_prefix = {}
    for locant, prefix in found:
        by_prefix.setdefault(prefix, []).append(locant)
    parts = []
    for prefix in sorted(by_prefix):
        locants = sorted(by_prefix[prefix])
        multiplier = numerical_term(len(locants)) if len(locants) > 1 else ""
        parts.append(f"{','.join(map(str, locants))}-{multiplier}{prefix}")
    return "-".join(parts)


def chain_anion_name(chain):
    """'trihydrogen 2-thiodiphosphate' for the ester-anion word of the chain (P-106.2, P-106.3.2)."""
    hydrogen = ""
    if chain.hydrogens:
        hydrogen = (numerical_term(chain.hydrogens) if chain.hydrogens > 1 else "") + "hydrogen "
    if chain.length == 1:
        sulfur = chain.sulfur_per_phosphorus[0]
        if sulfur > 2:
            raise UnsupportedStructure("more than two sulfur atoms on one phosphorus")
        stem = {0: "phosphate", 1: "phosphorothioate", 2: "phosphorodithioate"}[sulfur]
        return hydrogen + stem
    prefix = _replacement_text(_replacements(chain))
    return f"{hydrogen}{prefix}{numerical_term(chain.length)}phosphate"
