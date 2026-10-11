"""A fullerene ring system as a component of a spiro union (P-27.6.3, P-24.5): a fullerene cannot be a spiro component
itself, but a homofullerene can, the atom inserted by 'homo' being the spiro atom, and so can a fullerene with a fused
ring, the spiro atom being an atom of that ring. The component is cited in square brackets when its name starts with
locants (`[1(9)a]homo(C60-Ih)[5,6]fullerene`) and the spiro atom has the locant of the homo atom (1a) or of the ring atom."""

from ._fullerene_attached import attached_reading
from ._fullerene_modified import modified_candidate

_MIN_ATOMS = 31
_MIN_RINGS = 11


def fullerene_component(mol, comp, spiro_atoms):
    """The spiro component dict of a modified fullerene ring system, with rings fused to it, whose spiro atoms are its
    homo atoms or atoms of its fused rings; else None."""
    if comp["rings"] < _MIN_RINGS or len(comp["atoms"]) < _MIN_ATOMS or not spiro_atoms:
        return None
    bare = frozenset(spiro_atoms)
    plain = modified_candidate(mol, comp["atoms"], bare=bare)
    if plain is not None and not (plain.edit.seco or plain.edit.cyclo):
        if plain.edit.saturated or not bare <= {x for x, _ in plain.edit.homo_atoms}:
            return None
        text = plain.text()
        locants = [candidate.atom_locants(n) for candidate in plain.ties for n in candidate.numberings]
    else:
        reading = attached_reading(mol, comp["atoms"], bare)
        if reading is None or reading.bridged is not None or reading.found.edit.saturated:
            return None
        ring_atoms = set().union(*(ring.atoms for ring in reading.rings)) if reading.rings else set()
        homo = {x for x, _ in reading.found.edit.homo_atoms}
        if not bare <= ring_atoms | homo:
            return None
        hydro, indicated, _ = reading.render(reading.found.numbering, reading.found.edit)
        if hydro or indicated:
            return None
        text = reading.text()
        locants = reading.locants()
    return {
        "name": _bracketed(text),
        "numberings": [(locant_of, ()) for locant_of in locants],
        "replacement": False,
        "fused": True,
    }


def _bracketed(text):
    """A name that starts with locants is cited with them in square brackets: 1(9)a-homo(...) -> [1(9)a]homo(...)."""
    if text[:1].isdigit() and "-" in text:
        locants, rest = text.split("-", 1)
        return f"[{locants}]{rest}"
    return text
