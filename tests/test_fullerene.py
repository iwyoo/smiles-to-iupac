from smiles_to_iupac import smiles_to_iupac
from smiles_to_iupac._fullerene import (
    _C84_D6H_ISOMER_24_SMILES,
    _CYCLOPROPA_C60_SMILES,
    _FULLERENE_C60_SMILES,
    _FULLERENE_C70_SMILES,
    _FULLERENE_C76_SMILES,
    _SILA_C60_SMILES,
)


def test_buckminsterfullerene():
    # C60, cross-checked against PubChem CID 123591's own connectivity
    # SMILES: 60 all-carbon atoms, every atom degree 3, ring perception of
    # exactly 12 five-membered and 20 six-membered rings.
    assert smiles_to_iupac(_FULLERENE_C60_SMILES) == "[60]fullerene"


def test_c70_fullerene():
    # C70-D5h(6), cross-checked against PubChem CID 16131935's own
    # connectivity SMILES: 70 all-carbon atoms, every atom degree 3, ring
    # perception of exactly 12 five-membered and 25 six-membered rings.
    assert smiles_to_iupac(_FULLERENE_C70_SMILES) == "(C70-D5h(6))[5,6]fullerene"


def test_c76_fullerene():
    # C76-D2, cross-checked against PubChem CID 56846604's own
    # connectivity SMILES: 76 all-carbon atoms, every atom degree 3, ring
    # perception of exactly 12 five-membered and 28 six-membered rings.
    assert smiles_to_iupac(_FULLERENE_C76_SMILES) == "(C76-D2)[5,6]fullerene"


def test_sila_c60_fullerene():
    # sila(C60-Ih), cross-checked against PubChem CID 101063510's own
    # connectivity SMILES: 60 skeletal atoms (59 carbon + 1 silicon),
    # every carbon degree 3, silicon degree 3 with 0 H and 0 radical
    # electrons, ring perception of exactly 12 five-membered and 20
    # six-membered rings -- P-27.5.1's own worked example cites no locant
    # since every C60-Ih vertex is symmetry-equivalent.
    assert smiles_to_iupac(_SILA_C60_SMILES) == "sila(C60-Ih)[5,6]fullerene"


def test_cyclopropa_c60_fullerene():
    # 3'H-cyclopropa[1,9](C60-Ih), cross-checked against PubChem CID
    # 11422743's own connectivity SMILES: a single methylene fused across
    # an intact 6,6-bond of the C60-Ih cage (formula C61H2) -- the two
    # bridgehead cage atoms remain directly bonded to each other, forming
    # a fused cyclopropane, matching P-27.6.1's own literal worked example
    # (not P-27.4.1's homofullerene, which has no such direct bond --
    # see _fullerene.py's module docstring). C60-Ih's 6,6-bonds are a
    # single symmetry orbit, so the locant is fixed regardless of which
    # bond was used.
    assert smiles_to_iupac(_CYCLOPROPA_C60_SMILES) == "3'H-cyclopropa[1,9](C60-Ih)[5,6]fullerene"


def test_c84_isomer_24_via_spiral_match():
    # (C84-D6h(24)), a real structure (PubChem CID 133108900, InChIKey
    # FQRWAZOLUJHNDT-UHFFFAOYSA-N) not in _fullerene.py's exact-match table:
    # identified instead via _fullerene_spiral.py's general ring-spiral
    # algorithm against the 24 known C84 IPR isomers -- the #997 "isomer
    # atlas" mechanism. Independently confirmed by a real paper's title,
    # "A minor isomer of C84 fullerene, D6h-C84(24)" (ScienceDirect), naming
    # this same isomer the same way.
    assert smiles_to_iupac(_C84_D6H_ISOMER_24_SMILES) == "(C84-D6h(24))[5,6]fullerene"


def test_benzene_still_resolves():
    # a sanity check that the fullerene shape check (60 atoms, all
    # degree-3 carbon) doesn't misfire on an ordinary small ring.
    assert smiles_to_iupac("c1ccccc1") == "benzene"


def test_smaller_cage_is_not_fullerene():
    # pentaprismane: a smaller all-carbon cage-like polycyclic (not 60
    # atoms) must not match the fullerene shape check; it's correctly
    # resolved via _polycyclic.py's hexacyclic support instead (see
    # tests/test_hexacyclic.py).
    assert smiles_to_iupac("C12C3C4C1C1C2C2C3C4C12") == "hexacyclo[4.4.0.0^2,5.0^3,9.0^4,8.0^7,10]decane"
