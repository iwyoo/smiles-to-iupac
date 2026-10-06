"""Table 10.1 stereoparents (P-101.2.7, Appendix 3): one configured SMILES per parent whose atom maps are the drawn
locants (n, +100 per prime, +1000 for a letter a-c, +10000 per superscript digit), the reference ring and drawing
sense of its α/β frame, and the steroid anchor. Atoms the drawing numbers only by element letter map to 90001 + the
position of their locant in UNNUMBERED. `scope` is "both", or the one engine that names the parent."""

from __future__ import annotations

from typing import NamedTuple

UNNUMBERED = (
    "A15c",
    "A15d",
    "A15m",
    "A15o",
    "A1c",
    "A1d",
    "A1m",
    "A1o",
    "A2c",
    "A2d",
    "A2m",
    "A2o",
    "A6c",
    "A6d",
    "A6m",
    "A6o",
    "A9c",
    "A9d",
    "A9m",
    "A9o",
    "Ac",
    "Ad",
    "Am",
    "Cm",
    "Oa",
    "Oa′",
    "Ob",
    "Oc",
    "Od",
    "Oe",
    "Oh",
    "Oh′",
    "Ok",
    "Oq",
    "Ox",
    "Ox′",
    "Oy",
    "Oy′",
    "Oz",
)
_SUPERSCRIPTS = "¹²³"
_PRIMES = "′″"


class Stereoparent(NamedTuple):
    smiles: str
    ref: str = ""
    anticlockwise: bool | None = None
    anchor: str = ""
    scope: str = "both"


def locant(map_number):
    if map_number >= 90000:
        return UNNUMBERED[map_number - 90001]
    number, letter = map_number % 100, (map_number // 1000) % 10
    superscript, prime = map_number // 10000, (map_number // 100) % 10
    return (
        f"{number}{chr(96 + letter) if letter else ''}"
        f"{_SUPERSCRIPTS[superscript - 1] if superscript else ''}{_PRIMES[prime - 1] if prime else ''}"
    )


STEREOPARENTS = {
    "18-oxayohimban": Stereoparent(
        "[cH:11]1[cH:10][cH:9][c:8]2[c:7]3[c:2]([nH:1][c:13]2[cH:12]1)[C@@H:3]1[CH2:14][C@@H:15]2[CH2:16][CH2:17][O:18][CH2:19][C@H:20]2[CH2:21][N:4]1[CH2:5][CH2:6]3",
        "1,2,7,8,13",
        True
    ),
    "2,2′-neolignane": Stereoparent(
        "[c:1]1([CH2:7][CH2:8][CH3:9])[c:2](-[c:102]2[c:101]([CH2:107][CH2:108][CH3:109])[cH:106][cH:105][cH:104][cH:103]2)[cH:3][cH:4][cH:5][cH:6]1",
        "",
        None,
        "",
        "table10.1"
    ),
    "2,3′-neolignane": Stereoparent(
        "[c:1]1([CH2:7][CH2:8][CH3:9])[c:2](-[c:103]2[cH:102][c:101]([CH2:107][CH2:108][CH3:109])[cH:106][cH:105][cH:104]2)[cH:3][cH:4][cH:5][cH:6]1",
        "",
        None,
        "",
        "table10.1"
    ),
    "2,4′-neolignane": Stereoparent(
        "[c:1]1([CH2:7][CH2:8][CH3:9])[c:2](-[c:104]2[cH:103][cH:102][c:101]([CH2:107][CH2:108][CH3:109])[cH:106][cH:105]2)[cH:3][cH:4][cH:5][cH:6]1",
        "",
        None,
        "",
        "table10.1"
    ),
    "2,7′-neolignane": Stereoparent(
        "[c:1]1([CH2:7][CH2:8][CH3:9])[c:2]([CH:107]([c:101]2[cH:102][cH:103][cH:104][cH:105][cH:106]2)[CH2:108][CH3:109])[cH:3][cH:4][cH:5][cH:6]1",
        "",
        None,
        "",
        "table10.1"
    ),
    "2,8′-neolignane": Stereoparent(
        "[c:1]1([CH2:7][CH2:8][CH3:9])[c:2]([CH:108]([CH2:107][c:101]2[cH:102][cH:103][cH:104][cH:105][cH:106]2)[CH3:109])[cH:3][cH:4][cH:5][cH:6]1",
        "",
        None,
        "",
        "table10.1"
    ),
    "2,9′-neolignane": Stereoparent(
        "[c:1]1([CH2:7][CH2:8][CH3:9])[c:2]([CH2:109][CH2:108][CH2:107][c:101]2[cH:102][cH:103][cH:104][cH:105][cH:106]2)[cH:3][cH:4][cH:5][cH:6]1",
        "",
        None,
        "",
        "table10.1"
    ),
    "21H,23H-porphyrin": Stereoparent(
        "[c:1]12[cH:2][cH:3][c:4]([cH:5][c:6]3[n:22][c:9]([cH:10][c:11]4[cH:12][cH:13][c:14]([cH:15][c:16]5[n:24][c:19]([cH:20]1)[CH:18]=[CH:17]5)[nH:23]4)[CH:8]=[CH:7]3)[nH:21]2",
        "",
        None,
        "",
        "appendix3"
    ),
    "21H-biline": Stereoparent(
        "[cH:1]1[cH:2][cH:3][c:4]([CH:5]=[C:6]2[CH:7]=[CH:8][C:9]([CH:10]=[C:11]3[CH:12]=[CH:13][C:14]([CH:15]=[C:16]4[CH:17]=[CH:18][CH:19]=[N:24]4)=[N:23]3)=[N:22]2)[nH:21]1"
    ),
    "3,3′-neolignane": Stereoparent(
        "[c:1]1([CH2:7][CH2:8][CH3:9])[cH:2][c:3](-[c:103]2[cH:102][c:101]([CH2:107][CH2:108][CH3:109])[cH:106][cH:105][cH:104]2)[cH:4][cH:5][cH:6]1"
    ),
    "3,4′-neolignane": Stereoparent(
        "[c:1]1([CH2:7][CH2:8][CH3:9])[cH:2][c:3](-[c:104]2[cH:103][cH:102][c:101]([CH2:107][CH2:108][CH3:109])[cH:106][cH:105]2)[cH:4][cH:5][cH:6]1",
        "",
        None,
        "",
        "table10.1"
    ),
    "3,7′-neolignane": Stereoparent(
        "[c:1]1([CH2:7][CH2:8][CH3:9])[cH:2][c:3]([CH:107]([c:101]2[cH:102][cH:103][cH:104][cH:105][cH:106]2)[CH2:108][CH3:109])[cH:4][cH:5][cH:6]1",
        "",
        None,
        "",
        "table10.1"
    ),
    "3,8′-neolignane": Stereoparent(
        "[c:1]1([CH2:7][CH2:8][CH3:9])[cH:2][c:3]([CH:108]([CH2:107][c:101]2[cH:102][cH:103][cH:104][cH:105][cH:106]2)[CH3:109])[cH:4][cH:5][cH:6]1",
        "",
        None,
        "",
        "table10.1"
    ),
    "3,9′-neolignane": Stereoparent(
        "[c:1]1([CH2:7][CH2:8][CH3:9])[cH:2][c:3]([CH2:109][CH2:108][CH2:107][c:101]2[cH:102][cH:103][cH:104][cH:105][cH:106]2)[cH:4][cH:5][cH:6]1",
        "",
        None,
        "",
        "table10.1"
    ),
    "4,4′-neolignane": Stereoparent(
        "[c:1]1([CH2:7][CH2:8][CH3:9])[cH:2][cH:3][c:4](-[c:104]2[cH:103][cH:102][c:101]([CH2:107][CH2:108][CH3:109])[cH:106][cH:105]2)[cH:5][cH:6]1",
        "",
        None,
        "",
        "table10.1"
    ),
    "4,7′-neolignane": Stereoparent(
        "[c:1]1([CH2:7][CH2:8][CH3:9])[cH:2][cH:3][c:4]([CH:107]([c:101]2[cH:102][cH:103][cH:104][cH:105][cH:106]2)[CH2:108][CH3:109])[cH:5][cH:6]1",
        "",
        None,
        "",
        "table10.1"
    ),
    "4,8′-neolignane": Stereoparent(
        "[c:1]1([CH2:7][CH2:8][CH3:9])[cH:2][cH:3][c:4]([CH:108]([CH2:107][c:101]2[cH:102][cH:103][cH:104][cH:105][cH:106]2)[CH3:109])[cH:5][cH:6]1",
        "",
        None,
        "",
        "table10.1"
    ),
    "4,9′-neolignane": Stereoparent(
        "[c:1]1([CH2:7][CH2:8][CH3:9])[cH:2][cH:3][c:4]([CH2:109][CH2:108][CH2:107][c:101]2[cH:102][cH:103][cH:104][cH:105][cH:106]2)[cH:5][cH:6]1",
        "",
        None,
        "",
        "table10.1"
    ),
    "7,7′-neolignane": Stereoparent(
        "[c:1]1([CH:7]([CH2:8][CH3:9])[CH:107]([c:101]2[cH:102][cH:103][cH:104][cH:105][cH:106]2)[CH2:108][CH3:109])[cH:2][cH:3][cH:4][cH:5][cH:6]1",
        "",
        None,
        "",
        "table10.1"
    ),
    "7,8′-neolignane": Stereoparent(
        "[c:1]1([CH:7]([CH2:8][CH3:9])[CH:108]([CH2:107][c:101]2[cH:102][cH:103][cH:104][cH:105][cH:106]2)[CH3:109])[cH:2][cH:3][cH:4][cH:5][cH:6]1",
        "",
        None,
        "",
        "table10.1"
    ),
    "7,9′-neolignane": Stereoparent(
        "[c:1]1([CH:7]([CH2:8][CH3:9])[CH2:109][CH2:108][CH2:107][c:101]2[cH:102][cH:103][cH:104][cH:105][cH:106]2)[cH:2][cH:3][cH:4][cH:5][cH:6]1",
        "",
        None,
        "",
        "table10.1"
    ),
    "8,9′-neolignane": Stereoparent(
        "[c:1]1([CH2:7][CH:8]([CH3:9])[CH2:109][CH2:108][CH2:107][c:101]2[cH:102][cH:103][cH:104][cH:105][cH:106]2)[cH:2][cH:3][cH:4][cH:5][cH:6]1",
        "",
        None,
        "",
        "table10.1"
    ),
    "9,9′-neolignane": Stereoparent(
        "[c:1]1([CH2:7][CH2:8][CH2:9][CH2:109][CH2:108][CH2:107][c:101]2[cH:102][cH:103][cH:104][cH:105][cH:106]2)[cH:2][cH:3][cH:4][cH:5][cH:6]1",
        "",
        None,
        "",
        "table10.1"
    ),
    "abietane": Stereoparent(
        "[CH3:16][CH:15]([CH3:17])[C@H:13]1[CH2:12][CH2:11][C@H:9]2[C@@H:8]([CH2:7][CH2:6][C@H:5]3[C:4]([CH3:18])([CH3:19])[CH2:3][CH2:2][CH2:1][C@:10]23[CH3:20])[CH2:14]1",
        "1,2,3,4,5,10",
        True,
        "5:a"
    ),
    "aconitane": Stereoparent(
        "[CH3:18][C@:4]12[CH2:3][CH2:2][CH2:1][C@:11]34[CH:17]([NH:20][CH2:19]1)[C@H:7]([CH2:6][C@H:5]23)[C@H:8]1[CH2:15][CH2:16][C@@H:13]2[CH2:14][C@H:9]1[C@H:10]4[CH2:12]2"
    ),
    "ajmalan": Stereoparent(
        "[CH3:18][CH2:19][C@@H:20]1[CH2:21][N:4]2[C@H:5]3[CH2:6][C@:7]45[CH2:17][CH:16]3[C@H:15]1[CH2:14][C@H:3]2[C@@H:2]4[N:1]([CH3:22])[c:13]1[cH:12][cH:11][cH:10][cH:9][c:8]15"
    ),
    "akuammilan": Stereoparent(
        "[CH3:18]/[CH:19]=[C:20]1/[CH2:21][N:4]2[CH2:5][CH2:6][C@:7]34[C:2](=[N:1][c:13]5[cH:12][cH:11][cH:10][cH:9][c:8]53)[C@@H:3]2[CH2:14][C@@H:15]1[CH:16]4[CH3:17]"
    ),
    "alstophyllan": Stereoparent(
        "[CH3:18][CH2:19][C:20]1=[CH:21][O:24][CH2:17][C@@H:16]2[C@H:15]1[CH2:14][C@H:3]1[c:2]3[c:7]([c:8]4[cH:9][cH:10][cH:11][cH:12][c:13]4[n:1]3[CH3:22])[CH2:6][C@@H:5]2[N:4]1[CH3:23]"
    ),
    "ambrosane": Stereoparent(
        "[CH3:12][CH:11]([CH3:13])[C@@H:7]1[CH2:8][CH2:9][C@H:10]([CH3:14])[C@@H:1]2[CH2:2][CH2:3][CH2:4][C@@:5]2([CH3:15])[CH2:6]1",
        "1,2,3,4,5",
        True,
        "5:b"
    ),
    "androstane": Stereoparent(
        "[CH3:18][C@@:13]12[CH2:17][CH2:16][CH2:15][C@H:14]1[C@@H:8]1[CH2:7][CH2:6][CH:5]3[CH2:4][CH2:3][CH2:2][CH2:1][C@:10]3([CH3:19])[C@H:9]1[CH2:11][CH2:12]2",
        "",
        None,
        "8:b",
        "table10.1"
    ),
    "aporphine": Stereoparent(
        "[cH:1]1[cH:2][cH:3][c:1003]2[c:3011]3[c:2011]1-[c:1011]1[cH:11][cH:10][cH:9][cH:8][c:1007]1[CH2:7][C@@H:1006]3[N:6]([CH3:12])[CH2:5][CH2:4]2",
        "1,2,3,3a,11c,11b",
        False
    ),
    "aristolane": Stereoparent(
        "[CH3:14][C@@H:4]1[CH2:3][CH2:2][CH2:1][C@@H:10]2[CH2:9][CH2:8][C@@H:7]3[C@@H:6]([C:11]3([CH3:12])[CH3:13])[C@@:5]12[CH3:15]",
        "1,2,3,4,5,10",
        True,
        "5:a"
    ),
    "aspidofractinine": Stereoparent(
        "[cH:16]1[cH:15][cH:14][c:13]2[c:18]([cH:17]1)[NH:1][C:2]13[CH2:3][CH2:4][C:5]4([CH2:6][CH2:7][CH2:8][N:9]5[CH2:10][CH2:11][C@:12]21[C@H:19]45)[CH2:20][CH2:21]3"
    ),
    "aspidospermidine": Stereoparent(
        "[CH3:21][CH2:20][C@@:5]12[CH2:6][CH2:7][CH2:8][N:9]3[CH2:10][CH2:11][C@:12]4([c:13]5[cH:14][cH:15][cH:16][cH:17][c:18]5[NH:1][C@@H:2]4[CH2:3][CH2:4]1)[C@@H:19]23",
        "1,2,12,13,18",
        True,
        "19:b"
    ),
    "atidane": Stereoparent(
        "[CH3:17][C@H:16]1[CH2:15][C@:8]23[CH2:14][CH2:13][C@H:12]1[CH2:11][C@H:9]2[C@@:10]12[CH2:1][CH2:2][CH2:3][C@@:4]([CH3:18])([CH2:19][NH:21][CH2:20]1)[C@H:5]2[CH2:6][CH2:7]3"
    ),
    "atisane": Stereoparent(
        "[CH3:17][C@@H:16]1[CH2:15][C@@:8]23[CH2:14][CH2:13][C@@H:12]1[CH2:11][C@@H:9]2[C@@:10]1([CH3:20])[CH2:1][CH2:2][CH2:3][C:4]([CH3:18])([CH3:19])[C@@H:5]1[CH2:6][CH2:7]3",
        "",
        None,
        "5:a"
    ),
    "atisine": Stereoparent(
        "[CH2:17]=[C:16]1[C@H:12]2[CH2:13][CH2:14][C@@:8]3([CH2:7][CH2:6][C@@H:5]4[C@@:4]5([CH3:18])[CH2:3][CH2:2][CH2:1][C@@:10]4([C@@H:20]4[O:24][CH2:23][CH2:22][N:21]4[CH2:19]5)[C@@H:9]3[CH2:11]2)[C@@H:15]1[OH:90031]"
    ),
    "berbaman": Stereoparent(
        "[cH:13]1[cH:14][c:9]2[cH:10][c:11]([cH:12]1)[O:90027][c:112]1[cH:111][cH:110][c:109]([cH:114][cH:113]1)[CH2:115][C@@H:101]1[NH:102][CH2:103][CH2:104][c:1104]3[cH:105][cH:106][c:107]([cH:108][c:1108]31)[O:90025][c:8]1[cH:7][cH:6][cH:5][c:1004]3[c:1008]1[C@@H:1]([CH2:15]2)[NH:2][CH2:3][CH2:4]3"
    ),
    "berbine": Stereoparent(
        "[cH:1]1[cH:2][cH:3][cH:4][c:1004]2[c:2013]1[C@H:1013]1[N:7]([CH2:6][CH2:5]2)[CH2:8][c:1008]2[cH:9][cH:10][cH:11][cH:12][c:1012]2[CH2:13]1",
        "1,2,3,4,4a,13b",
        False
    ),
    "beyerane": Stereoparent(
        "[CH3:18][C:4]1([CH3:19])[CH2:3][CH2:2][CH2:1][C@:10]2([CH3:20])[C@@H:5]1[CH2:6][CH2:7][C@@:8]13[CH2:15][CH2:16][C@@:13]([CH3:17])([CH2:12][CH2:11][C@@H:9]21)[CH2:14]3",
        "",
        None,
        "5:b"
    ),
    "bisabolane": Stereoparent(
        "[CH3:12][CH:11]([CH3:13])[CH2:7][CH2:8][CH2:9][CH:10]([CH3:15])[CH:1]1[CH2:2][CH2:3][CH:4]([CH3:14])[CH2:5][CH2:6]1"
    ),
    "bornane": Stereoparent(
        "[CH3:10][C:1]12[CH2:2][CH2:3][CH:4]([CH2:5][CH2:6]1)[C:7]2([CH3:8])[CH3:9]",
        "",
        None,
        "",
        "table10.1"
    ),
    "bufanolide": Stereoparent(
        "[CH3:18][C@:13]12[CH2:12][CH2:11][C@H:9]3[C@@H:8]([CH2:7][CH2:6][CH:5]4[CH2:4][CH2:3][CH2:2][CH2:1][C@:10]34[CH3:19])[C@@H:14]1[CH2:15][CH2:16][C@@H:17]2[C@H:20]1[CH2:22][CH2:23][C:24](=[O:90027])[O:90025][CH2:21]1",
        "1,2,3,4,5,10",
        True,
        "8:b"
    ),
    "cadinane": Stereoparent(
        "[CH3:12][CH:11]([CH3:13])[C@@H:7]1[CH2:8][CH2:9][C@H:10]([CH3:15])[C@@H:1]2[CH2:2][CH2:3][C@H:4]([CH3:14])[CH2:5][C@@H:6]12",
        "1,2,3,4,5,6",
        True,
        "6:b"
    ),
    "campestane": Stereoparent(
        "[CH3:26][CH:25]([CH3:27])[C@H:24]([CH3:10024])[CH2:23][CH2:22][C@@H:20]([CH3:21])[C@H:17]1[CH2:16][CH2:15][C@H:14]2[C@@H:8]3[CH2:7][CH2:6][CH:5]4[CH2:4][CH2:3][CH2:2][CH2:1][C@:10]4([CH3:19])[C@H:9]3[CH2:11][CH2:12][C@:13]12[CH3:18]",
        "",
        None,
        "8:b",
        "table10.1"
    ),
    "carane": Stereoparent(
        "[CH3:10][CH:3]1[CH2:4][CH2:5][CH:6]2[CH:1]([CH2:2]1)[C:7]2([CH3:8])[CH3:9]",
        "",
        None,
        "",
        "table10.1"
    ),
    "cardanolide": Stereoparent(
        "[CH2:1]1[CH2:2][CH2:3][CH2:4][CH:5]2[CH2:6][CH2:7][C@@H:8]3[C@@H:9]([C@@:10]12[CH3:19])[CH2:11][CH2:12][C@@:13]1([CH3:18])[C@@H:14]3[CH2:15][CH2:16][C@@H:17]1[C@H:20]1[CH2:21][O:90025][C:23](=[O:90027])[CH2:22]1",
        "1,2,3,4,5,10",
        True,
        "8:b"
    ),
    "caryophyllane": Stereoparent(
        "[CH3:14][CH:3]1[CH2:4][CH2:5][CH2:6][CH:7]([CH3:15])[CH:8]2[CH2:10][C:11]([CH3:12])([CH3:13])[CH:9]2[CH2:1][CH2:2]1"
    ),
    "cedrane": Stereoparent(
        "[CH3:15][C@@H:8]1[CH2:9][CH2:10][C@@:1]23[CH2:11][C@@H:7]1[C:6]([CH3:13])([CH3:14])[C@@H:5]2[CH2:4][CH2:3][C@H:2]3[CH3:12]",
        "",
        None,
        "5:b"
    ),
    "cephalotaxine": Stereoparent(
        "[CH3:21][O:90025][C:2]1=[CH:1][C@:5]23[CH2:6][CH2:7][CH2:8][N:9]2[CH2:10][CH2:11][c:12]2[cH:17][c:16]4[c:15]([cH:14][c:13]2[C@@H:4]3[C@@H:3]1[OH:90031])[O:18][CH2:19][O:20]4"
    ),
    "cepham": Stereoparent(
        "[O:90025]=[C:8]1[CH2:7][C@H:6]2[S:1][CH2:2][CH2:3][CH2:4][N:5]12",
        "1,2,3,4,5,6",
        False,
        "6:a"
    ),
    "cevane": Stereoparent(
        "[CH3:27][C@H:25]1[CH2:24][CH2:23][C@H:22]2[C@H:20]([CH3:21])[C@H:17]3[CH2:16][CH2:15][C@@H:14]4[C@@H:12]([CH2:11][C@H:9]5[C@H:8]4[CH2:7][CH2:6][C@@H:5]4[CH2:4][CH2:3][CH2:2][CH2:1][C@:10]54[CH3:19])[C@@H:13]3[CH2:18][N:28]2[CH2:26]1",
        "1,2,3,4,5,10",
        True
    ),
    "chelidonine": Stereoparent(
        "[CH3:25][N:5]1[CH2:6][c:19]2[c:23]([cH:10][cH:9][c:8]3[c:7]2[O:20][CH2:21][O:22]3)[C@@H:13]2[C@H:14]1[c:18]1[cH:4][c:3]3[c:2]([cH:1][c:24]1[CH2:12][C@@H:11]2[OH:90031])[O:15][CH2:16][O:17]3",
        "1,2,3,4,18,24",
        False
    ),
    "cholane": Stereoparent(
        "[CH3:24][CH2:23][CH2:22][C@@H:20]([CH3:21])[C@H:17]1[CH2:16][CH2:15][C@H:14]2[C@@H:8]3[CH2:7][CH2:6][C@@H:5]4[CH2:4][CH2:3][CH2:2][CH2:1][C@:10]4([CH3:19])[C@H:9]3[CH2:11][CH2:12][C@:13]12[CH3:18]",
        "",
        None,
        "8:b",
        "table10.1"
    ),
    "cholestane": Stereoparent(
        "[CH3:26][CH:25]([CH3:27])[CH2:24][CH2:23][CH2:22][C@@H:20]([CH3:21])[C@H:17]1[CH2:16][CH2:15][C@H:14]2[C@@H:8]3[CH2:7][CH2:6][CH:5]4[CH2:4][CH2:3][CH2:2][CH2:1][C@:10]4([CH3:19])[C@H:9]3[CH2:11][CH2:12][C@:13]12[CH3:18]",
        "",
        None,
        "8:b",
        "table10.1"
    ),
    "cinchonan": Stereoparent(
        "[CH2:11]=[CH:10][C@H:3]1[CH2:2][N:1]2[CH2:6][CH2:5][C@H:4]1[CH2:7][C@@H:8]2[CH2:9][c:104]1[cH:103][cH:102][n:101][c:109]2[cH:108][cH:107][cH:106][cH:105][c:110]12"
    ),
    "conanine": Stereoparent(
        "[CH3:21][C@H:20]1[C@H:17]2[CH2:16][CH2:15][C@H:14]3[C@@H:8]4[CH2:7][CH2:6][CH:5]5[CH2:4][CH2:3][CH2:2][CH2:1][C@:10]5([CH3:19])[C@H:9]4[CH2:11][CH2:12][C@:13]23[CH2:18][N:22]1[CH3:23]",
        "",
        None,
        "10:b"
    ),
    "corrin": Stereoparent(
        "[CH:1]12[CH2:2][CH2:3][C:4](=[N:21]1)[CH:5]=[C:6]1[CH2:7][CH2:8][C:9](=[N:22]1)[CH:10]=[C:11]1[CH2:12][CH2:13][C:14](=[N:23]1)[CH:15]=[C:16]1[CH2:17][CH2:18][CH:19]2[NH:24]1"
    ),
    "corynan": Stereoparent(
        "[CH3:17][CH2:16][C@H:15]1[CH2:14][C@H:3]2[c:2]3[nH:1][c:13]4[cH:12][cH:11][cH:10][cH:9][c:8]4[c:7]3[CH2:6][CH2:5][N:4]2[CH2:21][C@@H:20]1[CH2:19][CH3:18]",
        "1,2,7,8,13",
        True
    ),
    "corynoxan": Stereoparent(
        "[CH3:18][CH2:19][C@@H:20]1[CH2:21][N:4]2[CH2:5][CH2:6][C@@:7]3([CH2:2][NH:1][c:13]4[cH:12][cH:11][cH:10][cH:9][c:8]43)[C@@H:3]2[CH2:14][C@@H:15]1[CH2:16][CH3:17]",
        "3,14,15,20,21,4",
        True
    ),
    "crinan": Stereoparent(
        "[cH:7]1[c:14]2[c:18]([cH:10][c:9]3[c:8]1[O:15][CH2:16][O:17]3)[C@:19]13[CH2:1][CH2:2][CH2:3][CH2:4][C@H:13]1[N:5]([CH2:12][CH2:11]3)[CH2:6]2"
    ),
    "curan": Stereoparent(
        "[CH3:18][CH2:19][C@@H:20]1[CH2:21][N:4]2[CH2:5][CH2:6][C@:7]34[c:8]5[cH:9][cH:10][cH:11][cH:12][c:13]5[NH:1][C@H:2]3[C@H:16]([CH3:17])[C@H:15]1[CH2:14][C@H:3]24"
    ),
    "dammarane": Stereoparent(
        "[CH3:26][CH:25]([CH3:27])[CH2:24][CH2:23][CH2:22][C@@H:20]([CH3:21])[C@H:17]1[CH2:16][CH2:15][C@:14]2([CH3:30])[C@@H:13]1[CH2:12][CH2:11][C@@H:9]1[C@@:10]3([CH3:19])[CH2:1][CH2:2][CH2:3][C:4]([CH3:28])([CH3:29])[C@@H:5]3[CH2:6][CH2:7][C@@:8]21[CH3:18]",
        "1,2,3,4,5,10",
        True,
        "5:a"
    ),
    "daphnane": Stereoparent(
        "[CH3:23][CH2:22][CH2:21][C@@:7]12[C@@H:2]3[C@@H:3]([CH:17]([CH3:18])[CH3:19])[CH2:4][CH2:5][C@@:6]1([CH3:20])[C@@H:15]1[CH2:14][CH2:13][C@@:12]4([CH2:11][CH2:10][CH2:9][C@@H:8]24)[N:1]3[CH2:16]1"
    ),
    "dendrobane": Stereoparent(
        "[CH3:16][CH:15]([CH3:17])[C@H:8]1[C@@H:7]2[CH2:12][O:13][C@H:9]1[C@H:10]1[N:1]([CH3:14])[CH2:2][C@H:3]3[CH2:4][CH2:5][C@@H:6]2[C@@:11]13[CH3:18]"
    ),
    "drimane": Stereoparent(
        "[CH3:12][C@H:8]1[CH2:7][CH2:6][C@H:5]2[C:4]([CH3:13])([CH3:14])[CH2:3][CH2:2][CH2:1][C@:10]2([CH3:15])[C@H:9]1[CH3:11]",
        "1,2,3,4,5,10",
        True,
        "5:a"
    ),
    "eburnamenine": Stereoparent(
        "[CH3:21][CH2:20][C@:16]12[CH:15]=[CH:14][n:1]3[c:2]4[c:7]([c:8]5[cH:9][cH:10][cH:11][cH:12][c:13]35)[CH2:6][CH2:5][N:4]([CH2:19][CH2:18][CH2:17]1)[C@H:3]24",
        "1,2,7,8,13",
        True,
        "3:b"
    ),
    "emetan": Stereoparent(
        "[CH3:17][CH2:16][C@H:3]1[CH2:4][N:5]2[CH2:6][CH2:7][c:8]3[cH:9][cH:10][cH:11][cH:12][c:13]3[C@@H:14]2[CH2:1][C@@H:2]1[CH2:15][C@H:101]1[NH:102][CH2:103][CH2:104][c:110]2[cH:105][cH:106][cH:107][cH:108][c:109]21",
        "1,2,3,4,5,14",
        True
    ),
    "eremophilane": Stereoparent(
        "[CH3:12][CH:11]([CH3:13])[C@@H:7]1[CH2:8][CH2:9][C@H:10]2[CH2:1][CH2:2][CH2:3][C@H:4]([CH3:14])[C@@:5]2([CH3:15])[CH2:6]1",
        "1,2,3,4,5,10",
        True,
        "5:b"
    ),
    "ergoline": Stereoparent(
        "[cH:13]1[cH:12][c:11]2[c:16]3[c:3]([cH:2][nH:1][c:15]3[cH:14]1)[CH2:4][C@H:5]1[NH:6][CH2:7][CH2:8][CH2:9][C@H:10]21",
        "1,2,3,16,15",
        True,
        "5:b"
    ),
    "ergostane": Stereoparent(
        "[CH3:26][CH:25]([CH3:27])[C@@H:24]([CH3:10024])[CH2:23][CH2:22][C@@H:20]([CH3:21])[C@H:17]1[CH2:16][CH2:15][C@H:14]2[C@@H:8]3[CH2:7][CH2:6][CH:5]4[CH2:4][CH2:3][CH2:2][CH2:1][C@:10]4([CH3:19])[C@H:9]3[CH2:11][CH2:12][C@:13]12[CH3:18]",
        "",
        None,
        "8:b",
        "table10.1"
    ),
    "ergotaman": Stereoparent(
        "[nH:1]1[cH:2][c:3]2[c:16]3[c:11]([cH:12][cH:13][cH:14][c:15]13)[C:10]1=[CH:9][C@@H:8]([CH2:18][NH:19][C@H:102]3[O:101][C@H:112]4[N:104]([CH2:103]3)[CH2:105][CH2:106][N:107]3[CH2:108][CH2:109][CH2:110][C@H:111]43)[CH2:7][N:6]([CH3:17])[C@@H:5]1[CH2:4]2",
        "1,2,3,16,15",
        True
    ),
    "erythrinan": Stereoparent(
        "[cH:16]1[cH:15][cH:14][c:13]2[c:12]([cH:17]1)[CH2:11][CH2:10][N:9]1[CH2:8][CH2:7][C@@H:6]3[CH2:1][CH2:2][CH2:3][CH2:4][C@:5]213",
        "1,2,3,4,5,6",
        True,
        "6:b"
    ),
    "estrane": Stereoparent(
        "[CH3:18][C@@:13]12[CH2:17][CH2:16][CH2:15][C@H:14]1[C@@H:8]1[CH2:7][CH2:6][CH:5]3[CH2:4][CH2:3][CH2:2][CH2:1][C@@H:10]3[C@H:9]1[CH2:11][CH2:12]2",
        "",
        None,
        "8:b",
        "table10.1"
    ),
    "eudesmane": Stereoparent(
        "[CH3:12][CH:11]([CH3:13])[C@@H:7]1[CH2:8][CH2:9][C@@:10]2([CH3:15])[CH2:1][CH2:2][CH2:3][C@H:4]([CH3:14])[C@@H:5]2[CH2:6]1",
        "1,2,3,4,5,10",
        True,
        "5:a"
    ),
    "evonimine": Stereoparent(
        "[CH3:90003][C:90001](=[O:90002])[O:90004][CH2:15][C@:10]12[C@H:9]([O:90020][C:90017]([CH3:90019])=[O:90018])[C:8](=[O:90029])[C@H:7]3[C@@H:6]([O:90016][C:90013]([CH3:90015])=[O:90014])[C@@:5]14[O:29][C@@:11]3([CH3:13])[CH2:12][O:16][C:17](=[O:90030])[c:18]1[cH:19][cH:20][cH:21][n:22][c:23]1[CH2:24][CH2:25][C@H:26]([CH3:30])[C:27](=[O:90033])[O:28][C@@H:3]([C@H:2]([O:90012][C:90009]([CH3:90011])=[O:90010])[C@@H:1]2[O:90008][C:90005]([CH3:90007])=[O:90006])[C@:4]4([CH3:14])[OH:90031]"
    ),
    "evonine": Stereoparent(
        "[CH3:90003][C:90001](=[O:90002])[O:90004][CH2:15][C:10]12[CH:9]([O:90020][C:90017]([CH3:90019])=[O:90018])[C:8](=[O:90029])[CH:7]3[CH:6]([O:90016][C:90013]([CH3:90015])=[O:90014])[C:5]14[O:28][C:11]3([CH3:13])[CH2:12][O:16][C:17](=[O:90030])[c:18]1[cH:19][cH:20][cH:21][n:22][c:23]1[CH:24]([CH3:29])[CH:25]([CH3:30])[C:26](=[O:90033])[O:27][CH:3]([CH:2]([O:90012][C:90009]([CH3:90011])=[O:90010])[CH:1]2[O:90008][C:90005]([CH3:90007])=[O:90006])[C:4]4([CH3:14])[OH:90031]"
    ),
    "fenchane": Stereoparent(
        "[CH3:10][C:1]12[CH2:2][CH2:3][CH:4]([CH2:7]1)[C:5]([CH3:8])([CH3:9])[CH2:6]2",
        "",
        None,
        "",
        "table10.1"
    ),
    "flavan": Stereoparent(
        "[cH:104]1[cH:103][cH:102][c:101]([CH:2]2[CH2:3][CH2:4][c:1004]3[cH:5][cH:6][cH:7][cH:8][c:1008]3[O:1]2)[cH:106][cH:105]1"
    ),
    "formosanan": Stereoparent(
        "[CH:17]1=[CH:16][C@H:15]2[CH2:14][C@@H:3]3[N:4]([CH2:5][CH2:6][C@:7]34[CH2:2][NH:1][c:13]3[cH:12][cH:11][cH:10][cH:9][c:8]34)[CH2:21][C@@H:20]2[CH2:19][O:18]1",
        "3,14,15,20,21,4",
        True
    ),
    "furostan": Stereoparent(
        "[CH3:27][CH:25]([CH3:26])[CH2:24][CH2:23][CH:22]1[O:90025][C@H:16]2[CH2:15][C@H:14]3[C@@H:8]4[CH2:7][CH2:6][CH:5]5[CH2:4][CH2:3][CH2:2][CH2:1][C@:10]5([CH3:19])[C@H:9]4[CH2:11][CH2:12][C@:13]3([CH3:18])[C@H:17]2[C@@H:20]1[CH3:21]",
        "",
        None,
        "8:b"
    ),
    "galanthamine": Stereoparent(
        "[CH3:13][O:90025][c:6]1[cH:7][cH:8][c:1008]2[c:2012]3[c:1005]1[O:5][C@H:1004]1[CH2:4][C@@H:3]([OH:90031])[CH:2]=[CH:1][C@@:1012]31[CH2:12][CH2:11][N:10]([CH3:14])[CH2:9]2",
        "1,2,3,4,4a,12a",
        True
    ),
    "galanthan": Stereoparent(
        "[cH:9]1[cH:10][cH:11][c:14]2[c:13]([cH:8]1)[CH2:7][N:6]1[CH2:5][CH2:4][C@@H:12]3[CH2:3][CH2:2][CH2:1][C@@H:15]2[C@H:16]13",
        "1,2,3,12,16,15",
        False
    ),
    "gammacerane": Stereoparent(
        "[CH3:29][C:22]1([CH3:30])[CH2:21][CH2:20][CH2:19][C@:18]2([CH3:28])[C@H:13]3[CH2:12][CH2:11][C@@H:9]4[C@@:10]5([CH3:25])[CH2:1][CH2:2][CH2:3][C:4]([CH3:23])([CH3:24])[C@@H:5]5[CH2:6][CH2:7][C@@:8]4([CH3:26])[C@:14]3([CH3:27])[CH2:15][CH2:16][C@@H:17]12",
        "1,2,3,4,5,10",
        True,
        "5:a"
    ),
    "germacrane": Stereoparent(
        "[CH2:1]1[CH2:2][CH2:3][C@H:4]([CH3:14])[CH2:5][CH2:6][C@H:7]([CH:11]([CH3:12])[CH3:13])[CH2:8][CH2:9][C@@H:10]1[CH3:15]",
        "1,2,3,4,5,6,7,8,9,10",
        True,
        "10:a"
    ),
    "gibbane": Stereoparent(
        "[CH2:3]1[CH2:2][CH2:1][C@H:1010]2[CH2:10][C@:1009]34[CH2:9][CH2:8][C@H:7]([CH2:6][CH2:5][CH:2004]3[CH:1004]2[CH2:4]1)[CH2:11]4",
        "",
        None,
        "10a:b"
    ),
    "gonane": Stereoparent(
        "[CH2:3]1[CH2:2][CH2:1][C@H:10]2[CH:5]([CH2:4]1)[CH2:6][CH2:7][C@H:8]1[C@@H:14]3[CH2:15][CH2:16][CH2:17][C@H:13]3[CH2:12][CH2:11][C@H:9]21",
        "",
        None,
        "8:b",
        "table10.1"
    ),
    "gorgostane": Stereoparent(
        "[CH3:26][CH:25]([CH3:27])[C@@H:24]([CH3:10024])[C@@:23]1([CH3:10023])[CH2:10022][C@@H:22]1[C@@H:20]([CH3:21])[C@H:17]1[CH2:16][CH2:15][C@H:14]2[C@@H:8]3[CH2:7][CH2:6][CH:5]4[CH2:4][CH2:3][CH2:2][CH2:1][C@:10]4([CH3:19])[C@H:9]3[CH2:11][CH2:12][C@:13]12[CH3:18]",
        "1,2,3,4,5,10",
        True,
        "8:b"
    ),
    "grayanotoxane": Stereoparent(
        "[CH3:20][C@@H:10]1[C@@H:1]2[CH2:2][CH2:3][C:4]([CH3:18])([CH3:19])[C@H:5]2[CH2:6][CH2:7][C@@:8]23[CH2:14][C@@H:13]([CH2:12][CH2:11][C@@H:9]12)[C@@H:15]([CH3:17])[CH2:16]3",
        "",
        None,
        "5:b"
    ),
    "guaiane": Stereoparent(
        "[CH3:12][CH:11]([CH3:13])[C@@H:7]1[CH2:8][CH2:9][C@H:10]([CH3:14])[C@@H:1]2[CH2:2][CH2:3][C@H:4]([CH3:15])[C@@H:5]2[CH2:6]1",
        "1,2,3,4,5",
        True,
        "5:a"
    ),
    "hasubanan": Stereoparent(
        "[cH:2]1[cH:3][cH:4][c:12]2[c:11]([cH:1]1)[CH2:10][CH2:9][C@@:14]13[CH2:8][CH2:7][CH2:6][CH2:5][C@@:13]21[CH2:15][CH2:16][NH:17]3"
    ),
    "hetisan": Stereoparent(
        "[CH2:17]=[C:16]1[CH2:15][C@:8]23[CH2:7][C@H:6]4[C@@H:5]5[C@@:4]6([CH3:18])[CH2:3][CH2:2][CH2:1][C@:10]57[CH:20]([CH:14]2[CH2:13][C@H:12]1[CH2:11][C@H:9]37)[N:21]4[CH2:19]6"
    ),
    "himachalane": Stereoparent(
        "[CH3:13][CH:4]1[CH2:3][CH2:2][C@H:1]2[CH:11]([CH3:12])[CH2:10][CH2:9][CH2:8][C:7]([CH3:14])([CH3:15])[C@H:6]2[CH2:5]1",
        "1,2,3,4,5,6",
        True,
        "1:a"
    ),
    "hopane": Stereoparent(
        "[CH3:29][CH:22]([CH3:30])[C@H:21]1[CH2:20][CH2:19][C@:18]2([CH3:28])[C@H:13]3[CH2:12][CH2:11][C@@H:9]4[C@@:10]5([CH3:25])[CH2:1][CH2:2][CH2:3][C:4]([CH3:23])([CH3:24])[C@@H:5]5[CH2:6][CH2:7][C@@:8]4([CH3:26])[C@:14]3([CH3:27])[CH2:15][CH2:16][C@@H:17]12",
        "1,2,3,4,5,10",
        True,
        "5:a"
    ),
    "humulane": Stereoparent(
        "[CH3:12][CH:11]1[CH2:1][CH2:2][CH2:3][CH:4]([CH3:13])[CH2:5][CH2:6][C:7]([CH3:14])([CH3:15])[CH2:8][CH2:9][CH2:10]1"
    ),
    "ibogamine": Stereoparent(
        "[CH3:21][CH2:20][C@H:4]1[CH2:3][C@@H:2]2[CH2:1][C@H:18]3[c:17]4[nH:16][c:15]5[cH:14][cH:13][cH:12][cH:11][c:10]5[c:9]4[CH2:8][CH2:7][N:6]([CH2:19]2)[C@@H:5]13"
    ),
    "isoflavan": Stereoparent(
        "[cH:104]1[cH:103][cH:102][c:101]([CH:3]2[CH2:2][O:1][c:1008]3[cH:8][cH:7][cH:6][cH:5][c:1004]3[CH2:4]2)[cH:106][cH:105]1"
    ),
    "kaurane": Stereoparent(
        "[CH3:17][C@@H:16]1[CH2:15][C@:8]23[CH2:7][CH2:6][C@H:5]4[C:4]([CH3:18])([CH3:19])[CH2:3][CH2:2][CH2:1][C@:10]4([CH3:20])[C@H:9]2[CH2:11][CH2:12][C@H:13]1[CH2:14]3",
        "",
        None,
        "5:a"
    ),
    "kopsan": Stereoparent(
        "[cH:16]1[cH:15][cH:14][c:13]2[c:18]([cH:17]1)[NH:1][C@:2]13[CH2:21][CH2:20][C@:5]45[CH2:6][CH2:7][CH2:8][N:9]6[CH2:10][CH:11]([CH2:22][C@@H:3]1[CH2:4]4)[C@:12]23[C@H:19]56"
    ),
    "labdane": Stereoparent(
        "[CH3:15][CH2:14][C@@H:13]([CH3:16])[CH2:12][CH2:11][C@H:9]1[C@@H:8]([CH3:20])[CH2:7][CH2:6][C@H:5]2[C:4]([CH3:18])([CH3:19])[CH2:3][CH2:2][CH2:1][C@:10]12[CH3:17]",
        "1,2,3,4,5,10",
        True,
        "5:a"
    ),
    "lanostane": Stereoparent(
        "[CH3:26][CH:25]([CH3:27])[CH2:24][CH2:23][CH2:22][C@@H:20]([CH3:21])[C@H:17]1[CH2:16][CH2:15][C@@:14]2([CH3:30])[C@@H:8]3[CH2:7][CH2:6][C@H:5]4[C:4]([CH3:28])([CH3:29])[CH2:3][CH2:2][CH2:1][C@:10]4([CH3:19])[C@H:9]3[CH2:11][CH2:12][C@:13]12[CH3:18]",
        "1,2,3,4,5,10",
        True,
        "5:a"
    ),
    "lignane": Stereoparent(
        "[CH3:9][CH:8]([CH2:7][c:1]1[cH:2][cH:3][cH:4][cH:5][cH:6]1)[CH:108]([CH3:109])[CH2:107][c:101]1[cH:102][cH:103][cH:104][cH:105][cH:106]1"
    ),
    "lunarine": Stereoparent(
        "[O:90025]=[C:3]1[CH2:4][CH2:5][C@:6]23/[CH:13]=[CH:14]/[C:15](=[O:90027])[NH:16][CH2:17][CH2:18][CH2:19][CH2:20][NH:21][CH2:22][CH2:23][CH2:24][NH:25][C:26](=[O:90028])/[CH:27]=[CH:28]/[c:9]4[cH:10][cH:11][c:12]([c:7]2[cH:8]4)[O:29][C@@H:1]3[CH2:2]1"
    ),
    "lupane": Stereoparent(
        "[CH3:29][CH:20]([CH3:30])[C@@H:19]1[CH2:21][CH2:22][C@:17]2([CH3:28])[CH2:16][CH2:15][C@:14]3([CH3:27])[C@H:13]([CH2:12][CH2:11][C@@H:9]4[C@@:10]5([CH3:25])[CH2:1][CH2:2][CH2:3][C:4]([CH3:23])([CH3:24])[C@@H:5]5[CH2:6][CH2:7][C@@:8]34[CH3:26])[C@@H:18]12",
        "1,2,3,4,5,10",
        True,
        "5:a"
    ),
    "lycopodane": Stereoparent(
        "[CH3:16][CH:15]1[CH2:8][C@H:7]2[CH2:6][CH2:5][C@H:4]3[CH2:3][CH2:2][CH2:1][N:17]4[CH2:9][CH2:10][CH2:11][C@H:12]2[C@:13]34[CH2:14]1"
    ),
    "lycorenan": Stereoparent(
        "[CH:4]1=[C:12]2[CH2:3][CH2:2][NH:1][C@H:17]2[C@@H:16]2[c:15]3[cH:11][cH:10][cH:9][cH:8][c:14]3[CH2:7][O:6][C@@H:13]2[CH2:5]1",
        "1,2,3,12,17",
        False
    ),
    "lythran": Stereoparent(
        "[CH:13]1=[CH:14]\\[c:105]2[cH:104][cH:103][cH:102][c:101]([cH:106]2)-[c:201]2[cH:206][cH:205][cH:204][cH:203][c:202]2[C@@H:4]2[CH2:3][C@H:2]([CH2:1][C@@H:10]3[CH2:9][CH2:8][CH2:7][CH2:6][N:5]23)[O:11][CH2:12]/1"
    ),
    "lythranidine": Stereoparent(
        "[CH3:207][O:90028][c:202]1[cH:203][cH:204][c:205]2[cH:206][c:201]1-[c:101]1[cH:102][c:103]([cH:104][cH:105][c:106]1[OH:90029])[CH2:6][CH2:7][CH:8]([OH:90025])[CH2:9][CH:10]1[CH2:1][CH2:2][CH2:3][CH:4]([CH2:11][CH:12]([OH:90027])[CH2:13][CH2:14]2)[NH:5]1"
    ),
    "matridine": Stereoparent(
        "[CH2:13]1[CH2:14][CH2:15][N:16]2[CH2:17][C@@H:5]3[CH2:4][CH2:3][CH2:2][N:1]4[CH2:10][CH2:9][CH2:8][C@@H:7]([C@H:6]34)[C@H:11]2[CH2:12]1",
        "1,2,3,4,5,6",
        False,
        "11:b"
    ),
    "morphinan": Stereoparent(
        "[cH:2]1[cH:3][cH:4][c:12]2[c:11]([cH:1]1)[CH2:10][C@H:9]1[NH:17][CH2:16][CH2:15][C@@:13]23[CH2:5][CH2:6][CH2:7][CH2:8][C@@H:14]13",
        "",
        None,
        "13:b"
    ),
    "neoflavan": Stereoparent(
        "[cH:104]1[cH:103][cH:102][c:101]([CH:4]2[CH2:3][CH2:2][O:1][c:1008]3[cH:8][cH:7][cH:6][cH:5][c:1004]32)[cH:106][cH:105]1"
    ),
    "nupharidine": Stereoparent(
        "[CH3:11][C@H:7]1[CH2:8][CH2:9][C@H:10]2[C@H:1]([CH3:12])[CH2:2][CH2:3][C@@H:4]([c:103]3[cH:104][cH:105][o:101][cH:102]3)[N@@+:5]2([O-:90025])[CH2:6]1",
        "1,2,3,4,5,10",
        True
    ),
    "oleanane": Stereoparent(
        "[CH3:29][C:20]1([CH3:30])[CH2:21][CH2:22][C@:17]2([CH3:28])[CH2:16][CH2:15][C@:14]3([CH3:27])[C@H:13]([CH2:12][CH2:11][C@@H:9]4[C@@:10]5([CH3:25])[CH2:1][CH2:2][CH2:3][C:4]([CH3:23])([CH3:24])[C@@H:5]5[CH2:6][CH2:7][C@@:8]34[CH3:26])[C@@H:18]2[CH2:19]1",
        "1,2,3,4,5,10",
        True,
        "5:a"
    ),
    "ophiobolane": Stereoparent(
        "[CH3:21][CH:20]([CH3:22])[CH2:19][CH2:18][CH2:17][C@H:15]([CH3:16])[C@H:14]1[CH2:13][CH2:12][C@:11]2([CH3:23])[CH2:1][C@H:2]3[C@H:6]([CH2:5][CH2:4][C@@H:3]3[CH3:24])[C@@H:7]([CH3:25])[CH2:8][CH2:9][C@@H:10]12",
        "1,2,6,7,8,9,10,11",
        False,
        "11:b"
    ),
    "ormosanine": Stereoparent(
        "[CH2:21]1[CH2:22][CH2:23][C@H:18]([C@:9]23[CH2:8][C@H:7]([CH2:17][C@@H:16]4[CH2:15][CH2:14][CH2:13][NH:12][C@@H:11]24)[C@H:6]2[CH2:5][CH2:4][CH2:3][CH2:2][N:1]2[CH2:10]3)[NH:19][CH2:20]1"
    ),
    "oxyacanthan": Stereoparent(
        "[cH:111]1[cH:110][c:109]2[cH:114][c:113]([cH:112]1)[O:18][c:12]1[cH:11][cH:10][c:9]([cH:14][cH:13]1)[CH2:15][C@@H:1]1[NH:2][CH2:3][CH2:4][c:16]3[cH:5][cH:6][cH:7][c:8]([c:17]31)[O:118][c:107]1[cH:106][cH:105][c:116]3[c:117]([cH:108]1)[C@@H:101]([CH2:115]2)[NH:102][CH2:103][CH2:104]3"
    ),
    "p-menthane": Stereoparent(
        "[CH3:7][CH:1]1[CH2:2][CH2:3][CH:4]([CH:8]([CH3:9])[CH3:10])[CH2:5][CH2:6]1",
        "",
        None,
        "",
        "table10.1"
    ),
    "pancracine": Stereoparent(
        "[OH:90025][C@H:2]1[CH:1]=[C:1011]2[C@H:11]3[CH2:12][N:5]([CH2:6][c:1006]4[cH:7][c:8]5[c:9]([cH:10][c:1010]43)[O:15][CH2:13][O:14]5)[C@H:1004]2[CH2:4][C@@H:3]1[OH:90027]"
    ),
    "penam": Stereoparent(
        "[O:90025]=[C:7]1[CH2:6][C@H:5]2[S:1][CH2:2][CH2:3][N:4]12",
        "1,2,3,4,5",
        False,
        "5:a"
    ),
    "picrasane": Stereoparent(
        "[CH3:18][C@@H:4]1[CH2:3][CH2:2][CH2:1][C@:10]2([CH3:19])[C@H:9]3[CH2:11][CH2:12][C@H:13]([CH3:21])[C@@H:14]4[CH2:15][CH2:16][O:17][C@H:7]([CH2:6][C@@H:5]12)[C@:8]34[CH3:20]",
        "1,2,3,4,5,10",
        True,
        "5:a"
    ),
    "pimarane": Stereoparent(
        "[CH3:16][CH2:15][C@:13]1([CH3:17])[CH2:12][CH2:11][C@H:9]2[C@@H:8]([CH2:7][CH2:6][C@H:5]3[C:4]([CH3:18])([CH3:19])[CH2:3][CH2:2][CH2:1][C@:10]23[CH3:20])[CH2:14]1",
        "1,2,3,4,5,10",
        True,
        "5:a"
    ),
    "pinane": Stereoparent(
        "[CH3:10][CH:2]1[CH2:3][CH2:4][CH:5]2[CH2:7][CH:1]1[C:6]2([CH3:8])[CH3:9]",
        "",
        None,
        "",
        "table10.1"
    ),
    "podocarpane": Stereoparent(
        "[CH3:15][C:4]1([CH3:16])[CH2:3][CH2:2][CH2:1][C@:10]2([CH3:17])[C@H:9]3[CH2:11][CH2:12][CH2:13][CH2:14][C@@H:8]3[CH2:7][CH2:6][C@@H:5]12",
        "1,2,3,4,5,10",
        True,
        "5:a"
    ),
    "poriferastane": Stereoparent(
        "[CH3:20024][CH2:10024][C@@H:24]([CH2:23][CH2:22][C@@H:20]([CH3:21])[C@H:17]1[CH2:16][CH2:15][C@H:14]2[C@@H:8]3[CH2:7][CH2:6][CH:5]4[CH2:4][CH2:3][CH2:2][CH2:1][C@:10]4([CH3:19])[C@H:9]3[CH2:11][CH2:12][C@:13]12[CH3:18])[CH:25]([CH3:26])[CH3:27]",
        "",
        None,
        "8:b",
        "table10.1"
    ),
    "porphyrin": Stereoparent(
        "[c:1]12[cH:2][cH:3][c:4]([cH:5][c:6]3[n:22][c:9]([cH:10][c:11]4[cH:12][cH:13][c:14]([cH:15][c:16]5[n:24][c:19]([cH:20]1)[CH:18]=[CH:17]5)[nH:23]4)[CH:8]=[CH:7]3)[nH:21]2",
        "",
        None,
        "",
        "table10.1"
    ),
    "pregnane": Stereoparent(
        "[CH3:21][CH2:20][C@H:17]1[CH2:16][CH2:15][C@H:14]2[C@@H:8]3[CH2:7][CH2:6][C@@H:5]4[CH2:4][CH2:3][CH2:2][CH2:1][C@:10]4([CH3:19])[C@H:9]3[CH2:11][CH2:12][C@:13]12[CH3:18]",
        "",
        None,
        "8:b",
        "table10.1"
    ),
    "prostane": Stereoparent(
        "[CH3:20][CH2:19][CH2:18][CH2:17][CH2:16][CH2:15][CH2:14][CH2:13][C@H:12]1[CH2:11][CH2:10][CH2:9][C@@H:8]1[CH2:7][CH2:6][CH2:5][CH2:4][CH2:3][CH2:2][CH3:1]",
        "8,9,10,11,12",
        True,
        "8:b"
    ),
    "protostane": Stereoparent(
        "[CH3:26][CH:25]([CH3:27])[CH2:24][CH2:23][CH2:22][C@@H:20]([CH3:21])[C@H:17]1[CH2:16][CH2:15][C@@:14]2([CH3:30])[C@H:13]1[CH2:12][CH2:11][C@H:9]1[C@@:10]3([CH3:19])[CH2:1][CH2:2][CH2:3][C:4]([CH3:28])([CH3:29])[C@@H:5]3[CH2:6][CH2:7][C@:8]21[CH3:18]",
        "1,2,3,4,5,10",
        True,
        "5:a"
    ),
    "retinal": Stereoparent(
        "[CH3:18][C:5]1=[C:6](/[CH:7]=[CH:8]/[C:9]([CH3:19])=[CH:10]/[CH:11]=[CH:12]/[C:13]([CH3:20])=[CH:14]/[CH:15]=[O:90025])[C:1]([CH3:16])([CH3:17])[CH2:2][CH2:3][CH2:4]1"
    ),
    "rheadan": Stereoparent(
        "[cH:11]1[cH:12][cH:13][c:14]2[c:9]([cH:10]1)[CH2:8][O:7][C@@H:6]1[c:5]3[cH:4][cH:3][cH:2][cH:1][c:19]3[CH2:18][CH2:17][NH:16][C@H:15]21",
        "1,2,3,4,5,19",
        True
    ),
    "rodiasine": Stereoparent(
        "[CH3:120][O:90029][c:106]1[cH:105][c:116]2[c:117]3[cH:108][c:107]1[O:119][c:8]1[c:7]([O:90027][CH3:20])[c:6]([O:90025][CH3:19])[cH:5][c:16]4[c:17]1[C@H:1]([CH2:15][c:9]1[cH:14][cH:13][c:12]([O:90028][CH3:21])[c:11]([cH:10]1)-[c:111]1[cH:110][c:109]([cH:114][cH:113][c:112]1[OH:90031])[CH2:115][C@H:101]3[N:102]([CH3:118])[CH2:103][CH2:104]2)[N:2]([CH3:18])[CH2:3][CH2:4]4"
    ),
    "rosane": Stereoparent(
        "[CH3:16][CH2:15][C@:13]1([CH3:17])[CH2:12][CH2:11][C@:9]2([CH3:20])[C@H:8]([CH2:7][CH2:6][C@@H:5]3[C@H:10]2[CH2:1][CH2:2][CH2:3][C:4]3([CH3:18])[CH3:19])[CH2:14]1",
        "1,2,3,4,5,10",
        True,
        "5:a"
    ),
    "samandarine": Stereoparent(
        "[CH3:18][C@:13]12[CH2:12][CH2:11][C@H:9]3[C@@H:8]([CH2:7][CH2:6][C@@H:5]4[CH2:1004][C@H:4]5[NH:3][CH2:2][C@H:1]([O:20]5)[C@:10]34[CH3:19])[C@@H:14]1[CH2:15][C@H:16]([OH:90031])[CH2:17]2"
    ),
    "sarpagan": Stereoparent(
        "[CH3:18]/[CH:19]=[C:20]1/[CH2:21][N:4]2[C@H:3]3[CH2:14][C@@H:15]1[C@@H:16]([CH3:17])[C@@H:5]2[CH2:6][c:7]1[c:2]3[nH:1][c:13]2[cH:12][cH:11][cH:10][cH:9][c:8]12"
    ),
    "senecionan": Stereoparent(
        "[C:1]12=[CH:2][CH2:3][N:4]3[CH2:5][CH2:6][C@H:7]([C@@H:8]13)[O:17][CH2:16]/[C:15](=[CH:20]\\[CH3:21])[CH2:14][C@@H:13]([CH3:19])[C@H:12]([CH3:18])[CH2:11][O:10][CH2:9]2",
        "1,2,3,4,8",
        False
    ),
    "solanidane": Stereoparent(
        "[CH3:27][C@H:25]1[CH2:24][CH2:23][C@@H:22]2[C@@H:20]([CH3:21])[C@H:17]3[C@H:16]([CH2:15][C@H:14]4[C@@H:8]5[CH2:7][CH2:6][CH:5]6[CH2:4][CH2:3][CH2:2][CH2:1][C@:10]6([CH3:19])[C@H:9]5[CH2:11][CH2:12][C@:13]34[CH3:18])[N:28]2[CH2:26]1",
        "1,2,3,4,5,10",
        True
    ),
    "sparteine": Stereoparent(
        "[CH2:4]1[CH2:3][CH2:2][N:1]2[CH2:10][C@@H:9]3[CH2:8][C@@H:7]([CH2:17][N:16]4[CH2:15][CH2:14][CH2:13][CH2:12][C@@H:11]34)[C@H:6]2[CH2:5]1"
    ),
    "spirosolane": Stereoparent(
        "[CH3:27][CH:25]1[CH2:24][CH2:23][C:22]2([NH:28][CH2:26]1)[O:90025][C@H:16]1[CH2:15][C@H:14]3[C@@H:8]4[CH2:7][CH2:6][CH:5]5[CH2:4][CH2:3][CH2:2][CH2:1][C@:10]5([CH3:19])[C@H:9]4[CH2:11][CH2:12][C@:13]3([CH3:18])[C@H:17]1[C@@H:20]2[CH3:21]",
        "1,2,3,4,5,10",
        True,
        "8:b"
    ),
    "spirostan": Stereoparent(
        "[CH2:1]1[CH2:2][CH2:3][CH2:4][CH:5]2[CH2:6][CH2:7][C@@H:8]3[C@@H:9]([C@@:10]12[CH3:19])[CH2:11][CH2:12][C@@:13]1([CH3:18])[C@H:14]3[CH2:15][C@H:16]2[C@@H:17]1[C@H:20]([CH3:21])[C:22]1([CH2:23][CH2:24][CH:25]([CH3:27])[CH2:26][O:90027]1)[O:90025]2",
        "1,2,3,4,5,10",
        True,
        "8:b"
    ),
    "stigmastane": Stereoparent(
        "[CH3:20024][CH2:10024][C@H:24]([CH2:23][CH2:22][C@@H:20]([CH3:21])[C@H:17]1[CH2:16][CH2:15][C@H:14]2[C@@H:8]3[CH2:7][CH2:6][CH:5]4[CH2:4][CH2:3][CH2:2][CH2:1][C@:10]4([CH3:19])[C@H:9]3[CH2:11][CH2:12][C@:13]12[CH3:18])[CH:25]([CH3:26])[CH3:27]",
        "",
        None,
        "8:b",
        "table10.1"
    ),
    "strychnidine": Stereoparent(
        "[CH:22]1=[C:21]2[CH2:20][N:19]3[CH2:18][CH2:17][C@:7]45[c:6]6[cH:1][cH:2][cH:3][cH:4][c:5]6[N:9]6[CH2:10][CH2:11][C@H:12]([O:24][CH2:23]1)[C@@H:13]([C@@H:8]46)[C@H:14]2[CH2:15][C@H:16]35"
    ),
    "taxane": Stereoparent(
        "[CH3:20][C@@H:4]1[CH2:5][CH2:6][CH2:7][C@@:8]2([CH3:19])[CH2:9][CH2:10][C@H:11]3[C@H:12]([CH3:18])[CH2:13][CH2:14][C@@H:1]([CH2:2][C@H:3]12)[C:15]3([CH3:16])[CH3:17]",
        "",
        None,
        "1:b"
    ),
    "tazettine": Stereoparent(
        "[CH3:16][O:90025][C@@H:3]1[CH:2]=[CH:1][C@@:1012]23[c:3012]4[cH:12][c:11]5[c:10]([cH:9][c:1008]4[CH2:8][O:7][C@:1006]2([OH:90031])[CH2:6][N:5]([CH3:17])[C@H:2012]3[CH2:4]1)[O:13][CH2:14][O:15]5",
        "1,2,3,4,12b,12a",
        False
    ),
    "thromboxane": Stereoparent(
        "[CH3:20][CH2:19][CH2:18][CH2:17][CH2:16][CH2:15][CH2:14][CH2:13][C@H:12]1[O:90025][CH2:11][CH2:10][CH2:9][C@@H:8]1[CH2:7][CH2:6][CH2:5][CH2:4][CH2:3][CH2:2][CH3:1]",
        "8,9,10,11,_25,12",
        True,
        "8:b"
    ),
    "thujane": Stereoparent(
        "[CH3:10][CH:4]1[CH2:3][CH2:2][C:1]2([CH:7]([CH3:8])[CH3:9])[CH2:6][CH:5]12",
        "",
        None,
        "",
        "table10.1"
    ),
    "trichothecane": Stereoparent(
        "[CH3:16][CH:9]1[CH2:8][CH2:7][C@@:6]2([CH3:15])[C@@H:11]([CH2:10]1)[O:1][C@@H:2]1[CH2:3][CH2:4][C@@:5]2([CH3:14])[C@@H:12]1[CH3:13]",
        "",
        None,
        "6:a"
    ),
    "tropane": Stereoparent(
        "[CH3:9][N:8]1[C@@H:1]2[CH2:2][CH2:3][CH2:4][C@H:5]1[CH2:6][CH2:7]2"
    ),
    "tubocuraran": Stereoparent(
        "[cH:113]1[cH:114][c:109]2[cH:110][c:111]([cH:112]1)[O:18][c:7]1[cH:6][cH:5][c:16]3[c:17]([cH:8]1)[C@H:1]([CH2:15][c:9]1[cH:10][cH:11][c:12]([cH:13][cH:14]1)[O:118][c:108]1[cH:107][cH:106][cH:105][c:116]4[c:117]1[C@@H:101]([CH2:115]2)[NH:102][CH2:103][CH2:104]4)[NH:2][CH2:3][CH2:4]3"
    ),
    "tubulosan": Stereoparent(
        "[CH3:17][CH2:16][C@H:3]1[CH2:4][N:5]2[CH2:6][CH2:7][c:8]3[cH:9][cH:10][cH:11][cH:12][c:13]3[C@@H:14]2[CH2:1][C@@H:2]1[CH2:15][C@H:101]1[NH:102][CH2:103][CH2:104][c:105]2[c:113]1[nH:112][c:111]1[cH:110][cH:109][cH:108][cH:107][c:106]21",
        "1,2,3,4,5,14",
        True
    ),
    "ursane": Stereoparent(
        "[CH3:29][C@@H:19]1[C@H:18]2[C@H:13]3[CH2:12][CH2:11][C@@H:9]4[C@@:10]5([CH3:25])[CH2:1][CH2:2][CH2:3][C:4]([CH3:23])([CH3:24])[C@@H:5]5[CH2:6][CH2:7][C@@:8]4([CH3:26])[C@:14]3([CH3:27])[CH2:15][CH2:16][C@@:17]2([CH3:28])[CH2:22][CH2:21][C@H:20]1[CH3:30]",
        "1,2,3,4,5,10",
        True,
        "5:a"
    ),
    "veratraman": Stereoparent(
        "[CH3:18][C:13]1=[C:12]2[CH2:11][C@H:9]3[C@@H:8]([CH2:7][CH:6]=[C:5]4[CH2:4][CH2:3][CH2:2][CH2:1][C@:10]34[CH3:19])[C@@H:14]2[CH2:15][CH2:16][C@H:17]1[C@H:20]([CH3:21])[C@H:22]1[CH2:23][CH2:24][C@H:25]([CH3:27])[CH2:26][NH:28]1",
        "1,2,3,4,5,10",
        True
    ),
    "vincaleukoblastine": Stereoparent(
        "[CH3:121][CH2:120][C@:104]1([OH:90032])[CH2:103][C@@H:102]2[CH2:119][N:106]([CH2:107][CH2:108][c:109]3[c:117]([nH:116][c:115]4[cH:114][cH:113][cH:112][cH:111][c:110]34)[C@@:118]([C:122](=[O:90036])[O:90038][CH3:123])([c:15]3[cH:14][c:13]4[c:18]([cH:17][c:16]3[O:90034][CH3:25])[N:1]([CH3:22])[C@H:2]3[C@@:3]([OH:90031])([C:23](=[O:90035])[O:90037][CH3:24])[C@H:4]([O:90039][C:90021]([CH3:90023])=[O:90022])[C@:5]5([CH2:20][CH3:21])[CH:6]=[CH:7][CH2:8][N:9]6[CH2:10][CH2:11][C@:12]43[C@H:19]56)[CH2:101]2)[CH2:105]1"
    ),
    "vincane": Stereoparent(
        "[CH3:21][CH2:20][C@:16]12[CH2:17][CH2:18][CH2:19][N:4]3[CH2:5][CH2:6][c:7]4[c:2]([n:1]([c:13]5[cH:12][cH:11][cH:10][cH:9][c:8]45)[CH2:14][CH2:15]1)[C@H:3]23",
        "1,2,7,8,13",
        True
    ),
    "vobasan": Stereoparent(
        "[CH3:18]/[CH:19]=[C:20]1/[CH2:21][N:4]([CH3:22])[C@@H:5]2[CH2:6][c:7]3[c:2]([nH:1][c:13]4[cH:12][cH:11][cH:10][cH:9][c:8]34)[CH2:3][CH2:14][C@@H:15]1[C@@H:16]2[CH3:17]"
    ),
    "vobtusine": Stereoparent(
        "[CH3:90024][O:90037][C:23](=[O:90035])[C:3]1=[C:2]2[NH:1][c:18]3[cH:17][cH:16][cH:15][cH:14][c:13]3[C@@:12]23[CH2:11][CH2:10][N:9]2[CH2:8][C@:7]4([CH2:122][C@@H:103]5[CH2:104][C@:105]67[CH2:120][CH2:121][O:124][C@H:106]6[CH2:107][CH2:108][N:109]6[CH2:110][CH2:111][C@@:112]8([c:113]9[cH:114][cH:115][cH:116][c:117]([O:90026][CH3:125])[c:118]9[N:101]([CH2:123]4)[C@@:102]58[OH:90032])[C@H:119]76)[C@@H:6]4[O:22][CH2:21][CH2:20][C@:5]4([CH2:4]1)[C@@H:19]32"
    ),
    "yohimban": Stereoparent(
        "[cH:11]1[cH:10][cH:9][c:8]2[c:7]3[c:2]([nH:1][c:13]2[cH:12]1)[C@@H:3]1[CH2:14][C@@H:15]2[CH2:16][CH2:17][CH2:18][CH2:19][C@H:20]2[CH2:21][N:4]1[CH2:5][CH2:6]3",
        "1,2,7,8,13",
        True,
        "3:a"
    ),
    "β,β-carotene": Stereoparent(
        "[CH3:18][C:5]1=[C:6](/[CH:7]=[CH:8]/[C:9]([CH3:19])=[CH:10]/[CH:11]=[CH:12]/[C:13]([CH3:20])=[CH:14]/[CH:15]=[CH:115]/[CH:114]=[C:113]([CH3:120])/[CH:112]=[CH:111]/[CH:110]=[C:109]([CH3:119])/[CH:108]=[CH:107]/[C:106]2=[C:105]([CH3:118])[CH2:104][CH2:103][CH2:102][C:101]2([CH3:116])[CH3:117])[C:1]([CH3:16])([CH3:17])[CH2:2][CH2:3][CH2:4]1"
    ),
    "β,γ-carotene": Stereoparent(
        "[CH2:118]=[C:105]1[CH2:104][CH2:103][CH2:102][C:101]([CH3:116])([CH3:117])[CH:106]1/[CH:107]=[CH:108]/[C:109]([CH3:119])=[CH:110]/[CH:111]=[CH:112]/[C:113]([CH3:120])=[CH:114]/[CH:115]=[CH:15]/[CH:14]=[C:13]([CH3:20])/[CH:12]=[CH:11]/[CH:10]=[C:9]([CH3:19])/[CH:8]=[CH:7]/[C:6]1=[C:5]([CH3:18])[CH2:4][CH2:3][CH2:2][C:1]1([CH3:16])[CH3:17]"
    ),
    "β,ε-carotene": Stereoparent(
        "[CH3:118][C:105]1=[CH:104][CH2:103][CH2:102][C:101]([CH3:116])([CH3:117])[CH:106]1/[CH:107]=[CH:108]/[C:109]([CH3:119])=[CH:110]/[CH:111]=[CH:112]/[C:113]([CH3:120])=[CH:114]/[CH:115]=[CH:15]/[CH:14]=[C:13]([CH3:20])/[CH:12]=[CH:11]/[CH:10]=[C:9]([CH3:19])/[CH:8]=[CH:7]/[C:6]1=[C:5]([CH3:18])[CH2:4][CH2:3][CH2:2][C:1]1([CH3:16])[CH3:17]"
    ),
    "β,κ-carotene": Stereoparent(
        "[CH3:18][C:5]1=[C:6](/[CH:7]=[CH:8]/[C:9]([CH3:19])=[CH:10]/[CH:11]=[CH:12]/[C:13]([CH3:20])=[CH:14]/[CH:15]=[CH:115]/[CH:114]=[C:113]([CH3:120])/[CH:112]=[CH:111]/[CH:110]=[C:109]([CH3:119])/[CH:108]=[CH:107]/[CH2:106][C:105]2([CH3:118])[CH2:104][CH2:103][CH2:102][C:101]2([CH3:116])[CH3:117])[C:1]([CH3:16])([CH3:17])[CH2:2][CH2:3][CH2:4]1"
    ),
    "β,φ-carotene": Stereoparent(
        "[CH3:18][C:5]1=[C:6](/[CH:7]=[CH:8]/[C:9]([CH3:19])=[CH:10]/[CH:11]=[CH:12]/[C:13]([CH3:20])=[CH:14]/[CH:15]=[CH:115]/[CH:114]=[C:113]([CH3:120])/[CH:112]=[CH:111]/[CH:110]=[C:109]([CH3:119])/[CH:108]=[CH:107]/[c:106]2[c:105]([CH3:118])[cH:104][cH:103][c:102]([CH3:117])[c:101]2[CH3:116])[C:1]([CH3:16])([CH3:17])[CH2:2][CH2:3][CH2:4]1"
    ),
    "β,χ-carotene": Stereoparent(
        "[CH3:18][C:5]1=[C:6](/[CH:7]=[CH:8]/[C:9]([CH3:19])=[CH:10]/[CH:11]=[CH:12]/[C:13]([CH3:20])=[CH:14]/[CH:15]=[CH:115]/[CH:114]=[C:113]([CH3:120])/[CH:112]=[CH:111]/[CH:110]=[C:109]([CH3:119])/[CH:108]=[CH:107]/[c:106]2[cH:105][cH:104][c:103]([CH3:118])[c:102]([CH3:117])[c:101]2[CH3:116])[C:1]([CH3:16])([CH3:17])[CH2:2][CH2:3][CH2:4]1"
    ),
    "β,ψ-carotene": Stereoparent(
        "[CH3:116][C:101]([CH3:117])=[CH:102][CH2:103][CH2:104][C:105]([CH3:118])=[CH:106]/[CH:107]=[CH:108]/[C:109]([CH3:119])=[CH:110]/[CH:111]=[CH:112]/[C:113]([CH3:120])=[CH:114]/[CH:115]=[CH:15]/[CH:14]=[C:13]([CH3:20])/[CH:12]=[CH:11]/[CH:10]=[C:9]([CH3:19])/[CH:8]=[CH:7]/[C:6]1=[C:5]([CH3:18])[CH2:4][CH2:3][CH2:2][C:1]1([CH3:16])[CH3:17]"
    ),
    "γ,γ-carotene": Stereoparent(
        "[CH2:18]=[C:5]1[CH2:4][CH2:3][CH2:2][C:1]([CH3:16])([CH3:17])[CH:6]1/[CH:7]=[CH:8]/[C:9]([CH3:19])=[CH:10]/[CH:11]=[CH:12]/[C:13]([CH3:20])=[CH:14]/[CH:15]=[CH:115]/[CH:114]=[C:113]([CH3:120])/[CH:112]=[CH:111]/[CH:110]=[C:109]([CH3:119])/[CH:108]=[CH:107]/[CH:106]1[C:105](=[CH2:118])[CH2:104][CH2:103][CH2:102][C:101]1([CH3:116])[CH3:117]"
    ),
    "γ,ε-carotene": Stereoparent(
        "[CH2:18]=[C:5]1[CH2:4][CH2:3][CH2:2][C:1]([CH3:16])([CH3:17])[CH:6]1/[CH:7]=[CH:8]/[C:9]([CH3:19])=[CH:10]/[CH:11]=[CH:12]/[C:13]([CH3:20])=[CH:14]/[CH:15]=[CH:115]/[CH:114]=[C:113]([CH3:120])/[CH:112]=[CH:111]/[CH:110]=[C:109]([CH3:119])/[CH:108]=[CH:107]/[CH:106]1[C:105]([CH3:118])=[CH:104][CH2:103][CH2:102][C:101]1([CH3:116])[CH3:117]"
    ),
    "γ,κ-carotene": Stereoparent(
        "[CH2:18]=[C:5]1[CH2:4][CH2:3][CH2:2][C:1]([CH3:16])([CH3:17])[CH:6]1/[CH:7]=[CH:8]/[C:9]([CH3:19])=[CH:10]/[CH:11]=[CH:12]/[C:13]([CH3:20])=[CH:14]/[CH:15]=[CH:115]/[CH:114]=[C:113]([CH3:120])/[CH:112]=[CH:111]/[CH:110]=[C:109]([CH3:119])/[CH:108]=[CH:107]/[CH2:106][C:105]1([CH3:118])[CH2:104][CH2:103][CH2:102][C:101]1([CH3:116])[CH3:117]"
    ),
    "γ,φ-carotene": Stereoparent(
        "[CH2:18]=[C:5]1[CH2:4][CH2:3][CH2:2][C:1]([CH3:16])([CH3:17])[CH:6]1/[CH:7]=[CH:8]/[C:9]([CH3:19])=[CH:10]/[CH:11]=[CH:12]/[C:13]([CH3:20])=[CH:14]/[CH:15]=[CH:115]/[CH:114]=[C:113]([CH3:120])/[CH:112]=[CH:111]/[CH:110]=[C:109]([CH3:119])/[CH:108]=[CH:107]/[c:106]1[c:105]([CH3:118])[cH:104][cH:103][c:102]([CH3:117])[c:101]1[CH3:116]"
    ),
    "γ,χ-carotene": Stereoparent(
        "[CH2:18]=[C:5]1[CH2:4][CH2:3][CH2:2][C:1]([CH3:16])([CH3:17])[CH:6]1/[CH:7]=[CH:8]/[C:9]([CH3:19])=[CH:10]/[CH:11]=[CH:12]/[C:13]([CH3:20])=[CH:14]/[CH:15]=[CH:115]/[CH:114]=[C:113]([CH3:120])/[CH:112]=[CH:111]/[CH:110]=[C:109]([CH3:119])/[CH:108]=[CH:107]/[c:106]1[cH:105][cH:104][c:103]([CH3:118])[c:102]([CH3:117])[c:101]1[CH3:116]"
    ),
    "γ,ψ-carotene": Stereoparent(
        "[CH2:18]=[C:5]1[CH2:4][CH2:3][CH2:2][C:1]([CH3:16])([CH3:17])[CH:6]1/[CH:7]=[CH:8]/[C:9]([CH3:19])=[CH:10]/[CH:11]=[CH:12]/[C:13]([CH3:20])=[CH:14]/[CH:15]=[CH:115]/[CH:114]=[C:113]([CH3:120])/[CH:112]=[CH:111]/[CH:110]=[C:109]([CH3:119])/[CH:108]=[CH:107]/[CH:106]=[C:105]([CH3:118])[CH2:104][CH2:103][CH:102]=[C:101]([CH3:116])[CH3:117]"
    ),
    "ε,ε-carotene": Stereoparent(
        "[CH3:18][C:5]1=[CH:4][CH2:3][CH2:2][C:1]([CH3:16])([CH3:17])[CH:6]1/[CH:7]=[CH:8]/[C:9]([CH3:19])=[CH:10]/[CH:11]=[CH:12]/[C:13]([CH3:20])=[CH:14]/[CH:15]=[CH:115]/[CH:114]=[C:113]([CH3:120])/[CH:112]=[CH:111]/[CH:110]=[C:109]([CH3:119])/[CH:108]=[CH:107]/[CH:106]1[C:105]([CH3:118])=[CH:104][CH2:103][CH2:102][C:101]1([CH3:116])[CH3:117]"
    ),
    "ε,κ-carotene": Stereoparent(
        "[CH3:18][C:5]1=[CH:4][CH2:3][CH2:2][C:1]([CH3:16])([CH3:17])[CH:6]1/[CH:7]=[CH:8]/[C:9]([CH3:19])=[CH:10]/[CH:11]=[CH:12]/[C:13]([CH3:20])=[CH:14]/[CH:15]=[CH:115]/[CH:114]=[C:113]([CH3:120])/[CH:112]=[CH:111]/[CH:110]=[C:109]([CH3:119])/[CH:108]=[CH:107]/[CH2:106][C:105]1([CH3:118])[CH2:104][CH2:103][CH2:102][C:101]1([CH3:116])[CH3:117]"
    ),
    "ε,φ-carotene": Stereoparent(
        "[CH3:18][C:5]1=[CH:4][CH2:3][CH2:2][C:1]([CH3:16])([CH3:17])[CH:6]1/[CH:7]=[CH:8]/[C:9]([CH3:19])=[CH:10]/[CH:11]=[CH:12]/[C:13]([CH3:20])=[CH:14]/[CH:15]=[CH:115]/[CH:114]=[C:113]([CH3:120])/[CH:112]=[CH:111]/[CH:110]=[C:109]([CH3:119])/[CH:108]=[CH:107]/[c:106]1[c:105]([CH3:118])[cH:104][cH:103][c:102]([CH3:117])[c:101]1[CH3:116]"
    ),
    "ε,χ-carotene": Stereoparent(
        "[CH3:18][C:5]1=[CH:4][CH2:3][CH2:2][C:1]([CH3:16])([CH3:17])[CH:6]1/[CH:7]=[CH:8]/[C:9]([CH3:19])=[CH:10]/[CH:11]=[CH:12]/[C:13]([CH3:20])=[CH:14]/[CH:15]=[CH:115]/[CH:114]=[C:113]([CH3:120])/[CH:112]=[CH:111]/[CH:110]=[C:109]([CH3:119])/[CH:108]=[CH:107]/[c:106]1[cH:105][cH:104][c:103]([CH3:118])[c:102]([CH3:117])[c:101]1[CH3:116]"
    ),
    "ε,ψ-carotene": Stereoparent(
        "[CH3:116][C:101]([CH3:117])=[CH:102][CH2:103][CH2:104][C:105]([CH3:118])=[CH:106]/[CH:107]=[CH:108]/[C:109]([CH3:119])=[CH:110]/[CH:111]=[CH:112]/[C:113]([CH3:120])=[CH:114]/[CH:115]=[CH:15]/[CH:14]=[C:13]([CH3:20])/[CH:12]=[CH:11]/[CH:10]=[C:9]([CH3:19])/[CH:8]=[CH:7]/[CH:6]1[C:5]([CH3:18])=[CH:4][CH2:3][CH2:2][C:1]1([CH3:16])[CH3:17]"
    ),
    "κ,κ-carotene": Stereoparent(
        "[CH3:20][C:13](/[CH:12]=[CH:11]/[CH:10]=[C:9]([CH3:19])/[CH:8]=[CH:7]/[CH2:6][C:5]1([CH3:18])[CH2:4][CH2:3][CH2:2][C:1]1([CH3:16])[CH3:17])=[CH:14]\\[CH:15]=[CH:115]\\[CH:114]=[C:113]([CH3:120])\\[CH:112]=[CH:111]\\[CH:110]=[C:109]([CH3:119])\\[CH:108]=[CH:107]\\[CH2:106][C:105]1([CH3:118])[CH2:104][CH2:103][CH2:102][C:101]1([CH3:116])[CH3:117]"
    ),
    "κ,φ-carotene": Stereoparent(
        "[CH3:20][C:13](/[CH:12]=[CH:11]/[CH:10]=[C:9]([CH3:19])/[CH:8]=[CH:7]/[CH2:6][C:5]1([CH3:18])[CH2:4][CH2:3][CH2:2][C:1]1([CH3:16])[CH3:17])=[CH:14]\\[CH:15]=[CH:115]\\[CH:114]=[C:113]([CH3:120])\\[CH:112]=[CH:111]\\[CH:110]=[C:109]([CH3:119])\\[CH:108]=[CH:107]\\[c:106]1[c:105]([CH3:118])[cH:104][cH:103][c:102]([CH3:117])[c:101]1[CH3:116]"
    ),
    "κ,χ-carotene": Stereoparent(
        "[CH3:20][C:13](/[CH:12]=[CH:11]/[CH:10]=[C:9]([CH3:19])/[CH:8]=[CH:7]/[CH2:6][C:5]1([CH3:18])[CH2:4][CH2:3][CH2:2][C:1]1([CH3:16])[CH3:17])=[CH:14]\\[CH:15]=[CH:115]\\[CH:114]=[C:113]([CH3:120])\\[CH:112]=[CH:111]\\[CH:110]=[C:109]([CH3:119])\\[CH:108]=[CH:107]\\[c:106]1[cH:105][cH:104][c:103]([CH3:118])[c:102]([CH3:117])[c:101]1[CH3:116]"
    ),
    "κ,ψ-carotene": Stereoparent(
        "[CH3:116][C:101]([CH3:117])=[CH:102][CH2:103][CH2:104][C:105]([CH3:118])=[CH:106]/[CH:107]=[CH:108]/[C:109]([CH3:119])=[CH:110]/[CH:111]=[CH:112]/[C:113]([CH3:120])=[CH:114]/[CH:115]=[CH:15]/[CH:14]=[C:13]([CH3:20])/[CH:12]=[CH:11]/[CH:10]=[C:9]([CH3:19])/[CH:8]=[CH:7]/[CH2:6][C:5]1([CH3:18])[CH2:4][CH2:3][CH2:2][C:1]1([CH3:16])[CH3:17]"
    ),
    "φ,φ-carotene": Stereoparent(
        "[CH3:20][C:13](/[CH:12]=[CH:11]/[CH:10]=[C:9]([CH3:19])/[CH:8]=[CH:7]/[c:6]1[c:5]([CH3:18])[cH:4][cH:3][c:2]([CH3:17])[c:1]1[CH3:16])=[CH:14]\\[CH:15]=[CH:115]\\[CH:114]=[C:113]([CH3:120])\\[CH:112]=[CH:111]\\[CH:110]=[C:109]([CH3:119])\\[CH:108]=[CH:107]\\[c:106]1[c:105]([CH3:118])[cH:104][cH:103][c:102]([CH3:117])[c:101]1[CH3:116]"
    ),
    "φ,χ-carotene": Stereoparent(
        "[CH3:120][C:113](/[CH:112]=[CH:111]/[CH:110]=[C:109]([CH3:119])/[CH:108]=[CH:107]/[c:106]1[cH:105][cH:104][c:103]([CH3:118])[c:102]([CH3:117])[c:101]1[CH3:116])=[CH:114]\\[CH:115]=[CH:15]\\[CH:14]=[C:13]([CH3:20])\\[CH:12]=[CH:11]\\[CH:10]=[C:9]([CH3:19])\\[CH:8]=[CH:7]\\[c:6]1[c:5]([CH3:18])[cH:4][cH:3][c:2]([CH3:17])[c:1]1[CH3:16]"
    ),
    "φ,ψ-carotene": Stereoparent(
        "[CH3:116][C:101]([CH3:117])=[CH:102][CH2:103][CH2:104][C:105]([CH3:118])=[CH:106]/[CH:107]=[CH:108]/[C:109]([CH3:119])=[CH:110]/[CH:111]=[CH:112]/[C:113]([CH3:120])=[CH:114]/[CH:115]=[CH:15]/[CH:14]=[C:13]([CH3:20])/[CH:12]=[CH:11]/[CH:10]=[C:9]([CH3:19])/[CH:8]=[CH:7]/[c:6]1[c:5]([CH3:18])[cH:4][cH:3][c:2]([CH3:17])[c:1]1[CH3:16]"
    ),
    "χ,χ-carotene": Stereoparent(
        "[CH3:20][C:13](/[CH:12]=[CH:11]/[CH:10]=[C:9]([CH3:19])/[CH:8]=[CH:7]/[c:6]1[cH:5][cH:4][c:3]([CH3:18])[c:2]([CH3:17])[c:1]1[CH3:16])=[CH:14]\\[CH:15]=[CH:115]\\[CH:114]=[C:113]([CH3:120])\\[CH:112]=[CH:111]\\[CH:110]=[C:109]([CH3:119])\\[CH:108]=[CH:107]\\[c:106]1[cH:105][cH:104][c:103]([CH3:118])[c:102]([CH3:117])[c:101]1[CH3:116]"
    ),
    "χ,ψ-carotene": Stereoparent(
        "[CH3:116][C:101]([CH3:117])=[CH:102][CH2:103][CH2:104][C:105]([CH3:118])=[CH:106]/[CH:107]=[CH:108]/[C:109]([CH3:119])=[CH:110]/[CH:111]=[CH:112]/[C:113]([CH3:120])=[CH:114]/[CH:115]=[CH:15]/[CH:14]=[C:13]([CH3:20])/[CH:12]=[CH:11]/[CH:10]=[C:9]([CH3:19])/[CH:8]=[CH:7]/[c:6]1[cH:5][cH:4][c:3]([CH3:18])[c:2]([CH3:17])[c:1]1[CH3:16]"
    ),
    "ψ,ψ-carotene": Stereoparent(
        "[CH3:16][C:1]([CH3:17])=[CH:2][CH2:3][CH2:4][C:5]([CH3:18])=[CH:6]/[CH:7]=[CH:8]/[C:9]([CH3:19])=[CH:10]/[CH:11]=[CH:12]/[C:13]([CH3:20])=[CH:14]/[CH:15]=[CH:115]/[CH:114]=[C:113]([CH3:120])/[CH:112]=[CH:111]/[CH:110]=[C:109]([CH3:119])/[CH:108]=[CH:107]/[CH:106]=[C:105]([CH3:118])[CH2:104][CH2:103][CH:102]=[C:101]([CH3:116])[CH3:117]"
    ),
}
