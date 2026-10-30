"""
Konténer lemezszabás kalkulátor

Fő változások az előző verzióhoz képest:
- Valódi MaxRects elhelyezés (minden szabad téglalapot karban tart), opcionálisan
  csak egyenes végigvágás (táblaolló / guillotine mód).
- Minden megnyitott lemez maradéka "nyitva marad": minden darab először a már
  megkezdett lemezek maradékába próbál bekerülni, csak ha sehol nem fér el,
  nyit új lemezt. Így a 2., 3., ... konténer darabjai is a korábbi lemezek
  maradékait használják.
- A hibás maradék-összevonás (merge_residuals) kikerült: az nem szomszédos
  téglalapokat is összevont, ami nem létező anyagot hozott létre.
- Lemezméret-optimalizálás: a végén minden lemezt megpróbál kisebb (olcsóbb)
  táblaméretre cserélni, ha a rajta lévő darabok elférnek rajta.
- Több rendezési sorrend és elhelyezési szabály kipróbálása, a legolcsóbb nyer.
- Túl nagy darabok darabolása (szükség esetén mindkét irányban), azonos
  darabtípusok azonos módon darabolva (nincs kombinatorikus robbanás).
- Vágási rés (fűrészrés / plazma rés) megadható.
- Cél: 1. a legjobb lemezkihasználtság (a felhasznált táblák összterülete a
  lehető legkisebb), 2. a lehető legkevesebb tábla, 3. a legkevesebb varrat.
  Az ár nem befolyásolja a választást.
- Táblaszám-csökkentés: a gyengén kihasznált lemezek darabjait megpróbálja
  más lemezek darabjaival egy táblára összevonni.
- Hegesztett toldás: ha egy darab egyben nem fér a maradékokba, egy egyenes
  varrattal KÉT részre bontva a megkezdett lemezek maradékaiba szabja
  (egy alkatrész több részből is állhat, de egyik rész sem lehet kisebb,
  mint a termék legkisebb alkatrészének területének 40%-a – a % állítható). Csak akkor toldat,
  ha ezzel táblát spórol; azonos táblaszámnál a kevesebb varrat nyer.
- Táblaméretek KIZÁRÓLAG az anyagolás.txt fájlból (opcionálisan darabárral).
- Konténertípusok: Sittes (téglalap oldalak), ABO_A (derékszögű trapéz oldalak:
  felső él, NETTÓ magasság (+55 mm perem), alsó él – a két oldal 180°-ban egymásba fordítva vagy
  külön szabva), ABO_B (orros oldalak: teljes hossz, felső és alsó egyenes él, magasság). ABO típusoknál a hátlap, homloklap
  és oldalak felső élén 30x25 mm-es ráhajlított perem (kiterítve +55 mm).
  Az alakos oldalak valós területével számol, a DXF a valódi kontúrt rajzolja.
"""

import sys
import os
import re
import time
import copy
import dataclasses
import json
import random
import math
import itertools
from dataclasses import dataclass, field
from typing import List, Tuple, Optional, Dict
import threading
import multiprocessing

# ---------- Geometria (alakos alkatrészekhez) ----------

Poly = List[Tuple[float, float]]
Shape = Tuple[str, Poly]          # (alkatrész neve, sokszög a darab saját koordinátáiban)
AREA_EPS = 1.0                    # mm² – ennél kisebb metszet = nincs metszet


HOLE_SUFFIX = " (furat)"           # az alkatrész belső kivágásának címkéje


def is_hole(label: str) -> bool:
    return label.endswith(HOLE_SUFFIX)


def poly_area(poly: Poly) -> float:
    a = 0.0
    for i in range(len(poly)):
        x1, y1 = poly[i]
        x2, y2 = poly[(i + 1) % len(poly)]
        a += x1 * y2 - x2 * y1
    return abs(a) / 2.0


def clip_poly_rect(poly: Poly, x0: float, y0: float, x1: float, y1: float) -> Poly:
    """Sutherland–Hodgman vágás tengelyirányú téglalapra."""
    def clip(pts, inside, inter):
        out = []
        for i in range(len(pts)):
            cur, prev = pts[i], pts[i - 1]
            if inside(cur):
                if not inside(prev):
                    out.append(inter(prev, cur))
                out.append(cur)
            elif inside(prev):
                out.append(inter(prev, cur))
        return out

    def ix(xc):
        return lambda p, q: (xc, p[1] + (q[1] - p[1]) * (xc - p[0]) / (q[0] - p[0]))

    def iy(yc):
        return lambda p, q: (p[0] + (q[0] - p[0]) * (yc - p[1]) / (q[1] - p[1]), yc)

    pts = list(poly)
    for inside, inter in ((lambda p: p[0] >= x0, ix(x0)), (lambda p: p[0] <= x1, ix(x1)),
                          (lambda p: p[1] >= y0, iy(y0)), (lambda p: p[1] <= y1, iy(y1))):
        if not pts:
            break
        pts = clip(pts, inside, inter)
    return pts


def cut_line_length(shapes: List[Shape], axis: str, pos: float,
                    lo: float, hi: float) -> float:
    """Egy egyenes vágás (varrat) anyagban futó hossza: axis 'H' = vízszintes
    vonal y=pos (x lo..hi), 'W' = függőleges vonal x=pos (y lo..hi)."""
    d = 0.5
    total = 0.0
    for _, poly in shapes:
        if axis == "H":
            c = clip_poly_rect(poly, lo, pos - d, hi, pos + d)
        else:
            c = clip_poly_rect(poly, pos - d, lo, pos + d, hi)
        if len(c) >= 3:
            total += poly_area(c) / (2 * d)
    return total


def clip_shapes(shapes: List[Shape], x0, y0, x1, y1) -> List[Shape]:
    """A darab alakjainak a (x0,y0)-(x1,y1) téglalapba eső része, áthelyezve (0,0)-ba."""
    out = []
    for label, poly in shapes:
        c = clip_poly_rect(poly, x0, y0, x1, y1)
        if len(c) >= 3 and poly_area(c) > AREA_EPS:
            out.append((label, [(x - x0, y - y0) for x, y in c]))
    return out


def rect_poly(w: float, h: float) -> Poly:
    return [(0, 0), (w, 0), (w, h), (0, h)]


# ---------- Adatszerkezetek ----------

@dataclass(frozen=True)
class SheetSize:
    """width/height: a VÁGHATÓ terület (erre szabunk); nom_w/nom_h: a megvásárolt
    tábla névleges mérete (0 = azonos), off: a darabok eltolása a tábla szélétől."""
    width: int
    height: int
    price: Optional[float] = None   # ha nincs ár, a terület a költség
    nom_w: int = 0
    nom_h: int = 0
    off: int = 0
    thick: float = 0.0              # lemezvastagság (mm) – a számítás végén kerül rá

    def weight(self) -> float:
        """A teljes tábla súlya (kg)."""
        return self.area / 1e6 * self.thick * STEEL_KG_PER_M2_MM

    @property
    def nw(self) -> int:
        return self.nom_w or self.width

    @property
    def nh(self) -> int:
        return self.nom_h or self.height

    @property
    def area(self) -> int:
        # a felhasznált anyag a teljes (névleges) tábla
        return self.nw * self.nh

    def label(self) -> str:
        t = f", {self.thick:g} mm vastag" if self.thick else ""
        if (self.nw, self.nh) == (self.width, self.height):
            return f"{self.nw} x {self.nh} mm{t}"
        return f"{self.nw} x {self.nh} mm (vágható {self.width} x {self.height}){t}"

    def cost(self) -> float:
        return self.price if self.price is not None else float(self.area)


@dataclass
class Piece:
    width: int
    height: int
    name: str
    unit: int = 1          # hányadik konténerhez tartozik (1..qty)
    base: str = ""         # eredeti alkatrész neve (toldott részeknél)
    part: int = 0          # 0 = egyben szabott, 1.. = hegesztett toldás része
    tag: str = ""          # megjelenített azonosító, pl. "K2" vagy "Sittes 2."
    min_area: int = 0      # toldásnál egy rész legkisebb megengedett területe (mm²)
    shapes: List[Shape] = field(default_factory=list)   # üres = teljes téglalap
    seam: float = 0.0      # a darab egyik élén lévő toldási varrat hossza (mm)
    no_weld: bool = False  # nem toldható (pl. DXF-ből betöltött egyéb alkatrész)


@dataclass
class PlacedPiece:
    sheet_index: int
    x: int
    y: int
    width: int
    height: int
    name: str
    unit: int = 1
    rotated: bool = False
    from_residual: bool = False   # korábbi konténer lemezének maradékából
    base: str = ""
    part: int = 0
    tag: str = ""
    shapes: List[Shape] = field(default_factory=list)   # természetes (forgatás előtti) állásban
    seam: float = 0.0      # toldási varrat hossza, amely ehhez a darabhoz tartozik (mm)

    @property
    def net_area(self) -> float:
        if not self.shapes:
            return float(self.width * self.height)
        return sum(-poly_area(p) if is_hole(lab) else poly_area(p) for lab, p in self.shapes)


@dataclass
class FreeRect:
    x: int
    y: int
    width: int
    height: int



# ---------- Anyagadatbázis (anyagolás.txt) ----------

def _read_text(path: str) -> str:
    # Windows Jegyzettömb gyakran ANSI (cp1250) kódolással ment
    for enc in ("utf-8-sig", "cp1250", "latin-1"):
        try:
            with open(path, encoding=enc) as f:
                return f.read()
        except UnicodeDecodeError:
            continue
    raise ValueError("Az anyagfájl kódolását nem sikerült felismerni.")


def _parse_number(tok: str) -> float:
    # "1.250" / "1 250" ezres tagolás → 1250; "1250,5" → 1250.5
    if re.fullmatch(r"\d{1,3}(\.\d{3})+", tok):
        return float(tok.replace(".", ""))
    return float(tok.replace(",", "."))


def load_sheet_sizes(path: str) -> List[SheetSize]:
    """
    Soronként egy táblaméret. CSAK az itt felsorolt méretek használhatók.
    Elfogadott formák (példák):
        1500x3000
        1500 x 3000
        1500;3000;62000        (harmadik, 100-nál nagyobb szám: darabár)
        3mm 1250x2500          (100 alatti szám, pl. vastagság, figyelmen kívül)
    '#' utáni rész megjegyzés.
    """
    sizes = []
    for raw in _read_text(path).splitlines():
        line = raw.split("#", 1)[0]
        line = re.sub(r",(?=\d{3})", " ", line)                  # "1250,2500" → két szám
        toks = re.findall(r"\d{1,3}(?:\.\d{3})+(?!\d)|\d+(?:[.,]\d+)?", line)
        nums = [_parse_number(t) for t in toks]
        big = [n for n in nums if n >= 100]                       # mm méretek / ár
        if len(big) < 2:
            continue
        w, h = int(round(big[0])), int(round(big[1]))
        price = big[2] if len(big) >= 3 else None
        sizes.append(SheetSize(min(w, h), max(w, h), price))
    if not sizes:
        raise ValueError("Az anyagfájlban nem találtam érvényes táblaméretet.")
    # ha csak néhány sorban van ár, a költség nem összevethető → területtel számolunk
    if any(s.price is None for s in sizes):
        sizes = [SheetSize(s.width, s.height, None) for s in sizes]
    # azonos méret csak egyszer
    uniq: Dict[Tuple[int, int], SheetSize] = {}
    for s in sizes:
        uniq.setdefault((s.width, s.height), s)
    return sorted(uniq.values(), key=lambda s: s.area)


DEFAULT_MAX_CUT = 2995     # mm – a lézer legnagyobb vágási hossza (lézeres vágásnál)
STEEL_KG_PER_M2_MM = 7.85  # acél: 7,85 kg/dm³ → 1 m² 1 mm vastag lemez 7,85 kg
DEFAULT_MARGIN = 0         # mm – széltávolság a tábla minden szélén


def usable_sheets(sheets: List[SheetSize], margin: float = 0,
                  max_cut: Optional[float] = None) -> List[SheetSize]:
    """A táblák vágható területe: minden szélen 'margin' mm marad, és egyik irányban
    sem lehet hosszabb a gép legnagyobb vágási hosszánál (pl. lézer: 2995 mm)."""
    m = int(round(margin or 0))
    out = []
    for s in sheets:
        w, h = s.nw, s.nh
        if max_cut:
            w, h = min(w, int(max_cut)), min(h, int(max_cut))
        w, h = w - 2 * m, h - 2 * m
        if w <= 0 or h <= 0:
            raise ValueError(f"A széltávolság ({m} mm) túl nagy a {s.nw}x{s.nh} táblához.")
        out.append(SheetSize(min(w, h), max(w, h), s.price, s.nw, s.nh, m))
    return out


def piece_shapes(p: Piece) -> List[Shape]:
    return p.shapes if p.shapes else [(p.base or p.name, rect_poly(p.width, p.height))]


def is_plain_rect(p: Piece) -> bool:
    sh = piece_shapes(p)
    return len(sh) == 1 and abs(poly_area(sh[0][1]) - p.width * p.height) < AREA_EPS


# ---------- Valódi alakos elhelyezés: a darab foglalt területe ----------

_NEST = [True]           # alakos elhelyezés be/ki (a keresés váltogatja: a gyorsabb
                         # befoglaló-téglalapos és a pontosabb alakos módot)


def set_nesting(on: bool):
    _NEST[0] = bool(on)


COVER_STEPS = 5          # a ferde élek mentén ennyi lépcsőre bontjuk a foglalt sávot
_COVER_CACHE: Dict[int, tuple] = {}
_COVER_CONTENT: Dict[tuple, Optional[list]] = {}


def _poly_cover(poly: Poly, w: int, h: int, horizontal: bool) -> List[Tuple[int, int, int, int]]:
    """A sokszöget lefedő lépcsős téglalapok (x0, y0, x1, y1): vízszintes (vagy
    függőleges) sávokra bontja, minden sávban a sokszög ott lévő szélességét foglalja.
    A ferde élek mentén COVER_STEPS lépcső, így a levágott háromszögek nagy része szabad."""
    ax = 1 if horizontal else 0
    other = 1 - ax
    lim = h if horizontal else w
    brk = sorted(set(min(max(p[ax], 0), lim) for p in poly))
    cuts = []
    for a, b in zip(brk, brk[1:]):
        if b - a < 1:
            continue
        n = COVER_STEPS if b - a >= 4 * COVER_STEPS else 1
        cuts.extend(a + (b - a) * i / n for i in range(n))
    if brk:
        cuts.append(brk[-1])
    strips = []
    for a, b in zip(cuts, cuts[1:]):
        lo_a, hi_a = math.floor(a + 1e-6), math.ceil(b - 1e-6)
        if hi_a <= lo_a:
            continue
        c = (clip_poly_rect(poly, -1e9, a, 1e9, b) if horizontal
             else clip_poly_rect(poly, a, -1e9, b, 1e9))
        if len(c) < 3:
            continue
        lo_o = max(0, math.floor(min(p[other] for p in c) + 1e-6))
        hi_o = min(w if horizontal else h, math.ceil(max(p[other] for p in c) - 1e-6))
        if hi_o <= lo_o:
            continue
        if strips and strips[-1][2:] == (lo_o, hi_o) and strips[-1][1] >= lo_a:
            strips[-1] = (strips[-1][0], hi_a, lo_o, hi_o)        # azonos szélesség: összevonás
        else:
            strips.append((lo_a, hi_a, lo_o, hi_o))
    if horizontal:
        return [(x0, y0, x1, y1) for y0, y1, x0, x1 in strips]
    return [(x0, y0, x1, y1) for x0, x1, y0, y1 in strips]


def piece_cover(p) -> Optional[List[Tuple[int, int, int, int]]]:
    """A darab (Piece/PlacedPiece, természetes állásban) ténylegesen foglalt területe
    téglalapok listájaként, vagy None, ha a befoglaló téglalap (nem éri meg bontani)."""
    if not p.shapes or p.part > 0 or not _NEST[0]:
        return None           # téglalap, hegesztett toldás csíkja, vagy kikapcsolva
    hit = _COVER_CACHE.get(id(p.shapes))
    if hit is not None and hit[0] is p.shapes and hit[1] == (p.width, p.height):
        return hit[2]
    w, h = (p.height, p.width) if getattr(p, "rotated", False) else (p.width, p.height)
    # tartalom szerinti gyorsítótár (a másolt darabok alakja új lista, de ugyanaz)
    ckey = (w, h, tuple(tuple(pt for pt in poly) for _, poly in p.shapes))
    res = _COVER_CONTENT.get(ckey, False)
    if res is not False:
        _COVER_CACHE[id(p.shapes)] = (p.shapes, (p.width, p.height), res)
        return res
    rects = []
    for _, poly in p.shapes:
        a = _poly_cover(poly, w, h, True)
        b = _poly_cover(poly, w, h, False)
        area = lambda rs: sum((x1 - x0) * (y1 - y0) for x0, y0, x1, y1 in rs)
        rects.extend(a if area(a) <= area(b) else b)
    # egymást teljesen tartalmazó téglalapok kiszűrése
    rects = [r for i, r in enumerate(rects)
             if not any(j != i and o[0] <= r[0] and o[1] <= r[1] and o[2] >= r[2] and o[3] >= r[3]
                        and (o != r or j < i) for j, o in enumerate(rects))]
    covered = sum((x1 - x0) * (y1 - y0) for x0, y0, x1, y1 in rects)
    res = None if (not rects or covered >= 0.97 * w * h) else rects
    if len(_COVER_CACHE) > 20000:
        _COVER_CACHE.clear()
    if len(_COVER_CONTENT) > 5000:
        _COVER_CONTENT.clear()
    _COVER_CACHE[id(p.shapes)] = (p.shapes, (p.width, p.height), res)
    _COVER_CONTENT[ckey] = res
    return res


_POCKET_CACHE: Dict[int, tuple] = {}


def piece_pockets(p) -> List[Tuple[int, int]]:
    """Az alakos darab befoglaló téglalapjában maradó legnagyobb szabad téglalapok
    méretei (rövidebb, hosszabb oldal) – ide férhet más darab."""
    cv = piece_cover(p)
    if cv is None:
        return []
    hit = _POCKET_CACHE.get(id(p.shapes))
    if hit is not None and hit[0] is p.shapes:
        return hit[1]
    w, h = (p.height, p.width) if getattr(p, "rotated", False) else (p.width, p.height)
    pk = SheetPacker(SheetSize(w, h), 0, 0, False)
    pk._split_maxrects([FreeRect(x0, y0, x1 - x0, y1 - y0) for x0, y0, x1, y1 in cv])
    res = [(min(r.width, r.height), max(r.width, r.height)) for r in pk.free_rects]
    if len(_POCKET_CACHE) > 20000:
        _POCKET_CACHE.clear()
    _POCKET_CACHE[id(p.shapes)] = (p.shapes, res)
    return res


def smallest_dims(pieces) -> Optional[Tuple[int, int]]:
    """A legkisebb darab (rövidebb, hosszabb oldal) – a zsebek hasznosságához.
    None, ha több alakos darab is van (azok egymás zsebeibe is beférhetnek)."""
    if not pieces:
        return (0, 0)
    if sum(1 for p in pieces if piece_cover(p) is not None) >= 2:
        return None
    return min(((min(p.width, p.height), max(p.width, p.height)) for p in pieces),
               key=lambda d: (d[0] * d[1], d))


def occ_area(p, smallest: Optional[Tuple[int, int]] = None) -> int:
    """A darab által ténylegesen elfoglalt (lépcsős) terület – becslésekhez. Ha
    'smallest' meg van adva, a lépcsős területet csak akkor számolja, ha a legkisebb
    darab befér a darab mellett maradó zsebek valamelyikébe (különben a zseb
    úgysem hasznosítható, és a befoglaló téglalap a reális becslés)."""
    c = piece_cover(p)
    if c is None:
        return p.width * p.height
    if smallest is not None:
        a, b = smallest
        if not any(a <= pw and b <= ph for pw, ph in piece_pockets(p)):
            return p.width * p.height
    return sum((x1 - x0) * (y1 - y0) for x0, y0, x1, y1 in c)


_ROT180_CACHE: Dict[int, tuple] = {}


def shapes_rot180(p: Piece) -> List[Shape]:
    """A darab alakjai 180°-kal elforgatva (a befoglaló téglalap ugyanaz)."""
    hit = _ROT180_CACHE.get(id(p.shapes))
    if hit is not None and hit[0] is p.shapes:
        return hit[1]
    w, h = p.width, p.height
    res = [(lab, [(w - x, h - y) for x, y in poly]) for lab, poly in p.shapes]
    if len(_ROT180_CACHE) > 20000:
        _ROT180_CACHE.clear()
    _ROT180_CACHE[id(p.shapes)] = (p.shapes, res)
    return res


# ---------- Egy lemez elhelyező (MaxRects / guillotine) ----------

class SheetPacker:
    """
    A vágási rést úgy kezeljük, hogy minden darabot 'kerf'-fel megnövelünk,
    a lemezt pedig szintén 'kerf'-fel (így a lemez szélén nem kell rés).
    """

    def __init__(self, sheet: SheetSize, sheet_index: int, kerf: int = 0,
                 guillotine: bool = False, pocket_min: Optional[Tuple[int, int]] = None):
        self.sheet = sheet
        # a legkisebb darab mérete: alakos darab mellé csak akkor hagyunk zsebet,
        # ha ez befér (különben fölösleges felaprózni a szabad területet)
        self.pocket_min = pocket_min
        self.sheet_index = sheet_index
        self.kerf = kerf
        self.guillotine = guillotine
        self.free_rects: List[FreeRect] = [
            FreeRect(0, 0, sheet.width + kerf, sheet.height + kerf)
        ]
        self.placed: List[PlacedPiece] = []
        self.used_area = 0
        self.has_pockets = False     # van-e már alakos darab (zsebekkel) a lemezen

    # --- pontozás: kisebb = jobb ---
    @staticmethod
    def _score(fr: FreeRect, w: int, h: int, heuristic: str):
        dw, dh = fr.width - w, fr.height - h
        if heuristic == "BSSF":     # legkisebb rövid oldali maradék
            return (min(dw, dh), max(dw, dh))
        if heuristic == "BAF":      # legkisebb terület-maradék
            return (fr.width * fr.height - w * h, min(dw, dh))
        # BL: legalsó, majd legbaloldalibb pozíció
        return (fr.y + h, fr.x)

    def find_position(self, piece: Piece, heuristic: str):
        """Visszaadja (score, free_rect, w, h, rotated) vagy None."""
        k = self.kerf
        best = None
        for fr in self.free_rects:
            for w, h, rot in ((piece.width, piece.height, False),
                              (piece.height, piece.width, True)):
                if rot and piece.width == piece.height:
                    continue
                if w + k <= fr.width and h + k <= fr.height:
                    s = self._score(fr, w + k, h + k, heuristic)
                    if best is None or s < best[0]:
                        best = (s, fr, w, h, rot)
        if best is None and not self.guillotine and self.has_pockets:
            best = self.find_nested(piece)
        return best

    def find_nested(self, piece: Piece):
        """Alakos darab egymásba illesztése: a darab befoglaló téglalapja nem fér
        el szabad helyen, de a valódi (lépcsős) alakja igen – pl. egy másik
        ferde/orros darab mellett maradt sarokban. 0°, 90°, 180° és 270°-os
        állásokat is kipróbál. Visszaad: (score, None, w, h, rot, x, y, shapes)."""
        if piece_cover(piece) is None:
            return None
        k = self.kerf
        SW, SH = self.sheet.width + k, self.sheet.height + k
        free = self.free_rects
        if occ_area(piece) > sum(o.width * o.height for o in free):
            return None
        # gyors kizárás: a legnagyobb foglalt sávnak valamelyik szabad helyre be kell férnie
        big = max(piece_cover(piece), key=lambda c: (c[2] - c[0]) * (c[3] - c[1]))
        bw, bh = big[2] - big[0] + k, big[3] - big[1] + k
        if not any((bw <= o.width and bh <= o.height) or (bh <= o.width and bw <= o.height)
                   for o in free):
            return None
        best = None
        variants = [(piece.shapes, piece_cover(piece))]
        r180 = shapes_rot180(piece)
        tmp = Piece(piece.width, piece.height, "", shapes=r180)
        variants.append((r180, piece_cover(tmp)))
        for shapes, cv in variants:
            for rot in (False, True):
                if rot:
                    w, h = piece.height, piece.width
                    rects = [(piece.height - y1, x0, piece.height - y0, x1) for x0, y0, x1, y1 in cv]
                else:
                    w, h = piece.width, piece.height
                    rects = cv
                if w + k > SW or h + k > SH:
                    continue
                for fr in free:
                    for c in rects:
                        x, y = fr.x - c[0], fr.y - c[1]
                        if x < 0 or y < 0 or x + w + k > SW or y + h + k > SH:
                            continue
                        sc = (y + h, x)
                        if best is not None and sc >= best[0]:
                            continue
                        ok = True
                        for x0, y0, x1, y1 in rects:
                            ax0, ay0, ax1, ay1 = x + x0, y + y0, x + x1 + k, y + y1 + k
                            if not any(o.x <= ax0 and o.y <= ay0 and ax1 <= o.x + o.width
                                       and ay1 <= o.y + o.height for o in free):
                                ok = False
                                break
                        if ok:
                            best = (sc, (x, y, w, h, rot, shapes, rects))
        if best is None:
            return None
        x, y, w, h, rot, shapes, rects = best[1]
        return ((10 ** 12, best[0]), None, w, h, rot, x, y, shapes)

    def place(self, piece: Piece, pos) -> PlacedPiece:
        _, fr, w, h, rot = pos[:5]
        k = self.kerf
        if len(pos) > 5:
            # egymásba illesztett alakos darab (find_nested)
            x, y, shapes = pos[5], pos[6], pos[7]
            piece = Piece(piece.width, piece.height, piece.name, piece.unit, piece.base,
                          piece.part, piece.tag, piece.min_area, shapes, piece.seam)
        else:
            x, y = fr.x, fr.y
        cover = None if self.guillotine else piece_cover(piece)
        if cover is not None and self.pocket_min is not None:
            a, b = self.pocket_min
            if not any(a + k <= pw and b + k <= ph for pw, ph in piece_pockets(piece)):
                cover = None
        if self.guillotine:
            self._split_guillotine(fr, w + k, h + k)
        elif cover is None:
            self._split_maxrects([FreeRect(x, y, w + k, h + k)])
        else:
            # alakos darab: csak a valóban foglalt (lépcsős) részt vesszük el, a
            # ferde élek melletti üres sarkokba más darab kerülhet
            nh = w            # elforgatott darab természetes magassága
            rs = []
            for x0, y0, x1, y1 in cover:
                if rot:
                    x0, y0, x1, y1 = nh - y1, x0, nh - y0, x1
                rs.append(FreeRect(x + x0, y + y0, x1 - x0 + k, y1 - y0 + k))
            self._split_maxrects(rs)
            self.has_pockets = True
        pp = PlacedPiece(self.sheet_index, x, y, w, h, piece.name,
                         unit=piece.unit, rotated=rot,
                         base=piece.base or piece.name, part=piece.part,
                         tag=piece.tag or f"K{piece.unit}",
                         shapes=piece_shapes(piece), seam=piece.seam)
        self.placed.append(pp)
        self.used_area += w * h
        return pp

    # --- MaxRects: minden metszett szabad téglalapból max. 4 új ---
    def _split_maxrects(self, used_list: List[FreeRect]):
        free = self.free_rects
        for used in used_list:
            new_free = []
            for fr in free:
                if (used.x >= fr.x + fr.width or used.x + used.width <= fr.x or
                        used.y >= fr.y + fr.height or used.y + used.height <= fr.y):
                    new_free.append(fr)
                    continue
                if used.x > fr.x:
                    new_free.append(FreeRect(fr.x, fr.y, used.x - fr.x, fr.height))
                if used.x + used.width < fr.x + fr.width:
                    nx = used.x + used.width
                    new_free.append(FreeRect(nx, fr.y, fr.x + fr.width - nx, fr.height))
                if used.y > fr.y:
                    new_free.append(FreeRect(fr.x, fr.y, fr.width, used.y - fr.y))
                if used.y + used.height < fr.y + fr.height:
                    ny = used.y + used.height
                    new_free.append(FreeRect(fr.x, ny, fr.width, fr.y + fr.height - ny))
            free = new_free if len(used_list) == 1 else self._prune(new_free)
        self.free_rects = self._prune(free) if len(used_list) == 1 else free

    # --- Guillotine: egyenes végigvágás, a rövidebb maradék mentén ---
    def _split_guillotine(self, fr: FreeRect, w: int, h: int):
        self.free_rects.remove(fr)
        dw, dh = fr.width - w, fr.height - h
        if dw < dh:   # vízszintes vágás a teljes szélességben
            right = FreeRect(fr.x + w, fr.y, dw, h)
            top = FreeRect(fr.x, fr.y + h, fr.width, dh)
        else:         # függőleges vágás a teljes magasságban
            right = FreeRect(fr.x + w, fr.y, dw, fr.height)
            top = FreeRect(fr.x, fr.y + h, w, dh)
        for r in (right, top):
            if r.width > 0 and r.height > 0:
                self.free_rects.append(r)

    @staticmethod
    def _prune(rects: List[FreeRect]) -> List[FreeRect]:
        """A más szabad téglalapban teljesen benne lévők elhagyása (az azonosakból
        az elsőt tartjuk meg); a sorrend az eredeti marad."""
        rects = [r for r in rects if r.width > 0 and r.height > 0]
        n = len(rects)
        order = sorted(range(n), key=lambda i: (-(rects[i].width * rects[i].height), i))
        kept: List[FreeRect] = []
        flag = [False] * n
        for i in order:
            r = rects[i]
            rx1, ry1 = r.x + r.width, r.y + r.height
            for o in kept:
                if o.x <= r.x and o.y <= r.y and rx1 <= o.x + o.width and ry1 <= o.y + o.height:
                    break
            else:
                kept.append(r)
                flag[i] = True
        return [r for i, r in enumerate(rects) if flag[i]]


# ---------- Túl nagy darabok darabolása ----------

def fits_any(w: int, h: int, sheets: List[SheetSize], kerf: int = 0) -> bool:
    return any((w <= s.width and h <= s.height) or (h <= s.width and w <= s.height)
               for s in sheets)


def _split_len(total: int, limit: int, mode: str) -> List[int]:
    n = math.ceil(total / limit)
    if n <= 1:
        return [total]
    if mode == "equal":
        base = total // n
        parts = [base] * n
        parts[-1] += total - base * n
        return parts
    # "max": maximális darabok + maradék
    parts = [limit] * (n - 1)
    parts.append(total - limit * (n - 1))
    return parts


MAX_PARTS = 8   # egy alkatrész ennyi részből állhat legfeljebb (a fő korlát a min. részterület)


class Cancelled(Exception):
    """A felhasználó megszakította a számítást."""


_CANCEL = threading.Event()     # a GUI állítja be a „Megszakítás” gombbal


def _check_cancel():
    if _CANCEL.is_set():
        raise Cancelled("A számítás megszakítva.")


def generate_cut_variants(w: int, h: int, sheets: List[SheetSize],
                          min_area: int = 0) -> List[List[Tuple[int, int]]]:
    """
    Egy túl nagy darabtípus kötelező kettévágási változatai (a többrészes
    felosztást a toldásos kitöltés végzi, a min. részterület betartásával).
    """
    if fits_any(w, h, sheets):
        return [[(w, h)]]
    variants = []
    seen = set()
    for total, other, along_w in ((w, h, True), (h, w, False)):
        for sh in sheets:
            for limit_cut, limit_other in ((sh.width, sh.height), (sh.height, sh.width)):
                if other > limit_other or total > 2 * limit_cut:
                    continue
                for s1 in (limit_cut, total - limit_cut, total // 2):
                    s2 = total - s1
                    if not (0 < s1 <= limit_cut and 0 < s2):
                        continue
                    if min(s1, s2) * other < min_area:
                        continue
                    parts = [(s1, other), (s2, other)] if along_w else [(other, s1), (other, s2)]
                    if not all(fits_any(a, b, sheets) for a, b in parts):
                        continue
                    key = tuple(sorted(parts))
                    if key not in seen:
                        seen.add(key)
                        variants.append(parts)
    return variants


# ---------- Több lemezes elhelyezés ----------

SORTS = {
    "terület": lambda p: (p.width * p.height, max(p.width, p.height)),
    "hosszabb oldal": lambda p: (max(p.width, p.height), p.width * p.height),
    "kerület": lambda p: (p.width + p.height, p.width * p.height),
}


def _pack_single(pieces: List[Piece], sheet: SheetSize, kerf: int, guillotine: bool,
                 sheet_index: int) -> Optional[List[PlacedPiece]]:
    """Megpróbálja az összes darabot egyetlen adott méretű lemezre tenni."""
    sm = smallest_dims(pieces)
    for sort_key in SORTS.values():
        for heur in ("BSSF", "BAF", "BL"):
            pk = SheetPacker(sheet, sheet_index, kerf, guillotine, sm)
            ok = True
            for p in sorted(pieces, key=sort_key, reverse=True):
                pos = pk.find_position(p, heur)
                if pos is None:
                    ok = False
                    break
                pk.place(p, pos)
            if ok:
                return pk.placed
    return None


def _to_pieces(placed: List[PlacedPiece]) -> List[Piece]:
    return [Piece(pp.height if pp.rotated else pp.width,
                  pp.width if pp.rotated else pp.height,
                  pp.name, pp.unit, pp.base, pp.part, pp.tag, 0, pp.shapes, pp.seam)
            for pp in placed]


def _merge_sheets(groups, sheets: List[SheetSize], kerf: int, guillotine: bool):
    """
    Táblaszám-csökkentés: a leggyengébben kihasznált lemez darabjait
    megpróbáljuk egy másik lemez darabjaival együtt egyetlen (legnagyobb)
    táblára szabni. Ha sikerül, eggyel kevesebb tábla kell. Ismétli, amíg talál.
    """
    big = max(sheets, key=lambda s: s.area)
    groups = [(s, _to_pieces(pl), pl) for s, pl in groups]

    sm = smallest_dims([p for _, ps, _ in groups for p in ps])

    def area(ps):
        return sum(occ_area(p, sm) for p in ps)

    improved = True
    while improved and len(groups) > 1:
        improved = False
        order = sorted(range(len(groups)), key=lambda i: area(groups[i][1]))
        for i in order:
            for j in sorted(range(len(groups)), key=lambda j: area(groups[j][1])):
                if i == j or area(groups[i][1]) + area(groups[j][1]) > big.area:
                    continue
                if big.area > groups[i][0].area + groups[j][0].area:
                    continue   # az összevonás nőne a felhasznált táblaterületen
                combined = groups[i][1] + groups[j][1]
                placed = _pack_single(combined, big, kerf, guillotine, 0)
                if placed is not None:
                    new = (big, combined, placed)
                    groups = [g for k, g in enumerate(groups) if k not in (i, j)] + [new]
                    improved = True
                    break
            if improved:
                break
    return [(s, pl) for s, _, pl in groups]


def count_welds(placed: List[PlacedPiece]) -> int:
    """
    Hegesztési varratok száma: alkatrészenként (hány lemezdarabból áll - 1).
    Alakos (pl. páros oldal) daraboknál az alkatrészeket külön számolja.
    """
    groups: Dict[Tuple[int, str], int] = {}
    for pp in placed:
        labels = {lab for lab, _ in pp.shapes} if pp.shapes else {pp.base or pp.name}
        for lab in labels:
            groups[(pp.unit, lab)] = groups.get((pp.unit, lab), 0) + 1
    return sum(n - 1 for n in groups.values())


def _cut_ok(shapes: List[Shape], region_a, region_b, min_area: int) -> bool:
    """Igaz, ha egy vágás után minden alkatrész-darab vagy üres, vagy >= min_area."""
    for _, poly in shapes:
        for (x0, y0, x1, y1) in (region_a, region_b):
            if x1 <= x0 or y1 <= y0:
                continue
            c = clip_poly_rect(poly, x0, y0, x1, y1)
            a = poly_area(c) if len(c) >= 3 else 0.0
            if AREA_EPS < a < min_area:
                return False
    return True


def _fill_with_parts(piece: Piece, packers: List[SheetPacker], weld_min: int,
                     max_parts: int, new_sheet: Optional[SheetSize], kerf: int,
                     guillotine: bool) -> Optional[List[SheetPacker]]:
    """
    A darabot hegesztett toldásokkal, csíkokra bontva a megkezdett lemezek
    maradékaiba szabja. A csíkok a darab teljes szélességében/magasságában
    futnak (egyenes varrat). Egyik alkatrész-rész sem lehet kisebb a darab
    min_area értékénél (alakos daraboknál a valós, ferde élekkel vágott
    területet számolja). Ha new_sheet meg van adva, szükség esetén új lemezt
    is nyithat. Siker esetén az új packer-listát adja vissza.
    """
    if piece.part > 0 or piece.no_weld:
        return None   # már egy toldott rész (vagy nem toldható darab) – tovább nem bontható
    pks = copy.deepcopy(packers)
    shapes = piece_shapes(piece)
    simple = is_plain_rect(piece)
    W, H = piece.width, piece.height
    x0, y0 = 0, 0          # a még elhelyezetlen maradék helye a darab saját koordinátáiban
    k = kerf
    part = 0
    just_opened = False

    def total_of(ax):
        return H if ax == "H" else W

    while W > 0 and H > 0:
        if part >= max_parts:
            return None
        best = None
        memo: Dict[tuple, Optional[int]] = {}    # (tengely, hely) → vágási hossz
        for pk in pks:
            for fr in pk.free_rects:
                aw, ah = fr.width - k, fr.height - k
                if aw <= 0 or ah <= 0:
                    continue
                opts = []
                if W <= aw:
                    opts.append(("H", ah, False))
                if W <= ah:
                    opts.append(("H", aw, True))
                if H <= ah:
                    opts.append(("W", aw, False))
                if H <= aw:
                    opts.append(("W", ah, True))
                for axis, avail, rot in opts:
                    total = H if axis == "H" else W
                    keep = W if axis == "H" else H
                    mk = (axis, avail)
                    if mk in memo:
                        s = memo[mk]
                        if s is None:
                            continue
                        sc = (s * keep, -pk.sheet_index)
                        if best is None or sc > best[0]:
                            best = (sc, pk, fr, axis, s, keep, rot)
                        continue
                    memo[mk] = None
                    s = min(total, avail)
                    if part == max_parts - 1 and s < total:
                        continue          # az utolsó résznek egyben kell elférnie
                    if simple:
                        lim = max(1, math.ceil(piece.min_area / keep))
                        rest = total - s
                        if 0 < rest < lim:
                            s = total - lim
                        if s <= 0 or (s < total and s < lim):
                            continue
                    else:
                        # alakos darab: a vágási helyet lépésenként csökkentjük,
                        # amíg minden alkatrész-darab eléri a min. területet
                        step = max(10, total // 60)
                        found = False
                        while s > 0:
                            if axis == "H":
                                ra = (x0, y0, x0 + W, y0 + s)
                                rb = (x0, y0 + s, x0 + W, y0 + H)
                            else:
                                ra = (x0, y0, x0 + s, y0 + H)
                                rb = (x0 + s, y0, x0 + W, y0 + H)
                            has_material = any(
                                poly_area(c) > AREA_EPS
                                for c in (clip_poly_rect(pl, *ra) for _, pl in shapes) if len(c) >= 3)
                            if has_material and (s == total or _cut_ok(shapes, ra, rb, piece.min_area)):
                                found = True
                                break
                            s -= step
                        if not found:
                            continue
                    memo[mk] = s
                    sc = (s * keep, -pk.sheet_index)
                    if best is None or sc > best[0]:
                        best = (sc, pk, fr, axis, s, keep, rot)
        if best is None:
            if new_sheet is None or just_opened:
                return None     # még egy üres táblára sem vágható szabályosan
            pks.append(SheetPacker(new_sheet, len(pks), kerf, guillotine,
                                   pks[0].pocket_min if pks else None))
            just_opened = True
            continue
        just_opened = False
        _, pk, fr, axis, s, keep, rot = best
        part += 1
        if axis == "H":
            sw, sh = keep, s       # természetes állás: (W x s)
            region = (x0, y0, x0 + W, y0 + s)
        else:
            sw, sh = s, keep       # (s x H)
            region = (x0, y0, x0 + s, y0 + H)
        content = clip_shapes(shapes, *region)
        seam = 0.0
        if s < total_of(axis):
            if axis == "H":
                seam = cut_line_length(shapes, "H", y0 + s, x0, x0 + W)
            else:
                seam = cut_line_length(shapes, "W", x0 + s, y0, y0 + H)
        if not simple and content:
            # a csík lemezét az anyagot tartalmazó részre szűkítjük
            bx0 = min(x for _, pl in content for x, _ in pl)
            by0 = min(y for _, pl in content for _, y in pl)
            bx0, by0 = max(0, math.floor(bx0 + 1e-6)), max(0, math.floor(by0 + 1e-6))
            content = [(lab, [(x - bx0, y - by0) for x, y in pl]) for lab, pl in content]
            sw = min(sw, math.ceil(max(x for _, pl in content for x, _ in pl) - 1e-6))
            sh = min(sh, math.ceil(max(y for _, pl in content for _, y in pl) - 1e-6))
        pw, ph = (sh, sw) if rot else (sw, sh)
        strip = Piece(sw, sh, f"{piece.name} #{part}", piece.unit,
                      piece.base or piece.name, part, piece.tag, 0, content, seam)
        pk.place(strip, (None, fr, pw, ph, rot))
        if axis == "H":
            H -= s
            y0 += s
        else:
            W -= s
            x0 += s
        if not simple and W > 0 and H > 0:
            # a maradék tartományt az anyagot tartalmazó befoglalóra szűkítjük
            rest = clip_shapes(shapes, x0, y0, x0 + W, y0 + H)
            if not rest:
                break
            rx0 = max(x0, math.floor(x0 + min(x for _, pl in rest for x, _ in pl) + 1e-6))
            ry0 = max(y0, math.floor(y0 + min(y for _, pl in rest for _, y in pl) + 1e-6))
            rx1 = min(x0 + W, math.ceil(x0 + max(x for _, pl in rest for x, _ in pl) - 1e-6))
            ry1 = min(y0 + H, math.ceil(y0 + max(y for _, pl in rest for _, y in pl) - 1e-6))
            x0, y0, W, H = rx0, ry0, rx1 - rx0, ry1 - ry0
    if part == 1:
        return None   # egyben ment volna – azt a normál út kezeli
    return pks


def pack_pieces_into_sheets(pieces: List[Piece], sheets: List[SheetSize], kerf: int,
                            guillotine: bool, sort_key, heuristic: str,
                            new_sheet: SheetSize, weld_min: int = 0,
                            max_parts: int = MAX_PARTS):
    """
    Minden darab a már megnyitott lemezek közül oda kerül, ahol a legjobban
    illeszkedik (ez a maradékfelhasználás). Ha egyben sehol nem fér el, és a
    toldás engedélyezett, hegesztett részekből próbálja a maradékokba szabni.
    Új lemez csak akkor nyílik, ha ez sem megy.
    """
    packers: List[SheetPacker] = []
    biggest = max(sheets, key=lambda s: s.area)
    sm = smallest_dims(pieces)
    for p in sorted(pieces, key=sort_key, reverse=True):
        best = None
        for pk in packers:
            pos = pk.find_position(p, heuristic)
            if pos is not None and (best is None or pos[0] < best[1][0]):
                best = (pk, pos)
        if best is None and weld_min > 0 and packers:
            filled = _fill_with_parts(p, packers, weld_min, max_parts, None, kerf, guillotine)
            if filled is not None:
                packers = filled
                continue
        if best is None:
            # új lemez: az alapértelmezett méret, vagy ha nem fér rá, a legkisebb ami jó
            candidates = [new_sheet] + sorted(sheets, key=lambda s: s.area)
            pk = None
            for s in candidates:
                trial = SheetPacker(s, len(packers), kerf, guillotine, sm)
                pos = trial.find_position(p, heuristic)
                if pos is not None:
                    pk = trial
                    break
            if pk is None:
                if weld_min > 0:
                    # túl nagy darab: toldással, szükség szerint új lemez(ek)kel
                    filled = _fill_with_parts(p, packers, weld_min, MAX_PARTS, biggest,
                                              kerf, guillotine)
                    if filled is not None:
                        packers = filled
                        continue
                raise ValueError(f"A(z) {p.name} ({p.width}x{p.height}) egyik táblára sem fér rá.")
            packers.append(pk)
            best = (pk, pos)
        best[0].place(p, best[1])
    packers = [pk for pk in packers if pk.placed]
    for i, pk in enumerate(packers):
        pk.sheet_index = i
        for pp in pk.placed:
            pp.sheet_index = i

    groups = [(pk.sheet, pk.placed) for pk in packers]
    groups = _merge_sheets(groups, sheets, kerf, guillotine)
    groups.sort(key=lambda g: (min(p.unit for p in g[1]), -sum(p.width * p.height for p in g[1])))

    # Lemezméret-csökkentés: a darabszám marad, de minden lemezt a lehető
    # legkisebb (legolcsóbb) táblára teszünk → jobb kihasználtság
    result_sheets: List[SheetSize] = []
    result_placed: List[PlacedPiece] = []
    for idx, (sheet, placed) in enumerate(groups):
        orig = _to_pieces(placed)
        sm = smallest_dims(orig)
        for s in sorted(sheets, key=lambda s: s.area):
            if s.area >= sheet.area or s.area < sum(occ_area(p, sm) for p in orig):
                continue
            repacked = _pack_single(orig, s, kerf, guillotine, idx)
            if repacked is not None:
                sheet, placed = s, repacked
                break
        for pp in placed:
            pp.sheet_index = idx
        result_sheets.append(sheet)
        result_placed.extend(placed)

    _mark_residual(result_placed)
    return result_placed, result_sheets


def _mark_residual(result_placed: List[PlacedPiece]):
    """A lemez annak a konténernek a táblája, amelyiknek a legnagyobb darabja
    rajta van; a többi konténer darabjai ennek maradékából készülnek."""
    owner: Dict[int, Tuple[int, int]] = {}
    for pp in result_placed:
        a = pp.width * pp.height
        cur = owner.get(pp.sheet_index)
        if cur is None or (a, -pp.unit) > (cur[0], -cur[1]):
            owner[pp.sheet_index] = (a, pp.unit)
    for pp in result_placed:
        pp.from_residual = pp.unit != owner[pp.sheet_index][1]


def _sheet_mixes(sheets: List[SheetSize], need: float, below: float, max_count: int = 12):
    """Táblaméret-kombinációk (multihalmazok), amelyek összterülete need ≤ T < below,
    növekvő összterület (majd darabszám) szerint."""
    sizes = sorted(sheets, key=lambda s: -s.area)
    out = []

    def rec(i, chosen, area):
        if area >= below or len(chosen) > max_count:
            return
        if area >= need:
            out.append(list(chosen))
        if i == len(sizes):
            return
        for j in range(i, len(sizes)):
            chosen.append(sizes[j])
            rec(j, chosen, area + sizes[j].area)
            chosen.pop()

    rec(0, [], 0)
    out.sort(key=lambda m: (sum(s.area for s in m), len(m)))
    return out


def exact_repack(pieces: List[Piece], sheets: List[SheetSize], kerf: int, guillotine: bool,
                 below_area: float, deadline: float):
    """
    Pontosabb keresés kevés darabnál: sorra veszi a táblaméret-kombinációkat
    a legkisebb összterülettől felfelé (csak a jelenlegi legjobbnál kisebbeket),
    és mélységi kereséssel szétosztja a darabokat a táblákra. Az első sikeres
    kombináció a legjobb kihasználtságú. Visszaad: (placed, used) vagy None.
    Gyorsítások: az egy táblára férés eredménye gyorsítótárban (a darabok méretei
    alapján), az azonos méretű darabok sorrendje nem számít (szimmetriatörés), és
    a keresés egyre nagyobb lépésszámmal többször fut végig a kombinációkon, így a
    reménytelen kombinációk nem fogják le az időt.
    """
    if not pieces or len(pieces) > 40:
        return None
    pcs = sorted(pieces, key=lambda p: (-(p.width * p.height), -max(p.width, p.height),
                                        min(p.width, p.height)))
    # darabtípusok: méret + (alakos darabnál) a foglalt terület alakja
    types: Dict[tuple, int] = {}
    reps: List[Piece] = []
    dims = []
    for p in pcs:
        cv = piece_cover(p)
        sig = (p.width, p.height, tuple(cv) if cv else None)
        if sig not in types:
            types[sig] = len(reps)
            reps.append(Piece(p.width, p.height, "", shapes=p.shapes))
        dims.append((min(p.width, p.height), max(p.width, p.height), types[sig]))
    sm = smallest_dims(pcs)
    occ = [occ_area(p, sm) for p in pcs]
    need = sum(occ)
    cache: Dict[tuple, bool] = {}

    def fits(group_dims, sh):
        key = (sh.width, sh.height, tuple(sorted(group_dims)))
        r = cache.get(key)
        if r is None:
            tmp = [reps[t] for _, _, t in group_dims]
            r = _pack_single(tmp, sh, kerf, guillotine, 0) is not None
            cache[key] = r
        return r

    mixes = [m for m in _sheet_mixes(sheets, need, below_area)
             if all(any((w <= s.width and h <= s.height) or (h <= s.width and w <= s.height)
                        for s in m) for w, h, _ in dims)]
    # előbb a befoglaló téglalapok szerinti alsó korlát fölötti kombinációk (ezek
    # többnyire megoldhatók), csak utána a kisebbek, amelyekhez egymásba illesztés kell
    bbox_need = sum(p.width * p.height for p in pcs)
    found = None
    all_mixes = mixes
    for phase in (0, 1):
        mixes = [m for m in all_mixes
                 if (sum(s.area for s in m) >= bbox_need) == (phase == 0)
                 and (found is None or sum(s.area for s in m) < found[0])]
        for node_limit in (2000, 20000, 200000, 2000000):
            if not mixes or time.time() > deadline:
                break
            exhausted = True          # minden kombináció végig lett-e keresve
            for mix in mixes:
                if time.time() > deadline:
                    exhausted = False
                    break
                if found is not None and sum(s.area for s in mix) >= found[0]:
                    break
                bins = mix
                groups: List[List[int]] = [[] for _ in bins]
                gdims: List[list] = [[] for _ in bins]
                load = [0] * len(bins)
                where = [0] * len(pcs)
                nodes = [0]
                aborted = [False]

                def dfs(i):
                    if i == len(pcs):
                        return True
                    nodes[0] += 1
                    if nodes[0] > node_limit or (nodes[0] & 255 == 0 and time.time() > deadline):
                        aborted[0] = True
                        return False
                    a = occ[i]
                    # azonos méretű darabok: csak nem csökkenő táblasorszámmal (szimmetria)
                    start = where[i - 1] if i > 0 and dims[i] == dims[i - 1] else 0
                    tried_empty = set()
                    for b in range(start, len(bins)):
                        sh = bins[b]
                        if load[b] + a > sh.width * sh.height:
                            continue
                        if not groups[b]:
                            key = (sh.width, sh.height)
                            if key in tried_empty:
                                continue
                            tried_empty.add(key)
                        gdims[b].append(dims[i])
                        if fits(gdims[b], sh):
                            groups[b].append(i)
                            load[b] += a
                            where[i] = b
                            if dfs(i + 1):
                                return True
                            load[b] -= a
                            groups[b].pop()
                        gdims[b].pop()
                        if aborted[0]:
                            return False
                    return False

                if dfs(0):
                    placed, used = [], []
                    for b, sh in enumerate(bins):
                        if not groups[b]:
                            continue
                        idx = len(used)
                        pl = _pack_single([pcs[i] for i in groups[b]], sh, kerf, guillotine, idx)
                        if pl is None:
                            break
                        for pp in pl:
                            pp.sheet_index = idx
                        placed.extend(pl)
                        used.append(sh)
                    else:
                        _mark_residual(placed)
                        found = (sum(s.area for s in used), placed, used)
                        break       # növekvő területsorrend: ennél jobb ebben a körben nincs
                if aborted[0]:
                    exhausted = False
            if found is not None:
                # a kisebb, még be nem fejezett kombinációkat nagyobb lépésszámmal újra
                mixes = [m for m in mixes if sum(s.area for s in m) < found[0]]
            if exhausted and found is None:
                break
    return (found[1], found[2]) if found else None


def solution_cost(sheets: List[SheetSize], placed: List[PlacedPiece]) -> Tuple[int, int, int, float]:
    """
    Rangsor (kisebb = jobb):
      1. felhasznált táblák összterülete – azonos beépülő anyagnál ez a
         LEGJOBB KIHASZNÁLTSÁG (fő szempont, az ár nem számít bele)
      2. táblák darabszáma
      3. hegesztési varratok száma
      4. ár (csak döntetlennél, ha az anyagfájlban van ár)
    """
    return (sum(s.area for s in sheets), len(sheets), count_welds(placed),
            sum(s.cost() for s in sheets))


def _search(pieces, sheets, kerf, guillotine, fast: bool, weld_min: int = 0, max_parts: int = MAX_PARTS,
            deadline: Optional[float] = None):
    best = None
    sort_items = list(SORTS.values())[:2] if fast else list(SORTS.values())
    heurs = ("BSSF", "BAF") if fast else ("BSSF", "BAF", "BL")
    new_sheet_opts = [max(sheets, key=lambda s: s.area)] if fast else sheets
    for sk in sort_items:
        for h in heurs:
            for ns in new_sheet_opts:
                for wm in ([0, weld_min] if weld_min else [0]):
                    _check_cancel()
                    if best is not None and deadline is not None and time.time() > deadline:
                        return best
                    try:
                        placed, used = pack_pieces_into_sheets(pieces, sheets, kerf, guillotine,
                                                               sk, h, ns, wm, max_parts)
                    except ValueError:
                        if wm == 0 and weld_min:
                            continue
                        raise
                    c = solution_cost(used, placed)
                    if best is None or c < best[0]:
                        best = (c, placed, used)
    if best is None:
        raise ValueError("Nincs megoldás.")
    return best


@dataclass
class PartDef:
    """Egy szabandó lemezdarab: befoglaló téglalap + a benne lévő alkatrész(ek) alakja."""
    name: str
    width: int
    height: int
    shapes: List[Shape]
    alts: list = field(default_factory=list)   # alternatív felbontások (PartDef-listák)
    pair_params: Optional[tuple] = None        # ABO_A oldalpár: (felső, nettó mag., alsó, perem)
    no_weld: bool = False                      # nem toldható / nem darabolható (egyéb alkatrész)

    @property
    def is_rect(self) -> bool:
        return len(self.shapes) == 1 and abs(poly_area(self.shapes[0][1]) - self.width * self.height) < AREA_EPS


def rect_part(name: str, w: float, h: float) -> PartDef:
    w, h = int(round(w)), int(round(h))
    return PartDef(name, w, h, [(name, rect_poly(w, h))])


def abo_a_side_pair(top: float, height: float, bottom: float, flange: float = 0) -> PartDef:
    """
    ABO_A oldal: derékszögű trapéz NETTÓ magassággal ('height'); a ferde él az
    alsó él végétől a NETTÓ magasság vonaláig fut. Fölötte a teljes felső él
    ('top') hosszában 'flange' mm-es peremsáv (kiterítve, alapból 55 mm) –
    a jobb felső sarokban így egy rövid függőleges él van.
    A két oldalt egymással 180°-ban elforgatva egy közös ferde vágással szabjuk,
    vagy külön-külön – a program a jobbat választja.
    """
    if bottom >= top:
        raise ValueError("ABO_A: az alsó él legyen rövidebb a felső élnél.")
    f = flange
    W = top + bottom
    Htot = height + 2 * f
    # páros elrendezés: az 1. oldal f-fel feljebb tolva, peremmel felül;
    # a 2. oldal 180°-kal elforgatva, peremmel alul – közös ferde vágás
    s1p = [(0, f), (bottom, f), (top, height + f), (top, height + 2 * f), (0, height + 2 * f)]
    s2p = [(W - x, Htot - y) for x, y in s1p]
    # külön szabott oldal
    s1 = [(0, 0), (bottom, 0), (top, height), (top, height + f), (0, height + f)]
    bw, bh = int(math.ceil(top)), int(math.ceil(height + f))
    separate = [PartDef("Oldal 1", bw, bh, [("Oldal 1", list(s1))]),
                PartDef("Oldal 2", bw, bh, [("Oldal 2", list(s1))])]
    pd = PartDef("Oldalpár", int(math.ceil(W)), int(math.ceil(Htot)),
                 [("Oldal 1", s1p), ("Oldal 2", s2p)], alts=[separate])
    pd.pair_params = (top, height, bottom, f)
    return pd


def abo_a_pair_shifted(top: float, height: float, bottom: float, f: float, a: float) -> PartDef:
    """
    Oldalpár, amelyben a 2. (elforgatott) oldalt a közös ferde él MENTÉN 'a' mm-rel
    balra (és arányosan lefelé) csúsztatjuk: a pár keskenyebb lesz (felső+alsó−a),
    cserébe magasabb. Így pl. egy 1550 mm széles pár 1500-as táblára is befér.
    """
    W = top + bottom
    Htot = height + 2 * f
    b = a * height / (top - bottom)
    s1p = [(0, f), (bottom, f), (top, height + f), (top, height + 2 * f), (0, height + 2 * f)]
    s2p = [(W - x, Htot - y) for x, y in s1p]
    s1p = [(x, y + b) for x, y in s1p]
    s2p = [(x - a, y) for x, y in s2p]
    return PartDef("Oldalpár", int(math.ceil(W - a)), int(math.ceil(Htot + b)),
                   [("Oldal 1", s1p), ("Oldal 2", s2p)])


def abo_a_info(top: float, height: float, bottom: float, flange: float) -> Tuple[float, float, float]:
    """(teljes magasság, ferde él hossza a nettó magasságig, belső szög fokban)."""
    H = height + flange
    slant = math.hypot(top - bottom, height)
    angle = 180.0 - math.degrees(math.atan2(height, top - bottom))
    return H, slant, angle


@dataclass
class Product:
    name: str
    qty: int
    parts: List[PartDef]
    assembly_weld: float = 0.0     # összeállítási varrat hossza konténerenként (mm)
    info: str = ""                 # tájékoztató sor (pl. ABO_A ferde él, szög)
    merge: str = ""                # összevonás kulcsa (MERGE_OPTIONS)
    ctype: str = ""                # konténertípus (Sittes, ABO_A, ABO_B)
    thick: float = 2.0             # lemezvastagság (mm)
    part_thick: Dict[str, float] = field(default_factory=dict)   # eltérő vastagságú alkatrészek
    sub_of: str = ""               # vastagság szerinti résztermék: "" = fő rész

    def __post_init__(self):
        # régi forma támogatása: (szélesség, magasság, név) hármasok
        self.parts = [p if isinstance(p, PartDef) else rect_part(p[2], p[0], p[1])
                      for p in self.parts]


OTHER_TYPE = "Egyéb (DXF)"   # rajzolt, DXF-ből betöltött alkatrész
CONTAINER_TYPES = ["Sittes", "ABO_A", "ABO_B", OTHER_TYPE]


DEFAULT_FLANGE = (30, 25)    # ABO peremek: 30 x 25 mm, kiterítve 55 mm


MERGE_OPTIONS = {
    "": "nincs összevonás",
    "hatlap": "Alja + Hátlap egy lemezből hajlítva",
    "homlok": "Alja + Homloklap egy lemezből hajlítva",
    "mind": "Hátlap + Alja + Homloklap egy lemezből hajlítva (U)",
}
EDGE_TOL = 1.0     # mm – ekkora eltérésnél még közös élnek tekintjük


def _shared_edge(alja, width) -> Optional[Tuple[float, float]]:
    """Ha 'width' egyezik az alja valamelyik élével: (közös él, alja mélysége)."""
    aw, ah = alja
    if abs(width - aw) <= EDGE_TOL:
        return aw, ah
    if abs(width - ah) <= EDGE_TOL:
        return ah, aw
    return None


def _convex_clip(subject: Poly, clip: Poly) -> Poly:
    """Konvex sokszög vágása konvex sokszöggel (Sutherland–Hodgman), mindkettő CCW."""
    def inside(p, a, b):
        return (b[0] - a[0]) * (p[1] - a[1]) - (b[1] - a[1]) * (p[0] - a[0]) >= -1e-9

    def inter(p, q, a, b):
        x1, y1, x2, y2 = p[0], p[1], q[0], q[1]
        x3, y3, x4, y4 = a[0], a[1], b[0], b[1]
        den = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4)
        if abs(den) < 1e-12:
            return q
        t = ((x1 - x3) * (y3 - y4) - (y1 - y3) * (x3 - x4)) / den
        return (x1 + t * (x2 - x1), y1 + t * (y2 - y1))

    out = list(subject)
    for i in range(len(clip)):
        a, b = clip[i], clip[(i + 1) % len(clip)]
        inp, out = out, []
        if not inp:
            break
        for j in range(len(inp)):
            cur, prev = inp[j], inp[j - 1]
            if inside(cur, a, b):
                if not inside(prev, a, b):
                    out.append(inter(prev, cur, a, b))
                out.append(cur)
            elif inside(prev, a, b):
                out.append(inter(prev, cur, a, b))
    return out


def _overlap(p1: Poly, p2: Poly) -> float:
    c = _convex_clip(p1, p2)
    return poly_area(c) if len(c) >= 3 else 0.0


def abo_b_side_poly(A: float, H: float, C: float, ty: float, E: Optional[float] = None) -> Poly:
    """ABO_B oldal kontúrja (CCW): alul E hosszú él, orr csúcsa (A, ty), felül C hosszú él, magasság H."""
    E = E or C
    return [(0, 0), (E, 0), (A, ty), (C, H), (0, H)]


def abo_b_side_parts(A: float, B: float, C: float, f: float,
                     D: Optional[float] = None, E: Optional[float] = None) -> List[PartDef]:
    """
    ABO_B oldalpár valós (orros) alakkal. A = teljes hossz az orr csúcsáig,
    C = felső egyenes él, E = alsó egyenes él (üresen = C), B = nettó magasság,
    f = perem (a teljes magassághoz adódik), D = az orr csúcsának magassága
    alulról (alapból középen).
    Visszaad: [egyesével szabott oldalak] és – ha anyagot spórol – egy 180°-ban
    egymásba fordított pár alternatívát.
    """
    if C >= A:
        raise ValueError("ABO_B: a felső egyenes él (C) legyen rövidebb a teljes hossznál (A).")
    E = E or C
    if E >= A:
        raise ValueError("ABO_B: az alsó egyenes él legyen rövidebb a teljes hossznál (A).")
    H = B + f
    ty = D if D else H / 2
    if not (0 < ty < H):
        raise ValueError("ABO_B: az orr csúcsának magassága (D) a lap magasságán belül legyen.")
    p = abo_b_side_poly(A, H, C, ty, E)
    single = [PartDef("Oldal 1", int(math.ceil(A)), int(math.ceil(H)), [("Oldal 1", p)]),
              PartDef("Oldal 2", int(math.ceil(A)), int(math.ceil(H)), [("Oldal 2", list(p))])]
    # egymásba fordított pár: a 2. oldal 180°-kal elforgatva, függőlegesen eltolva,
    # vízszintesen a lehető legközelebb csúsztatva – a legkisebb befoglaló téglalap nyer
    best = None
    rot = [(A - x, H - y) for x, y in p]
    for k in range(0, 21):
        dy = H * k / 20.0
        lo, hi = 0.0, A          # bináris keresés a legkisebb átfedésmentes eltolásra
        q_hi = [(x + hi, y + dy) for x, y in rot]
        if _overlap(p, q_hi) > 1.0:
            continue
        for _ in range(30):
            mid = (lo + hi) / 2
            q = [(x + mid, y + dy) for x, y in rot]
            if _overlap(p, q) > 1.0:
                lo = mid
            else:
                hi = mid
        dx = hi
        W = max(A, dx + A)
        Ht = max(H, dy + H)
        area = W * Ht
        if best is None or area < best[0]:
            best = (area, dx, dy, W, Ht)
    alts = []
    if best and best[0] < 0.97 * 2 * A * H:
        _, dx, dy, W, Ht = best
        q = [(x + dx, y + dy) for x, y in rot]
        pair = PartDef("Oldalpár", int(math.ceil(W - 1e-3)), int(math.ceil(Ht - 1e-3)),
                       [("Oldal 1", p), ("Oldal 2", q)])
        alts.append(pair)
    return single, alts


def _panel_poly(w_bottom: float, w_top: float, h: float, W: float) -> Poly:
    """Szimmetrikus trapéz (vagy téglalap) lap: alsó él w_bottom, felső w_top,
    magasság h, vízszintesen egy W széles sáv közepére igazítva, alul y=0."""
    return [((W - w_bottom) / 2, 0), ((W + w_bottom) / 2, 0),
            ((W + w_top) / 2, h), ((W - w_top) / 2, h)]


def trapezoid_part(name: str, w_bottom: float, w_top: float, h: float) -> PartDef:
    """Sittes homlok/hátlap: alul keskenyebb, felül szélesebb (vagy fordítva) trapéz."""
    if abs(w_top - w_bottom) < 0.5:
        return rect_part(name, w_bottom, h)
    W = max(w_bottom, w_top)
    return PartDef(name, int(math.ceil(W)), int(math.ceil(h)),
                   [(name, _panel_poly(w_bottom, w_top, h, W))])


def trapezoid_pair(n1: str, n2: str, w_bottom: float, w_top: float, h: float) -> PartDef:
    """Két egyforma trapéz egymáshoz képest 180°-kal elforgatva, közös ferde vágással."""
    Wm = max(w_bottom, w_top)
    t1 = _panel_poly(w_bottom, w_top, h, Wm)
    sh = (w_bottom + w_top) / 2
    t2 = [(Wm - x + sh, h - y) for x, y in t1]
    W = Wm + sh
    pd = PartDef(f"{n1}+{n2} pár", int(math.ceil(W)), int(math.ceil(h)),
                 [(n1, t1), (n2, t2)],
                 alts=[[trapezoid_part(n1, w_bottom, w_top, h),
                        trapezoid_part(n2, w_bottom, w_top, h)]])
    return pd


def merged_bottom_part(merge: str, alja, kw, k_sheet_h, fw, f_sheet_h,
                       k_top: Optional[float] = None,
                       f_top: Optional[float] = None) -> Optional[Tuple[PartDef, float, int]]:
    """
    Összevont (egy lemezből hajlított) aljas darab. Visszaad: (PartDef, közös él
    hossza, hajlítások száma) vagy None, ha az élek nem egyeznek.
    Elrendezés alulról felfelé: Homloklap | Alja | Hátlap. A hátlap és a homloklap
    lehet trapéz is (sittes): az alján lévő közös él egyezik az alja élével,
    a túlsó éle k_top / f_top széles.
    """
    k_top = kw if k_top is None else k_top
    f_top = fw if f_top is None else f_top
    use_k = merge in ("hatlap", "mind")
    use_f = merge in ("homlok", "mind")
    e1 = _shared_edge(alja, kw) if use_k else None
    e2 = _shared_edge(alja, fw) if use_f else None
    if (use_k and not e1) or (use_f and not e2):
        return None
    if e1 and e2 and abs(e1[0] - e2[0]) > EDGE_TOL:
        return None
    E, D = e1 or e2
    W = max([E] + ([k_top] if use_k else []) + ([f_top] if use_f else []))
    shapes = []
    y = 0.0
    if use_f:   # homloklap alul, fejjel lefelé (a közös él felül), valódi trapéz alakkal
        poly = _panel_poly(E, f_top, f_sheet_h, W)
        shapes.append(("Homloklap", [(x, f_sheet_h - yy) for x, yy in poly]))
        y += f_sheet_h
    shapes.append(("Alja", [(x + (W - E) / 2, yy + y) for x, yy in rect_poly(E, D)]))
    y += D
    if use_k:
        poly = _panel_poly(E, k_top, k_sheet_h, W)
        shapes.append(("Hátlap", [(x, yy + y) for x, yy in poly]))
        y += k_sheet_h
    name = {"hatlap": "Alja+Hátlap", "homlok": "Alja+Homloklap",
            "mind": "Hátlap+Alja+Homloklap"}[merge]
    return PartDef(name, int(math.ceil(W)), int(math.ceil(y)), shapes), E, (2 if merge == "mind" else 1)


def make_product_variants(name: str, ctype: str, qty: int, alja, hatlap, homlok, oldal,
                          flange=DEFAULT_FLANGE, weld_sides: int = 1,
                          merges=("", "hatlap", "homlok", "mind"),
                          extras: Optional[Dict[str, float]] = None) -> List["Product"]:
    """A termék összes megvalósítható összevonási változata (az élhossz-egyezés alapján)."""
    out = []
    for m in merges:
        try:
            out.append(make_product(name, ctype, qty, alja, hatlap, homlok, oldal,
                                    flange, weld_sides, merge=m, extras=extras))
        except ValueError:
            if m == "":
                raise
    return out


def make_product(name: str, ctype: str, qty: int, alja, hatlap, homlok, oldal,
                 flange=DEFAULT_FLANGE, weld_sides: int = 1, merge: str = "",
                 extras: Optional[Dict[str, float]] = None) -> Product:
    """oldal: Sittes és ABO_B → (szélesség, magasság); ABO_A → (felső hossz, magasság, alsó hossz).
    ABO típusoknál a hátlap, a homloklap és az oldalak FELSŐ élére ráhajlított
    perem kerül: a lemez ott a perem kiterített szélességével (30+25 mm) nagyobb.
    A magasság-mezőkbe a kész (hajlított) méretet kell írni."""
    f = sum(flange) if ctype in ("ABO_A", "ABO_B") else 0
    if ctype == "ABO_A":
        # ABO_A homloklap: a ferde oldalélre kerül, nettó magassága = a ferde él
        # hossza (a nettó magasságig), ehhez jön a saját peremje (+55 mm).
        # A felhasználó csak a szélességét adja meg.
        _, slant, _ = abo_a_info(oldal[0], oldal[1], oldal[2], f)
        homlok_sheet_h = slant + f
        homlok = (homlok[0], slant)
    else:
        homlok_sheet_h = homlok[1] + f
    ex = extras or {}
    # Sittes: a homlok- és hátlap trapéz lehet: a táblázat szélessége a felső él (C),
    # a rajzon megadott D az alsó él (az aljához csatlakozik); ha nincs D, téglalap
    k_top, f_top = hatlap[0], homlok[0]
    if ctype == "Sittes":
        if ex.get("hat_D"):
            hatlap = (float(ex["hat_D"]), hatlap[1])
        if ex.get("hom_D"):
            homlok = (float(ex["hom_D"]), homlok[1])
    saved_weld = 0.0
    if merge:
        mp = merged_bottom_part(merge, alja, hatlap[0], hatlap[1] + f, homlok[0], homlok_sheet_h,
                                k_top, f_top)
        if mp is None:
            raise ValueError(f"{MERGE_OPTIONS[merge]}: a közös élek hossza nem egyezik "
                             f"(alja {alja[0]:g}x{alja[1]:g}, hátlap szél. {hatlap[0]:g}, "
                             f"homloklap szél. {homlok[0]:g}).")
        part, edge, bends = mp
        saved_weld = edge * bends
        parts = [part]
        if merge == "hatlap":
            parts.append(trapezoid_part("Homloklap", homlok[0], f_top, homlok_sheet_h))
        elif merge == "homlok":
            parts.append(trapezoid_part("Hátlap", hatlap[0], k_top, hatlap[1] + f))
    else:
        # trapéz homlok-/hátlap valódi alakkal (a ferde élek melletti sarkok
        # felhasználhatók); két egyforma trapéz egymásba fordítva is szabható
        parts = [rect_part("Alja", *alja)]
        same = (abs(hatlap[0] - homlok[0]) < 0.5 and abs(k_top - f_top) < 0.5
                and abs(hatlap[1] + f - homlok_sheet_h) < 0.5)
        if same and abs(hatlap[0] - k_top) >= 0.5:
            # egymásba fordított pár; alternatívája a külön szabás
            parts.append(trapezoid_pair("Hátlap", "Homloklap", hatlap[0], k_top, homlok_sheet_h))
        else:
            parts += [trapezoid_part("Hátlap", hatlap[0], k_top, hatlap[1] + f),
                      trapezoid_part("Homloklap", homlok[0], f_top, homlok_sheet_h)]
    if ctype == "Sittes":
        parts += [rect_part("Oldal 1", oldal[0], oldal[1]), rect_part("Oldal 2", oldal[0], oldal[1])]
    elif ctype == "ABO_A":
        parts.append(abo_a_side_pair(*oldal[:3], flange=f))
    elif ctype == "ABO_B":
        # orros oldallap (A teljes hossz, C egyenes él, B nettó magasság + perem)
        if len(oldal) >= 3 and oldal[2]:
            D = float(ex.get("side_D") or 0) or None
            E = float(ex.get("side_E") or 0) or None
            single, alts = abo_b_side_parts(oldal[0], oldal[1], oldal[2], f, D, E)
            if alts:
                # a pár az alapeset, az egyenként szabott oldalak az alternatíva
                pair = alts[0]
                pair.alts = [single]
                parts.append(pair)
            else:
                parts += single
        else:
            parts += [rect_part("Oldal 1", oldal[0], oldal[1] + f),
                      rect_part("Oldal 2", oldal[0], oldal[1] + f)]
    else:
        raise ValueError(f"Ismeretlen konténertípus: {ctype}")
    info = ""
    if ctype == "ABO_A":
        Ht, sl, ang = abo_a_info(oldal[0], oldal[1], oldal[2], f)
        info = (f"ABO_A oldal: teljes magasság {Ht:.1f} mm (nettó {oldal[1]:g} + {f:g} perem), "
                f"ferde él {sl:.1f} mm, szög {ang:.1f}° → homloklap {homlok[0]:g} x {sl:.1f} mm "
                f"(+{f:g} perem = {sl + f:.1f} mm)")
    if ctype == "ABO_B" and len(oldal) >= 3 and oldal[2]:
        H = oldal[1] + f
        ty = float(ex.get("side_D") or 0) or H / 2
        E = float(ex.get("side_E") or 0) or oldal[2]
        up = math.hypot(oldal[0] - oldal[2], H - ty)
        lo = math.hypot(oldal[0] - E, ty)
        info = (f"ABO_B oldal: A {oldal[0]:g}, felső él C {oldal[2]:g}, alsó él {E:g}, teljes magasság {H:g} "
                f"(nettó {oldal[1]:g} + {f:g} perem), orr csúcsa {ty:g} mm-en; "
                f"ferde élek: felső {up:.1f}, alsó {lo:.1f} mm")
    if ctype == "Sittes" and (ex.get("hom_A") or ex.get("hat_A") or ex.get("hom_D") or ex.get("hat_D")):
        def btxt(pre, wb, wt, h):
            a = float(ex.get(pre + "_A") or 0)
            w = f"{wt:g} (alul {wb:g})" if abs(wb - wt) >= 0.5 else f"{wt:g}"
            return (f"{w} x {h:g} (hajlítás {a:g} mm-rel az alja felett)" if a
                    else f"{w} x {h:g}")
        info = (f"Sittes szabás: homloklap {btxt('hom', homlok[0], f_top, homlok[1])}; "
                f"hátlap {btxt('hat', hatlap[0], k_top, hatlap[1])}")
    if merge:
        info = (info + " | " if info else "") + \
            f"{MERGE_OPTIONS[merge]} (−{saved_weld / 1000:.2f} m varrat/db)"
    prod = Product(name, qty, parts,
                   (assembly_weld_length(ctype, alja, hatlap, homlok, oldal, k_top, f_top,
                                         side_E=float(ex.get("side_E") or 0) or None,
                                         side_D=float(ex.get("side_D") or 0) or None)
                    - saved_weld) * weld_sides,
                   info)
    prod.merge = merge
    prod.ctype = ctype
    return prod


def assembly_weld_length(ctype: str, alja, hatlap, homlok, oldal,
                         k_top: Optional[float] = None, f_top: Optional[float] = None,
                         side_E: Optional[float] = None, side_D: Optional[float] = None) -> float:
    """
    Összeállítási varrat egy konténerre (mm), egyoldali, végig hegesztve.
    Nyitott tetejű dobozként számol: minden alkatrész-él, amely nem a felső
    peremen van, két alkatrész közös éle → varrat = (kerületek összege −
    szabad felső élek) / 2. A kész (perem nélküli) méretekből számol.
    """
    k_top = hatlap[0] if k_top is None else k_top
    f_top = homlok[0] if f_top is None else f_top

    def trap_per(wb, wt, h):
        return wb + wt + 2 * math.hypot((wt - wb) / 2, h)

    per = (2 * (alja[0] + alja[1]) + trap_per(hatlap[0], k_top, hatlap[1])
           + trap_per(homlok[0], f_top, homlok[1]))
    if ctype == "ABO_A":
        top, h, bottom = oldal[:3]
        # a ferde él hegesztett része a hajlítási vonalig (nettó magasságig) tart
        slant_net = math.hypot(top - bottom, h) if h > 0 else 0.0
        side_per = top + bottom + h + slant_net
        side_top = top
    elif ctype == "ABO_B" and len(oldal) >= 3 and oldal[2]:
        A_, h, C_ = oldal[:3]
        E_ = side_E or C_
        ty = min(side_D, h) if side_D else h / 2
        side_per = C_ + E_ + h + math.hypot(A_ - C_, h - ty) + math.hypot(A_ - E_, ty)
        side_top = C_ + math.hypot(A_ - C_, h - ty)      # felső él + felső ferde él (perem)
    else:
        side_per = 2 * (oldal[0] + oldal[1])
        side_top = oldal[0]
    per += 2 * side_per
    free = k_top + f_top + 2 * side_top
    return max(0.0, (per - free) / 2)


def optimize_products(products: List[Product], sheets: List[SheetSize], kerf: int = 0,
                      guillotine: bool = False, weld: bool = True,
                      min_part_pct: float = 40.0, max_parts: int = MAX_PARTS,
                      time_limit: float = 15.0, deadline: Optional[float] = None,
                      seed: int = 0):
    """
    time_limit: a TELJES számítás időkerete másodpercben (deadline: ugyanez
    abszolút időpontként). A keret a körök között oszlik meg; minden kör
    legalább egy megoldást ad, utána a határidőnél abbahagyja.
    seed: a véletlen keresés kiindulópontja (párhuzamos futásnál eltérő).
    weld: maradékok felhasználása hegesztett toldással.
    min_part_pct: egy kettévágott alkatrész egyik része sem lehet kisebb, mint a
    termék legkisebb területű alkatrészének ennyi %-a (a kötelező, túlméret miatti
    vágásokra is érvényes).

    Több termék (különböző méretekkel és darabszámmal) közös szabása:
    minden termék minden darabja ugyanabból a lemezkészletből, közös
    maradékkezeléssel készül. Az azonos méretű alkatrészek ugyanúgy darabolódnak.
    """
    products = [p for p in products if p.qty > 0]
    if not products:
        raise ValueError("Nincs gyártandó termék.")
    t0 = time.time()
    if deadline is None:
        deadline = t0 + max(0.0, time_limit)
    budget = max(0.0, deadline - t0)
    n_pieces = sum(len(pr_.parts) * pr_.qty for pr_ in products)
    # kevés darabnál a pontos újraosztás (4. kör) kapja az idő nagyobb részét
    d1, d2 = t0 + 0.2 * budget, t0 + 0.35 * budget
    d3 = t0 + (0.55 if n_pieces <= 40 else 0.8) * budget

    weld_min = 1 if weld else 0     # belső kapcsoló a toldásos kitöltéshez
    min_area_of = {}
    for prod in products:
        # a legkisebb ALKATRÉSZ valós területe (alakos oldalaknál a sokszög területe)
        smallest = min(poly_area(poly) for pd in prod.parts for lab, poly in pd.shapes
                       if not is_hole(lab))
        min_area_of[prod.name] = int(smallest * min_part_pct / 100)

    # egyedi darabtípusok (+ min. részterület) és darabolási változataik
    def sig(pd: PartDef, ma: int):
        return (pd.width, pd.height, ma, pd.is_rect,
                tuple((lab, tuple(poly)) for lab, poly in pd.shapes))

    dims: List[tuple] = []
    for prod in products:
        for pd in prod.parts:
            key = sig(pd, min_area_of[prod.name])
            if key not in dims:
                dims.append(key)
    # változat = [(PartDef, eredeti név a varrathoz, rész sorszáma), ...]
    pd_of = {}
    for prod in products:
        for pd in prod.parts:
            pd_of.setdefault(sig(pd, min_area_of[prod.name]), pd)
    per_dim = []
    for key in dims:
        w, h, ma, is_rect, _ = key
        pd = pd_of[key]
        fits = fits_any(w, h, sheets)
        variants = []
        if pd.no_weld and not fits:
            raise ValueError(f"A(z) {pd.name} ({w}x{h} mm) egyik táblára sem fér rá "
                             f"(egyéb alkatrész, nem toldható).")
        if fits or weld_min:
            variants.append([(pd, "", 0, 0.0)])     # egyben (vagy toldással kitöltve)
        if is_rect and not pd.no_weld:
            for cut in generate_cut_variants(w, h, sheets, ma):
                if len(cut) > 1:
                    same_h = all(b == cut[0][1] for _, b in cut)
                    seam_len = cut[0][1] if same_h else cut[0][0]
                    variants.append([(rect_part(f"{pd.name}_{chr(65 + i)}", a, b), pd.name, i + 1,
                                      seam_len if i < len(cut) - 1 else 0.0)
                                     for i, (a, b) in enumerate(cut)])
        for alt in pd.alts:
            variants.append([(sub, "", 0, 0.0) for sub in alt])
        if pd.pair_params:
            # eltolt (keskenyebb) párok: pontosan a táblaméretekhez igazítva
            top_, h_, bot_, f_ = pd.pair_params
            dims_ = sorted({d for sh in sheets for d in (sh.width, sh.height)})
            for d in dims_:
                a = (top_ + bot_) - d
                if 0 < a <= bot_:
                    sp = abo_a_pair_shifted(top_, h_, bot_, f_, a)
                    if fits_any(sp.width, sp.height, sheets):
                        variants.append([(sp, "", 0, 0.0)])
        if not variants:
            raise ValueError(f"A(z) {pd.name} ({w}x{h} mm) toldás nélkül nem készíthető el a "
                             f"táblákból úgy, hogy egyik rész se legyen kisebb {ma / 1e6:.2f} m²-nél "
                             f"({min_part_pct:g}% szabály). Kapcsold be a hegesztett toldást.")
        per_dim.append(variants)

    # kombinációk: alapeset, egyenkénti eltérések, majd véletlen minták (max. 64)
    combos = [tuple(0 for _ in per_dim)]
    for i, vs in enumerate(per_dim):
        for j in range(1, len(vs)):
            c = [0] * len(per_dim)
            c[i] = j
            combos.append(tuple(c))
    total = 1
    for vs in per_dim:
        total *= len(vs)
    rnd = random.Random(42 + seed)
    tries = 0
    while len(combos) < min(64, total) and tries < 1000:
        tries += 1
        c = tuple(rnd.randrange(len(vs)) for vs in per_dim)
        if c not in combos:
            combos.append(c)
    combos = combos[:64]

    multi = len(products) > 1
    labels = []
    for pi, prod in enumerate(products):
        for k in range(1, prod.qty + 1):
            labels.append((f"{prod.name} {k}." if multi else f"K{k}", prod))

    def build(combo):
        pieces = []
        for unit, (tag, prod) in enumerate(labels, start=1):
            ma = min_area_of[prod.name]
            for pd in prod.parts:
                di = dims.index(sig(pd, ma))
                for sub, base, idx, seam in per_dim[di][combo[di]]:
                    pieces.append(Piece(sub.width, sub.height, sub.name, unit, base, idx,
                                        tag=tag, min_area=ma,
                                        shapes=[] if sub.is_rect else list(sub.shapes),
                                        seam=seam, no_weld=sub.no_weld))
        return pieces

    # alakos (egymásba illesztő) mód csak akkor, ha van nem téglalap darab
    has_shapes = any(not sub.is_rect for vs in per_dim for v in vs for sub, *_ in v)
    modes = (False, True) if has_shapes else (False,)

    # 1. kör: gyors értékelés minden darabolási kombinációra (előbb a gyors
    # befoglaló-téglalapos, majd az alakos elhelyezéssel)
    scored = []
    for nest in modes:
        set_nesting(nest)
        for combo in combos:
            if scored and time.time() > (d1 if nest else t0 + 0.6 * (d1 - t0)) :
                break
            try:
                r = _search(build(combo), sheets, kerf, guillotine, True, weld_min, max_parts,
                            deadline=d1)
            except ValueError:
                continue
            scored.append((r[0], combo, r))
    scored.sort(key=lambda t: t[0])

    # 2. kör: alapos keresés a legígéretesebb kombinációkon – a legjobb 3
    # kihasználtság ÉS a legjobb 3 táblaszám szerint (a kettő más-más
    # irányba viszi a keresést, a végső döntés a kihasználtság szerint megy)
    if not scored:
        raise ValueError("Nem találtam olyan szabást, amivel minden darab elkészíthető.")
    top = []
    for c in [t[1] for t in scored[:3]] + \
             [t[1] for t in sorted(scored, key=lambda t: (t[0][1], t[0][0]))[:3]]:
        if c not in top:
            top.append(c)
    best = scored[0][2]
    for nest in modes:
        set_nesting(nest)
        for combo in top:
            if time.time() > (d2 if nest or len(modes) == 1 else d1 + 0.5 * (d2 - d1)):
                break
            try:
                r = _search(build(combo), sheets, kerf, guillotine, False, weld_min, max_parts,
                            deadline=d2)
            except ValueError:
                continue
            if r[0] < best[0]:
                best = r
    if best is None:
        raise ValueError("Nem találtam olyan szabást, amivel minden darab elkészíthető.")

    # 3. kör: véletlen újraindítások a megadott időkereten belül
    rnd2 = random.Random(7 + 1000 * seed)
    biggest = max(sheets, key=lambda s: s.area)
    while time.time() < d3:
        _check_cancel()
        set_nesting(rnd2.random() < 0.5 if has_shapes else False)
        pieces = build(rnd2.choice(top))
        noise = {id(p): 0.7 + 0.6 * rnd2.random() for p in pieces}
        sk = lambda p, noise=noise: p.width * p.height * noise[id(p)]
        heur = rnd2.choice(("BSSF", "BAF", "BL"))
        ns = biggest if rnd2.random() < 0.6 else rnd2.choice(sheets)
        try:
            placed, used = pack_pieces_into_sheets(pieces, sheets, kerf, guillotine,
                                                   sk, heur, ns, weld_min, max_parts)
        except ValueError:
            continue
        c = solution_cost(used, placed)
        if c < best[0]:
            best = (c, placed, used)

    # 4. kör: pontos újraosztás kevés darabnál – vegyes táblaméretekkel is
    set_nesting(has_shapes)
    ex_deadline = deadline
    for combo in top:
        if time.time() > ex_deadline:
            break
        r = exact_repack(build(combo), sheets, kerf, guillotine, best[0][0],
                         min(ex_deadline, time.time() + (ex_deadline - time.time()) / 2)
                         if combo is not top[-1] else ex_deadline)
        if r is not None:
            c = solution_cost(r[1], r[0])
            if c < best[0]:
                best = (c, r[0], r[1])

    set_nesting(True)
    _, placed, used = best
    allowed = {(s.width, s.height) for s in sheets}
    for s in used:
        if (s.width, s.height) not in allowed:
            raise RuntimeError(f"Belső hiba: nem engedélyezett táblaméret ({s.width}x{s.height}).")
    return placed, used


def variant_key(products, placed, used):
    """Változatok rangsora: 1. táblák összterülete (= kihasználtság), 2. táblaszám,
    3. hegesztés hossza, 4. ár."""
    return (sum(s.area for s in used), len(used), round(total_weld(products, placed)),
            sum(s.cost() for s in used))


def total_weld(products: List[Product], placed: List[PlacedPiece]) -> float:
    return sum(p.assembly_weld * p.qty for p in products) + sum(x.seam for x in placed)


def _opt_worker(args):
    """Párhuzamos folyamatban futó egység: egy változat-kombináció optimalizálása."""
    combo, sheets, kerf, guillotine, weld, min_pct, deadline, seed = args
    try:
        placed, used = optimize_products(list(combo), sheets, kerf, guillotine, weld, min_pct,
                                         deadline=deadline, seed=seed)
    except ValueError:
        return None
    return placed, used


def worker_count() -> int:
    return max(1, min(8, (os.cpu_count() or 2) - 1))


def _run_tasks(pool, tasks, progress=None, label=""):
    """Feladatok futtatása – párhuzamosan, ha van folyamatkészlet. A megszakítást
    (Cancelled) a folyamatkészlet leállításával adja tovább."""
    if pool is None:
        out = []
        for i, t in enumerate(tasks):
            _check_cancel()
            if progress and label:
                progress(f"{label}: {i + 1}/{len(tasks)}")
            out.append(_opt_worker(t))
        return out
    pending = [pool.apply_async(_opt_worker, (t,)) for t in tasks]
    while True:
        if _CANCEL.is_set():
            pool.terminate()
            raise Cancelled("A számítás megszakítva.")
        done = sum(1 for r in pending if r.ready())
        if progress and label:
            progress(f"{label}: {done}/{len(tasks)}")
        if done == len(pending):
            return [r.get() for r in pending]
        pending[[r.ready() for r in pending].index(False)].wait(0.1)


def optimize_with_variants(variant_lists: List[List[Product]], sheets: List[SheetSize],
                           kerf: int = 0, guillotine: bool = False, weld: bool = True,
                           min_part_pct: float = 40.0, time_limit: float = 15.0,
                           progress=None, pool=None, workers: int = 1):
    """
    Termékenként több változat (pl. összevont alja+hátlap) közül a legjobb
    kombináció: 1. táblák összterülete, 2. táblaszám, 3. hegesztés (m), 4. költség.
    A teljes futás a time_limit-en belül marad. pool: multiprocessing.Pool –
    ha meg van adva, a változatokat és a véletlen kereséseket párhuzamosan futtatja.
    Visszaad: (placed, used, products, összehasonlító táblázat).
    """
    t0 = time.time()
    end = t0 + max(0.5, time_limit)
    workers = max(1, workers if pool is not None else 1)
    combos = list(itertools.product(*variant_lists))
    args = lambda combo, dl, seed: (combo, sheets, kerf, guillotine, weld, min_part_pct, dl, seed)

    if len(combos) == 1:
        # egy változat: párhuzamosan több, eltérő kiindulású keresés, a legjobb nyer
        res = _run_tasks(pool, [args(combos[0], end, sd) for sd in range(workers)],
                         progress, "Keresés")
        res = [r for r in res if r is not None]
        if not res:
            raise ValueError("Nem találtam olyan szabást, amivel minden darab elkészíthető.")
        placed, used = min(res, key=lambda r: variant_key(list(combos[0]), r[0], r[1]))
        return placed, used, list(combos[0]), []

    # 1. kör: gyors értékelés – a keret kb. 40%-a. Kevés kombinációnál mindet
    # kipróbálja; sok terméknél termékenként, a többit a legjobbon tartva
    # (így a próbák száma nem nő hatványozottan a termékek számával).
    p1_end = t0 + 0.4 * (end - t0)
    rows = []
    seen = {}

    def evaluate(batch, label):
        batch = [c for c in batch if tuple(map(id, c)) not in seen]
        if not batch:
            return
        rounds_ = math.ceil(len(batch) / workers)
        per = max(0.05, min(2.0, (p1_end - time.time()) / rounds_ / max(1, left_batches[0])))
        res = _run_tasks(pool, [args(c, time.time() + per * (1 + i // workers), 0)
                                for i, c in enumerate(batch)], progress, label)
        for combo, r in zip(batch, res):
            seen[tuple(map(id, combo))] = r
            if r is not None:
                rows.append((variant_key(list(combo), r[0], r[1]), combo, r[0], r[1]))

    if len(combos) <= max(16, 4 * workers):
        left_batches = [1]
        evaluate(combos, "Változatok összevetése")
    else:
        n = len(variant_lists)
        left_batches = [n + 1]
        cur = tuple(vl[0] for vl in variant_lists)
        evaluate([cur], "Változatok összevetése")
        for j in range(n):
            left_batches[0] = n - j
            evaluate([cur[:j] + (v,) + cur[j + 1:] for v in variant_lists[j]],
                     f"Változatok összevetése ({j + 1}/{n}. termék)")
            if rows:
                cur = min(rows, key=lambda r: r[0])[1]
    if not rows:
        raise ValueError("Egyik összevonási változattal sem készíthető el a termék.")
    rows.sort(key=lambda r: r[0])
    best_area = rows[0][0][0]
    # 2. kör: a maradék idő a legjobb változatokra (legfeljebb 1 táblával rosszabbak);
    # több magon minden jelölt több, eltérő kiindulású keresést kap
    biggest = max(s.area for s in sheets)
    cand = [r for r in rows if r[0][0] <= best_area + biggest][:max(3, min(workers, 6))]
    final = {id(r[1]): r for r in rows}
    if end - time.time() > 0.3:
        if pool is None:
            tasks = []
            left = end - time.time()
            for j, (_, combo, _, _) in enumerate(cand):
                tasks.append(args(combo, None, 1))
            res = []
            for j, t in enumerate(tasks):
                _check_cancel()
                if progress:
                    progress(f"A legjobb változatok finomítása: {j + 1}/{len(tasks)}")
                dl = time.time() + (end - time.time()) / (len(tasks) - j)
                res.append(_opt_worker(t[:6] + (dl, 1)))
        else:
            per = max(1, workers // len(cand))
            tasks = [args(combo, end, 1 + sd) for _, combo, _, _ in cand for sd in range(per)]
            res = _run_tasks(pool, tasks, progress, "A legjobb változatok finomítása")
        for t, r in zip(tasks, res):
            if r is None:
                continue
            combo = t[0]
            k2 = variant_key(list(combo), r[0], r[1])
            if k2 < final[id(combo)][0]:
                final[id(combo)] = (k2, combo, r[0], r[1])
    rows = sorted(final.values(), key=lambda r: r[0])
    best = rows[0]
    table = [(" / ".join(f"{p.name}: {MERGE_OPTIONS[p.merge]}" for p in combo), key[1], key[2],
              100.0 * sum(x.net_area for x in pl) / key[0] if key[0] else 0.0)
             for key, combo, pl, _ in rows]
    return best[2], best[3], list(best[1]), table


def optimize_by_thickness(variant_lists: List[List[Product]], sheets: List[SheetSize],
                          kerf: int = 0, guillotine: bool = False, weld: bool = True,
                          min_part_pct: float = 40.0, time_limit: float = 15.0,
                          progress=None, pool=None, workers: int = 1):
    """
    Különböző vastagságú termékek nem kerülhetnek közös táblára: vastagságonként
    külön szabás, a keresési idő a darabszám arányában oszlik meg. A felhasznált
    táblákra rákerül a vastagság (súlyszámításhoz), a lemezek sorszáma folytatólagos.
    """
    groups: Dict[float, List[List[Product]]] = {}
    for vl in variant_lists:
        pt = {n: t for n, t in vl[0].part_thick.items() if abs(t - vl[0].thick) > 1e-9}
        if not pt:
            groups.setdefault(float(vl[0].thick), []).append(vl)
            continue
        # eltérő vastagságú alkatrész (pl. Sittes alja): külön szabás a saját vastagságában;
        # a fő részhez tartozik az összeállítási varrat
        by_t: Dict[float, List[Product]] = {}
        for v in vl:
            for t in {float(v.thick)} | {float(x) for x in pt.values()}:
                parts = [p for p in v.parts if float(pt.get(p.name, v.thick)) == t]
                if not parts:
                    continue
                main = t == float(v.thick)
                by_t.setdefault(t, []).append(dataclasses.replace(
                    v, parts=parts, thick=t, part_thick={},
                    assembly_weld=v.assembly_weld if main else 0.0,
                    sub_of="" if main else ", ".join(p.name for p in parts)))
        for t, sub in by_t.items():
            groups.setdefault(t, []).append(sub)
    if len(groups) == 1:
        (t, vls), = groups.items()
        placed, used, prods, table = optimize_with_variants(
            vls, sheets, kerf, guillotine, weld, min_part_pct, time_limit=time_limit,
            progress=progress, pool=pool, workers=workers)
        return placed, [dataclasses.replace(u, thick=t) for u in used], prods, table
    t0 = time.time()
    end = t0 + max(0.5, time_limit)
    size = {t: sum(len(vl[0].parts) * vl[0].qty for vl in vls) for t, vls in groups.items()}
    order = sorted(groups)
    all_placed, all_used, all_prods, all_table = [], [], [], []
    for gi, t in enumerate(order):
        _check_cancel()
        rest = sum(size[x] for x in order[gi:])
        tl = max(0.3, (end - time.time()) * size[t] / rest)
        pre = f"{t:g} mm ({gi + 1}/{len(order)}): "
        prog = (lambda m, pre=pre: progress(pre + m)) if progress else None
        placed, used, prods, table = optimize_with_variants(
            groups[t], sheets, kerf, guillotine, weld, min_part_pct, time_limit=tl,
            progress=prog, pool=pool, workers=workers)
        off = len(all_used)
        unit_off = sum(p.qty for p in all_prods)
        single = prods[0].name if len(prods) == 1 else None
        for pp in placed:
            pp.sheet_index += off
            pp.unit += unit_off
            m = re.fullmatch(r"K(\d+)", pp.tag)
            if single and m:
                # egy termékes csoport: a több termékes jelölés (név + sorszám) kell
                pp.tag = f"{single} {m.group(1)}."
        all_placed += placed
        all_used += [dataclasses.replace(u, thick=t) for u in used]
        all_prods += prods
        all_table += [(f"{t:g} mm: {d}", n, w, u) for d, n, w, u in table]
    # egy termék egyszer szerepeljen (a vastagság szerinti részek közül a fő rész)
    merged: Dict[str, Product] = {}
    for p in all_prods:
        cur = merged.get(p.name)
        if cur is None or (cur.sub_of and not p.sub_of):
            merged[p.name] = p
    extra: Dict[str, Dict[str, float]] = {}
    for p in all_prods:
        if p.sub_of:
            extra.setdefault(p.name, {})[p.sub_of] = p.thick
    out_prods = []
    for name, p in merged.items():
        if name in extra:
            p = dataclasses.replace(
                p, part_thick=extra[name],
                info=(p.info + " | " if p.info else "") + "eltérő lemezvastagság – " +
                "; ".join(f"{n}: {t:g} mm" for n, t in extra[name].items()))
        out_prods.append(p)
    return all_placed, all_used, out_prods, all_table


def weight_summary(placed: List[PlacedPiece], sheets: List[SheetSize],
                   products: Optional[List["Product"]] = None) -> List[str]:
    """A felhasznált táblák súlya vastagságonként és összesen (bruttó = teljes tábla,
    nettó = a beépülő alkatrészek)."""
    lines = [f"Lemezek súlya (acél, {STEEL_KG_PER_M2_MM:g} kg / m² / mm):"]
    by_t: Dict[float, List[int]] = {}
    for i, s in enumerate(sheets):
        by_t.setdefault(s.thick, []).append(i)
    tot_b = tot_n = 0.0
    for t in sorted(by_t):
        idx = set(by_t[t])
        area = sum(sheets[i].area for i in idx) / 1e6
        net = sum(p.net_area for p in placed if p.sheet_index in idx) / 1e6
        wb, wn = area * t * STEEL_KG_PER_M2_MM, net * t * STEEL_KG_PER_M2_MM
        tot_b += wb
        tot_n += wn
        lines.append(f"  {t:g} mm: {len(idx)} tábla, {area:.2f} m² → bruttó {wb:.1f} kg "
                     f"(beépülő {net:.2f} m² → nettó {wn:.1f} kg)")
    lines.append(f"  ÖSSZES SÚLY: bruttó {tot_b:.1f} kg, nettó {tot_n:.1f} kg "
                 f"(hulladék {tot_b - tot_n:.1f} kg)")
    if products:
        for pr in products:
            th = {i: sheets[i].thick for i in range(len(sheets))}
            net = sum(p.net_area / 1e6 * th[p.sheet_index] * STEEL_KG_PER_M2_MM for p in placed
                      if p.tag.rsplit(" ", 1)[0] == pr.name or (len(products) == 1))
            tt = f"{pr.thick:g} mm" + ("".join(f", {n.lower()} {t:g} mm" for n, t in pr.part_thick.items())
                                       if pr.part_thick else "")
            lines.append(f"  {pr.name} ({tt}): nettó {net:.1f} kg összesen, "
                         f"{net / pr.qty:.1f} kg/db" if pr.qty else f"  {pr.name}: –")
    return lines


def optimize_pieces(base_parts: List[Tuple[int, int, str]], qty: int,
                    sheets: List[SheetSize], kerf: int = 0, guillotine: bool = False,
                    weld: bool = True, min_part_pct: float = 40.0, time_limit: float = 15.0):
    """Egy termék (egyszerűsített hívási forma)."""
    return optimize_products([Product("K", qty, base_parts)], sheets, kerf,
                             guillotine, weld, min_part_pct, time_limit=time_limit)


# ---------- DXF export ----------

def placed_polys(p: PlacedPiece) -> List[Shape]:
    """A darab alkatrész-alakjai a lemez koordinátáiban (forgatással együtt)."""
    shapes = p.shapes or [(p.base or p.name, rect_poly(p.width, p.height))]
    out = []
    for lab, poly in shapes:
        if p.rotated:
            nat_h = p.width            # természetes magasság = elhelyezett szélesség
            pts = [(p.x + (nat_h - v), p.y + u) for u, v in poly]
        else:
            pts = [(p.x + u, p.y + v) for u, v in poly]
        out.append((lab, pts))
    return out


def export_dxf(filename: str, placed: List[PlacedPiece], sheets: List[SheetSize]):
    import ezdxf
    doc = ezdxf.new("R2010")
    msp = doc.modelspace()
    gap = 200

    if "RESIDUAL" not in doc.layers:
        doc.layers.new(name="RESIDUAL", dxfattribs={"color": 1})

    base_y = 0
    for i, sheet in enumerate(sheets):
        tt = f" {sheet.thick:g} mm" if sheet.thick else ""
        layer = f"SHEET_{i + 1}" + (f"_{sheet.thick:g}mm".replace(".", ",") if sheet.thick else "")
        if layer not in doc.layers:
            doc.layers.new(name=layer)
        msp.add_lwpolyline([(0, base_y), (sheet.nw, base_y),
                            (sheet.nw, base_y + sheet.nh),
                            (0, base_y + sheet.nh)],
                           close=True, dxfattribs={"layer": layer})
        if (sheet.width, sheet.height) != (sheet.nw, sheet.nh):
            # a vágható terület határa (széltávolság / max. vágási hossz)
            o = sheet.off
            msp.add_lwpolyline([(o, base_y + o), (o + sheet.width, base_y + o),
                                (o + sheet.width, base_y + o + sheet.height),
                                (o, base_y + o + sheet.height)],
                               close=True, dxfattribs={"layer": layer, "linetype": "DASHED"}
                               if "DASHED" in doc.linetypes else {"layer": layer})
        title = msp.add_text(f"{i + 1}. lemez {sheet.nw}x{sheet.nh}{tt}",
                             dxfattribs={"height": 60, "layer": layer})
        title.dxf.insert = (sheet.nw + 50, base_y + sheet.nh - 60)

        for p in placed:
            if p.sheet_index != i:
                continue
            lay = "RESIDUAL" if p.from_residual else layer
            th = max(10, min(50, p.height / 5, p.width / 12))
            o = sheet.off
            for lab, pts in placed_polys(p):
                msp.add_lwpolyline([(x + o, y + base_y + o) for x, y in pts], close=True,
                                   dxfattribs={"layer": lay})
                cx = sum(x for x, _ in pts) / len(pts) + o
                cy = sum(y for _, y in pts) / len(pts) + o
                t = msp.add_text(f"{p.tag} {lab}", dxfattribs={"height": th, "layer": lay})
                t.dxf.insert = (cx - th * len(lab) / 3, cy + base_y)
        base_y += sheet.nh + gap

    doc.saveas(filename)


# ---------- Lézer export: minden alkatrész külön DXF-be ----------

def _norm_poly(poly: Poly) -> Poly:
    mx = min(x for x, _ in poly)
    my = min(y for _, y in poly)
    out = []
    for x, y in poly:
        p = (round(x - mx, 2), round(y - my, 2))
        if not out or (abs(out[-1][0] - p[0]) > 0.01 or abs(out[-1][1] - p[1]) > 0.01):
            out.append(p)
    if len(out) > 1 and abs(out[0][0] - out[-1][0]) < 0.01 and abs(out[0][1] - out[-1][1]) < 0.01:
        out.pop()
    # kollineáris pontok elhagyása
    clean = []
    n = len(out)
    for i in range(n):
        a, b, c = out[i - 1], out[i], out[(i + 1) % n]
        if abs((b[0] - a[0]) * (c[1] - b[1]) - (b[1] - a[1]) * (c[0] - b[0])) > 1e-3:
            clean.append(b)
    return clean if len(clean) >= 3 else out


def _stack_union(shapes: List[Shape]):
    """Egymás fölé illesztett (alul-felül vízszintes élű) lapok közös kontúrja
    + a köztük lévő hajlítási vonalak. Hajlított, egy lemezből készülő darabhoz."""
    parts = []
    for _, poly in shapes:
        ys = [y for _, y in poly]
        y0, y1 = min(ys), max(ys)
        bot = sorted(x for x, y in poly if abs(y - y0) < 0.01)
        top = sorted(x for x, y in poly if abs(y - y1) < 0.01)
        if not bot or not top:
            return None
        parts.append((y0, y1, bot[0], bot[-1], top[0], top[-1]))
    parts.sort()
    right = []
    left = []
    bends = []
    for i, (y0, y1, bl, br, tl, tr) in enumerate(parts):
        right += [(br, y0), (tr, y1)]
        left += [(bl, y0), (tl, y1)]
        if i + 1 < len(parts):
            n = parts[i + 1]
            if abs(n[0] - y1) > 0.5:
                return None              # nem érintkeznek – nem egy darab
            xa, xb = max(tl, n[2]), min(tr, n[3])
            if xb > xa:
                bends.append(((xa, y1), (xb, y1)))
    outline = right + list(reversed(left))
    return outline, bends


def _canonical(poly: Poly, bends, holes=()):
    """Az alkatrész egységes állása: a 8 lehetséges forgatás/tükrözés közül a
    lexikografikusan legkisebb – így a 180°-kal elforgatott vagy tükrözött (bal/jobb)
    azonos alkatrészek egy tételbe kerülnek. A szélesebb oldal vízszintes."""
    best = None
    for rot in range(4):
        for mir in (False, True):
            def T(p, rot=rot, mir=mir):
                x, y = p
                if mir:
                    x = -x
                for _ in range(rot):
                    x, y = -y, x
                return (x, y)
            pts = [T(p) for p in poly]
            mx = min(x for x, _ in pts)
            my = min(y for _, y in pts)
            pts = [(round(x - mx, 1), round(y - my, 1)) for x, y in pts]
            w = max(x for x, _ in pts)
            h = max(y for _, y in pts)
            if w + 0.05 < h:
                continue
            a = sum(pts[i][0] * pts[(i + 1) % len(pts)][1] - pts[(i + 1) % len(pts)][0] * pts[i][1]
                    for i in range(len(pts)))
            if a < 0:
                pts = list(reversed(pts))
            k0 = min(range(len(pts)), key=lambda i: pts[i])
            pts = pts[k0:] + pts[:k0]
            bb = []
            for a_, b_ in bends:
                pa, pb = T(a_), T(b_)
                pa = (round(pa[0] - mx, 1), round(pa[1] - my, 1))
                pb = (round(pb[0] - mx, 1), round(pb[1] - my, 1))
                bb.append(tuple(sorted((pa, pb))))
            hh = tuple(sorted(tuple((round(T(q)[0] - mx, 1), round(T(q)[1] - my, 1)) for q in h)
                              for h in holes))
            key = (tuple(pts), tuple(sorted(bb)), hh)
            if best is None or key < best:
                best = key
    pts, bb, hh = best
    if holes:
        return list(pts), [(a, b) for a, b in bb], [list(h) for h in hh]
    return list(pts), [(a, b) for a, b in bb]


def laser_parts(placed: List[PlacedPiece], sheets: Optional[List[SheetSize]] = None):
    """
    A kiszabandó alkatrészek egyenként (a lemezen elfoglalt helyüktől függetlenül):
    [(név, kontúr, hajlítási vonalak, darabszám, konténerek)] – az azonos alakúak
    összevonva, darabszámmal.
      • egymásba fordított pár (pl. ABO_A oldalpár) → külön-külön alkatrész
      • egy lemezből hajlított darab (pl. Alja+Hátlap) → egy kontúr + hajlítási vonal
      • toldott részek → a résznek megfelelő kontúr
    """
    items: Dict[tuple, dict] = {}
    for pp in placed:
        shapes = pp.shapes or [(pp.base or pp.name, rect_poly(pp.width, pp.height))]
        bent = "+" in pp.name and "pár" not in pp.name and len(shapes) > 1
        groups = []
        if bent:
            u = _stack_union(shapes)
            if u:
                groups.append((pp.name.split(" #")[0].replace("_", " "), u[0], u[1]))
            else:
                bent = False
        if not bent:
            for lab, poly in shapes:
                if not is_hole(lab):          # a furatok a saját alkatrészükhöz tartoznak
                    groups.append((lab, poly, []))
        for lab, poly, bends in groups:
            mx = min(x for x, _ in poly)
            my = min(y for _, y in poly)
            sb = [((p[0] - mx, p[1] - my), (q[0] - mx, q[1] - my)) for p, q in bends]
            holes = [[(x - mx, y - my) for x, y in hp] for hl, hp in shapes if hl == lab + HOLE_SUFFIX]
            if holes:
                # _norm_poly a saját minimumához tol el – a furatokat ugyanígy kell
                npoly0 = _norm_poly(poly)
                npoly, nb, nholes = _canonical(npoly0, sb, holes)
            else:
                npoly, nb = _canonical(_norm_poly(poly), sb)
                nholes = []
            base = re.sub(r"\s+\d+$", "", lab)          # "Oldal 1" és "Oldal 2" egy tétel
            if pp.part > 0:
                base = f"{base} toldás"
            th = sheets[pp.sheet_index].thick if sheets else 0
            if th:
                base = f"{th:g}mm {base}"
            key = (base, tuple(npoly), tuple(nb), tuple(tuple(h) for h in nholes))
            it = items.setdefault(key, {"name": base, "poly": npoly, "bends": nb, "holes": nholes,
                                        "count": 0, "tags": set()})
            it["count"] += 1
            it["tags"].add(pp.tag.rsplit(" ", 1)[0] if " " in pp.tag else pp.tag)
    out = list(items.values())
    out.sort(key=lambda d: (d["name"], -max(x for x, _ in d["poly"]) * max(y for _, y in d["poly"])))
    return out


def _safe_name(s: str) -> str:
    tr = str.maketrans("áéíóöőúüűÁÉÍÓÖŐÚÜŰ", "aeiooouuuAEIOOOUUU")
    s = s.translate(tr)
    s = re.sub(r"[^A-Za-z0-9_+-]+", "_", s).strip("_")
    return s or "alkatresz"


def export_laser(folder: str, placed: List[PlacedPiece], prefix: str = "",
                 sheets: Optional[List[SheetSize]] = None) -> List[str]:
    """Minden alkatrész-típus külön DXF-be (mm): a teljes kiterített kontúr a VAGAS
    rétegen (a perem/hajlítás mérete benne van, hajlítási vonal nélkül), felirat a
    FELIRAT rétegen; az azonos (forgatott/tükrözött) alkatrészek egy fájlba kerülnek.
    + darabjegyzek.txt."""
    import ezdxf
    os.makedirs(folder, exist_ok=True)
    parts = laser_parts(placed, sheets)
    written, lines = [], []
    for i, p in enumerate(parts, start=1):
        w = max(x for x, _ in p["poly"])
        h = max(y for _, y in p["poly"])
        fname = f"{i:02d}_{_safe_name(prefix + p['name'])}_{w:.0f}x{h:.0f}_{p['count']}db.dxf"
        doc = ezdxf.new("R2010", setup=True)   # setup: szabványos vonaltípusok (DASHED)
        doc.header["$INSUNITS"] = 4          # mm
        doc.layers.new(name="VAGAS", dxfattribs={"color": 7})
        doc.layers.new(name="FELIRAT", dxfattribs={"color": 3})
        msp = doc.modelspace()
        # csak a teljes (kiterített) kontúr – a perem és a hajlítások mérete benne van,
        # hajlítási vonalat nem rajzolunk
        msp.add_lwpolyline(p["poly"], close=True, dxfattribs={"layer": "VAGAS"})
        for hp in p.get("holes", []):
            msp.add_lwpolyline(hp, close=True, dxfattribs={"layer": "VAGAS"})
        th = max(10.0, min(40.0, h / 10))
        t = msp.add_text(f"{p['name']} {p['count']} db", dxfattribs={"height": th, "layer": "FELIRAT"})
        t.dxf.insert = (w * 0.1, h * 0.5)
        path = os.path.join(folder, fname)
        doc.saveas(path)
        written.append(path)
        lines.append(f"{fname}\t{p['name']}\t{w:.1f} x {h:.1f} mm\t{p['count']} db")
    with open(os.path.join(folder, "darabjegyzek.txt"), "w", encoding="utf-8") as f:
        f.write("Fájl\tAlkatrész\tBefoglaló méret\tDarab\n")
        f.write("\n".join(lines) + "\n")
        f.write(f"\nÖsszesen {sum(p['count'] for p in parts)} db alkatrész, {len(parts)} fájl.\n")
    return written



# ---------- Egyéb alkatrész DXF-ből ----------

_DXF_UNITS = {0: 1.0, 1: 25.4, 2: 304.8, 4: 1.0, 5: 10.0, 6: 1000.0}   # $INSUNITS → mm


def _point_in_poly(pt, poly) -> bool:
    x, y = pt
    inside = False
    n = len(poly)
    for i in range(n):
        x1, y1 = poly[i]
        x2, y2 = poly[(i + 1) % n]
        if (y1 > y) != (y2 > y) and x < x1 + (y - y1) * (x2 - x1) / (y2 - y1):
            inside = not inside
    return inside


def _simplify(poly: Poly, tol: float = 0.2) -> Poly:
    """Közel egy egyenesbe eső és egybeeső pontok elhagyása (ívek közelítése után)."""
    out = []
    for p in poly:
        if not out or math.hypot(p[0] - out[-1][0], p[1] - out[-1][1]) > 1e-3:
            out.append(p)
    if len(out) > 2 and math.hypot(out[0][0] - out[-1][0], out[0][1] - out[-1][1]) <= 1e-3:
        out.pop()
    changed = True
    while changed and len(out) > 3:
        changed = False
        keep = []
        n = len(out)
        prev_removed = False
        for i in range(n):
            a, b, c = out[i - 1], out[i], out[(i + 1) % n]
            L = math.hypot(c[0] - a[0], c[1] - a[1])
            d = abs((c[0] - a[0]) * (a[1] - b[1]) - (a[0] - b[0]) * (c[1] - a[1])) / L if L else 0
            if d < tol and not prev_removed and i > 0:
                changed = True          # szomszédos pontokat egy körben nem hagyunk el
                prev_removed = True
                continue
            prev_removed = False
            keep.append(b)
        out = keep if len(keep) >= 3 else out
    return out


def _min_box_angle(pts) -> float:
    """A legkisebb területű befoglaló téglalap szöge (a konvex burok élei mentén)."""
    P = sorted(set((round(x, 3), round(y, 3)) for x, y in pts))
    if len(P) < 3:
        return 0.0

    def cross(o, a, b):
        return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])
    lo, hi = [], []
    for p in P:
        while len(lo) >= 2 and cross(lo[-2], lo[-1], p) <= 0:
            lo.pop()
        lo.append(p)
    for p in reversed(P):
        while len(hi) >= 2 and cross(hi[-2], hi[-1], p) <= 0:
            hi.pop()
        hi.append(p)
    hull = lo[:-1] + hi[:-1]
    best = (None, 0.0)
    for i in range(len(hull)):
        a, b = hull[i], hull[(i + 1) % len(hull)]
        ang = math.atan2(b[1] - a[1], b[0] - a[0])
        c, s_ = math.cos(-ang), math.sin(-ang)
        xs = [x * c - y * s_ for x, y in hull]
        ys = [x * s_ + y * c for x, y in hull]
        area = (max(xs) - min(xs)) * (max(ys) - min(ys))
        if best[0] is None or area < best[0] - 1e-6:
            best = (area, ang)
    return best[1]


def load_dxf_parts(path: str, name: str) -> List[PartDef]:
    """
    Rajzolt alkatrész(ek) beolvasása DXF-ből. A zárt kontúrokat (vonalláncok, vonalak +
    ívek összefűzve, körök, spline-ok) keresi; a legkülső kontúrok külön alkatrészek, a
    bennük lévő kontúrok furatok. Az alkatrészt a legkisebb befoglaló téglalap szerint
    elforgatja (ferdén rajzolt darab is jól szabható). Méretek mm-ben ($INSUNITS szerint).
    """
    import ezdxf
    from ezdxf import path as ezpath
    try:
        doc = ezdxf.readfile(path)
    except Exception as e:
        raise ValueError(f"A DXF fájl nem olvasható: {e}")
    scale = _DXF_UNITS.get(int(doc.header.get("$INSUNITS", 0) or 0), 1.0)
    kinds = {"LINE", "ARC", "CIRCLE", "LWPOLYLINE", "POLYLINE", "SPLINE", "ELLIPSE"}

    def walk(entities, depth=0):
        for e in entities:
            t = e.dxftype()
            if t == "INSERT" and depth < 5:
                try:
                    yield from walk(e.virtual_entities(), depth + 1)
                except Exception:
                    continue
            elif t in kinds:
                yield e

    loops, opens = [], []
    for e in walk(doc.modelspace()):
        try:
            p = ezpath.make_path(e)
            pts = [(v.x * scale, v.y * scale) for v in p.flattening(0.2 / scale)]
        except Exception:
            continue
        if len(pts) < 2:
            continue
        closed = math.hypot(pts[0][0] - pts[-1][0], pts[0][1] - pts[-1][1]) < 0.05
        if e.dxftype() == "CIRCLE" or closed:
            loops.append(pts)
        else:
            opens.append(pts)
    # mértékegység-ellenőrzés: ha a megadott egységgel 20 m-nél nagyobb lenne a rajz,
    # valójában mm-ben rajzolták (sok program alapból méter egységet ír a fejlécbe)
    allp = [p for l in loops + opens for p in l] if (loops or opens) else []
    if allp and scale > 1:
        ext = max(max(x for x, _ in allp) - min(x for x, _ in allp),
                  max(y for _, y in allp) - min(y for _, y in allp))
        if ext > 20000:
            loops = [[(x / scale, y / scale) for x, y in l] for l in loops]
            opens = [[(x / scale, y / scale) for x, y in l] for l in opens]
    # nyitott elemek (vonalak, ívek) összefűzése zárt kontúrokká
    tol = 0.5
    while opens:
        chain = opens.pop()
        grew = True
        while grew:
            grew = False
            for i, seg in enumerate(opens):
                for rev in (False, True):
                    s_ = list(reversed(seg)) if rev else seg
                    if math.hypot(chain[-1][0] - s_[0][0], chain[-1][1] - s_[0][1]) < tol:
                        chain = chain + s_[1:]
                    elif math.hypot(chain[0][0] - s_[-1][0], chain[0][1] - s_[-1][1]) < tol:
                        chain = s_[:-1] + chain
                    else:
                        continue
                    opens.pop(i)
                    grew = True
                    break
                if grew:
                    break
            if math.hypot(chain[0][0] - chain[-1][0], chain[0][1] - chain[-1][1]) < tol and len(chain) > 3:
                break
        if math.hypot(chain[0][0] - chain[-1][0], chain[0][1] - chain[-1][1]) < tol and len(chain) > 3:
            loops.append(chain)
    loops = [_simplify(l) for l in loops]
    loops = [l for l in loops if len(l) >= 3 and poly_area(l) > 1.0]
    if not loops:
        raise ValueError("A DXF-ben nem találtam zárt kontúrt (vonallánc, kör vagy összefüggő "
                         "vonalak/ívek). A méretvonalakat és szövegeket a program figyelmen kívül hagyja.")
    loops.sort(key=poly_area, reverse=True)
    depth = [sum(1 for j in range(i) if _point_in_poly(loops[i][0], loops[j])) for i in range(len(loops))]
    outers = [i for i in range(len(loops)) if depth[i] % 2 == 0]
    parts = []
    for k, oi in enumerate(outers, start=1):
        holes = [loops[i] for i in range(len(loops)) if depth[i] == depth[oi] + 1
                 and _point_in_poly(loops[i][0], loops[oi])
                 and not any(depth[j] == depth[oi] + 1 and j != i and
                             _point_in_poly(loops[i][0], loops[j]) for j in range(len(loops)))]
        ang = _min_box_angle(loops[oi])
        c, s_ = math.cos(-ang), math.sin(-ang)
        rot = lambda P: [(x * c - y * s_, x * s_ + y * c) for x, y in P]
        o = rot(loops[oi])
        hs = [rot(h) for h in holes]
        mx, my = min(x for x, _ in o), min(y for _, y in o)
        o = [(x - mx, y - my) for x, y in o]
        hs = [[(x - mx, y - my) for x, y in h] for h in hs]
        w, h = max(x for x, _ in o), max(y for _, y in o)
        if h > w:                      # fekvő állás (a szabás úgyis forgat)
            o = [(y, w - x) for x, y in o]
            hs = [[(y, w - x) for x, y in hh] for hh in hs]
            w, h = h, w
        nm = name if len(outers) == 1 else f"{name} {k}"
        shapes = [(nm, o)] + [(nm + HOLE_SUFFIX, hh) for hh in hs]
        pd = PartDef(nm, int(math.ceil(round(w, 2))), int(math.ceil(round(h, 2))), shapes, no_weld=True)
        parts.append(pd)
    return parts


def dxf_parts_summary(parts: List[PartDef]) -> str:
    lines = []
    for pd in parts:
        net = sum(-poly_area(p) if is_hole(l) else poly_area(p) for l, p in pd.shapes)
        nh = sum(1 for l, _ in pd.shapes if is_hole(l))
        lines.append(f"{pd.name}: {pd.width} x {pd.height} mm, nettó {net / 1e6:.3f} m²"
                     + (f", {nh} furat" if nh else ""))
    return "\n".join(lines)


def other_sheet_report(products: List["Product"], sheets: List[SheetSize], kerf: int,
                       guillotine: bool) -> List[str]:
    """Egyéb (DXF) alkatrészek: táblaméretenként hány tábla kellene, ha csak abból
    szabnánk – vastagságonként. Segít eldönteni, melyik táblából érdemes kiszabni."""
    others = [p for p in products if p.ctype == OTHER_TYPE]
    if not others:
        return []
    lines = ["", "Egyéb (DXF) alkatrészek – ha csak egyféle táblából szabnánk:"]
    by_t: Dict[float, List[Product]] = {}
    for p in others:
        by_t.setdefault(float(p.thick), []).append(p)
    for t, prods in sorted(by_t.items()):
        pieces = []
        for pr in prods:
            for u in range(1, pr.qty + 1):
                for pd in pr.parts:
                    pieces.append(Piece(pd.width, pd.height, pd.name, u, tag=f"{pr.name} {u}.",
                                        shapes=list(pd.shapes), no_weld=True))
        net = sum(sum(-poly_area(p) if is_hole(l) else poly_area(p) for l, p in pc.shapes)
                  for pc in pieces)
        names = ", ".join(f"{p.name} ({p.qty} db)" for p in prods)
        lines.append(f"  {t:g} mm – {names} (gyors becslés):")
        rows = []
        for s in sorted(sheets, key=lambda s: s.area):
            if not all(fits_any(pc.width, pc.height, [s]) for pc in pieces):
                lines.append(f"    {s.nw} x {s.nh}: nem fér rá minden darab")
                continue
            try:
                _, placed, used = _search(pieces, [s], kerf, guillotine, True,
                                          deadline=time.time() + 1.5)
            except ValueError:
                continue
            area = sum(u.area for u in used)
            rows.append((area, len(used), s))
            lines.append(f"    {s.nw} x {s.nh}: {len(used)} tábla, {area / 1e6:.2f} m², "
                         f"kihasználtság {100 * net / area:.1f} %, súly {area / 1e6 * t * STEEL_KG_PER_M2_MM:.1f} kg")
        if rows:
            best = min(rows, key=lambda r: (r[0], r[1]))
            lines.append(f"    → legkevesebb anyag: {best[2].nw} x {best[2].nh} ({best[1]} tábla)")
    return lines

# ---------- Export árajánlathoz (Acél árkalkulátor, SZABAS_EXPORT.md) ----------

DEFAULT_THICKNESS = "2"          # mm – alapértelmezett lemezvastagság (kézzel átírható)

MATERIAL_NAMES = [("hideg", "Hidegen hengerelt"), ("meleg", "Melegen hengerelt"),
                  ("horgany", "Horganyzott"), ("rozsdamentes", "Rozsdamentes"),
                  ("inox", "Rozsdamentes")]


def load_material_info(path: str) -> Tuple[Optional[float], Optional[str]]:
    """
    Az anyagfájlból a lemezvastagság és az anyag, ha meg van adva, pl.
        vastagság: 2 mm          vagy egy táblasorban: 2mm 1250x2500
        anyag: hidegen hengerelt (hideg, meleg, horganyzott, rozsdamentes/inox)
    Vastagságot csak akkor ad vissza, ha a fájlban egyetlen érték szerepel.
    """
    thick, mat = set(), None
    for raw in _read_text(path).splitlines():
        line = raw.split("#", 1)[0].lower()
        for t in re.findall(r"(\d+(?:[.,]\d+)?)\s*mm\b", line):
            v = float(t.replace(",", "."))
            if 0 < v < 100:
                thick.add(v)
        m = re.search(r"vastags[aá]g\w*\s*[:=]?\s*(\d+(?:[.,]\d+)?)", line)
        if m:
            thick.add(float(m.group(1).replace(",", ".")))
        if mat is None:
            for key, name in MATERIAL_NAMES:
                if key in line:
                    mat = name
                    break
    return (thick.pop() if len(thick) == 1 else None), mat


def export_arajanlathoz(path, projekt, termekek, tablak, beepulo_m2, tabla_m2,
                        vastagsag_mm=None, anyag=None, suly_kg=None):
    """termekek: [(nev, tipus, darab)], tablak: [(szel_mm, hossz_mm, db[, vastagság_mm])].
    Több vastagságnál a vastagság táblánként szerepel, a közös vastagsag_mm kimarad."""
    adat = {
        "forras": "Konténer lemezszabás kalkulátor",
        "verzio": 1,
        "projekt": projekt,
        "termekek": [{"nev": n, "tipus": t, "darab": int(d)} for n, t, d in termekek],
        "tablak": [dict({"szel": int(w), "hossz": int(h), "db": int(n)},
                        **({"vastagsag_mm": float(t[0])} if t and t[0] else {}))
                   for w, h, n, *t in tablak],
        "tabla_terulet_m2": round(tabla_m2, 3),
        "beepulo_terulet_m2": round(beepulo_m2, 3),
    }
    if suly_kg:
        adat["suly_kg"] = round(suly_kg, 1)
    if vastagsag_mm:
        adat["vastagsag_mm"] = float(vastagsag_mm)
    if anyag:
        adat["anyag"] = anyag
    with open(path, "w", encoding="utf-8") as f:
        f.write(json.dumps(adat, ensure_ascii=False, indent=2))


def quote_export_data(placed: List[PlacedPiece], sheets: List[SheetSize],
                      products: List["Product"]):
    """Ugyanazokból az adatokból, mint a „Felhasznált lemezek száma” blokk."""
    summary: Dict[Tuple[float, int, int], int] = {}
    for s in sheets:
        summary[(s.thick, s.nw, s.nh)] = summary.get((s.thick, s.nw, s.nh), 0) + 1
    tablak = [(w, h, n, t) for (t, w, h), n in sorted(summary.items())]
    termekek = [(p.name, p.ctype, p.qty) for p in products]
    types = list(dict.fromkeys(p.ctype for p in products if p.ctype))
    projekt = "+".join(types) if types else (products[0].name if products else "szabas")
    return (projekt, termekek, tablak,
            sum(p.net_area for p in placed) / 1e6, sum(s.area for s in sheets) / 1e6)


# ---------- Szöveges eredmény ----------

def format_result(placed: List[PlacedPiece], sheets: List[SheetSize], kerf: int,
                  products: Optional[List["Product"]] = None) -> str:
    lines = []
    total_area = sum(s.area for s in sheets)
    used_area = sum(p.net_area for p in placed)
    has_price = all(s.price is not None for s in sheets)

    lines.append(f"Felhasznált lemezek száma: {len(sheets)}")
    summary: Dict[Tuple[float, int, int], int] = {}
    for s in sheets:
        summary[(s.thick, s.nw, s.nh)] = summary.get((s.thick, s.nw, s.nh), 0) + 1
    for (t, w, h), n in sorted(summary.items()):
        tt = f"  ({t:g} mm)" if t else ""
        lines.append(f"  {n} db  {w} x {h} mm{tt}")
    lines.append(f"Összes táblaterület: {total_area / 1e6:.2f} m²")
    lines.append(f"Beépülő terület:     {used_area / 1e6:.2f} m²")
    lines.append(f"Kihasználtság:       {100 * used_area / total_area:.1f} %")
    if has_price:
        lines.append(f"Anyagköltség:        {sum(s.price for s in sheets):,.0f}".replace(",", " "))
    if kerf:
        lines.append(f"Vágási rés:          {kerf} mm")
    if any(s.thick for s in sheets):
        lines.append("")
        lines += weight_summary(placed, sheets, products)

    if products and len(products) > 1:
        lines.append("")
        lines.append("Termékenkénti anyagfelhasználás (a közös táblák területarányosan elosztva):")
        for prod in products:
            net = sum(p.net_area for p in placed
                      if p.tag.rsplit(" ", 1)[0] == prod.name)
            same = [i for i, s in enumerate(sheets) if s.thick == getattr(prod, "thick", s.thick)]
            g_used = sum(p.net_area for p in placed if p.sheet_index in set(same))
            share = net / g_used if g_used else 0
            lines.append(f"  {prod.name}: {prod.qty} db, beépülő {net / 1e6:.2f} m², "
                         f"arányos táblaigény {share * len(same):.1f} tábla "
                         f"({share * len(same) / prod.qty:.2f} tábla/db)")
        lines.append("")

    n_res = sum(1 for p in placed if p.from_residual)
    lines.append(f"Másik konténer táblájának maradékából szabott darab: {n_res} db")
    welds = count_welds(placed)
    seam_total = sum(p.seam for p in placed)
    if products and any(getattr(pr, "info", "") for pr in products):
        lines.append("")
        for pr in products:
            if getattr(pr, "info", ""):
                lines.append(f"{pr.name}: {pr.info}")
    if products and any(pr.assembly_weld for pr in products):
        asm_total = sum(pr.assembly_weld * pr.qty for pr in products)
        lines.append("")
        lines.append("Hegesztés összesen:")
        for pr in products:
            if pr.assembly_weld:
                lines.append(f"  {pr.name}: összeállítás {pr.assembly_weld / 1000:.2f} m/db × {pr.qty} db "
                             f"= {pr.assembly_weld * pr.qty / 1000:.2f} m")
        lines.append(f"  Toldási varratok (lemezek összetoldása): {seam_total / 1000:.2f} m")
        lines.append(f"  ÖSSZES HEGESZTÉS: {(asm_total + seam_total) / 1000:.2f} m")
        lines.append("")
    lines.append(f"Hegesztési toldások (varratok) száma: {welds} db, "
                 f"összhossza: {seam_total / 1000:.2f} m")
    if welds:
        groups: Dict[Tuple[int, str], List[Tuple[PlacedPiece, float]]] = {}
        for p in placed:
            for lab, poly in (p.shapes or [(p.base or p.name, rect_poly(p.width, p.height))]):
                groups.setdefault((p.unit, lab), []).append((p, poly_area(poly)))
        for (unit, lab), ps in sorted(groups.items()):
            if len(ps) < 2:
                continue
            desc = " + ".join(f"{a / 1e6:.2f} m² ({p.sheet_index + 1}. lemez)"
                              for p, a in sorted(ps, key=lambda t: (t[0].part, t[0].sheet_index)))
            seam_lab = sum(p.seam for p, _ in ps if len(p.shapes) <= 1)
            seam_txt = f", varrat {seam_lab / 1000:.2f} m" if seam_lab else ""
            lines.append(f"  {ps[0][0].tag} {lab}: {len(ps)} részből{seam_txt}: {desc}")
    lines.append("")

    for i, s in enumerate(sheets):
        on = [p for p in placed if p.sheet_index == i]
        ua = sum(p.net_area for p in on)
        lines.append(f"--- {i + 1}. lemez  {s.label()}  "
                     f"(kihasználtság {100 * ua / s.area:.1f} %) ---")
        for p in sorted(on, key=lambda p: (p.y, p.x)):
            tags = []
            if p.rotated:
                tags.append("forgatva")
            if p.from_residual:
                tags.append("maradékból")
            tag = f"  [{', '.join(tags)}]" if tags else ""
            inner = ""
            if p.shapes and not (len(p.shapes) == 1 and abs(p.net_area - p.width * p.height) < 1):
                nh = sum(1 for lab, _ in p.shapes if is_hole(lab))
                inner = "  {" + ", ".join(f"{lab} {poly_area(pl) / 1e6:.2f} m²" for lab, pl in p.shapes
                                          if not is_hole(lab)) + \
                    (f", {nh} furat, nettó {p.net_area / 1e6:.2f} m²" if nh else "") + "}"
            lines.append(f"  {p.tag:<10} {p.name:<14} {p.width} x {p.height} mm"
                         f"  pozíció: ({p.x}, {p.y}){tag}{inner}")
        lines.append("")
    return "\n".join(lines)


def format_variant_table(table) -> str:
    if not table:
        return ""
    lines = ["", "=== Összevonási változatok összevetése (kihasználtság, táblaszám, majd hegesztés szerint) ==="]
    for i, (desc, n, weld_mm, util) in enumerate(table):
        mark = "→ " if i == 0 else "  "
        lines.append(f"{mark}kihasználtság {util:.1f} %, {n} tábla, {weld_mm / 1000:.2f} m hegesztés – {desc}")
    lines.append("(A megjelölt változat került kiválasztásra. Más változatot az "
                 "'Alja összehajlítása' listában lehet kényszeríteni.)")
    return "\n".join(lines)



# ---------- ABO_A rajzos méretmegadás – a vázlat geometriája ----------
# Tisztán adat (nincs benne Qt), így a rajz külön is ellenőrizhető.

SKETCH_W, SKETCH_H = 780, 650
INPUT_W, INPUT_H = 84, 26

# mezőnév → (táblázat oszlop, felirat a súgóhoz)
ABO_A_FIELDS = {
    "side_C": (9,  "Oldallap C – felső él hossza"),
    "side_B": (10, "Oldallap B – nettó magasság (perem nélkül)"),
    "side_A": (11, "Oldallap A – alsó él hossza"),
    "hat_A":  (6,  "Hátlap A – nettó magasság"),
    "hat_B":  (5,  "Hátlap B – szélesség"),
    "hom_B":  (7,  "Homloklap B – szélesség"),
    "alja_A": (4,  "Alja A – hossz"),
    "alja_B": (3,  "Alja B – szélesség"),
}


WIDTH_SOURCE = "hat_B"                 # ezt írja be a felhasználó (konténer szélesség)
WIDTH_TARGETS = ("hom_B", "alja_A")    # ide automatikusan átkerül, de átírható


def abo_a_sketch() -> dict:
    """
    A rajz elemei vázlat-koordinátákban (bal felső sarok = 0,0):
      outline: [(pontlista, zárt?)] – fehér kontúr
      bend:    [(p1, p2)]          – piros szaggatott hajlítási vonal
      dims:    [(p1, p2, e1, e2)]  – zöld méretvonal p1–p2, segédvonal e1→p1, e2→p2
      texts:   [(x, y, szöveg, méret)]
      inputs:  {mező: (x_közép, y_közép)}
      calc:    (x, y) – a számolt értékek (ferde él, szög…) helye
      calc_front: (x_közép, y_közép) – a homloklap számolt magassága
    """
    o = {"outline": [], "bend": [], "dims": [], "texts": [], "inputs": {}}
    # --- oldallap (felül) ---
    x0, x1, xa = 150, 520, 330          # bal él, felső él vége, alsó él vége
    yt, yb, yn = 70, 260, 92            # felső él, alsó él, nettó magasság vonala (perem alja)
    o["outline"].append(([(x0, yt), (x1, yt), (x1, yn), (xa, yb), (x0, yb)], True))
    o["bend"] += [((x0, yt + 9), (x1 - 4, yt + 9)), ((x0, yn), (x1, yn))]
    o["texts"].append((x0 + 40, 170, "Oldallap", 13))
    o["dims"].append(((x0, 42), (x1, 42), (x0, yt), (x1, yt)))           # C
    o["inputs"]["side_C"] = ((x0 + x1) / 2, 24)
    o["dims"].append(((112, yn), (112, yb), (x0, yn), (x0, yb)))         # B (nettó)
    o["inputs"]["side_B"] = (58, (yn + yb) / 2)
    o["dims"].append(((x0, 290), (xa, 290), (x0, yb), (xa, yb)))         # A
    o["inputs"]["side_A"] = ((x0 + xa) / 2, 312)
    o["texts"].append((x1 + 12, yt + 14, "+ perem", 10))
    o["calc"] = (560, 120)
    # --- három téglalap (alul) ---
    yt2, yb2 = 420, 600
    blocks = [("Hátlap", 20, "hat_A", "hat_B", True),
              ("Homloklap", 275, None, "hom_B", True),
              ("Alja", 530, "alja_A", "alja_B", False)]
    for title, bx, fa, fb, flange in blocks:
        rx0, rx1 = bx + 95, bx + 235
        o["outline"].append(([(rx0, yt2), (rx1, yt2), (rx1, yb2), (rx0, yb2)], True))
        ynet = yt2 + 14 if flange else yt2
        if flange:
            o["bend"].append(((rx0, ynet), (rx1, ynet)))
        o["texts"].append((rx0 + 18, (yt2 + yb2) / 2 + 5, title, 13))
        o["dims"].append(((rx0, 392), (rx1, 392), (rx0, yt2), (rx1, yt2)))   # B (szélesség)
        o["inputs"][fb] = ((rx0 + rx1) / 2, 372)
        o["dims"].append(((rx0 - 22, ynet), (rx0 - 22, yb2), (rx0, ynet), (rx0, yb2)))  # A
        if fa:
            o["inputs"][fa] = (bx + 42, (ynet + yb2) / 2)
        else:
            o["calc_front"] = (bx + 42, (ynet + yb2) / 2)
            o["texts"].append((bx + 6, (ynet + yb2) / 2 - 20, "ferde él", 10))
    o["texts"].append((20, 636, "Minden magasság NETTÓ (perem nélkül) – a program hozzáadja a "
                                "peremet. A homloklap magassága a ferde él, a program számolja.", 10))
    return o


# ---------- Sittes rajzos méretmegadás – vázlat ----------

SITTES_W, SITTES_H = 960, 705

# mező → (táblázat oszlop vagy rejtett kulcs, felirat)
SITTES_FIELDS = {
    "side_A": (9,  "Oldal A – szélesség"),
    "side_B": (10, "Oldal B – magasság"),
    "hom_C":  (7,  "Homloklap C – szélesség"),
    "hom_B":  (8,  "Homloklap B – teljes magasság"),
    "hat_C":  (5,  "Hátlap C – szélesség"),
    "hat_B":  (6,  "Hátlap B – teljes magasság"),
    "hom_D":  ("hom_D", "Homloklap D – alsó él szélessége"),
    "hat_D":  ("hat_D", "Hátlap D – alsó él szélessége"),
    "alja_A": (4,  "Alja A – hossz"),
    "alja_B": (3,  "Alja B – szélesség"),
}
SITTES_OPTIONAL = ()
SITTES_LINKS = {"hom_C": ("hat_C", "hom_D"),
                "hom_D": ("hat_D", "alja_B"),
                "hom_B": ("hat_B",)}
# a táblázat lap-magasság oszlopa = B + A (teljes szabási magasság)
SITTES_DERIVED = {}     # a homlok- és hátlapnál nincs hajtás: B = teljes magasság


def sittes_sketch() -> dict:
    o = {"outline": [], "bend": [], "dims": [], "texts": [], "inputs": {},
         "ghost": [], "size": (SITTES_W, SITTES_H)}
    # --- oldal (bal felül): A felül, B jobb oldalon ---
    x0, x1, yt, yb = 120, 450, 80, 268
    o["ghost"].append(([(x0, yt), (x1, yt), (x1, yb), (x0, yb)], True))      # szabási téglalap
    o["outline"].append(([(x0, yt), (x1, yt), (x1, yt + 26), (x0, yt + 26)], True))   # felső keret
    o["outline"].append(([(x0 + 4, yt + 26), (x1 - 4, yt + 26), (x1 - 38, 210), (x1 - 88, yb),
                          (x0 + 58, yb)], True))
    o["texts"].append((x0 + 125, 190, "Oldal", 14))
    o["dims"].append(((x0, 55), (x1, 55), (x0, yt), (x1, yt)))
    o["inputs"]["side_A"] = ((x0 + x1) / 2, 33)
    o["dims"].append(((480, yt), (480, yb), (x1, yt), (x1 - 88, yb)))
    o["inputs"]["side_B"] = (527, 222)
    o["texts"].append((x0 + 6, yb + 18, "szabás: A × B téglalap", 9))
    # --- homloklap (jobb felül), hátlap (jobb alul): C szélesség, B teljes magasság (nincs hajtás) ---
    for title, dy, pre in (("Homloklap", 0, "hom"), ("Hátlap", 340, "hat")):
        x0, x1 = 640, 840
        yt, yb = 80 + dy, 268 + dy
        o["outline"].append(([(x0, yt), (x1, yt), (x1 - 18, yb), (x0 + 18, yb)], True))   # trapéz
        o["outline"].append(([(x0 - 6, yt + 6), (x1 + 6, yt + 6), (x1 + 6, yt + 14), (x0 - 6, yt + 14)], True))
        o["texts"].append((x0 + 45, yt + 85, title, 14))
        o["dims"].append(((x0, yt - 25), (x1, yt - 25), (x0, yt), (x1, yt)))         # C
        o["inputs"][pre + "_C"] = ((x0 + x1) / 2, yt - 47)
        o["dims"].append(((612, yt), (612, yb), (x0, yt), (x0 + 18, yb)))            # B
        o["inputs"][pre + "_B"] = (565, (yt + yb) / 2)
        o["dims"].append(((x0 + 18, yb + 22), (x1 - 18, yb + 22), (x0 + 18, yb), (x1 - 18, yb)))   # D
        o["inputs"][pre + "_D"] = ((x0 + x1) / 2, yb + 44)
    # --- alja (bal alul) ---
    x0, x1, yt, yb = 120, 450, 420, 608
    o["outline"].append(([(x0, yt), (x1, yt), (x1, yb), (x0, yb)], True))
    o["texts"].append((x0 + 140, 520, "Alja", 14))
    o["dims"].append(((x0, 395), (x1, 395), (x0, yt), (x1, yt)))
    o["inputs"]["alja_A"] = ((x0 + x1) / 2, 373)
    o["dims"].append(((95, yt), (95, yb), (x0, yt), (x0, yb)))
    o["inputs"]["alja_B"] = (48, (yt + yb) / 2)
    o["texts"].append((20, 694, "Homloklap/hátlap: trapéz, felül C, alul D széles (D = C → téglalap), "
                                "B = teljes magasság, nincs hajtás. A felső keret nem számít.", 9))
    return o


def sittes_calc(num, flange):
    lines, err = [], None
    for pre, title in (("hom", "Homloklap"), ("hat", "Hátlap")):
        c, b, d = num(pre + "_C"), num(pre + "_B"), num(pre + "_D")
        if c and b:
            shape = f" (trapéz, alul {d:g})" if d and abs(d - c) >= 0.5 else ""
            lines.append(f"{title}: szabás {c:g} × {b:g} mm{shape}")
    sa, sb = num("side_A"), num("side_B")
    if sa and sb:
        lines.append(f"Oldal: szabás {sa:g} × {sb:g} mm (2 db)")
    aa, ab = num("alja_A"), num("alja_B")
    if aa and ab:
        lines.append(f"Alja: szabás {aa:g} × {ab:g} mm")
    return ("\n".join(lines) or "Add meg a méreteket.").replace(".", ","), None, err


# ---------- Sittes rajzos méretmegadás – HAJTÁSSAL (kapcsolóval választható) ----------


# mező → (táblázat oszlop vagy rejtett kulcs, felirat)
SITTES_FOLD_FIELDS = {
    "side_A": (9,  "Oldal A – szélesség"),
    "side_B": (10, "Oldal B – magasság"),
    "hom_C":  (7,  "Homloklap C – szélesség"),
    "hom_B":  ("hom_B", "Homloklap B – magasság a hajlítási vonalig"),
    "hom_A":  ("hom_A", "Homloklap A – hajlítás alatti rész (nem kötelező)"),
    "hat_C":  (5,  "Hátlap C – szélesség"),
    "hat_B":  ("hat_B", "Hátlap B – magasság a hajlítási vonalig"),
    "hat_A":  ("hat_A", "Hátlap A – hajlítás alatti rész (nem kötelező)"),
    "hom_D":  ("hom_D", "Homloklap D – alsó él szélessége"),
    "hat_D":  ("hat_D", "Hátlap D – alsó él szélessége"),
    "alja_A": (4,  "Alja A – hossz"),
    "alja_B": (3,  "Alja B – szélesség"),
}
SITTES_FOLD_OPTIONAL = ("hom_A", "hat_A")
SITTES_FOLD_LINKS = {"hom_C": ("hat_C", "hom_D"),
                "hom_D": ("hat_D", "alja_B"),
                "hom_B": ("hat_B",),
                "hom_A": ("hat_A",)}
# a táblázat lap-magasság oszlopa = B + A (teljes szabási magasság)
SITTES_FOLD_DERIVED = {8: ("hom_B", "hom_A"), 6: ("hat_B", "hat_A")}


def sittes_fold_sketch() -> dict:
    o = {"outline": [], "bend": [], "dims": [], "texts": [], "inputs": {},
         "ghost": [], "size": (SITTES_W, SITTES_H)}
    # --- oldal (bal felül): A felül, B jobb oldalon ---
    x0, x1, yt, yb = 120, 450, 80, 268
    o["ghost"].append(([(x0, yt), (x1, yt), (x1, yb), (x0, yb)], True))      # szabási téglalap
    o["outline"].append(([(x0, yt), (x1, yt), (x1, yt + 26), (x0, yt + 26)], True))   # felső keret
    o["outline"].append(([(x0 + 4, yt + 26), (x1 - 4, yt + 26), (x1 - 38, 210), (x1 - 88, yb),
                          (x0 + 58, yb)], True))
    o["texts"].append((x0 + 125, 190, "Oldal", 14))
    o["dims"].append(((x0, 55), (x1, 55), (x0, yt), (x1, yt)))
    o["inputs"]["side_A"] = ((x0 + x1) / 2, 33)
    o["dims"].append(((480, yt), (480, yb), (x1, yt), (x1 - 88, yb)))
    o["inputs"]["side_B"] = (527, 222)
    o["texts"].append((x0 + 6, yb + 18, "szabás: A × B téglalap", 9))
    # --- homloklap (jobb felül), hátlap (jobb alul): C szélesség, B felső rész, A hajlítás alatti rész ---
    for title, dy, pre in (("Homloklap", 0, "hom"), ("Hátlap", 340, "hat")):
        x0, x1 = 640, 840
        yt, yb = 80 + dy, 268 + dy
        yd = yb - 60
        o["outline"].append(([(x0, yt), (x1, yt), (x1 - 18, yb), (x0 + 18, yb)], True))   # trapéz
        o["outline"].append(([(x0 - 6, yt + 6), (x1 + 6, yt + 6), (x1 + 6, yt + 14), (x0 - 6, yt + 14)], True))
        o["bend"].append(((x0 + 13, yd), (x1 - 13, yd)))
        o["dims"].append(((x0 + 18, yb + 22), (x1 - 18, yb + 22), (x0 + 18, yb), (x1 - 18, yb)))   # D
        o["inputs"][pre + "_D"] = ((x0 + x1) / 2, yb + 44)
        o["texts"].append((x0 + 45, yt + 85, title, 14))
        o["dims"].append(((x0, yt - 25), (x1, yt - 25), (x0, yt), (x1, yt)))         # C
        o["inputs"][pre + "_C"] = ((x0 + x1) / 2, yt - 47)
        o["dims"].append(((612, yt), (612, yd), (x0, yt), (x0, yd)))                 # B
        o["inputs"][pre + "_B"] = (565, yt + 38)
        o["dims"].append(((868, yd), (868, yb), (x1 - 13, yd), (x1 - 18, yb)))       # A
        o["inputs"][pre + "_A"] = (913, (yd + yb) / 2)
    # --- alja (bal alul) ---
    x0, x1, yt, yb = 120, 450, 420, 608
    o["outline"].append(([(x0, yt), (x1, yt), (x1, yb), (x0, yb)], True))
    o["texts"].append((x0 + 140, 520, "Alja", 14))
    o["dims"].append(((x0, 395), (x1, 395), (x0, yt), (x1, yt)))
    o["inputs"]["alja_A"] = ((x0 + x1) / 2, 373)
    o["dims"].append(((95, yt), (95, yb), (x0, yt), (x0, yb)))
    o["inputs"]["alja_B"] = (48, (yt + yb) / 2)
    o["texts"].append((20, 694, "Homloklap/hátlap: trapéz, felül C, alul D széles, magassága B + A; "
                                "a piros vonal a hajlítás. A felső keret nem számít.", 9))
    return o


def sittes_fold_calc(num, flange):
    lines, err = [], None
    for pre, title in (("hom", "Homloklap"), ("hat", "Hátlap")):
        c, b, a = num(pre + "_C"), num(pre + "_B"), num(pre + "_A") or 0
        d = num(pre + "_D")
        if c and b:
            extra = f" (felső rész {b:g} + hajlítás alatt {a:g})" if a else ""
            shape = f", trapéz, alul {d:g}" if d and abs(d - c) >= 0.5 else ""
            lines.append(f"{title}: szabás {c:g} × {b + a:g} mm{extra}{shape}")
    sa, sb = num("side_A"), num("side_B")
    if sa and sb:
        lines.append(f"Oldal: szabás {sa:g} × {sb:g} mm (2 db)")
    aa, ab = num("alja_A"), num("alja_B")
    if aa and ab:
        lines.append(f"Alja: szabás {aa:g} × {ab:g} mm")
    return ("\n".join(lines) or "Add meg a méreteket.").replace(".", ","), None, err


def abo_a_calc(num, flange):
    """(info szöveg, homloklap-cella szöveg, hibaüzenet OK-nál vagy None)"""
    a, b, c = num("side_A"), num("side_B"), num("side_C")
    if a and b and c and a < c:
        H, slant, ang = abo_a_info(c, b, a, flange)
        info = (f"Ferde él: {slant:.1f} mm\nSzög: {ang:.1f}°\n"
                f"Oldallap teljes magasság: {H:.0f} mm\n(nettó {b:g} + {flange:g} perem)\n"
                f"Homloklap: {slant:.1f} + {flange:g} = {slant + flange:.1f} mm").replace(".", ",")
        return info, f"{slant:.1f}".replace(".", ","), None
    if a and c and a >= c:
        msg = "Az oldallap alsó éle (A) legyen rövidebb a felső élénél (C)!"
        return msg, "–", msg
    return "Add meg az oldallap\nA, B, C méretét.", "–", None


# ---------- ABO_B rajzos méretmegadás – vázlat ----------

ABO_B_W, ABO_B_H = 900, 600

ABO_B_FIELDS = {
    "side_A": (9,  "Oldal A – teljes hossz (az orr csúcsáig)"),
    "side_C": (11, "Oldal C – felső egyenes él hossza"),
    "side_E": ("side_E", "Oldal – alsó egyenes él hossza"),
    "side_B": (10, "Oldal B – NETTÓ magasság (perem nélkül)"),
    "side_D": ("side_D", "Oldal D – az orr csúcsának magassága alulról (üresen = középen)"),
    "hom_A":  (7,  "Homloklap A – szélesség"),
    "hom_B":  (8,  "Homloklap B – NETTÓ magasság (perem nélkül)"),
    "hat_A":  (5,  "Hátlap A – szélesség"),
    "hat_B":  (6,  "Hátlap B – NETTÓ magasság (perem nélkül)"),
    "alja_A": (3,  "Alja A – szélesség"),
    "alja_B": (4,  "Alja B – hossz"),
}
ABO_B_LINKS = {"hat_A": ("hom_A", "alja_A"),
               "side_B": ("hat_B",)}


def abo_b_sketch() -> dict:
    o = {"outline": [], "bend": [], "dims": [], "texts": [], "inputs": {},
         "ghost": [], "size": (ABO_B_W, ABO_B_H)}
    # --- oldal (bal felül): orros alak, felül és a felső ferde élen perem ---
    x0, x1, xe, xt, yt, yb = 130, 330, 290, 420, 100, 240
    yn = yt + 12
    ty = (yt + yb) / 2
    o["outline"].append(([(x0, yt), (x1, yt), (xt, ty), (xe, yb), (x0, yb)], True))
    o["bend"].append(((x0, yn), (x1 - 5, yn)))
    o["bend"].append(((x1 - 5, yn), (xt - 13, ty)))
    o["texts"].append((x0 + 18, 185, "Oldalai", 13))
    o["dims"].append(((x0, 36), (xt, 36), (x0, yt), (xt, ty)))                       # A
    o["inputs"]["side_A"] = ((x0 + xt) / 2 + 40, 18)
    o["dims"].append(((x0, 74), (x1, 74), (x0, yt), (x1, yt)))                       # C
    o["inputs"]["side_C"] = ((x0 + x1) / 2, 56)
    o["dims"].append(((100, yn), (100, yb), (x0, yn), (x0, yb)))                     # B
    o["inputs"]["side_B"] = (52, (yn + yb) / 2)
    o["dims"].append(((xt + 22, ty), (xt + 22, yb), (xt, ty), (xe, yb)))             # D
    o["inputs"]["side_D"] = (xt + 70, (ty + yb) / 2)
    o["dims"].append(((x0, yb + 24), (xe, yb + 24), (x0, yb), (xe, yb)))             # alsó él
    o["inputs"]["side_E"] = ((x0 + xe) / 2, yb + 46)
    o["texts"].append((x0 + 18, 207, "a valós (orros) alakkal szabja", 9))
    # --- homloklap (jobb felül): A alul, B jobbra ---
    x0, x1, yt, yb = 560, 760, 60, 240
    yn = yt + 18
    o["outline"].append(([(x0, yt), (x1, yt), (x1, yb), (x0, yb)], True))
    o["bend"].append(((x0, yn), (x1, yn)))
    o["texts"].append((x0 + 40, 150, "Homloklap", 13))
    o["dims"].append(((790, yn), (790, yb), (x1, yn), (x1, yb)))
    o["inputs"]["hom_B"] = (835, (yn + yb) / 2)
    o["dims"].append(((x0, 265), (x1, 265), (x0, yb), (x1, yb)))
    o["inputs"]["hom_A"] = ((x0 + x1) / 2, 287)
    # --- hátlap (bal alul): A felül, B balra ---
    x0, x1, yt, yb = 150, 350, 380, 530
    yn = yt + 14
    o["outline"].append(([(x0, yt), (x1, yt), (x1, yb), (x0, yb)], True))
    o["bend"].append(((x0, yn), (x1, yn)))
    o["texts"].append((x0 + 50, 470, "Hátlap", 13))
    o["dims"].append(((x0, 355), (x1, 355), (x0, yt), (x1, yt)))
    o["inputs"]["hat_A"] = ((x0 + x1) / 2, 333)
    o["dims"].append(((120, yn), (120, yb), (x0, yn), (x0, yb)))
    o["inputs"]["hat_B"] = (66, (yn + yb) / 2)
    # --- alja (jobb alul): A felül, B jobbra ---
    x0, x1, yt, yb = 560, 760, 380, 530
    o["outline"].append(([(x0, yt), (x1, yt), (x1, yb), (x0, yb)], True))
    o["texts"].append((x0 + 70, 465, "Alja", 13))
    o["dims"].append(((x0, 355), (x1, 355), (x0, yt), (x1, yt)))
    o["inputs"]["alja_A"] = ((x0 + x1) / 2, 333)
    o["dims"].append(((790, yt), (790, yb), (x1, yt), (x1, yb)))
    o["inputs"]["alja_B"] = (835, (yt + yb) / 2)
    o["texts"].append((20, 588, "A magasságok NETTÓ értékek – a program hozzáadja a peremet (30+25). "
                                "A piros vonal a perem hajlítása.", 9))
    return o


def abo_b_calc(num, flange):
    lines, err = [], None
    a, b, c, d, e = num("side_A"), num("side_B"), num("side_C"), num("side_D"), num("side_E")
    if a and b and c and e:
        H = b + flange
        ty = d or H / 2
        if c >= a:
            err = "Az oldal felső egyenes éle (C) legyen rövidebb a teljes hossznál (A)."
        elif e >= a:
            err = "Az oldal alsó egyenes éle legyen rövidebb a teljes hossznál (A)."
        elif not (0 < ty < H):
            err = "Az orr csúcsa (D) a lap magasságán belül legyen."
        else:
            up, lo = math.hypot(a - c, H - ty), math.hypot(a - e, ty)
            area = poly_area(abo_b_side_poly(a, H, c, ty, e))
            lines.append(f"Oldal: {a:g} × {H:g} mm (nettó {b:g} + {flange:g} perem), felső él {c:g}, "
                         f"alsó él {e:g} mm, orr csúcsa {ty:g} mm-en, "
                         f"ferde élek {up:.1f} / {lo:.1f} mm, {area / 1e6:.3f} m² (2 db)")
    for key, title, cnt in (("hat", "Hátlap", 1), ("hom", "Homloklap", 1)):
        a, b = num(key + "_A"), num(key + "_B")
        if a and b:
            db = f" ({cnt} db)" if cnt > 1 else ""
            lines.append(f"{title}: szabás {a:g} × {b + flange:g} mm (nettó {b:g} + {flange:g} perem){db}")
    a, b = num("alja_A"), num("alja_B")
    if a and b:
        lines.append(f"Alja: szabás {a:g} × {b:g} mm")
    if err:
        lines.append("⚠ " + err)
    return ("\n".join(lines) or "Add meg a méreteket.").replace(".", ","), None, err


SKETCH_CONFIGS = {
    "ABO_A": dict(title="ABO_A méretek megadása (mm)", sketch=lambda: dict(abo_a_sketch(), size=(SKETCH_W, SKETCH_H)),
                  fields=ABO_A_FIELDS, optional=(), links={WIDTH_SOURCE: WIDTH_TARGETS},
                  calc=abo_a_calc),
    "ABO_B": dict(title="ABO_B méretek megadása (mm)", sketch=abo_b_sketch,
                  fields=ABO_B_FIELDS, optional=("side_D",), links=ABO_B_LINKS, calc=abo_b_calc),
    "Sittes": dict(title="Sittes konténer méretei (mm)", sketch=sittes_sketch,
                   fields=SITTES_FIELDS, optional=SITTES_OPTIONAL, links=SITTES_LINKS,
                   calc=sittes_calc, derived=SITTES_DERIVED, toggle=("Sittes_fold", False)),
    "Sittes_fold": dict(title="Sittes konténer méretei (mm) – hajtással", sketch=sittes_fold_sketch,
                        fields=SITTES_FOLD_FIELDS, optional=SITTES_FOLD_OPTIONAL,
                        links=SITTES_FOLD_LINKS, calc=sittes_fold_calc,
                        derived=SITTES_FOLD_DERIVED, toggle=("Sittes", True)),
}


# ---------- GUI ----------

PRODUCT_COLUMNS = ["Termék neve", "Típus", "Darab",
                   "Alja szél.", "Alja hossz",
                   "Hátlap szél.", "Hátlap mag.",
                   "Homlok szél.", "Homlok mag.",
                   "Oldal A", "Oldal B", "Oldal C", "Rajz", "Vastagság (mm)"]
COL_TYPE, COL_QTY, COL_FIRST_DIM, COL_SIDE = 1, 2, 3, 9
COL_DRAW = 12         # ABO_A: rajzos méretmegadás gombja
COL_THICK = 13        # lemezvastagság (mm), termékenként
COL_HOMLOK_H = 8      # ABO_A-nál számolt érték (= ferde él hossza), nem kell megadni
MERGE_MODES = [
    ("Automatikus (a legjobbat választja)", ("", "hatlap", "homlok", "mind")),
    ("Nincs összevonás", ("",)),
    ("Alja + Hátlap egy lemezből", ("hatlap",)),
    ("Alja + Homloklap egy lemezből", ("homlok",)),
    ("Hátlap + Alja + Homloklap (U)", ("mind",)),
]
PROGRAM_VERSION = "v2026.10.30 (egyéb alkatrész DXF-ből)"

SIDE_HELP = {
    OTHER_TYPE: ("– (a DXF rajzból)", "– (a DXF rajzból)", "– (a DXF rajzból)"),
    "Sittes": ("Oldal A = szélesség (hossz)", "Oldal B = magasság", "Oldal C = – (nem kell)"),
    "ABO_A": ("Oldal A = felső él hossza (pl. 1500)", "Oldal B = NETTÓ magasság, perem nélkül (pl. 775)",
              "Oldal C = alsó él hossza (pl. 500); a homloklap nettó magassága = ferde él (+55 perem)"),
    "ABO_B": ("Oldal A = teljes hossz (az orr csúcsáig)", "Oldal B = NETTÓ magasság (perem nélkül)",
              "Oldal C = felső egyenes él hossza; az alsó egyenes él a „Rajz...” ablakban adható meg"),
}


def run_gui():
    from PyQt5.QtWidgets import (
        QApplication, QWidget, QGridLayout, QLabel, QLineEdit, QCheckBox,
        QPushButton, QTextEdit, QFileDialog, QMessageBox, QTableWidget,
        QTableWidgetItem, QHeaderView, QHBoxLayout, QComboBox, QDialog,
        QVBoxLayout, QDialogButtonBox, QInputDialog
    )
    from PyQt5.QtGui import QPainter, QPen, QColor, QPolygonF, QBrush
    from PyQt5.QtCore import QPointF

    class SketchCanvas(QWidget):
        """Vázlat kirajzolása CAD-szerű sötét háttéren."""
        def __init__(self, sk, parent=None):
            super().__init__(parent)
            self.sk = sk
            self.setFixedSize(*sk["size"])

        def paintEvent(self, _ev):
            p = QPainter(self)
            p.setRenderHint(QPainter.Antialiasing)
            p.fillRect(self.rect(), QColor("#222831"))
            white, green, red, grey = (QColor("#e8e8e8"), QColor("#39d353"),
                                       QColor("#ff4d4d"), QColor("#6b7380"))
            pen = QPen(grey, 1.0)
            pen.setStyle(Qt.DashLine)
            p.setPen(pen)
            for pts, _closed in self.sk.get("ghost", []):
                p.drawPolygon(QPolygonF([QPointF(x, y) for x, y in pts]))
            p.setPen(QPen(white, 1.6))
            for pts, closed in self.sk["outline"]:
                poly = QPolygonF([QPointF(x, y) for x, y in pts])
                if closed:
                    p.drawPolygon(poly)
                else:
                    p.drawPolyline(poly)
            pen = QPen(red, 1.4)
            pen.setStyle(Qt.DashDotLine)
            p.setPen(pen)
            for (a, b) in self.sk["bend"]:
                p.drawLine(QPointF(*a), QPointF(*b))
            p.setPen(QPen(green, 1.2))
            p.setBrush(QBrush(green))
            for p1, p2, e1, e2 in self.sk["dims"]:
                p.drawLine(QPointF(*e1), QPointF(*p1))
                p.drawLine(QPointF(*e2), QPointF(*p2))
                p.drawLine(QPointF(*p1), QPointF(*p2))
                for tip, tail in ((p1, p2), (p2, p1)):
                    dx, dy = tail[0] - tip[0], tail[1] - tip[1]
                    ln = math.hypot(dx, dy) or 1
                    ux, uy = dx / ln, dy / ln
                    a1 = QPointF(tip[0] + 9 * ux - 3.5 * uy, tip[1] + 9 * uy + 3.5 * ux)
                    a2 = QPointF(tip[0] + 9 * ux + 3.5 * uy, tip[1] + 9 * uy - 3.5 * ux)
                    p.drawPolygon(QPolygonF([QPointF(*tip), a1, a2]))
            p.setBrush(Qt.NoBrush)
            p.setPen(QPen(white))
            for x, y, txt, size in self.sk["texts"]:
                f = p.font()
                f.setPointSize(size)
                p.setFont(f)
                p.drawText(QPointF(x, y), txt)
            p.end()

    class SketchDialog(QDialog):
        """Méretek megadása a műhelyrajz szerinti vázlaton (ABO_A, Sittes)."""
        def __init__(self, ctype: str, values: Dict[str, str], flange: float, parent=None):
            super().__init__(parent)
            cfg = SKETCH_CONFIGS[ctype]
            self.cfg = cfg
            self.setWindowTitle(cfg["title"])
            self.flange = flange
            lay = QVBoxLayout(self)
            if "toggle" in cfg:
                # Sittes: hajtás nélkül (teljes magasság) vagy hajtással (B + A)
                self.chk_fold = QCheckBox("Hajtás megadása a homlok- és hátlapnál "
                                          "(B a hajlítási vonalig + A a hajtás)")
                self.chk_fold.setChecked(cfg["toggle"][1])
                self.chk_fold.toggled.connect(lambda _on: self.done(2))
                lay.addWidget(self.chk_fold)
            self.canvas = SketchCanvas(cfg["sketch"](), self)
            lay.addWidget(self.canvas)
            sk = self.canvas.sk
            style = ("QLineEdit{background:#11151a;color:#ffffff;border:1px solid #39d353;"
                     "border-radius:3px;padding:1px 4px;font-weight:bold}")
            self.style = style
            self.style_auto = style.replace("#39d353", "#5b8f66").replace(
                "font-weight:bold", "font-weight:bold;font-style:italic;color:#bfe9c6")
            self.style_opt = style.replace("#39d353", "#8a8f98")
            self.edits = {}
            for name, (cx, cy) in sk["inputs"].items():
                e = QLineEdit(values.get(name, ""), self.canvas)
                e.setGeometry(int(cx - INPUT_W / 2), int(cy - INPUT_H / 2), INPUT_W, INPUT_H)
                e.setAlignment(Qt.AlignCenter)
                e.setToolTip(cfg["fields"][name][1])
                e.setStyleSheet(self.style_opt if name in cfg["optional"] else style)
                e.textChanged.connect(lambda _t, n=name: self.on_change(n))
                self.edits[name] = e
            # összekapcsolt mezők: a forrás értéke átkerül, amíg kézzel át nem írják
            self.src_of = {t: src for src, ts in cfg["links"].items() for t in ts}
            self.linked = {}
            for t, src in self.src_of.items():
                txt = self.edits[t].text().strip()
                self.linked[t] = (not txt) or txt == self.edits[src].text().strip()
                self._style_linked(t)
            # üres, összekapcsolt mező: a forrás értékével indul (a lánc mentén is)
            for _ in range(3):
                for t, src in self.src_of.items():
                    if self.linked[t] and not self.edits[t].text().strip():
                        self.edits[t].blockSignals(True)
                        self.edits[t].setText(self.edits[src].text())
                        self.edits[t].blockSignals(False)
            self.front = None
            if "calc_front" in sk:
                cx, cy = sk["calc_front"]
                self.front = QLabel("–", self.canvas)
                self.front.setGeometry(int(cx - INPUT_W / 2), int(cy - INPUT_H / 2), INPUT_W, INPUT_H)
                self.front.setAlignment(Qt.AlignCenter)
                self.front.setStyleSheet("QLabel{color:#ffd166;border:1px dashed #ffd166;border-radius:3px}")
                self.front.setToolTip("Homloklap nettó magassága = a ferde él hossza (számolt)")
            if "calc" in sk:
                x, y = sk["calc"]
                self.info = QLabel("", self.canvas)
                self.info.setGeometry(int(x), int(y), sk["size"][0] - int(x) - 10, 150)
            else:
                self.info = QLabel("")
                lay.addWidget(self.info)
            self.info.setStyleSheet("QLabel{color:#ffd166;background:#222831;font-size:10pt;padding:4px}")
            self.info.setAlignment(Qt.AlignLeft | Qt.AlignTop)
            self.info.setWordWrap(True)
            bb = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
            bb.accepted.connect(self.try_accept)
            bb.rejected.connect(self.reject)
            lay.addWidget(bb)
            self.recalc()

        def num(self, name) -> Optional[float]:
            try:
                v = float(self.edits[name].text().strip().replace(",", "."))
                return v if v > 0 else None
            except (ValueError, KeyError):
                return None

        def _style_linked(self, t):
            e = self.edits[t]
            base = self.style_opt if t in self.cfg["optional"] else self.style
            e.setStyleSheet(self.style_auto if self.linked[t] else base)
            e.setToolTip(self.cfg["fields"][t][1] + (" – automatikusan átvett érték, átírható"
                                                     if self.linked[t] else ""))

        def on_change(self, name):
            for t in self.cfg["links"].get(name, ()):
                if self.linked.get(t):
                    # a jel nincs letiltva: a láncolt mezők is követik (pl. C → D → alja)
                    self.edits[t].setText(self.edits[name].text())
            if name in self.src_of:
                # kézi átírás: leválik; ha üresre vagy a forrással egyezőre írják, újra követi
                txt = self.edits[name].text().strip()
                self.linked[name] = (not txt) or txt == self.edits[self.src_of[name]].text().strip()
                self._style_linked(name)
            self.recalc()

        def recalc(self):
            info, front, _err = self.cfg["calc"](self.num, self.flange)
            self.info.setText(info)
            if self.front is not None and front is not None:
                self.front.setText(front)

        def try_accept(self):
            missing = [self.cfg["fields"][n][1] for n in self.edits
                       if n not in self.cfg["optional"] and self.num(n) is None]
            bad_opt = [self.cfg["fields"][n][1] for n in self.edits
                       if n in self.cfg["optional"] and self.edits[n].text().strip()
                       and self.num(n) is None]
            if missing or bad_opt:
                QMessageBox.warning(self, "Hiányzó vagy hibás méret",
                                    "Ellenőrizd ezeket:\n• " + "\n• ".join(missing + bad_opt))
                return
            _info, _front, err = self.cfg["calc"](self.num, self.flange)
            if err:
                QMessageBox.warning(self, "Hibás méret", err)
                return
            self.accept()

        def values(self) -> Dict[str, str]:
            return {n: e.text().strip() for n, e in self.edits.items()}
    from PyQt5.QtGui import QFont
    from PyQt5.QtCore import Qt, QThread, pyqtSignal, QTimer

    class CalcThread(QThread):
        """A számítás külön szálon fut, így az ablak közben nem akad meg."""
        progress = pyqtSignal(str)
        done = pyqtSignal(object)
        failed = pyqtSignal(str)

        def __init__(self, owner, args, kwargs):
            super().__init__(owner)
            self.owner, self.args, self.kwargs = owner, args, kwargs

        def run(self):
            o = self.owner
            try:
                if o.pool is None and o.workers > 1:
                    self.progress.emit("Párhuzamos számítás indítása…")
                    try:
                        o.pool = multiprocessing.get_context("spawn").Pool(o.workers)
                    except Exception:
                        o.pool, o.workers = None, 1
                res = optimize_by_thickness(*self.args, progress=self.progress.emit,
                                             pool=o.pool, workers=o.workers, **self.kwargs)
                # egyéb (DXF) alkatrészek: táblaméretenkénti összehasonlítás
                prods0 = [vl[0] for vl in self.args[0]]
                if any(p.ctype == OTHER_TYPE for p in prods0):
                    self.progress.emit("Egyéb alkatrészek: táblaméretek összevetése…")
                    set_nesting(True)
                    rep = other_sheet_report(prods0, self.args[1], self.args[2], self.args[3])
                    res = tuple(res) + (rep,)
                self.done.emit(res)
            except Cancelled:
                if o.pool is not None:
                    o.pool.terminate()
                    o.pool = None
                self.failed.emit("")
            except Exception as e:
                self.failed.emit(str(e) or e.__class__.__name__)

    class ContainerApp(QWidget):
        def __init__(self):
            super().__init__()
            self.setWindowTitle(f"Konténer lemezszabás kalkulátor – {PROGRAM_VERSION}")
            self.resize(1100, 800)
            self.sheets: List[SheetSize] = []
            self.material_path = ""
            self.last = None
            self.pool = None
            self.workers = worker_count()
            self.thread = None
            self.calc_t0 = 0.0
            self.calc_limit = 0.0
            self.calc_msg = ""
            self.tick = QTimer(self)
            self.tick.setInterval(500)
            self.tick.timeout.connect(self.show_progress)
            self.last_products = []
            self.material_full_path = ""

            lay = QGridLayout(self)
            row = 0

            # --- termékek táblázata ---
            lay.addWidget(QLabel("Gyártandó termékek (méretek mm-ben):"), row, 0, 1, 3)
            row += 1
            self.tbl = QTableWidget(0, len(PRODUCT_COLUMNS))
            self.extras: Dict[int, Dict[str, str]] = {}
            self._updating = False
            self.tbl.itemChanged.connect(self.on_item_changed)
            self.tbl.setHorizontalHeaderLabels(PRODUCT_COLUMNS)
            self.tbl.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
            self.tbl.setMinimumHeight(140)
            lay.addWidget(self.tbl, row, 0, 1, 3)
            row += 1
            help_txt = "<b>Oldalméretek típusonként:</b><br>" + "<br>".join(
                f"<b>{t}:</b> " + "; ".join(v) for t, v in SIDE_HELP.items())
            help_txt += ("<br><i>ABO_A: a két trapéz alakú oldalt a program 180°-ban egymásba "
                         "fordítva is megpróbálja kiszabni (közös ferde vágás). "
                         "ABO típusoknál a hátlap, a homloklap és az oldalak felső élére a program "
                         "automatikusan rászámolja a perem kiterített szélességét (alapból 30+25 mm) – "
                         "a magasságokat kész méretben add meg.</i>")
            lbl_help = QLabel(help_txt)
            lbl_help.setWordWrap(True)
            lay.addWidget(lbl_help, row, 0, 1, 3)
            row += 1
            btns = QHBoxLayout()
            b_add = QPushButton("+ Termék hozzáadása")
            b_add.clicked.connect(lambda: self.add_product())
            b_dxf = QPushButton("+ Egyéb alkatrész (DXF)…")
            b_dxf.setToolTip("Megrajzolt alkatrész betöltése DXF-ből: a program kiszámolja, "
                             "melyik táblaméretből szabható ki a legjobban")
            b_dxf.clicked.connect(self.add_other_dxf)
            b_del = QPushButton("Kijelölt termék törlése")
            b_del.clicked.connect(self.delete_product)
            btns.addWidget(b_add)
            btns.addWidget(b_dxf)
            btns.addWidget(b_del)
            btns.addStretch()
            lay.addLayout(btns, row, 0, 1, 3)
            row += 1
            self.add_product()

            self.field_labels = {}

            def field_(label, default=""):
                nonlocal row
                lbl = QLabel(label)
                lay.addWidget(lbl, row, 0)
                e = QLineEdit(default)
                self.field_labels[id(e)] = lbl
                lay.addWidget(e, row, 1)
                row += 1
                return e

            self.txt_kerf = field_("Vágási rés (mm):", "0")
            self.txt_min_pct = field_("Toldott rész min. területe (% a legkisebb alkatrészből):", "40")
            self.txt_time = field_("Keresési idő (másodperc, ennyi ideig számol; több = jobb eredmény):", "30")
            self.txt_flange = field_("ABO perem (mm, két szár, pl. 30x25):", "30x25")
            self.txt_flange.textChanged.connect(lambda _t: self.refresh_all_rows())
            self.txt_margin = field_("Széltávolság a tábla szélén (mm, minden oldalon):", str(DEFAULT_MARGIN))
            self.update_flange_enabled()

            self.chk_guillotine = QCheckBox("Csak egyenes végigvágás (táblaolló)")
            lay.addWidget(self.chk_guillotine, row, 0, 1, 2)
            row += 1

            self.chk_lasercut = QCheckBox(f"Lézeres vágás (a táblán legfeljebb {DEFAULT_MAX_CUT} mm "
                                          f"hosszan vág; kikapcsolva a teljes, 3000 mm-es hossz használható)")
            lay.addWidget(self.chk_lasercut, row, 0, 1, 3)
            row += 1

            self.chk_weld = QCheckBox("Hegesztett toldás engedélyezése (maradékok felhasználása)")
            self.chk_weld.setChecked(True)
            lay.addWidget(self.chk_weld, row, 0, 1, 2)
            row += 1

            self.chk_double = QCheckBox("Összeállítási varrat kétoldali (belül + kívül)")
            lay.addWidget(self.chk_double, row, 0, 1, 2)
            row += 1

            lay.addWidget(QLabel("Alja összehajlítása (kevesebb varrat):"), row, 0)
            self.cmb_merge = QComboBox()
            self.cmb_merge.addItems([m[0] for m in MERGE_MODES])
            lay.addWidget(self.cmb_merge, row, 1)
            row += 1

            self.lbl_sheets = QLabel()
            lay.addWidget(self.lbl_sheets, row, 0, 1, 2)
            btn_load = QPushButton("Anyagfájl betöltése…")
            btn_load.clicked.connect(self.choose_material_file)
            lay.addWidget(btn_load, row, 2)
            row += 1

            self.btn_calc = QPushButton("Optimalizált számítás")
            self.btn_calc.clicked.connect(self.calculate)
            lay.addWidget(self.btn_calc, row, 0)
            self.btn_cancel = QPushButton("Megszakítás")
            self.btn_cancel.setEnabled(False)
            self.btn_cancel.clicked.connect(self.cancel_calc)
            lay.addWidget(self.btn_cancel, row, 1)
            self.btn_dxf = QPushButton("DXF export")
            self.btn_dxf.setEnabled(False)
            self.btn_dxf.clicked.connect(self.export)
            lay.addWidget(self.btn_dxf, row, 2)
            row += 1

            self.btn_quote = QPushButton("Export árajánlathoz")
            self.btn_quote.setEnabled(False)
            self.btn_quote.setToolTip("A számítás eredményét JSON fájlba menti az Acél árkalkulátornak")
            self.btn_quote.clicked.connect(self.export_quote)
            lay.addWidget(self.btn_quote, row, 2)
            row += 1

            self.chk_laser = QCheckBox("Export lézerre: minden alkatrész külön DXF fájlba "
                                       "(azonos alkatrészek egy fájlban, darabszámmal)")
            lay.addWidget(self.chk_laser, row, 0, 1, 3)
            row += 1

            self.txt_result = QTextEdit()
            self.txt_result.setReadOnly(True)
            self.txt_result.setFont(QFont("Consolas", 9))
            lay.addWidget(self.txt_result, row, 0, 1, 3)
            lay.setRowStretch(row, 1)

            here = os.path.dirname(os.path.abspath(sys.argv[0]))
            found = None
            for d in (here, os.getcwd()):
                for name in ("anyagolás.txt", "anyagolas.txt", "Anyagolás.txt"):
                    p = os.path.join(d, name)
                    if os.path.exists(p):
                        found = p
                        break
                if found:
                    break
            if found:
                self.load_material(found)
            else:
                self.lbl_sheets.setText("Nincs betöltött anyagfájl – töltsd be az anyagolás.txt-t!")

        # --- termékek ---
        def add_product(self, values=None):
            r = self.tbl.rowCount()
            self.tbl.insertRow(r)
            if values is None:
                values = [f"Termék {r + 1}", "Sittes", "1"] + [""] * (len(PRODUCT_COLUMNS) - 3)
            values = list(values) + [""] * (len(PRODUCT_COLUMNS) - len(values))
            if not str(values[COL_THICK]).strip():
                values[COL_THICK] = DEFAULT_THICKNESS
            for c, v in enumerate(values):
                if c == COL_TYPE:
                    cb = QComboBox()
                    cb.addItems(CONTAINER_TYPES)
                    cb.setCurrentText(str(v) if str(v) in CONTAINER_TYPES else "Sittes")
                    cb.currentTextChanged.connect(lambda _t, cb=cb: self.update_row_type(cb))
                    self.tbl.setCellWidget(r, c, cb)
                elif c != COL_DRAW:
                    self.tbl.setItem(r, c, QTableWidgetItem(str(v)))
            btn = QPushButton("Rajz…")
            btn.setToolTip("Méretek megadása a rajz szerint (ABO_A, Sittes)")
            btn.clicked.connect(lambda _c=False, b=btn: self.open_sketch_for(b))
            self.tbl.setCellWidget(r, COL_DRAW, btn)
            self.update_row(r)

        def row_of(self, widget):
            for r in range(self.tbl.rowCount()):
                if self.tbl.cellWidget(r, COL_TYPE) is widget:
                    return r
            return -1

        def update_row_type(self, cb):
            r = self.row_of(cb)
            if r >= 0:
                self.update_row(r)
                if self.row_type(r) == OTHER_TYPE:
                    if not self.row_extras(r).get("dxf"):
                        self.load_other_dxf(r)
                    return
                # ABO_A-ra váltáskor, ha még nincsenek oldalméretek, rögtön nyílik a rajz
                if self.row_type(r) in ("ABO_A", "ABO_B") and not any(
                        (self.tbl.item(r, c).text().strip() if self.tbl.item(r, c) else "")
                        for c in range(COL_SIDE, COL_SIDE + 3)):
                    self.open_sketch(r)

        def open_sketch_for(self, btn):
            for r in range(self.tbl.rowCount()):
                if self.tbl.cellWidget(r, COL_DRAW) is btn:
                    self.open_sketch(r)
                    return

        def row_extras(self, r) -> Dict[str, str]:
            """A táblázatban nem látszó méretek (pl. sittes C, D) – a sor Rajz gombjához kötve."""
            btn = self.tbl.cellWidget(r, COL_DRAW)
            return self.extras.setdefault(id(btn), {}) if btn is not None else {}

        def add_other_dxf(self):
            r = self.tbl.rowCount()
            self.add_product([f"Egyéb {r + 1}", OTHER_TYPE, "1"])
            if not self.load_other_dxf(r):
                self.tbl.removeRow(r)

        def load_other_dxf(self, r) -> bool:
            """Egyéb alkatrész: DXF kiválasztása és beolvasása (ellenőrzéssel)."""
            extras = self.row_extras(r)
            start = os.path.dirname(extras.get("dxf", "")) if extras.get("dxf") else ""
            path, _ = QFileDialog.getOpenFileName(self, "Alkatrész rajza (DXF)", start,
                                                  "DXF rajz (*.dxf)")
            if not path:
                return False
            it = self.tbl.item(r, 0)
            name = it.text().strip() if it else ""
            if not name or re.fullmatch(r"(Termék|Egyéb) \d+", name):
                name = os.path.splitext(os.path.basename(path))[0]
            try:
                parts = load_dxf_parts(path, name)
            except Exception as e:
                QMessageBox.warning(self, "DXF betöltése", str(e))
                return False
            extras["dxf"] = path
            if it is None:
                self.tbl.setItem(r, 0, QTableWidgetItem(name))
            else:
                it.setText(name)
            self.update_row(r)
            QMessageBox.information(self, "DXF betöltve",
                                    f"{os.path.basename(path)}\n\n{dxf_parts_summary(parts)}\n\n"
                                    "A darabszámot és a lemezvastagságot a táblázatban add meg.")
            return True

        def open_sketch(self, r):
            ctype = self.row_type(r)
            if ctype == OTHER_TYPE:
                self.load_other_dxf(r)
                return
            if ctype not in SKETCH_CONFIGS:
                return
            extras = self.row_extras(r)
            key = "Sittes_fold" if ctype == "Sittes" and extras.get("sittes_fold") == "1" else ctype
            carry = None
            while True:
                cfg = SKETCH_CONFIGS[key]
                fields = cfg["fields"]
                vals = {}
                for name, (target, _) in fields.items():
                    if isinstance(target, str):
                        vals[name] = extras.get(target, "")
                    else:
                        it = self.tbl.item(r, target)
                        vals[name] = it.text().strip() if it else ""
                # ha a rejtett felbontás (pl. B + A) még nincs meg, a táblázat teljes értékéből indul
                for col, (fb, fa) in cfg.get("derived", {}).items():
                    if not vals.get(fb) and not vals.get(fa):
                        it = self.tbl.item(r, col)
                        vals[fb] = it.text().strip() if it else ""
                if carry is not None:
                    vals.update({k_: v_ for k_, v_ in carry.items() if k_ in fields})
                dlg = SketchDialog(key, vals, self.current_flange(), self)
                res = dlg.exec_()
                if res == 2:
                    # hajtás kapcsoló: a beírt értékek átkerülnek a másik rajzra
                    carry = dict(dlg.values())
                    other = cfg["toggle"][0]
                    for pre in ("hom", "hat"):
                        b_, a_ = carry.get(pre + "_B", ""), carry.get(pre + "_A", "")
                        if other == "Sittes":          # hajtásról teljes magasságra: B + A
                            try:
                                tot = float(b_.replace(",", ".")) + (float(a_.replace(",", ".")) if a_ else 0)
                                carry[pre + "_B"] = f"{tot:g}"
                            except ValueError:
                                pass
                        else:                           # teljes magasságról hajtásra: B = teljes, A üres
                            carry[pre + "_A"] = ""
                    key = other
                    continue
                if res != QDialog.Accepted:
                    return
                break
            self._updating = True
            try:
                for name, v in dlg.values().items():
                    target = fields[name][0]
                    if isinstance(target, str):
                        extras[target] = v
                        continue
                    it = self.tbl.item(r, target)
                    if it is None:
                        self.tbl.setItem(r, target, QTableWidgetItem(v))
                    else:
                        it.setText(v)
                vals2 = dlg.values()
                for col, (fb, fa) in cfg.get("derived", {}).items():
                    try:
                        tot = float(vals2[fb].replace(",", ".")) + \
                            (float(vals2[fa].replace(",", ".")) if vals2[fa] else 0.0)
                        txt = f"{tot:g}"
                    except (ValueError, KeyError):
                        continue
                    it = self.tbl.item(r, col)
                    if it is None:
                        self.tbl.setItem(r, col, QTableWidgetItem(txt))
                    else:
                        it.setText(txt)
                if ctype == "Sittes":
                    self.ask_alja_thick(r, extras)
                    extras["sittes_fold"] = "1" if key == "Sittes_fold" else ""
                    if key == "Sittes":
                        # a hajtásos felbontás nem érvényes többé (a teljes magasság a táblázatban van)
                        for k_ in ("hom_A", "hat_A", "hom_B", "hat_B"):
                            extras[k_] = ""
            finally:
                self._updating = False
            self.fill_abo_a_front(r)

        def row_thick(self, r) -> float:
            it = self.tbl.item(r, COL_THICK)
            try:
                txt = (it.text().strip() if it else "") or str(DEFAULT_THICKNESS)
                return float(txt.replace(",", "."))
            except ValueError:
                return float(DEFAULT_THICKNESS)

        def ask_alja_thick(self, r, extras):
            """Sittes: az alja lemezvastagsága egyezik-e a többi lemezével; ha nem, bekéri."""
            t = self.row_thick(r)
            ans = QMessageBox.question(
                self, "Alja lemezvastagsága",
                f"Az alja lemezvastagsága megegyezik a többi oldal lemezvastagságával ({t:g} mm)?",
                QMessageBox.Yes | QMessageBox.No, QMessageBox.Yes)
            if ans == QMessageBox.Yes:
                extras["alja_thick"] = ""
            else:
                try:
                    cur = float(extras.get("alja_thick", "").replace(",", ".") or t)
                except ValueError:
                    cur = t
                v, ok = QInputDialog.getDouble(self, "Alja lemezvastagsága",
                                               "Az alja lemezvastagsága (mm):", cur, 0.1, 99.0, 2)
                if ok:
                    extras["alja_thick"] = f"{v:g}" if abs(v - t) > 1e-9 else ""
            self.show_alja_thick(r)

        def show_alja_thick(self, r):
            btn = self.tbl.cellWidget(r, COL_DRAW)
            if btn is None:
                return
            if self.row_type(r) == OTHER_TYPE:
                p = self.row_extras(r).get("dxf", "")
                btn.setText(f"DXF: {os.path.basename(p)}" if p else "DXF betöltése…")
                btn.setToolTip(p or "Az alkatrész rajzának (DXF) betöltése")
                return
            at = self.row_extras(r).get("alja_thick", "") if self.row_type(r) == "Sittes" else ""
            btn.setText(f"Rajz… (alja {at} mm)" if at else "Rajz…")
            btn.setToolTip("Méretek megadása a rajz szerint (ABO_A, Sittes)" +
                           (f"\nAz alja lemezvastagsága: {at} mm (eltér a többitől)" if at else ""))

        def update_row(self, r):
            """Típus szerint: tooltipek az oldal-oszlopokon, Sittesnél az 'Oldal C' tiltva."""
            ctype = self.row_type(r)
            btn = self.tbl.cellWidget(r, COL_DRAW)
            if btn is not None and hasattr(btn, "setEnabled"):
                btn.setEnabled(ctype in SKETCH_CONFIGS or ctype == OTHER_TYPE)
                self.show_alja_thick(r)
            self.update_flange_enabled()
            for c in range(COL_FIRST_DIM, COL_SIDE):
                it = self.tbl.item(r, c)
                if it is None:
                    it = QTableWidgetItem("")
                    self.tbl.setItem(r, c, it)
                if ctype == OTHER_TYPE:
                    it.setFlags(it.flags() & ~Qt.ItemIsEditable & ~Qt.ItemIsEnabled)
                    it.setToolTip("Egyéb alkatrésznél a méretek a DXF rajzból jönnek")
                elif c != COL_HOMLOK_H:
                    it.setFlags(it.flags() | Qt.ItemIsEditable | Qt.ItemIsEnabled)
            for i, c in enumerate(range(COL_SIDE, COL_SIDE + 3)):
                it = self.tbl.item(r, c)
                if it is None:
                    it = QTableWidgetItem("")
                    self.tbl.setItem(r, c, it)
                it.setToolTip(SIDE_HELP[ctype][i])
                if ctype == OTHER_TYPE or (ctype not in ("ABO_A", "ABO_B") and i == 2):
                    if ctype != OTHER_TYPE:
                        it.setText("")
                    it.setFlags(it.flags() & ~Qt.ItemIsEditable & ~Qt.ItemIsEnabled)
                else:
                    it.setFlags(it.flags() | Qt.ItemIsEditable | Qt.ItemIsEnabled)
            it = self.tbl.item(r, COL_HOMLOK_H)
            if it is None:
                it = QTableWidgetItem("")
                self.tbl.setItem(r, COL_HOMLOK_H, it)
            if ctype == OTHER_TYPE:
                it.setFlags(it.flags() & ~Qt.ItemIsEditable & ~Qt.ItemIsEnabled)
            elif ctype == "ABO_A":
                it.setToolTip("ABO_A-nál a homloklap NETTÓ magassága = a ferde oldalél hossza "
                              "(a nettó magasságig); a program számolja ki, és hozzáadja a peremet")
                it.setFlags(it.flags() & ~Qt.ItemIsEditable)
                self.fill_abo_a_front(r)
            else:
                if it.text().startswith("="):
                    it.setText("")
                it.setToolTip("Homloklap nettó magassága")
                it.setFlags(it.flags() | Qt.ItemIsEditable | Qt.ItemIsEnabled)

        def on_item_changed(self, item):
            if self._updating or item is None:
                return
            if COL_SIDE <= item.column() < COL_SIDE + 3:
                self.fill_abo_a_front(item.row())

        def current_flange(self) -> float:
            try:
                fl = [float(t.replace(",", ".")) for t in
                      re.split(r"[xX×*; ]+", self.txt_flange.text().strip()) if t]
                return sum(fl) if len(fl) == 2 else sum(DEFAULT_FLANGE)
            except (ValueError, AttributeError):
                return sum(DEFAULT_FLANGE)

        def fill_abo_a_front(self, r):
            """ABO_A: a homloklap magasságát (= ferde él) automatikusan beírja."""
            if self.row_type(r) != "ABO_A":
                return
            it = self.tbl.item(r, COL_HOMLOK_H)
            if it is None:
                return
            try:
                a, b, c = (float(self.tbl.item(r, col).text().replace(",", "."))
                           for col in range(COL_SIDE, COL_SIDE + 3))
                if not (a > 0 and b > 0 and c > 0 and c < a):
                    raise ValueError
                _, slant, _ = abo_a_info(a, b, c, self.current_flange())
                txt = f"= {slant:.1f}".replace(".", ",")
            except (ValueError, AttributeError):
                txt = "= (ferde él)"
            self._updating = True
            try:
                it.setText(txt)
            finally:
                self._updating = False

        def refresh_all_rows(self):
            for r in range(self.tbl.rowCount()):
                self.fill_abo_a_front(r)

        def row_type(self, r) -> str:
            cb = self.tbl.cellWidget(r, COL_TYPE)
            return cb.currentText() if cb is not None else "Sittes"

        def delete_product(self):
            rows = sorted({i.row() for i in self.tbl.selectedIndexes()}, reverse=True)
            if not rows and self.tbl.rowCount() > 0:
                rows = [self.tbl.rowCount() - 1]
            for r in rows:
                self.tbl.removeRow(r)
            self.update_flange_enabled()

        def update_flange_enabled(self):
            """Az ABO perem mező csak akkor aktív, ha van ABO típusú termék."""
            if not hasattr(self, "txt_flange"):
                return
            on = any(self.row_type(r) in ("ABO_A", "ABO_B") for r in range(self.tbl.rowCount()))
            self.txt_flange.setEnabled(on)
            lbl = self.field_labels.get(id(self.txt_flange))
            if lbl is not None:
                lbl.setEnabled(on)
            self.txt_flange.setToolTip("" if on else "Csak ABO_A / ABO_B konténernél van perem")

        def read_products(self) -> List[List[Product]]:
            prods = []
            names = set()
            for r in range(self.tbl.rowCount()):
                def cell(c):
                    it = self.tbl.item(r, c)
                    return it.text().strip() if it else ""
                name = cell(0) or f"Termék {r + 1}"
                ctype = self.row_type(r)
                if ctype == OTHER_TYPE:
                    path = self.row_extras(r).get("dxf", "")
                    if not path:
                        continue
                    if name in names:
                        raise ValueError(f"A(z) „{name}” név kétszer szerepel – adj egyedi neveket.")
                    names.add(name)
                    try:
                        qty = int(cell(COL_QTY))
                        thick = float((cell(COL_THICK) or DEFAULT_THICKNESS).replace(",", "."))
                        if qty < 0 or not (0 < thick < 100):
                            raise ValueError
                    except ValueError:
                        raise ValueError(f"„{name}”: a darabszám egész szám, a vastagság pozitív szám legyen.")
                    if qty == 0:
                        continue
                    try:
                        parts = load_dxf_parts(path, name)
                    except Exception as e:
                        raise ValueError(f"„{name}”: {e}")
                    pr = Product(name, qty, parts, 0.0, "DXF: " + os.path.basename(path) + " – " +
                                 dxf_parts_summary(parts).replace("\n", "; "),
                                 ctype=OTHER_TYPE, thick=thick)
                    prods.append([pr])
                    continue
                n_side = 3 if ctype in ("ABO_A", "ABO_B") else 2
                dim_cols = [c for c in range(COL_FIRST_DIM, COL_SIDE + n_side)
                            if not (ctype == "ABO_A" and c == COL_HOMLOK_H)]
                if not any(cell(c) for c in dim_cols):
                    continue   # teljesen üres sor
                if name in names:
                    raise ValueError(f"A(z) „{name}” név kétszer szerepel – adj egyedi neveket.")
                names.add(name)
                try:
                    qty = int(cell(COL_QTY))
                    v = [float(cell(c).replace(",", ".")) for c in dim_cols]
                    if qty < 0 or min(v) <= 0:
                        raise ValueError
                except ValueError:
                    raise ValueError(f"„{name}” ({ctype}): a darabszám és minden szükséges méret "
                                     f"pozitív szám legyen (mm). Oldalméretek: "
                                     + "; ".join(SIDE_HELP[ctype][:n_side]))
                if qty == 0:
                    continue
                try:
                    thick = float((cell(COL_THICK) or DEFAULT_THICKNESS).replace(",", "."))
                    if not (0 < thick < 100):
                        raise ValueError
                except ValueError:
                    raise ValueError(f"„{name}”: a lemezvastagság pozitív szám legyen (mm), pl. 2 vagy 1,5.")
                if ctype == "ABO_A":
                    aw, ah, kw, kh, fw = v[:5]
                    fh = 0.0          # a ferde élből számolja
                    side = v[5:]
                else:
                    aw, ah, kw, kh, fw, fh = v[:6]
                    side = v[6:]
                try:
                    mode = self.cmb_merge.currentIndex()
                    merges = MERGE_MODES[mode][1]
                    ex = {}
                    if ctype in ("Sittes", "ABO_B"):
                        for k_, v_ in self.row_extras(r).items():
                            try:
                                if v_.strip():
                                    ex[k_] = float(v_.replace(",", "."))
                            except ValueError:
                                raise ValueError(f"„{name}”: hibás méret a rajzon ({k_}: {v_}).")
                    vs = make_product_variants(name, ctype, qty, (aw, ah), (kw, kh), (fw, fh), side,
                                               extras=ex,
                                               flange=getattr(self, "flange", DEFAULT_FLANGE),
                                               weld_sides=2 if self.chk_double.isChecked() else 1,
                                               merges=merges)
                    if mode > 1 and not any(v.merge for v in vs):
                        raise ValueError(f"{MERGE_MODES[mode][0]} nem lehetséges: a közös élek "
                                         f"hossza nem egyezik (alja {aw:g}x{ah:g}, hátlap szél. "
                                         f"{kw:g}, homloklap szél. {fw:g}).")
                    if mode > 1:
                        vs = [v for v in vs if v.merge] or vs
                    at = ex.get("alja_thick") if ctype == "Sittes" else None
                    if at and abs(at - thick) > 1e-9:
                        # eltérő vastagságú alja: nem hajlítható egy lemezből a hátlappal/homloklappal
                        if mode > 1:
                            raise ValueError(f"{MERGE_MODES[mode][0]} nem lehetséges, mert az alja "
                                             f"lemezvastagsága ({at:g} mm) eltér a többitől ({thick:g} mm).")
                        vs = [v for v in vs if not v.merge]
                    for v_ in vs:
                        v_.thick = thick
                        if at and abs(at - thick) > 1e-9:
                            v_.part_thick = {"Alja": at}
                    prods.append(vs)
                except ValueError as e:
                    raise ValueError(f"„{name}”: {e}")
            if not prods:
                raise ValueError("Adj meg legalább egy terméket méretekkel és darabszámmal.")
            return prods

        # --- anyagfájl ---
        def update_sheet_label(self, src):
            txt = ", ".join(f"{s.width}x{s.height}" for s in self.sheets)
            self.lbl_sheets.setText(f"Táblaméretek ({src}): {txt}")
            self.material_path = src

        def load_material(self, path):
            try:
                self.sheets = load_sheet_sizes(path)
                self.material_full_path = path
                try:
                    th, _mat = load_material_info(path)
                except Exception:
                    th = None
                self.update_sheet_label(os.path.basename(path))
            except Exception as e:
                QMessageBox.warning(self, "Anyagfájl", f"Nem sikerült beolvasni:\n{e}")

        def choose_material_file(self):
            p, _ = QFileDialog.getOpenFileName(self, "Anyagfájl", "", "Szövegfájl (*.txt)")
            if p:
                self.load_material(p)

        # --- számítás ---
        def calculate(self):
            if not self.sheets:
                QMessageBox.warning(self, "Hiba", "Előbb töltsd be az anyagolás.txt fájlt.")
                return
            try:
                kerf = int(self.txt_kerf.text() or 0)
                min_pct = float((self.txt_min_pct.text() or "0").replace(",", "."))
                time_limit = float((self.txt_time.text() or "0").replace(",", "."))
                fl = [float(t.replace(",", ".")) for t in
                      re.split(r"[xX×*; ]+", self.txt_flange.text().strip()) if t]
                if len(fl) != 2 or min(fl) < 0:
                    raise ValueError
                self.flange = tuple(fl)
                if time_limit < 0:
                    raise ValueError
                if kerf < 0 or not (0 <= min_pct <= 50):
                    raise ValueError
            except ValueError:
                QMessageBox.warning(self, "Hiba", "Ellenőrizd a vágási rést (egész mm), a toldás %-át (0–50), a keresési időt és a perem méretét (pl. 30x25).")
                return
            try:
                products = self.read_products()
            except ValueError as e:
                QMessageBox.warning(self, "Hiba", str(e))
                return
            try:
                margin = float((self.txt_margin.text() or "0").replace(",", "."))
                max_cut = DEFAULT_MAX_CUT if self.chk_lasercut.isChecked() else None
                if margin < 0:
                    raise ValueError("A széltávolság nem lehet negatív (mm).")
                opt_sheets = usable_sheets(self.sheets, margin, max_cut)
            except ValueError as e:
                QMessageBox.warning(self, "Hiba", str(e) or "Ellenőrizd a széltávolságot és a max. vágási hosszt.")
                return
            self.calc_kerf = kerf
            _CANCEL.clear()
            self.btn_calc.setEnabled(False)
            self.btn_cancel.setEnabled(True)
            self.btn_dxf.setEnabled(False)
            self.btn_quote.setEnabled(False)
            self.calc_t0, self.calc_limit = time.time(), time_limit
            self.calc_msg = "Számolás folyamatban…"
            self.show_progress()
            self.tick.start()
            self.thread = CalcThread(self, (products, opt_sheets, kerf, self.chk_guillotine.isChecked(),
                                            self.chk_weld.isChecked(), min_pct),
                                     dict(time_limit=time_limit))
            self.thread.progress.connect(self.on_progress)
            self.thread.done.connect(self.on_calc_done)
            self.thread.failed.connect(self.on_calc_failed)
            self.thread.start()

        def on_progress(self, msg):
            self.calc_msg = msg
            self.show_progress()

        def show_progress(self):
            el = time.time() - self.calc_t0
            self.txt_result.setPlainText(
                f"{self.calc_msg}\n\nEltelt idő: {el:.0f} / {self.calc_limit:g} mp"
                f"  ({self.workers} szálon számol)")

        def _calc_finished(self):
            self.tick.stop()
            self.btn_calc.setEnabled(True)
            self.btn_cancel.setEnabled(False)

        def cancel_calc(self):
            _CANCEL.set()
            self.btn_cancel.setEnabled(False)
            self.calc_msg = "Megszakítás…"
            self.show_progress()

        def on_calc_failed(self, msg):
            self._calc_finished()
            if not msg:
                self.txt_result.setPlainText("A számítás megszakítva.")
                return
            self.txt_result.setPlainText("")
            QMessageBox.critical(self, "Hiba", msg)

        def on_calc_done(self, res):
            self._calc_finished()
            placed, used, products, table = res[:4]
            other_rep = res[4] if len(res) > 4 else []
            el = time.time() - self.calc_t0
            self.last = (placed, used)
            self.last_products = products
            allowed = ", ".join(f"{s.nw}x{s.nh}" for s in self.sheets)
            usable = ""
            if used and any((s.width, s.height) != (s.nw, s.nh) for s in used):
                m = used[0].off
                usable = (f"Vágható terület: széltávolság {m} mm"
                          + (f", lézeres vágás: max. {DEFAULT_MAX_CUT} mm hossz"
                             if self.chk_lasercut.isChecked() else "") + "\n")
            self.txt_result.setPlainText(
                f"Felhasználható táblák ({self.material_path}): {allowed}\n{usable}"
                f"Számítási idő: {el:.1f} mp ({self.workers} szálon)\n\n"
                + format_result(placed, used, self.calc_kerf, products)
                + ("\n".join(other_rep) + "\n" if other_rep else "")
                + format_variant_table(table))
            self.btn_dxf.setEnabled(True)
            self.btn_quote.setEnabled(True)

        def closeEvent(self, ev):
            _CANCEL.set()
            if self.thread is not None and self.thread.isRunning():
                self.thread.wait(3000)
            if self.pool is not None:
                self.pool.terminate()
                self.pool = None
            super().closeEvent(ev)

        def export_quote(self):
            if not self.last:
                return
            projekt, termekek, tablak, beepulo, tabla = quote_export_data(
                self.last[0], self.last[1], self.last_products)
            anyag = None
            if self.material_full_path:
                try:
                    _t, anyag = load_material_info(self.material_full_path)
                except Exception:
                    pass
            ths = sorted({s.thick for s in self.last[1]})
            vastag = ths[0] if len(ths) == 1 else None
            suly = sum(s.weight() for s in self.last[1])
            safe = re.sub(r'[\\/:*?"<>|]+', "_", projekt)
            path, _ = QFileDialog.getSaveFileName(self, "Export árajánlathoz", f"{safe}_szabas.json",
                                                  "JSON fájl (*.json)")
            if not path:
                return
            if not path.lower().endswith(".json"):
                path += ".json"
            try:
                export_arajanlathoz(path, projekt, termekek, tablak, beepulo, tabla, vastag, anyag, suly)
            except Exception as e:
                QMessageBox.critical(self, "Hiba", f"Export sikertelen:\n{e}")
                return
            extra = []
            extra.append(f"vastagság: {vastag:g} mm" if vastag else
                         "vastagság táblánként: " + ", ".join(f"{t:g} mm" for t in ths))
            extra.append(f"súly: {suly:.1f} kg")
            extra.append(f"anyag: {anyag}" if anyag else "anyag: nincs megadva")
            QMessageBox.information(self, "Kész", f"Árajánlat-export elkészült ({', '.join(extra)}):\n{path}")

        def export(self):
            if not self.last:
                return
            if self.chk_laser.isChecked():
                folder = QFileDialog.getExistingDirectory(self, "Mappa a lézeres DXF fájloknak")
                if folder:
                    try:
                        files = export_laser(folder, self.last[0], sheets=self.last[1])
                        QMessageBox.information(
                            self, "Kész", f"{len(files)} alkatrész-DXF + darabjegyzek.txt elkészült:\n{folder}")
                    except Exception as e:
                        QMessageBox.critical(self, "Hiba", f"Lézer export sikertelen:\n{e}")
                return
            path, _ = QFileDialog.getSaveFileName(self, "DXF export", "", "DXF Files (*.dxf)")
            if path:
                try:
                    export_dxf(path, *self.last)
                    QMessageBox.information(self, "Kész", f"DXF export sikeres:\n{path}")
                except Exception as e:
                    QMessageBox.critical(self, "Hiba", f"DXF export sikertelen:\n{e}")

    app = QApplication(sys.argv)
    win = ContainerApp()
    win.show()
    sys.exit(app.exec_())


if __name__ == "__main__":
    multiprocessing.freeze_support()     # exe-ként (PyInstaller) a párhuzamos számításhoz kell
    run_gui()
