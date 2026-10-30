"""
Webes (iPhone / böngésző) felület híd a kontener_szabas.py számításához.

A kontener_szabas.py változatlanul ugyanaz, mint az asztali program; ez a fájl
csak a PyQt5-ös ablak (ContainerApp) adatkezelését veszi át: a termékek
beolvasását, a számítás indítását és az exportokat. A böngészőben a Pyodide
futtatja, egy Web Workerben (worker.js).
"""

import io
import json
import os
import re
import time
import zipfile

import kontener_szabas as ks

WORK = "/work"
MATERIAL_PATH = WORK + "/anyagolas.txt"
os.makedirs(WORK, exist_ok=True)

_LAST = {}          # az utolsó sikeres számítás: placed, used, products, material


def info() -> str:
    """A felület felépítéséhez szükséges állandók (típusok, oszlopok, vázlatok)."""
    sketches = {}
    for key, cfg in ks.SKETCH_CONFIGS.items():
        sketches[key] = {
            "title": cfg["title"],
            "sketch": cfg["sketch"](),
            "fields": {n: [t, lab] for n, (t, lab) in cfg["fields"].items()},
            "optional": list(cfg["optional"]),
            "links": {k: list(v) for k, v in cfg["links"].items()},
            "derived": {str(k): list(v) for k, v in cfg.get("derived", {}).items()},
            "toggle": list(cfg["toggle"]) if "toggle" in cfg else None,
        }
    return json.dumps({
        "version": ks.PROGRAM_VERSION,
        "types": ks.CONTAINER_TYPES,
        "other": ks.OTHER_TYPE,
        "columns": ks.PRODUCT_COLUMNS,
        "side_help": ks.SIDE_HELP,
        "merge_modes": [m[0] for m in ks.MERGE_MODES],
        "default_thickness": ks.DEFAULT_THICKNESS,
        "default_margin": ks.DEFAULT_MARGIN,
        "max_cut": ks.DEFAULT_MAX_CUT,
        "default_flange": list(ks.DEFAULT_FLANGE),
        "sketches": sketches,
    }, ensure_ascii=False)


def _num_from(values):
    def num(name):
        try:
            v = float(str(values.get(name, "")).strip().replace(",", "."))
            return v if v > 0 else None
        except ValueError:
            return None
    return num


def sketch_calc(key: str, values_json: str, flange: float) -> str:
    """A rajzos ablak számolt értékei: (info, homloklap-cella, hiba)."""
    info_txt, front, err = ks.SKETCH_CONFIGS[key]["calc"](_num_from(json.loads(values_json)), flange)
    return json.dumps([info_txt, front, err], ensure_ascii=False)


def material_info(text: str) -> str:
    """Az anyagfájl (anyagolás.txt) ellenőrzése: táblaméretek, vastagság, anyag."""
    with open(MATERIAL_PATH, "w", encoding="utf-8") as f:
        f.write(text)
    try:
        sheets = ks.load_sheet_sizes(MATERIAL_PATH)
    except Exception as e:
        return json.dumps({"error": str(e)}, ensure_ascii=False)
    try:
        th, mat = ks.load_material_info(MATERIAL_PATH)
    except Exception:
        th, mat = None, None
    return json.dumps({"sheets": [f"{s.width}x{s.height}" for s in sheets],
                       "thick": th, "material": mat}, ensure_ascii=False)


def dxf_summary(path: str, name: str) -> str:
    """Egyéb alkatrész DXF-ének beolvasása (ellenőrzés + összefoglaló)."""
    try:
        parts = ks.load_dxf_parts(path, name)
    except Exception as e:
        return json.dumps({"error": str(e)}, ensure_ascii=False)
    return json.dumps({"summary": ks.dxf_parts_summary(parts)}, ensure_ascii=False)


def _flange(text: str):
    fl = [float(t.replace(",", ".")) for t in re.split(r"[xX×*; ]+", (text or "").strip()) if t]
    if len(fl) != 2 or min(fl) < 0:
        raise ValueError
    return tuple(fl)


def read_products(rows, merge_mode: int, flange, double_weld: bool):
    """Az asztali program ContainerApp.read_products megfelelője: a sorok a felület
    termékkártyáiból jönnek (cells: oszlopindex → szöveg, extras: rejtett méretek)."""
    COL_QTY, COL_FIRST_DIM, COL_SIDE = ks.COL_QTY, ks.COL_FIRST_DIM, ks.COL_SIDE
    prods = []
    names = set()
    for r, row in enumerate(rows):
        cells = row.get("cells", {})

        def cell(c):
            return str(cells.get(str(c), "")).strip()

        name = str(row.get("name", "")).strip() or f"Termék {r + 1}"
        ctype = row.get("type") or "Sittes"
        qty_txt = str(row.get("qty", "")).strip()
        thick_txt = str(row.get("thick", "")).strip()
        extras = row.get("extras", {}) or {}
        if ctype == ks.OTHER_TYPE:
            path = row.get("dxf_path", "")
            if not path:
                continue
            if name in names:
                raise ValueError(f"A(z) „{name}” név kétszer szerepel – adj egyedi neveket.")
            names.add(name)
            try:
                qty = int(qty_txt)
                thick = float((thick_txt or ks.DEFAULT_THICKNESS).replace(",", "."))
                if qty < 0 or not (0 < thick < 100):
                    raise ValueError
            except ValueError:
                raise ValueError(f"„{name}”: a darabszám egész szám, a vastagság pozitív szám legyen.")
            if qty == 0:
                continue
            try:
                parts = ks.load_dxf_parts(path, name)
            except Exception as e:
                raise ValueError(f"„{name}”: {e}")
            pr = ks.Product(name, qty, parts, 0.0, "DXF: " + str(row.get("dxf_name", "")) + " – " +
                            ks.dxf_parts_summary(parts).replace("\n", "; "),
                            ctype=ks.OTHER_TYPE, thick=thick)
            prods.append([pr])
            continue
        n_side = 3 if ctype in ("ABO_A", "ABO_B") else 2
        dim_cols = [c for c in range(COL_FIRST_DIM, COL_SIDE + n_side)
                    if not (ctype == "ABO_A" and c == ks.COL_HOMLOK_H)]
        if not any(cell(c) for c in dim_cols):
            continue   # teljesen üres sor
        if name in names:
            raise ValueError(f"A(z) „{name}” név kétszer szerepel – adj egyedi neveket.")
        names.add(name)
        try:
            qty = int(qty_txt)
            v = [float(cell(c).replace(",", ".")) for c in dim_cols]
            if qty < 0 or min(v) <= 0:
                raise ValueError
        except ValueError:
            raise ValueError(f"„{name}” ({ctype}): a darabszám és minden szükséges méret "
                             f"pozitív szám legyen (mm). Oldalméretek: "
                             + "; ".join(ks.SIDE_HELP[ctype][:n_side]))
        if qty == 0:
            continue
        try:
            thick = float((thick_txt or ks.DEFAULT_THICKNESS).replace(",", "."))
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
            merges = ks.MERGE_MODES[merge_mode][1]
            ex = {}
            if ctype in ("Sittes", "ABO_B"):
                for k_, v_ in extras.items():
                    v_ = str(v_)
                    try:
                        if v_.strip():
                            ex[k_] = float(v_.replace(",", "."))
                    except ValueError:
                        raise ValueError(f"„{name}”: hibás méret a rajzon ({k_}: {v_}).")
            vs = ks.make_product_variants(name, ctype, qty, (aw, ah), (kw, kh), (fw, fh), side,
                                          extras=ex, flange=flange,
                                          weld_sides=2 if double_weld else 1,
                                          merges=merges)
            if merge_mode > 1 and not any(v.merge for v in vs):
                raise ValueError(f"{ks.MERGE_MODES[merge_mode][0]} nem lehetséges: a közös élek "
                                 f"hossza nem egyezik (alja {aw:g}x{ah:g}, hátlap szél. "
                                 f"{kw:g}, homloklap szél. {fw:g}).")
            if merge_mode > 1:
                vs = [v for v in vs if v.merge] or vs
            at = ex.get("alja_thick") if ctype == "Sittes" else None
            if at and abs(at - thick) > 1e-9:
                # eltérő vastagságú alja: nem hajlítható egy lemezből a hátlappal/homloklappal
                if merge_mode > 1:
                    raise ValueError(f"{ks.MERGE_MODES[merge_mode][0]} nem lehetséges, mert az alja "
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


def calculate(params_json: str, progress=None) -> str:
    """Az „Optimalizált számítás” gomb (ContainerApp.calculate + on_calc_done).
    Hibánál {"error": ...}, sikernél a szöveges eredmény és a lemezek rajza."""
    p = json.loads(params_json)
    _LAST.clear()
    with open(MATERIAL_PATH, "w", encoding="utf-8") as f:
        f.write(p.get("material", ""))
    try:
        sheets = ks.load_sheet_sizes(MATERIAL_PATH)
    except Exception as e:
        return json.dumps({"error": f"Anyagfájl: {e}"}, ensure_ascii=False)
    try:
        kerf = int(p.get("kerf") or 0)
        min_pct = float(str(p.get("min_pct") or "0").replace(",", "."))
        time_limit = float(str(p.get("time") or "0").replace(",", "."))
        flange = _flange(p.get("flange", ""))
        if time_limit < 0 or kerf < 0 or not (0 <= min_pct <= 50):
            raise ValueError
    except ValueError:
        return json.dumps({"error": "Ellenőrizd a vágási rést (egész mm), a toldás %-át (0–50), "
                                    "a keresési időt és a perem méretét (pl. 30x25)."},
                          ensure_ascii=False)
    try:
        products = read_products(p.get("rows", []), int(p.get("merge", 0)), flange,
                                 bool(p.get("double")))
    except ValueError as e:
        return json.dumps({"error": str(e)}, ensure_ascii=False)
    lasercut = bool(p.get("lasercut"))
    try:
        margin = float(str(p.get("margin") or "0").replace(",", "."))
        max_cut = ks.DEFAULT_MAX_CUT if lasercut else None
        if margin < 0:
            raise ValueError("A széltávolság nem lehet negatív (mm).")
        opt_sheets = ks.usable_sheets(sheets, margin, max_cut)
    except ValueError as e:
        return json.dumps({"error": str(e) or "Ellenőrizd a széltávolságot és a max. vágási hosszt."},
                          ensure_ascii=False)
    guillotine = bool(p.get("guillotine"))
    t0 = time.time()
    try:
        res = ks.optimize_by_thickness(products, opt_sheets, kerf, guillotine,
                                       bool(p.get("weld")), min_pct, time_limit=time_limit,
                                       progress=progress, pool=None, workers=1)
        prods0 = [vl[0] for vl in products]
        if any(pr.ctype == ks.OTHER_TYPE for pr in prods0):
            if progress:
                progress("Egyéb alkatrészek: táblaméretek összevetése…")
            ks.set_nesting(True)
            rep = ks.other_sheet_report(prods0, opt_sheets, kerf, guillotine)
            res = tuple(res) + (rep,)
    except ks.Cancelled:
        return json.dumps({"error": "A számítás megszakítva."}, ensure_ascii=False)
    except Exception as e:
        return json.dumps({"error": str(e) or e.__class__.__name__}, ensure_ascii=False)

    placed, used, products, table = res[:4]
    other_rep = res[4] if len(res) > 4 else []
    el = time.time() - t0
    _LAST.update(placed=placed, used=used, products=products)
    allowed = ", ".join(f"{s.nw}x{s.nh}" for s in sheets)
    usable = ""
    if used and any((s.width, s.height) != (s.nw, s.nh) for s in used):
        m = used[0].off
        usable = (f"Vágható terület: széltávolság {m} mm"
                  + (f", lézeres vágás: max. {ks.DEFAULT_MAX_CUT} mm hossz" if lasercut else "") + "\n")
    text = (f"Felhasználható táblák (anyagfájl): {allowed}\n{usable}"
            f"Számítási idő: {el:.1f} mp\n\n"
            + ks.format_result(placed, used, kerf, products)
            + ("\n".join(other_rep) + "\n" if other_rep else "")
            + ks.format_variant_table(table))
    return json.dumps({"text": text, "drawing": _drawing(placed, used),
                       "quote": json.loads(quote())}, ensure_ascii=False)


def _drawing(placed, used):
    """A lemezek kiosztása a képernyős rajzhoz (mm, a tábla bal alsó sarkától)."""
    out = []
    for i, s in enumerate(used):
        pcs = []
        for pp in placed:
            if pp.sheet_index != i:
                continue
            pcs.append({
                "tag": pp.tag, "res": bool(pp.from_residual), "part": pp.part,
                "shapes": [[lab, [[round(x + s.off, 1), round(y + s.off, 1)] for x, y in pts]]
                           for lab, pts in ks.placed_polys(pp)],
            })
        out.append({"w": s.nw, "h": s.nh, "cw": s.width, "ch": s.height, "off": s.off,
                    "thick": s.thick, "pieces": pcs})
    return out


def quote() -> str:
    """Az „Export árajánlathoz” fájl tartalma (SZABAS_EXPORT.md)."""
    placed, used, products = _LAST["placed"], _LAST["used"], _LAST["products"]
    projekt, termekek, tablak, beepulo, tabla = ks.quote_export_data(placed, used, products)
    anyag = None
    try:
        _t, anyag = ks.load_material_info(MATERIAL_PATH)
    except Exception:
        pass
    ths = sorted({s.thick for s in used})
    vastag = ths[0] if len(ths) == 1 else None
    suly = sum(s.weight() for s in used)
    path = WORK + "/arajanlat.json"
    ks.export_arajanlathoz(path, projekt, termekek, tablak, beepulo, tabla, vastag, anyag, suly)
    with open(path, encoding="utf-8") as f:
        adat = json.load(f)
    safe = re.sub(r'[\\/:*?"<>|]+', "_", projekt)
    return json.dumps({"data": adat, "filename": f"{safe}_szabas.json"}, ensure_ascii=False)


def export_dxf() -> str:
    """DXF export: a fájl útvonala a Pyodide fájlrendszerében."""
    path = WORK + "/szabas.dxf"
    ks.export_dxf(path, _LAST["placed"], _LAST["used"])
    return path


def export_laser() -> str:
    """Lézer export: alkatrészenként külön DXF + darabjegyzék, egy ZIP-be csomagolva."""
    folder = WORK + "/lezer"
    if os.path.isdir(folder):
        for n in os.listdir(folder):
            os.remove(os.path.join(folder, n))
    files = ks.export_laser(folder, _LAST["placed"], sheets=_LAST["used"])
    path = WORK + "/lezer.zip"
    with zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED) as z:
        for n in sorted(os.listdir(folder)):
            z.write(os.path.join(folder, n), n)
    return json.dumps({"path": path, "count": len(files)})
