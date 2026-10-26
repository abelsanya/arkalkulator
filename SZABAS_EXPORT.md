# Szabásprogram → árajánlat-kalkulátor export

Ez a leírás a **Konténer lemezszabás kalkulátor** (asztali program) fejlesztéséhez készült. A program egy új gombbal (**„Export árajánlathoz”**) JSON-fájlt ment. Az Árkalkulátor ezt a fájlt be tudja tölteni: a projekt Anyagok részén a „Lemez betöltése a szabás-számításból” sávban kell az „Exportfájl betöltése…” gombra kattintani.

## A fájl formája (UTF-8 JSON, kiterjesztése .json)

```json
{
  "forras": "Konténer lemezszabás kalkulátor",
  "verzio": 1,
  "projekt": "ABO_A",
  "termekek": [
    { "nev": "Termék 1", "tipus": "ABO_A", "darab": 4 }
  ],
  "vastagsag_mm": 2.0,
  "anyag": "Hidegen hengerelt",
  "tablak": [
    { "szel": 1000, "hossz": 2000, "db": 4 }
  ],
  "tabla_terulet_m2": 8.00,
  "beepulo_terulet_m2": 7.15
}
```

| Mező | Kötelező | Jelentés |
|---|---|---|
| `tablak` | igen | A felhasznált táblák méretenként: szélesség és hossz mm-ben, valamint darabszám. Ugyanaz, mint a „Felhasznált lemezek száma” blokk. |
| `beepulo_terulet_m2` | nem, de ajánlott | Ebből számolja a kalkulátor a nettó (beépülő) súlyt. |
| `termekek[].darab` | nem, de ajánlott | A darabok összege lesz a kalkulátorban a projekt darabszáma, így az egységár helyes. |
| `vastagsag_mm` | nem | Ha meg van adva, a kalkulátor ilyen vastag lemezt választ. Ha nincs, a kalkulátorban előbb ki kell választani a lemezt. |
| `anyag` | nem | A lemez fajtája az árlistában, pl. „Hidegen hengerelt”, „Melegen hengerelt”, „Horganyzott”, „Rozsdamentes”. Alapértelmezés: hidegen hengerelt. |
| `projekt`, `forras`, `verzio`, `tabla_terulet_m2` | nem | Tájékoztató adatok. |

A táblaméretnek szerepelnie kell az árlistában az adott vastagságban (pl. 1000x2000, 1250x2500, 1500x3000), különben a kalkulátor hibaüzenetet ad, és nem változtat semmit.

## Python minta az asztali programhoz

```python
import json
from pathlib import Path

def export_arajanlathoz(path, projekt, termekek, tablak, beepulo_m2, tabla_m2,
                        vastagsag_mm=None, anyag=None):
    """termekek: [(nev, tipus, darab)], tablak: [(szel_mm, hossz_mm, db)]"""
    adat = {
        "forras": "Konténer lemezszabás kalkulátor",
        "verzio": 1,
        "projekt": projekt,
        "termekek": [{"nev": n, "tipus": t, "darab": int(d)} for n, t, d in termekek],
        "tablak": [{"szel": int(w), "hossz": int(h), "db": int(n)} for w, h, n in tablak],
        "tabla_terulet_m2": round(tabla_m2, 3),
        "beepulo_terulet_m2": round(beepulo_m2, 3),
    }
    if vastagsag_mm:
        adat["vastagsag_mm"] = float(vastagsag_mm)
    if anyag:
        adat["anyag"] = anyag
    Path(path).write_text(json.dumps(adat, ensure_ascii=False, indent=2), encoding="utf-8")
```

A gombot az „Optimalizált számítás” után érdemes engedélyezni. Ugyanazokból az adatokból kell dolgoznia, amelyekből a program a „Felhasznált lemezek száma” blokkot írja ki. A fájlnév javaslata: `<projekt>_szabas.json`.

A sima szöveges eredmény bemásolása továbbra is működik, ha nincs kéznél exportfájl.
