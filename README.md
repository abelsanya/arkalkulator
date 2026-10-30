# Árkalkulátor

Árajánlat-kalkulátor a Tímár Vasker, tetolemez.hu és metallonova.hu árai alapján.

- `index.html`: a kalkulátor
- `prices.json`: árlista (hetente frissül)
- `termektipusok.json`: alap terméktípusok (Billenő konténer, Dohányzó), első megnyitáskor töltődnek be

A projektek a böngészőben tárolódnak.

## Konténer lemezszabás (szabas.html)

Az asztali szabásprogram (PyQt5) webes, telefonon is használható változata. A Python
számítás változatlanul fut a böngészőben (Pyodide, WebAssembly), csak a felület új.

- `szabas.html`: a felület. iPhone-on Safariban: Megosztás → „Hozzáadás a főképernyőhöz”, így külön alkalmazásként indul.
- `szabas/kontener_szabas.py`: az asztali program forrása, változtatás nélkül. Új verziónál elég ezt a fájlt lecserélni (és a `szabas.html`-ben az `ASSET_VER` értékét növelni, hogy a telefon ne a régit használja).
- `szabas/web.py`: híd a felület és a számítás között (a PyQt-s ablak adatkezelése helyett).
- `szabas/worker.js`: a számítás külön szálon fut, így az oldal közben nem akad meg.

Az „Átvitel az Árkalkulátorba” gomb az eredményt átadja az `index.html`-nek: ott a kezdőlapon és a projektben megjelenik, és egy gombbal betölthető a projekt Lemez soraiba. Ugyanaz az adat, mint az „Export árajánlathoz” JSON fájlban (`SZABAS_EXPORT.md`).

Eltérések az asztali programtól: a böngészőben egy szálon számol (több keresési időt érdemes adni), a megszakítás újraindítja a számítómotort, a DXF funkciók első használatkor letöltik az ezdxf modult. A lézeres export egy ZIP-fájlt ad.
