// Acél árkalkulátor – közös projekttár Google táblázatban.
// A táblázatban: Bővítmények > Apps Script, ide bemásolni, majd Telepítés > Új telepítés > Webes alkalmazás
// (Végrehajtás: én; Hozzáférés: bárki). A kapott /exec URL kell a kalkulátornak.
const SHEETS = { proj: 'Projektek', tpl: 'Termektipusok' };
const CHUNK = 45000; // egy cella max. 50 000 karakter

function sheet_(kind) {
  if (!SHEETS[kind]) throw new Error('ismeretlen kind: ' + kind);
  const ss = SpreadsheetApp.getActive();
  let sh = ss.getSheetByName(SHEETS[kind]);
  if (!sh) { sh = ss.insertSheet(SHEETS[kind]); sh.appendRow(['id', 'név', 'módosítva', 'adat (JSON)']); sh.setFrozenRows(1); }
  return sh;
}
function out_(o) { return ContentService.createTextOutput(JSON.stringify(o)).setMimeType(ContentService.MimeType.JSON); }
// a darabok elé "~" kerül, hogy a táblázat ne értelmezze képletnek
function readAll_(kind) {
  const v = sheet_(kind).getDataRange().getValues();
  return v.slice(1).filter((r) => r[0] !== '').map((r) => JSON.parse(r.slice(3).filter((c) => c !== '').map((c) => String(c).slice(1)).join('')));
}

function doGet(e) {
  try { return out_({ ok: true, items: readAll_(e.parameter.kind) }); }
  catch (err) { return out_({ ok: false, error: String(err) }); }
}

function doPost(e) {
  const lock = LockService.getScriptLock();
  lock.waitLock(20000);
  try {
    const b = JSON.parse(e.postData.contents);
    const sh = sheet_(b.kind);
    const last = sh.getLastRow();
    const ids = last > 1 ? sh.getRange(2, 1, last - 1, 1).getValues().map((r) => String(r[0])) : [];
    let row = ids.indexOf(String(b.id)) + 2; // 1 = nincs meg
    if (b.action === 'remove') { if (row > 1) sh.deleteRow(row); return out_({ ok: true }); }
    if (b.action !== 'save') throw new Error('ismeretlen action');
    const obj = Object.assign({}, b.data, { id: b.id });
    const s = JSON.stringify(obj), parts = [];
    for (let i = 0; i < s.length; i += CHUNK) parts.push('~' + s.slice(i, i + CHUNK));
    const vals = ["'" + b.id, "'" + (obj.nev || ''), "'" + (obj.modositva || ''), ...parts];
    if (row < 2) row = sh.getLastRow() + 1;
    else if (sh.getLastColumn() > 3) sh.getRange(row, 4, 1, sh.getLastColumn() - 3).clearContent();
    sh.getRange(row, 1, 1, vals.length).setValues([vals]);
    return out_({ ok: true });
  } catch (err) {
    return out_({ ok: false, error: String(err) });
  } finally {
    lock.releaseLock();
  }
}
