// A szabásszámítás a böngészőben: Pyodide (Python WebAssemblyben) egy Web Workerben,
// így a hosszú számítás alatt sem akad meg az oldal. A Python kód ugyanaz, mint az
// asztali programé (kontener_szabas.py), a web.py csak az adatokat adja át.
const PYODIDE = 'https://cdn.jsdelivr.net/pyodide/v0.27.7/full/';
importScripts(PYODIDE + 'pyodide.js');

const ver = new URL(self.location.href).searchParams.get('v') || '';
let py = null, web = null, ezdxf = null;

const status = (msg) => postMessage({ type: 'status', msg });

const ready = (async () => {
  status('Python betöltése…');
  py = await loadPyodide({ indexURL: PYODIDE });
  status('A szabásprogram betöltése…');
  for (const f of ['kontener_szabas.py', 'web.py']) {
    const r = await fetch(f + (ver ? '?v=' + ver : ''), { cache: 'no-cache' });
    if (!r.ok) throw new Error(f + ' nem tölthető le (' + r.status + ')');
    py.FS.writeFile('/home/pyodide/' + f, await r.text());
  }
  web = py.pyimport('web');
  status('');
})();

// az ezdxf csak a DXF-es funkciókhoz kell, ezért csak akkor töltjük le
function needEzdxf() {
  if (!ezdxf) {
    ezdxf = (async () => {
      status('DXF-modul (ezdxf) letöltése…');
      await py.loadPackage('micropip');
      await py.runPythonAsync("import micropip\nawait micropip.install('ezdxf==1.4.5')");
      status('');
    })();
    ezdxf.catch(() => { ezdxf = null; });
  }
  return ezdxf;
}

// a DXF rajzok (egyéb alkatrészek) a Pyodide fájlrendszerébe
function writeDxfs(dxfs) {
  py.FS.mkdirTree('/work/dxf');
  for (const [path, bytes] of Object.entries(dxfs || {})) py.FS.writeFile(path, bytes);
}

const ops = {
  info: () => JSON.parse(web.info()),
  sketch_calc: ({ key, values, flange }) => JSON.parse(web.sketch_calc(key, JSON.stringify(values), flange)),
  material: ({ text }) => JSON.parse(web.material_info(text)),
  async dxf_load({ path, bytes, name }) {
    await needEzdxf();
    writeDxfs({ [path]: bytes });
    return JSON.parse(web.dxf_summary(path, name));
  },
  async calc({ params, dxfs }) {
    if (dxfs && Object.keys(dxfs).length) await needEzdxf();
    writeDxfs(dxfs);
    const progress = (msg) => postMessage({ type: 'progress', msg });
    return JSON.parse(web.calculate(JSON.stringify(params), progress));
  },
  async export_dxf() {
    await needEzdxf();
    status('DXF készítése…');
    const path = web.export_dxf();
    status('');
    return py.FS.readFile(path);
  },
  async export_laser() {
    await needEzdxf();
    status('Lézeres DXF fájlok készítése…');
    const r = JSON.parse(web.export_laser());
    status('');
    return { count: r.count, bytes: py.FS.readFile(r.path) };
  },
};

onmessage = async (e) => {
  const { id, op, args } = e.data;
  try {
    await ready;
    const result = await ops[op](args || {});
    postMessage({ type: 'result', id, ok: true, result });
  } catch (err) {
    status('');
    // a Python hibaüzenet utolsó sora a lényeg
    const msg = String(err && err.message || err).trim().split('\n').filter(Boolean).pop() || 'Ismeretlen hiba';
    postMessage({ type: 'result', id, ok: false, error: msg });
  }
};

ready.then(() => postMessage({ type: 'ready' }),
  (err) => postMessage({ type: 'fatal', error: String(err && err.message || err) }));
