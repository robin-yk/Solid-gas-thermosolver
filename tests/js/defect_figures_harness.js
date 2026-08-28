/* Build the workspace model from the JS engines and render every defect
   figure. Usage: node defect_figures_harness.js <repo root> [outdir]
   Prints {id: svg, ...} so the gate can measure a figure without a
   browser, and writes the SVG files when an outdir is given. */
'use strict';
const fs = require('fs');
const path = require('path');

const root = process.argv[2] || path.resolve(__dirname, '..', '..');
const outdir = process.argv[3] || null;
const P = JSON.parse(fs.readFileSync(
  path.join(root, 'data', 'rutile_dft.json'), 'utf8'));
const SM = require(path.join(root, 'web', 'statmech.js'));
const CG = require(path.join(root, 'web', 'cgmc.js'));
const FIG = require(path.join(root, 'web', 'figures_defect.js'));

const T_C = P.defaults.T_C;
const vo = P.defaults.VO_total_umol_g;
const geom = SM.geometry(P);
const out = SM.distribute(P, T_C, vo, geom);
const bounds = SM.phaseBoundaries(P, T_C, geom);
const acc = SM.accessibility(P, out.umol_g, null, vo);

const xmax = Math.max(vo * 2.2, bounds.VO_cs_umol_g * 1.15);
const rows = [];
for (let i = 1; i <= 90; i++) {
  const v = xmax * i / 90;
  const d = SM.distribute(P, T_C, v, geom);
  const a = SM.accessibility(P, d.umol_g, null, v);
  rows.push({ VO_total_umol_g: v, surface: d.umol_g.surface,
    subsurface: d.umol_g.subsurface, bulk: d.umol_g.bulk,
    extended: d.extended_defects_umol_g,
    accessible: a.accessible_umol_g });
}

const tEnd = P.co2_kinetics.defaults.exposure_s;
const sys = CG.build(P, T_C, vo, { d_um: P.defaults.particle_diameter_um });
const r = CG.run(sys, tEnd, {});

const M = {
  T_C: T_C, vo: vo, geom: geom, out: out, bounds: bounds, acc: acc,
  sweep: { rows: rows, xmax: xmax },
  cg: { r: r, sys: sys, tEnd: tEnd, E_m_eV: P.kinetics_cgmc.E_m_bulk_eV,
        target_s: tEnd, target_pct: r.recovered_fraction * 100 }
};

const res = { ids: [], svg: {}, errors: {} };
Object.keys(FIG.FIGURES).forEach(id => {
  try {
    const svg = FIG.FIGURES[id](M);
    res.svg[id] = svg;
    res.ids.push(id);
    if (outdir) fs.writeFileSync(path.join(outdir, 'ws_' + id + '.svg'), svg);
  } catch (e) { res.errors[id] = e.message + '\n' + e.stack; }
});
/* export arithmetic, checked without a browser: the pixel size a figure
   declares, and the pHYs chunk the kit splices into a canvas PNG */
const K = require(path.join(root, 'web', 'figkit.js'));
res.px = {};
res.ids.forEach(id => { res.px[id] = K.pxAt600(res.svg[id]); });
res.px_wide = K.pxAt600(K.square({ wide: true }).done());

const ihdr = [0, 0, 0, 13, 73, 72, 68, 82, 0, 0, 8, 52, 0, 0, 8, 52,
              8, 6, 0, 0, 0, 0, 0, 0, 0];
const fake = new Uint8Array([137, 80, 78, 71, 13, 10, 26, 10]
  .concat(ihdr).concat([0, 0, 0, 0, 73, 69, 78, 68, 0, 0, 0, 0]));
const withDpi = K.pngWithDPI(fake.buffer, K.PX_PER_M);
res.phys = {
  grew_by: withDpi.length - fake.length,
  tag: String.fromCharCode(withDpi[37], withDpi[38], withDpi[39],
                           withDpi[40]),
  px_per_m: (withDpi[41] << 24 | withDpi[42] << 16 | withDpi[43] << 8
             | withDpi[44]) >>> 0,
  unit: withDpi[49],
  ihdr_intact: String.fromCharCode(withDpi[12], withDpi[13], withDpi[14],
                                   withDpi[15])
};
process.stdout.write(JSON.stringify(res));
