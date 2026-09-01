/* Build every figure the page draws, exactly the way the page builds it,
   and report what came out. No solver runs here: the payloads come from
   the committed site_data.json, which is what the page is handed too. */
'use strict';
const fs = require('fs');
const ROOT = process.argv[2];
const TF = require(ROOT + '/web/figures_thermo.js');
const SF = require(ROOT + '/web/figures_surface.js');
const D = JSON.parse(fs.readFileSync(ROOT + '/paper_outputs/site_data.json', 'utf8'));

const T = Number(process.argv[3] || 600);
const rows = D.equilibrium.rows;
const at = D.equilibrium.boundary.reduce((a, b) =>
  Math.abs(b.T_C - T) < Math.abs(a.T_C - T) ? b : a);
const r = rows.reduce((a, b) =>
  Math.abs(b.T_C - T) < Math.abs(a.T_C - T) ? b : a);

const conv = TF.conversionData(rows, T, 'TiO2', r.conversion_CO2_pct);
const marg = TF.marginData(rows, T);
const bnd = TF.boundaryData(D.equilibrium.boundary, null, T, at);
const cap = SF.capacityData(D.surface.rows, 0.9, 95);
const bd = SF.boundData(D.surface.rows, 0.9, 95);

const svg = {
  conversion: TF.conversion({ conversion: conv }),
  margin: TF.margin({ margin: marg }),
  boundary: TF.boundary({ boundary: bnd }),
  capacity: SF.capacity(cap),
  bound: SF.bound(bd),
};

console.log(JSON.stringify({
  T_C: T,
  payload: { conversion: conv, margin: marg, boundary: bnd,
             capacity: cap, bound: bd },
  bytes: Object.fromEntries(Object.entries(svg).map(([k, v]) => [k, v.length])),
  viewBoxes: Object.fromEntries(Object.entries(svg).map(([k, v]) =>
    [k, (v.match(/viewBox="([^"]+)"/) || [])[1]])),
  texts: Object.fromEntries(Object.entries(svg).map(([k, v]) =>
    [k, (v.match(/<text[^>]*>([^<]*)<\/text>/g) || [])
        .map(t => t.replace(/<[^>]+>/g, ''))])),
}));
