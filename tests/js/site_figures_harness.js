/* Build every figure the page draws, exactly the way the page builds it,
   and report what came out. No solver runs here: the payloads come from
   the committed site_data.json, which is what the page is handed too. */
'use strict';
const fs = require('fs');
const ROOT = process.argv[2];
const TF = require(ROOT + '/web/figures_thermo.js');
const PF = require(ROOT + '/web/figures_population.js');
const PM = require(ROOT + '/web/population.js');
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
const R = D.population.reported;
const prows = PM.partition(D.population.series, R.ns_umol_g, R.E_loss_eV,
                           R.nu_eff_s1);
const pop = PF.partitionData(prows);
const rat = PF.ratesData(prows);

const svg = {
  conversion: TF.conversion({ conversion: conv }),
  margin: TF.margin({ margin: marg }),
  boundary: TF.boundary({ boundary: bnd }),
  population: PF.population(pop),
  rates: PF.rates(rat),
};

console.log(JSON.stringify({
  T_C: T,
  payload: { conversion: conv, margin: marg, boundary: bnd,
             population: pop, rates: rat },
  bytes: Object.fromEntries(Object.entries(svg).map(([k, v]) => [k, v.length])),
  viewBoxes: Object.fromEntries(Object.entries(svg).map(([k, v]) =>
    [k, (v.match(/viewBox="([^"]+)"/) || [])[1]])),
  texts: Object.fromEntries(Object.entries(svg).map(([k, v]) =>
    [k, (v.match(/<text[^>]*>([^<]*)<\/text>/g) || [])
        .map(t => t.replace(/<[^>]+>/g, ''))])),
  /* one entry per drawn marker, so a gate can count the scatter rather
     than trust that it drew */
  marks: Object.fromEntries(Object.entries(svg).map(([k, v]) => [k, {
    circle: (v.match(/<circle[^>]*stroke="#0072B2"/g) || []).length,
    square: (v.match(/<rect[^>]*stroke="#D55E00"/g) || []).length,
    triangle: (v.match(/class="marker"[^>]*stroke="#777777"/g) || []).length,
  }])),
}));
