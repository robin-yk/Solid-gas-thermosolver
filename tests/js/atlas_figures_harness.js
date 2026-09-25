'use strict';

const fs = require('fs');
const path = require('path');
const root = process.argv[2];
const T_C = Number(process.argv[3] || 900);
const D = JSON.parse(fs.readFileSync(path.join(root, 'data',
  'activeset_data.json'), 'utf8'));
const AS = require(path.join(root, 'web', 'activeset.js'));
const FIG = require(path.join(root, 'web', 'figures_thermo.js'));
const solver = new AS.Solver(D);

function logPoint(lo, hi, i, n, centred) {
  const f = centred ? (i + 0.5) / n : i / (n - 1);
  return Math.exp(Math.log(lo) + (Math.log(hi) - Math.log(lo)) * f);
}

const T_LO = 400, T_HI = 1500, R_LO = 0.1, R_HI = 1e7;
const NT = 44, NR = 36;
const cells = [], boundary = [];
for (let i = 0; i < NT; i++) {
  const tc = T_LO + (T_HI - T_LO) * (i + 0.5) / NT;
  for (let j = 0; j < NR; j++) {
    const ratio = logPoint(R_LO, R_HI, j, NR, true);
    cells.push({ T_C: tc, ratio: ratio,
      conv_pct: solver.rwgsGasOnly(tc + 273.15, ratio).CO2_conversion_pct });
  }
}
for (let i = 0; i < 181; i++) {
  const tc = T_LO + (T_HI - T_LO) * i / 180;
  const b = solver.reductionBoundary(tc + 273.15, 'TiO2');
  boundary.push({ T_C: tc, ratio: (1 - b.y_CO2) / b.y_CO2,
                  phase: b.phase });
}
const exact = solver.rwgsGasOnly(T_C + 273.15, 1);
const map = { T_lo: T_LO, T_hi: T_HI, ratio_lo: R_LO, ratio_hi: R_HI,
  nT: NT, nR: NR, cells: cells, boundary: boundary, host: 'TiO2',
  current: { T_C: T_C, ratio: 1, conv_pct: exact.CO2_conversion_pct,
             literal: true } };

const trade = [];
for (let i = 0; i < 121; i++) {
  const ratio = logPoint(0.1, 10, i, 121, false);
  const q = solver.rwgsGasOnly(T_C + 273.15, ratio);
  trade.push({ ratio: ratio, co2: q.CO2_conversion_pct,
               h2: q.H2_utilization_pct,
               co: q.CO_yield_per_mol_feed_pct });
}
const coupling = [];
for (let i = 0; i < 49; i++) {
  const ratio = logPoint(R_LO, R_HI, i, 49, false);
  const gas = solver.rwgsGasOnly(T_C + 273.15, ratio);
  const both = solver.solve({ feed: { CO2: 1, H2: ratio },
                              T_K: T_C + 273.15 });
  coupling.push({ ratio: ratio, gas_only: gas.CO2_conversion_pct,
                   gas_solid: both.conversion_CO2_pct,
                   reduced_pct: both.reduced_pct,
                   phases: both.active_condensed_phases });
}
const b = solver.reductionBoundary(T_C + 273.15, 'TiO2');
const ratioBoundary = (1 - b.y_CO2) / b.y_CO2;
const payload = {
  map: map,
  trade: { T_C: T_C, ratio_lo: 0.1, ratio_hi: 10, rows: trade,
           current: { ratio: 1, co2: exact.CO2_conversion_pct } },
  coupling: { T_C: T_C, ratio_lo: R_LO, ratio_hi: R_HI,
              rows: coupling, boundary_ratio: ratioBoundary, host: 'TiO2' }
};
const svg = {
  map: FIG.operatingMap({ operatingMap: payload.map }),
  trade: FIG.ratioTradeoff({ ratioTradeoff: payload.trade }),
  coupling: FIG.solidCoupling({ solidCoupling: payload.coupling })
};

process.stdout.write(JSON.stringify({
  payload: payload,
  bytes: Object.fromEntries(Object.entries(svg).map(([k, v]) => [k, v.length])),
  viewBoxes: Object.fromEntries(Object.entries(svg).map(([k, v]) =>
    [k, (v.match(/viewBox="([^"]+)"/) || [])[1]])),
  text: Object.fromEntries(Object.entries(svg).map(([k, v]) => [k,
    (v.match(/<text[^>]*>([^<]*)<\/text>/g) || [])
      .map(t => t.replace(/<[^>]+>/g, ''))])),
}));
