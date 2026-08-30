/* Build the two equilibrium-workspace figure payloads in the browser
   modules over a grid of conditions and print them as JSON, so the
   Python gate can hold them against solidgas/activeset.py.

   The payload builders only reshape a solver result, so this harness is
   the whole of what the page adds between the engine and the drawing.

   Usage: node thermo_figures_harness.js <repo root> */
'use strict';
const fs = require('fs');
const path = require('path');

const root = process.argv[2] || path.resolve(__dirname, '..', '..');
const K = require(path.join(root, 'web', 'figkit.js'));
global.self = { FigKit: K };
const F = require(path.join(root, 'web', 'figures_thermo.js'));
const AS = require(path.join(root, 'web', 'activeset.js'));
const D = JSON.parse(fs.readFileSync(
  path.join(root, 'data', 'activeset_data.json'), 'utf8'));
const solver = new AS.Solver(D);

/* the panel's presets, plus a host that is already reduced and a feed
   with no oxygen in it at all */
const CASES = [
  { name: 'rwgs', feed: { CO2: 1, H2: 1 }, mass_g: 0.100, T_C: 900,
    initial_solid: 'TiO2' },
  { name: 'magneli', feed: { H2: 1 }, mass_g: 0.002, T_C: 1200,
    initial_solid: 'TiO2' },
  { name: 'reox', feed: { CO2: 1 }, mass_g: 0.100, T_C: 600,
    initial_solid: 'Ti4O7' },
  { name: 'h2co2', feed: { H2: 99, CO2: 1 }, mass_g: 0.002, T_C: 1000,
    initial_solid: 'TiO2' },
  { name: 'n2', feed: { N2: 1 }, mass_g: 0.100, T_C: 1400,
    initial_solid: 'Ti3O5' }
];
const SWEEP_N = 41;                 /* the page runs 121; 41 is the gate */
const T_LO = 300, T_HI = 1650;

function opts(c, T_C) {
  return { feed: c.feed, T_K: T_C + 273.15, P_atm: 1, V_cm3: 25,
           mass_g: c.mass_g, initial_solid: c.initial_solid };
}

const out = { sweep_n: SWEEP_N, T_lo: T_LO, T_hi: T_HI, cases: {} };
for (const c of CASES) {
  const one = solver.solve(opts(c, c.T_C));
  const rows = [];
  for (let i = 0; i < SWEEP_N; i++) {
    rows.push(solver.solve(opts(c, T_LO + (T_HI - T_LO) * i / (SWEEP_N - 1))));
  }
  out.cases[c.name] = {
    T_C: c.T_C,
    optimality: F.optimalityData(one),
    margin: F.marginData(rows, c.T_C)
  };
}
process.stdout.write(JSON.stringify(out));
