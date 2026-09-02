/* Run the browser side of the reduction line over the cases the Python
   side supplies, and print JSON. The spec is written by
   tests/test_activeset_boundary.py so neither side can be edited alone.

   Usage: node boundary_harness.js <repo root> <spec.json> */

'use strict';

const fs = require('fs');
const path = require('path');

const root = process.argv[2] || path.resolve(__dirname, '..', '..');
const spec = JSON.parse(fs.readFileSync(process.argv[3], 'utf8'));
const D = JSON.parse(fs.readFileSync(
  path.join(root, 'data', 'activeset_data.json'), 'utf8'));
const AS = require(path.join(root, 'web', 'activeset.js'));
const s = new AS.Solver(D);

global.self = global;
const KIT = require(path.join(root, 'web', 'figkit.js'));
global.FigKit = KIT; self.FigKit = KIT;
const FIG = require(path.join(root, 'web', 'figures_thermo.js'));

const out = { boundary: [], gas: [], equivalent: [], figures: {} };

out.axis = [];
spec.boundary.forEach(function (c) {
  out.boundary.push(s.reductionBoundary(c.T_C + 273.15, c.host, null,
                                        c.P_atm || 1));
});
spec.gas.forEach(function (c) {
  out.gas.push(s.gasOxygenPotential(c.feed, c.T_C + 273.15, c.P_atm || 1));
});
spec.equivalent.forEach(function (c) {
  out.equivalent.push(s.equivalentCO2Fraction(c.mu, c.T_C + 273.15,
                                              c.P_atm || 1));
});
(spec.axis || []).forEach(function (c) {
  out.axis.push(s.feedOnBoundaryAxis(c.feed, c.T_C + 273.15, c.P_atm || 1));
});

/* the two payloads the page hands its figures, so the gate sees the
   numbers that reach the header and not just the engine behind them */
if (spec.figure) {
  const f = spec.figure;
  const rows = [], curve = [];
  for (let i = 0; i < f.n; i++) {
    const T = f.T_lo + (f.T_hi - f.T_lo) * i / (f.n - 1);
    rows.push(s.solve({ feed: f.feed, T_K: T + 273.15, mass_g: f.mass_g }));
    const b = s.reductionBoundary(T + 273.15, 'TiO2');
    curve.push({ T_C: T, y_CO2: b.y_CO2, phase: b.phase });
  }
  const here = s.solve({ feed: f.feed, T_K: f.T_C + 273.15,
                         mass_g: f.mass_g });
  /* exactly what the page hands the figures, through the same helper */
  const pt = s.feedOnBoundaryAxis(f.feed, f.T_C + 273.15, 1);
  pt.T_C = f.T_C;
  pt.host = 'TiO2';
  out.figures.conversion = FIG.conversionData(rows, f.T_C, 'TiO2',
                                              here.conversion_CO2_pct);
  out.figures.boundary = FIG.boundaryData(curve, pt, f.T_C,
    s.reductionBoundary(f.T_C + 273.15, 'TiO2', null, 1));
  out.figures.svg = {
    conversion: FIG.conversion({ conversion: out.figures.conversion }).length,
    boundary: FIG.boundary({ boundary: out.figures.boundary }).length,
  };
}

process.stdout.write(JSON.stringify(out));
