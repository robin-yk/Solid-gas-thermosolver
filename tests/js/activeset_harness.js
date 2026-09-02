/* Run the JS active-set engine over every reference point and print JSON.
   Usage: node activeset_harness.js <repo root>
   The case list and the boundary H2O/H2 ratio come from the high-precision
   reference file, so the harness cannot drift from the oracle. */

'use strict';

const fs = require('fs');
const path = require('path');

const root = process.argv[2] || path.resolve(__dirname, '..', '..');
const D = JSON.parse(fs.readFileSync(path.join(root, 'data', 'activeset_data.json'), 'utf8'));
const REF = JSON.parse(fs.readFileSync(
  path.join(root, 'data', 'reference_results_high_precision.json'), 'utf8'));
const AS = require(path.join(root, 'web', 'activeset.js'));

const solver = new AS.Solver(D);

const CASES = {
  rwgs_1_1: { CO2: 1, H2: 1 },
  h2_rich: { CO2: 1, H2: 3 },
  product_mix: { CO2: 3, H2: 3, CO: 2, H2O: 2 },
  pure_h2: { H2: 1 },
  pure_n2: { N2: 1 },
  n2_10ppm_o2: { N2: 999990, O2: 10 },
};

const rows = [];
for (const [name, feed] of Object.entries(CASES)) {
  for (const T_C of REF.temperatures_C) {
    const r = solver.solve({ feed: feed, T_K: T_C + 273.15 });
    r.case = name;
    r.T_C_label = T_C;
    rows.push(r);
  }
}

/* Boundary flip points: the ratio comes from the oracle's boundary block. */
const bb = REF.boundary_validation.h2_h2o_boundary;
const ratio = Number(bb.H2O_over_H2_critical_ratio);
for (const [tag, fac] of [['h2_h2o_below_boundary', 1 - 1e-3],
                          ['h2_h2o_above_boundary', 1 + 1e-3]]) {
  const r = solver.solve({ feed: { H2: 1, H2O: ratio * fac }, T_K: 1173.15 });
  r.case = tag;
  r.T_C_label = 900;
  rows.push(r);
}

/* Permutation-invariance probe: same condition, shuffled solid list. */
const shuffled = ['Ti4O7', 'TiO2', 'Ti2O3', 'Ti10O19', 'Ti6O11',
                  'Ti8O15', 'Ti3O5', 'Ti5O9'];
const a = solver.solve({ feed: CASES.pure_n2, T_K: 1173.15 });
const b = solver.solve({ feed: CASES.pure_n2, T_K: 1173.15, solids: shuffled });
rows.push({
  case: 'permutation_probe', T_C_label: 900,
  identical: JSON.stringify(a.solid_moles) === JSON.stringify(b.solid_moles)
    && JSON.stringify(a.active_condensed_phases)
       === JSON.stringify(b.active_condensed_phases)
    && a.G_total_kJ === b.G_total_kJ,
});

process.stdout.write(JSON.stringify(rows));
