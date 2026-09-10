'use strict';

const fs = require('fs');
const path = require('path');
const root = process.argv[2];
const D = JSON.parse(fs.readFileSync(path.join(root, 'data',
  'activeset_data.json'), 'utf8'));
const REF = JSON.parse(fs.readFileSync(path.join(root, 'data',
  'reference_rwgs_gas_only.json'), 'utf8'));
const AS = require(path.join(root, 'web', 'activeset.js'));
const solver = new AS.Solver(D);

const rows = REF.rows.map(function (q) {
  return solver.rwgsGasOnly(Number(q.T_C) + 273.15,
                            Number(q.H2_CO2_ratio));
});
process.stdout.write(JSON.stringify(rows));
