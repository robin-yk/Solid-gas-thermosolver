/* Run web/population.js over the parameter grid the parity gate asks for. */
'use strict';
const fs = require('fs');
const ROOT = process.argv[2];
const spec = JSON.parse(process.argv[3]);
const P = require(ROOT + '/web/population.js');
const D = JSON.parse(fs.readFileSync(ROOT + '/paper_outputs/site_data.json', 'utf8'));
const series = D.population.series;

const out = { reported: P.REPORTED, cases: [] };
for (const c of spec) {
  const rows = P.partition(series, c.ns, c.E_loss, c.nu_eff);
  out.cases.push({
    label: c.label, ns: c.ns, E_loss: c.E_loss, nu_eff: c.nu_eff,
    rows: rows.map(q => ({
      sample: q.sample, n_iso: q.n_iso_umol_g, n_assoc: q.n_assoc_umol_g,
      n_below: q.n_below_umol_g, closure: q.closure_umol_g,
      r_CO_fitted: q.r_CO_fitted })),
    objective: {
      log: P.objective(series, c.ns, c.E_loss, c.nu_eff, 'log'),
      relative: P.objective(series, c.ns, c.E_loss, c.nu_eff, 'relative'),
      raw: P.objective(series, c.ns, c.E_loss, c.nu_eff, 'raw') },
  });
}
out.fit = P.fit(series, 'log');
const threw = f => { try { f(); return false; } catch (e) { return true; } };
out.rejects_zero_inventory = threw(() => P.integrate(0, 1800, 873.15, 15.3, 2.18, 5.5e6));
out.rejects_zero_capacity = threw(() => P.integrate(94, 1800, 873.15, 0, 2.18, 5.5e6));
console.log(JSON.stringify(out));
