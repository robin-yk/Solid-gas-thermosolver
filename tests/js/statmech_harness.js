/* Run the JS stat-mech engine over the cases supplied by the Python side
   and print JSON. Usage: node statmech_harness.js <repo root> <spec.json>
   The spec (lattice sizes, seeds, sweep counts) is written by
   tests/test_statmech_port.py so the two sides cannot drift. */

'use strict';

const fs = require('fs');
const path = require('path');

const root = process.argv[2] || path.resolve(__dirname, '..', '..');
const spec = JSON.parse(fs.readFileSync(process.argv[3], 'utf8'));
const P = JSON.parse(fs.readFileSync(path.join(root, 'rutile_dft.json'), 'utf8'));
const SM = require(path.join(root, 'statmech.js'));

const out = { rng: [], mc_cases: [], pipelines: [] };

const rng = SM.makeRng(12345);
for (let i = 0; i < 8; i++) out.rng.push(rng());

for (const c of spec.mc_cases) {
  const lat = SM.buildLattice(c.rows, c.row_sites, c.ss_layers,
                              c.bulk_layers, c.ordering);
  const r = SM.makeRng(c.seed);
  const occ = new Uint8Array(lat.n);
  const vac = [];
  SM.adjustFilling(lat, occ, vac, c.n_vac, r);
  const kt = SM.KB_EV * (c.T_C + 273.15);
  const res = SM.mcRun(lat, c.eps, kt, occ, vac, c.sweeps_eq,
                       c.sweeps_sample, r,
                       { widom_every: c.widom_every,
                         hist_every: c.hist_every, blocks: c.blocks });
  res.name = c.name;
  out.mc_cases.push(res);
}

for (const c of spec.pipelines) {
  const res = SM.runFull(P, { T_C: c.T_C, VO_total: c.VO_total,
                              zero_pairs: !!c.zero_pairs, mc: c.mc });
  out.pipelines.push({
    name: c.name,
    mu_V_eV: res.mu_V_eV,
    fractions: res.fractions,
    umol_g: res.umol_g,
    theta: res.theta,
    surface_cap_bound: res.surface_cap_bound,
    warn_kinds: res.warnings.map(w => w.kind),
    spacing: res.spacing,
    isotherm: res.isotherm.map(q => ({
      filling: q.filling, theta_s: q.theta_s, theta_ss: q.theta_ss,
      theta_b: q.theta_b, mu_eV: q.mu_eV, E_mean_eV: q.E_mean_eV,
      E_drift_eV: q.E_drift_eV, hist: q.hist })),
  });
}

process.stdout.write(JSON.stringify(out));
