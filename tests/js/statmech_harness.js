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
  const kt = SM.KB_EV * (c.T_C + 273.15);
  const eps = SM.energetics(P).eps;
  const fn = SM.thetaSInterp(res.isotherm, eps[0], kt);
  const sw = SM.sweep(P, c.T_C, res.geometry, spec.sweep_vo,
                      { theta_s_fn: fn, eps: eps });
  out.pipelines.push({
    name: c.name,
    mu_V_eV: res.mu_V_eV,
    fractions: res.fractions,
    umol_g: res.umol_g,
    theta: res.theta,
    surface_reconstruction_regime: res.surface_reconstruction_regime,
    warn_kinds: res.warnings.map(w => w.kind),
    spacing: res.spacing,
    accessibility: res.accessibility,
    sweep: sw,
    isotherm: res.isotherm.map(q => ({
      filling: q.filling, theta_s: q.theta_s, theta_ss: q.theta_ss,
      theta_b: q.theta_b, mu_eV: q.mu_eV, E_mean_eV: q.E_mean_eV,
      E_drift_eV: q.E_drift_eV, hist: q.hist,
      surf_vac: q.surf_vac, rows: q.rows, row_sites: q.row_sites })),
  });
}

if (spec.validity_cases && spec.validity_cases.length) {
  const AS = require(path.join(root, 'activeset.js'));
  const AD = JSON.parse(fs.readFileSync(
    path.join(root, 'activeset_data.json'), 'utf8'));
  const thermo = new AS.Solver(AD);
  const isRed = {};
  AD.solids.forEach(s => { isRed[s] = AD.ti3[s] > 0; });
  out.validity = [];
  for (const c of spec.validity_cases) {
    const R = thermo.solve({ feed: c.feed, T_K: c.T_K, mass_g: c.mass_g });
    const v = SM.phaseValidity(R.active_condensed_phases,
                               R.inactive_phase_reduced_costs, isRed);
    out.validity.push({ name: c.name, tier: v.tier,
      nearest_phase: v.nearest_phase,
      margin_per_mol_O_kJ: v.margin_per_mol_O_kJ,
      reduced_active: v.reduced_active,
      mu_O_kJ_mol: R.lambda_kJ_per_mol.O,
      active: R.active_condensed_phases });
  }
}

out.dielectric = [95, 10, 1000, 0].map(v => SM.dielectric(P, v));

process.stdout.write(JSON.stringify(out));
