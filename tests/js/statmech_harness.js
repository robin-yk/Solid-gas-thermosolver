/* Run the JS stat-mech engine over the cases supplied by the Python side
   and print JSON. Usage: node statmech_harness.js <repo root> <spec.json>
   The spec (lattice sizes, seeds, sweep counts) is written by
   tests/test_statmech_port.py so the two sides cannot drift. */

'use strict';

const fs = require('fs');
const path = require('path');

const root = process.argv[2] || path.resolve(__dirname, '..', '..');
const spec = JSON.parse(fs.readFileSync(process.argv[3], 'utf8'));
const P = JSON.parse(fs.readFileSync(path.join(root, 'data', 'rutile_dft.json'), 'utf8'));
const SM = require(path.join(root, 'web', 'statmech.js'));

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
  const eps = SM.energetics(P).eps;
  const sw = SM.sweep(P, c.T_C, res.geometry, spec.sweep_vo, { eps: eps });
  out.pipelines.push({
    name: c.name,
    mu_V_eV: res.mu_V_eV,
    fractions: res.fractions,
    umol_g: res.umol_g,
    theta: res.theta,
    regime: res.regime,
    surface_phase: res.surface_phase,
    extended_defects_umol_g: res.extended_defects_umol_g,
    phase_boundaries: res.phase_boundaries,
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
  const AS = require(path.join(root, 'web', 'activeset.js'));
  const AD = JSON.parse(fs.readFileSync(
    path.join(root, 'data', 'activeset_data.json'), 'utf8'));
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

out.analytic_flagship = (() => {
  const g = SM.geometry(P, {});
  const d = SM.distribute(P, 600.0, 95.0, g, {});
  return { mu_V_eV: d.mu_V_eV, fractions: d.fractions, umol_g: d.umol_g,
           regime: d.regime, surface_phase: d.surface_phase,
           extended_defects_umol_g: d.extended_defects_umol_g,
           phase_boundaries: d.phase_boundaries,
           warn_kinds: d.warnings.map(w => w.kind),
           access: SM.accessibility(P, d.umol_g, null, 95.0) };
})();

if (spec.ladder_vo) {
  const g = SM.geometry(P, {});
  out.ladder = spec.ladder_vo.map(vo => {
    const d = SM.distribute(P, 600.0, vo, g, {});
    return { vo: vo, mu_V_eV: d.mu_V_eV, regime: d.regime,
             umol_g: d.umol_g, extended: d.extended_defects_umol_g,
             phi: d.surface_phase.phi_reconstructed,
             warn_kinds: d.warnings.map(w => w.kind) };
  });
}

if (spec.shell_layers) {
  out.shell = spec.shell_layers.map(n => {
    const g = SM.geometry(P, { ss_layers: n });
    const d = SM.distribute(P, 600.0, 95.0, g, {});
    return { ss_layers: n, geometry: g, mu_V_eV: d.mu_V_eV,
             umol_g: d.umol_g, theta: d.theta, regime: d.regime,
             fractions: d.fractions,
             warn_kinds: d.warnings.map(w => w.kind) };
  });
}

/* Stage 1: the layer-resolved free energy, in every combination of
   resolution and interaction, so the port is compared on the branches the
   dataset ships switched off as well as the one it ships on. */
out.layers = [];
[1, 2, 3].forEach(n => ['off', 'on'].forEach(cp => {
  const g = SM.geometry(P, { resolved_layers: n });
  [1.0, 2.5, 5.0, 95.0, 400.0].forEach(vo => {
    const d = SM.distribute(P, 600.0, vo, g, { omega: cp });
    out.layers.push({
      n: n, omega: cp, vo: vo, mu_V_eV: d.mu_V_eV, regime: d.regime,
      N_layers: g.N_layers, shell_nm: g.shell_nm,
      theta: d.theta, umol_g: d.umol_g,
      extended: d.extended_defects_umol_g,
      layer_theta: d.layers.map(l => l.theta),
      layer_umol: d.layers.map(l => l.umol_g),
      layer_eps: d.layers.map(l => l.eps_eV),
      layer_omega: d.layers.map(l => l.omega_eV),
      boundaries: d.phase_boundaries
    });
  });
}));

/* the implicit occupancy on its own, away from any partition */
out.theta_of_mu = [];
[0.0, 0.15, 0.6, 1.0].forEach(w => {
  [-0.4, 0.0, 0.5, 0.85, 1.2].forEach(mu => {
    [0.0, 0.59, 1.31].forEach(eps => {
      out.theta_of_mu.push({ w: w, mu: mu, eps: eps,
        theta: SM.thetaOfMu(mu, eps, w, 8.617333262e-5 * 873.15) });
    });
  });
});

if (spec.cgmc) {
  const CG = require(path.join(root, 'web', 'cgmc.js'));
  const c = spec.cgmc;
  out.cgmc = {
    linear: c.times.map(t => CG.run(c.sys, t, { n_out: 8 }).final_eta),
    ladder: c.barriers.map(em => {
      const r = CG.run(CG.build(P, { cg: { E_m_bulk_eV: em } }), 300.0, {});
      return { recovery: r.recovered_fraction, steps: r.steps,
               mass: r.mass_check };
    }),
    closed: (() => {
      const sys = CG.build(P, { kin: c.no_fill, dist_fractions: c.all_bulk });
      const r = CG.run(sys, 300.0, {});
      return { final_eta: r.final_eta, steps: r.steps };
    })(),
    stoch: (() => {
      const sys = CG.build(P, { kin: c.no_fill, cg: c.stoch_cg });
      return CG.run(sys, 60.0, { mode: 'stoch', seed: 902 }).final_eta;
    })(),
  };
}

process.stdout.write(JSON.stringify(out));
