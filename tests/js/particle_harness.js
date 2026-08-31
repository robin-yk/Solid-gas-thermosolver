/* Run the browser particle engine over the cases the Python side supplies
   and print JSON.  Usage: node particle_harness.js <repo root> <spec.json>
   The spec is written by tests/test_particle_port.py so the two sides
   cannot drift apart by editing only one of them. */

'use strict';

const fs = require('fs');
const path = require('path');

const root = process.argv[2] || path.resolve(__dirname, '..', '..');
const spec = JSON.parse(fs.readFileSync(process.argv[3], 'utf8'));
const P = JSON.parse(fs.readFileSync(
  path.join(root, 'data', 'rutile_dft.json'), 'utf8'));
const PA = require(path.join(root, 'web', 'particle.js'));

const out = { phi: [], inverse: [], segregation: {}, equilibrium: [],
              crank: [], shells: [], spacecharge: [], transport: [] };

for (const c of spec.phi) out.phi.push(PA.phi(c.theta, c.ratio));
for (const c of spec.inverse) {
  out.inverse.push(PA.thetaOfPhi(c.x, c.compensated, c.ratio));
}
for (const key of spec.presets) {
  const s = PA.segregation(P, key);
  out.segregation[key] = {
    dE_seg_eV: s.dE_seg_eV, xi_nm: s.xi_nm, d110_nm: s.d110_nm,
    E: spec.depth_layers.map((n) => PA.formationEnergy(n * s.d110_nm, s))
  };
}
for (const c of spec.equilibrium) {
  const part = PA.particle(P, { d_um: c.d_um, preset: c.preset });
  const r = PA.equilibrate(c.vo, c.T_C, part, {
    ratio: c.ratio, compensated: c.compensated,
    reconstruction: c.reconstruction, shear_limit: c.shear_limit });
  out.equilibrium.push({
    mu_eV: r.mu_eV, theta_surface: r.surface.theta,
    surface_phase: r.surface.phase,
    reconstructed_fraction: r.surface.reconstructed_fraction,
    theta_interior_mean: r.theta_interior_mean,
    theta_bulk_plateau: r.theta_bulk_plateau,
    theta_first_layer: r.theta_first_layer,
    umol_surface: r.umol_surface, umol_interior: r.umol_interior,
    umol_extended: r.umol_extended, limited_by: r.limited_by,
    point_defect_ceiling_umol_g: r.point_defect_ceiling_umol_g,
    N_surface_umol_g: part.N_surface_umol_g,
    N_interior_umol_g: part.N_interior_umol_g,
    radius_nm: part.radius_nm
  });
  out.shells.push(PA.shellSplit(r, part, spec.shell_depths_nm)
    .map((q) => q.inventory_fraction));
  out.spacecharge.push(PA.assess(r, part, P).ratio_debye_over_xi);
  out.transport.push(PA.verdict(c.T_C, part, spec.exposure_s, P)
    .E_m_that_would_break_it_eV);
}
for (const f of spec.crank) out.crank.push(PA.tauForUptake(f));

process.stdout.write(JSON.stringify(out));
