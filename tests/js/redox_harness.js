/* Evaluate the browser redox and cycling engines over a grid and print
   the rows as JSON, so the Python gate can hold them against
   solidgas/redox.py and solidgas/cycling.py.

   Usage: node redox_harness.js <repo root> */
'use strict';
const fs = require('fs');
const path = require('path');

const root = process.argv[2] || path.resolve(__dirname, '..', '..');
const P = JSON.parse(fs.readFileSync(
  path.join(root, 'data', 'redox.json'), 'utf8'));
const R = require(path.join(root, 'web', 'redox.js'));
const CY = require(path.join(root, 'web', 'cycling.js'));

const st = P.structure;
const thSol = st.x_solubility / st.oxygen_per_formula;
const thMax = st.theta_surface_max != null ? st.theta_surface_max : 0.5;
const bands = P.regime_thresholds;

/* K over eight decades, and the enrichment axis the uniform model pins
   at 1 - including the two values the phase map is labelled with */
const KS = [];
for (let i = 0; i <= 64; i++) KS.push(Math.pow(10, -4 + 8 * i / 64));
const ES = [1, 10, 124, 125, 126, 1844, 1e4];

const steady = [];
for (const E of ES) {
  for (const K of KS) {
    const s = R.steady(K, { theta_sol: thSol, theta_surface_max: thMax,
                            enrichment: E, bands: bands });
    steady.push({ K: K, E: E, theta: s.theta,
                  theta_surface: s.theta_surface,
                  rate_over_k_ox: s.rate_over_k_ox, regime: s.regime,
                  beyond_solubility: s.beyond_solubility,
                  no_steady_state: s.no_steady_state, usable: s.usable });
  }
}

const windows = ES.map(function (E) {
  const w = R.usableWindow(E, thSol, thMax);
  return { E: E, K_max: w.K_max, K_second_phase: w.K_second_phase,
           K_no_steady_state: w.K_no_steady_state, binding: w.binding };
});

/* the two figure payloads, so the gate can check the page draws the same
   rows the exporter froze */
const crossing = R.crossingData(0.5, 201);
const phase = R.phaseData({ theta_sol: thSol, theta_surface_max: thMax });

/* cycling: the measured pair, a synthetic five-cycle round trip, and the
   discriminating experiment */
const measured = CY.fit([{ created: 94.8, recovered: 81.5 },
                         { created: 61.0, recovered: 51.6 }]);
const synth = CY.run([100, 80, 60, 50, 40], 0.7, 0.35);
const back = CY.fit(synth.map(function (r) {
  return { created: r.created_umol_g, recovered: r.recovered_umol_g }; }));
const two = CY.twoSamplePrediction(155.8, measured.f_rec, measured.lock, 2);
const two4 = CY.twoSamplePrediction(155.8, 0.7, 0.35, 4);

/* JSON has no infinity and K_second_phase is genuinely unbounded once
   the enrichment carries the face past the solubility, so it travels
   as a string rather than as the null JSON.stringify would write. */
function inf(key, v) {
  return typeof v === 'number' && !isFinite(v)
    ? (v > 0 ? 'inf' : '-inf') : v;
}

process.stdout.write(JSON.stringify({
  theta_sol: thSol, theta_surface_max: thMax,
  crossover_E: R.crossoverEnrichment(thSol, thMax),
  steady: steady, windows: windows,
  crossing: crossing, phase: phase,
  cycling: { measured: measured, synthetic_rows: synth,
             synthetic_fit: back, two_sample: two, two_sample_4: two4 }
}, inf));
