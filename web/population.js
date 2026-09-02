/* The isolated-surface-vacancy population model, in the browser.

   A mirror of solidgas/vacancy_population.py, reduction included. Writing
   q = n_iso + n_assoc, the loss term cancels out of the sum, so

       q(t)       = n_s [1 - exp(-g t / n_s)]
       n_below(t) = g t - q(t)

   are closed form and only one scalar equation is left,

       dn_iso/dt = g exp(-g t / n_s) - (2 k_loss / n_s) n_iso^2

   with n_iso(0) = 0 and n_assoc = q - n_iso. The closure is then exact
   rather than something to check.

   Python hands that equation to Radau. There is no implicit solver here,
   so this uses an exponentially-fitted step: over a step the equation is
   treated as constant generation against a linearised quadratic loss, both
   integrated in closed form, which is stable for any k_loss the page can
   be asked for. A parity gate holds the two to a relative tolerance rather
   than to the ulp, because they are different integrators of one equation
   and pretending otherwise would be a fiction. */

(function (root, factory) {
  'use strict';
  var api = factory();
  if (typeof module !== 'undefined' && module.exports) module.exports = api;
  root.Population = api;
})(typeof self !== 'undefined' ? self : this, function () {
  'use strict';

  var KB_EV = 8.617333262e-5;
  var STEPS = 4000;
  var NORMALISE_TO = 'R600';
  var REPORTED = { ns: 15.33381126, E_loss: 2.18450850, log10_nu: 6.74333249 };

  function rateConstant(T_K, E_loss, nu_eff) {
    var x = -E_loss / (KB_EV * T_K);
    return x > -700 ? nu_eff * Math.exp(x) : 0;
  }

  function integrate(nV, t_red, T_K, ns, E_loss, nu_eff, steps) {
    steps = steps || STEPS;
    if (!(nV > 0 && t_red > 0 && ns > 0)) {
      throw new Error('inventory, reduction time and capacity must be positive');
    }
    var g = nV / t_red;
    var k = rateConstant(T_K, E_loss, nu_eff);
    var lam = g / ns;
    var a = 2 * k / ns;
    var dt = t_red / steps;
    var iso = 0;

    /* Strang: half a step of loss, a whole step of generation, half a step
       of loss. Both sub-flows are exact, so the scheme is second order and
       stable for any k_loss. */
    function decay(h) {
      if (a > 0 && iso > 0) iso = iso / (1 + a * iso * h);
    }
    for (var i = 0; i < steps; i++) {
      var t = i * dt;
      decay(0.5 * dt);
      iso += g * (Math.exp(-lam * t) - Math.exp(-lam * (t + dt))) / lam;
      decay(0.5 * dt);
    }
    var q = ns * (1 - Math.exp(-lam * t_red));
    return { n_iso: iso, n_assoc: q - iso, n_below: g * t_red - q };
  }

  function partition(series, ns, E_loss, nu_eff, steps, normaliseTo) {
    normaliseTo = normaliseTo || NORMALISE_TO;
    var out = series.map(function (r) {
      var T_K = r.T_red_K != null ? r.T_red_K : r.T_red_C + 273.15;
      var y = integrate(r.nV_total_umol_g, r.t_red_s, T_K,
                        ns, E_loss, nu_eff, steps);
      return {
        sample: r.sample, treatment: r.treatment,
        T_red_C: r.T_red_C, T_red_K: T_K, t_red_s: r.t_red_s,
        nV_total_umol_g: r.nV_total_umol_g,
        r_CO_measured: r.r_CO_umol_g_s,
        n_iso_umol_g: y.n_iso, n_assoc_umol_g: y.n_assoc,
        n_below_umol_g: y.n_below,
        closure_umol_g: y.n_iso + y.n_assoc + y.n_below - r.nV_total_umol_g,
        below_fraction: y.n_below / r.nV_total_umol_g
      };
    });
    var ref = out.filter(function (q) { return q.sample === normaliseTo; })[0];
    if (!ref) throw new Error('no sample named ' + normaliseTo);
    var scale = ref.n_iso_umol_g > 0 ? ref.r_CO_measured / ref.n_iso_umol_g : 0;
    out.forEach(function (q) {
      q.r_CO_fitted = scale * q.n_iso_umol_g;
      q.residual = q.r_CO_fitted - q.r_CO_measured;
      q.residual_relative = q.residual / q.r_CO_measured;
      q.residual_log = q.r_CO_fitted > 0
        ? Math.log(q.r_CO_fitted / q.r_CO_measured) : Infinity;
    });
    return out;
  }

  function residuals(rows, kind) {
    var key = kind === 'raw' ? 'residual'
      : kind === 'relative' ? 'residual_relative' : 'residual_log';
    return rows.map(function (q) { return q[key]; });
  }

  function objective(series, ns, E_loss, nu_eff, kind, steps) {
    var rows;
    try { rows = partition(series, ns, E_loss, nu_eff, steps); }
    catch (e) { return 1e30; }
    var r = residuals(rows, kind || 'log');
    for (var i = 0; i < r.length; i++) if (!isFinite(r[i])) return 1e30;
    return r.reduce(function (a, v) { return a + v * v; }, 0);
  }

  function rateRatio(series, ns, E_loss, nu_eff, hi, lo, steps) {
    var by = {};
    partition(series, ns, E_loss, nu_eff, steps).forEach(function (q) {
      by[q.sample] = q;
    });
    return by[hi || 'R800'].r_CO_fitted / by[lo || 'R600'].r_CO_fitted;
  }

  function measuredRatio(series, hi, lo) {
    var by = {};
    series.forEach(function (r) { by[r.sample] = r.r_CO_umol_g_s; });
    return by[hi || 'R800'] / by[lo || 'R600'];
  }

  /* Nelder-Mead on (log n_s, E_loss, log10 nu_eff). */
  function fit(series, kind, nsBounds, start, steps, iters) {
    var b = nsBounds || [13.6, 67.8];
    kind = kind || 'log';
    iters = iters || 400;
    var searchSteps = steps || 1500;
    function unpack(p) {
      return [Math.min(Math.max(Math.exp(p[0]), b[0]), b[1]),
              Math.max(p[1], 0), Math.pow(10, p[2])];
    }
    function cost(p) {
      var u = unpack(p);
      return objective(series, u[0], u[1], u[2], kind, searchSteps);
    }
    var s0 = start || [Math.log(20), 2.2, 6.7];
    var pts = [s0.slice(), [s0[0] + 0.35, s0[1], s0[2]],
               [s0[0], s0[1] + 0.25, s0[2]], [s0[0], s0[1], s0[2] + 0.6]];
    var vals = pts.map(cost);
    for (var it = 0; it < iters; it++) {
      var ord = [0, 1, 2, 3].sort(function (x, y) { return vals[x] - vals[y]; });
      pts = ord.map(function (i) { return pts[i]; });
      vals = ord.map(function (i) { return vals[i]; });
      var cen = [0, 1, 2].map(function (j) {
        return (pts[0][j] + pts[1][j] + pts[2][j]) / 3;
      });
      var refl = [0, 1, 2].map(function (j) { return cen[j] + (cen[j] - pts[3][j]); });
      var fr = cost(refl);
      if (fr < vals[0]) {
        var ex = [0, 1, 2].map(function (j) { return cen[j] + 2 * (cen[j] - pts[3][j]); });
        var fe = cost(ex);
        if (fe < fr) { pts[3] = ex; vals[3] = fe; } else { pts[3] = refl; vals[3] = fr; }
      } else if (fr < vals[2]) { pts[3] = refl; vals[3] = fr; }
      else {
        var con = [0, 1, 2].map(function (j) {
          return cen[j] + 0.5 * (pts[3][j] - cen[j]);
        });
        var fc = cost(con);
        if (fc < vals[3]) { pts[3] = con; vals[3] = fc; }
        else {
          for (var i = 1; i < 4; i++) {
            pts[i] = [0, 1, 2].map(function (j) {
              return pts[0][j] + 0.5 * (pts[i][j] - pts[0][j]);
            });
            vals[i] = cost(pts[i]);
          }
        }
      }
    }
    var u = unpack(pts[0]);
    return { residual_kind: kind, ns: u[0], E_loss: u[1], nu_eff: u[2],
             log10_nu: Math.log10(u[2]),
             objective: objective(series, u[0], u[1], u[2], kind) };
  }

  return { KB_EV: KB_EV, STEPS: STEPS, REPORTED: REPORTED,
           rateConstant: rateConstant, integrate: integrate,
           partition: partition, residuals: residuals, objective: objective,
           rateRatio: rateRatio, measuredRatio: measuredRatio, fit: fit };
});
