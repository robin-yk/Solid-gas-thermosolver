/* Coarse-grained KMC of reoxidation transport, after Katsoulakis &
   Vlachos, J. Chem. Phys. 119, 9412 (2003). Line-for-line port of
   solidgas/cgmc.py: same cell construction, same rates, same adaptive
   pair-relaxation stepper with the exact surface-pair sink, same seeded
   tau-leap - verified against the Python reference and the 50-digit
   matrix-exponential oracle by tests/test_statmech_port.py. */

(function (root, factory) {
  'use strict';
  if (typeof module !== 'undefined' && module.exports) {
    module.exports = factory(require('./statmech.js'));
  } else {
    root.CGMC = factory(root.StatMech);
  }
})(typeof self !== 'undefined' ? self : this, function (SM) {
  'use strict';

  var KB_EV = SM.KB_EV;

  function kfillOf(p, kin) {
    var d = SM.co2KineticsParams(p, kin);
    var kt = KB_EV * (d.T_reox_C + 273.15);
    var pf = Math.pow(d.p_CO2_atm, d.p_order_m);
    return [d.A_eff_per_s_atm.surface * pf
            * Math.exp(-d.Ea_eV.surface / kt), d];
  }

  function build(p, opts) {
    opts = opts || {};
    var cgp = {};
    var k;
    for (k in p.kinetics_cgmc) cgp[k] = p.kinetics_cgmc[k];
    if (opts.cg) for (k in opts.cg) cgp[k] = opts.cg[k];
    var vo = opts.vo_total != null ? opts.vo_total
      : p.defaults.VO_total_umol_g;
    var dUm = opts.d_um != null ? opts.d_um
      : p.defaults.particle_diameter_um;
    var eps = opts.eps != null ? opts.eps : SM.energetics(p).eps;
    var aAng = p.lattice_constants_A.a;
    var cAng = p.lattice_constants_A.c;
    var aNm = aAng / Math.sqrt(2.0) / 10.0;
    var nO = 4.0 / (aAng * aAng * cAng) * 1000.0;
    var vO = 1.0 / nO;
    var sigB = 1.0 / (cAng * aAng * Math.sqrt(2.0)) * 100.0;
    var sigL = 4.0 * sigB;
    var rNm = dUm * 1000.0 / 2.0;
    var kfPair = kfillOf(p, opts.kin);
    var kf = kfPair[0];
    var kinUsed = kfPair[1];
    var kt = KB_EV * (kinUsed.T_reox_C + 273.15);
    var beta = 1.0 / kt;

    var m = Math.max(6, Math.floor(cgp.n_cells));
    var nBulk = m - 2;
    var rCore = rNm - 2.0 * aNm;
    var w = rCore / nBulk;

    var q = [4.0 * Math.PI * rNm * rNm * sigB,
             4.0 * Math.PI * (rNm - aNm) * (rNm - aNm) * sigL];
    var centers = [rNm - 0.5 * aNm, rNm - 1.5 * aNm];
    var epsC = [eps[0], eps[1]];
    var i;
    for (i = 0; i < nBulk; i++) {
      var rOut = rCore - i * w;
      var rIn = rOut - w;
      q.push(nO * 4.0 * Math.PI / 3.0
             * (rOut * rOut * rOut - rIn * rIn * rIn));
      centers.push(0.5 * (rOut + rIn));
      epsC.push(eps[2]);
    }
    var ifaces = [rNm - aNm, rCore];
    for (i = 1; i < nBulk; i++) ifaces.push(rCore - i * w);

    var eTs = cgp.E_m_bulk_eV + eps[2];
    var nu0 = cgp.nu0_per_s;
    var lam = [];
    for (var j = 0; j < m - 1; j++) {
      var area = 4.0 * Math.PI * ifaces[j] * ifaces[j];
      var dist = centers[j] - centers[j + 1];
      lam.push(nu0 * aNm * aNm * area * Math.exp(-beta * eTs)
               / (dist * vO * q[j] * q[j + 1]));
    }

    var nOParticle = nO * 4.0 * Math.PI / 3.0 * rNm * rNm * rNm;
    var nOPerG = 2.0 / p.molar_mass_g_mol * 1e6;
    var vTotal = vo / nOPerG * nOParticle;

    var df = opts.dist_fractions;
    if (df == null) {
      var geom = SM.geometry(p, { d_um: dUm });
      df = SM.distribute(p, kinUsed.T_reox_C, vo, geom,
                         { eps: eps }).fractions;
    }
    var eta = new Array(m);
    for (i = 0; i < m; i++) eta[i] = 0.0;
    var wantS = df.surface * vTotal;
    eta[0] = Math.min(wantS, 0.999999 * q[0]);
    var wantSS = df.subsurface * vTotal + (wantS - eta[0]);
    eta[1] = Math.min(wantSS, 0.999999 * q[1]);
    var rest = vTotal - eta[0] - eta[1];
    var qBulk = 0.0;
    for (i = 2; i < m; i++) qBulk += q[i];
    for (i = 2; i < m; i++) eta[i] = rest * q[i] / qBulk;

    return { m: m, q: q, eps_c: epsC, lam: lam, beta: beta, k_fill: kf,
             eta0: eta, V0: vTotal, R_nm: rNm, a_nm: aNm,
             centers: centers, T_reox_C: kinUsed.T_reox_C,
             E_m_eV: cgp.E_m_bulk_eV, eps_ctrl: cgp.eps_ctrl,
             column_counts: cgp.stochastic_column_counts };
  }

  function expm2(a, b, c, d, dt, x0, y0) {
    var mu = 0.5 * (a + d);
    var disc = 0.25 * (a - d) * (a - d) + b * c;
    var rt = Math.sqrt(disc > 0.0 ? disc : 0.0);
    var lp = (mu + rt) * dt;
    var lm = (mu - rt) * dt;
    var ep = Math.exp(lp < 700.0 ? lp : 700.0);
    var em = lm > -700.0 ? Math.exp(lm) : 0.0;
    var halfSum = 0.5 * (ep + em);
    var shOver = rt * dt < 1e-8 ? dt * ep : (ep - em) / (2.0 * rt);
    var x = halfSum * x0 + shOver * ((a - mu) * x0 + b * y0);
    var y = halfSum * y0 + shOver * (c * x0 + (d - mu) * y0);
    return [x > 0.0 ? x : 0.0, y > 0.0 ? y : 0.0];
  }

  function sweep(sys, q, eta, dt, forward) {
    var m = sys.m;
    var lam = sys.lam;
    var epsC = sys.eps_c;
    var beta = sys.beta;
    var kf = sys.k_fill;
    var removed = 0.0;
    var neu = eta.slice();
    var js = [];
    var j;
    if (forward) {
      for (j = 0; j < m - 1; j++) js.push(j);
    } else {
      for (j = m - 2; j >= 0; j--) js.push(j);
    }
    for (var ji = 0; ji < js.length; ji++) {
      j = js[ji];
      var k = j;
      var l = j + 1;
      var a1 = Math.exp(beta * epsC[k]);
      var a2 = Math.exp(beta * epsC[l]);
      if (j === 0 && kf > 0.0) {
        var r01 = lam[0] * a1 * (q[1] - neu[1]);
        var r10 = lam[0] * a2 * (q[0] - neu[0]);
        var tot0 = neu[0] + neu[1];
        var xy = expm2(-(r01 + kf), r10, r01, -r10, dt, neu[0], neu[1]);
        var x = xy[0];
        var y = xy[1];
        if (x > q[0]) x = q[0];
        if (y > q[1]) y = q[1];
        if (x + y > tot0) {
          var s2 = tot0 / (x + y);
          x *= s2;
          y *= s2;
        }
        removed += tot0 - x - y;
        neu[0] = x;
        neu[1] = y;
        continue;
      }
      var rk = lam[j] * a1 * (q[l] - neu[l]);
      var rl = lam[j] * a2 * (q[k] - neu[k]);
      var s = rk + rl;
      if (s <= 0.0) continue;
      var tot = neu[k] + neu[l];
      var aa = a1 - a2;
      var bb = a1 * (q[l] - tot) + a2 * (q[k] + tot);
      var cc = -a2 * tot * q[k];
      var lo = tot - q[l];
      if (lo < 0.0) lo = 0.0;
      var hi = tot < q[k] ? tot : q[k];
      var xstar;
      if (Math.abs(aa) * (Math.abs(lo) + Math.abs(hi))
          < 1e-14 * Math.abs(bb)) {
        xstar = -cc / bb;
      } else {
        var disc = bb * bb - 4.0 * aa * cc;
        var sq = Math.sqrt(disc > 0.0 ? disc : 0.0);
        var qq = -0.5 * (bb + (bb >= 0.0 ? sq : -sq));
        xstar = cc / qq;
        if (!(xstar >= lo && xstar <= hi)) xstar = qq / aa;
      }
      if (xstar < lo) xstar = lo;
      if (xstar > hi) xstar = hi;
      var xx = s * dt;
      var wfrac = 1.0 - (xx < 700.0 ? Math.exp(-xx) : 0.0);
      var xk = neu[k] + (xstar - neu[k]) * wfrac;
      if (xk < lo) xk = lo;
      if (xk > hi) xk = hi;
      neu[k] = xk;
      neu[l] = tot - xk;
    }
    return [neu, removed];
  }

  function poisson(lam_, rng) {
    if (lam_ <= 0.0) return 0.0;
    if (lam_ < 30.0) {
      var limit = Math.exp(-lam_);
      var k = 0;
      var prod = rng();
      while (prod > limit) {
        k += 1;
        prod *= rng();
      }
      return k;
    }
    var u1 = rng();
    var u2 = rng();
    var z = Math.sqrt(-2.0 * Math.log(1.0 - u1))
      * Math.cos(2.0 * Math.PI * u2);
    var n = Math.round(lam_ + Math.sqrt(lam_) * z);
    return n > 0 ? n : 0.0;
  }

  function leapCap(sys, q, eta, kf) {
    var m = sys.m;
    var lam = sys.lam;
    var epsC = sys.eps_c;
    var beta = sys.beta;
    var cap = 1e30;
    for (var j = 0; j < m - 1; j++) {
      var k = j;
      var l = j + 1;
      var pairs = [
        [lam[j] * Math.exp(beta * epsC[k]) * eta[k] * (q[l] - eta[l]),
         eta[k]],
        [lam[j] * Math.exp(beta * epsC[l]) * eta[l] * (q[k] - eta[k]),
         eta[l]]];
      for (var i = 0; i < 2; i++) {
        if (pairs[i][0] > 0.0) {
          var c = 0.2 * Math.max(pairs[i][1], 30.0) / pairs[i][0];
          if (c < cap) cap = c;
        }
      }
    }
    if (kf > 0.0 && eta[0] > 0.0) {
      var c2 = 0.2 * Math.max(eta[0], 30.0) / (kf * eta[0]);
      if (c2 < cap) cap = c2;
    }
    return cap;
  }

  function leap(sys, q, eta, dt, rng) {
    var m = sys.m;
    var lam = sys.lam;
    var epsC = sys.eps_c;
    var beta = sys.beta;
    var neu = eta.slice();
    for (var j = 0; j < m - 1; j++) {
      var k = j;
      var l = j + 1;
      var fkl = lam[j] * Math.exp(beta * epsC[k]) * eta[k]
        * (q[l] - eta[l]);
      var flk = lam[j] * Math.exp(beta * epsC[l]) * eta[l]
        * (q[k] - eta[k]);
      var moves = [[fkl, k, l], [flk, l, k]];
      for (var i = 0; i < 2; i++) {
        var x = moves[i][0] * dt;
        var n = x > 1e5 ? x : poisson(x, rng);
        if (n > neu[moves[i][1]]) n = neu[moves[i][1]];
        var room = q[moves[i][2]] - neu[moves[i][2]];
        if (n > room) n = room > 0.0 ? room : 0.0;
        neu[moves[i][1]] -= n;
        neu[moves[i][2]] += n;
      }
    }
    return neu;
  }

  function run(sys, tEnd, opts) {
    opts = opts || {};
    var mode = opts.mode || 'det';
    var m = sys.m;
    var eta = sys.eta0.slice();
    var scale = 1.0;
    var rng = null;
    if (mode === 'stoch') {
      scale = Math.max(1.0, sys.V0 / sys.column_counts);
      for (var i0 = 0; i0 < m; i0++) eta[i0] /= scale;
      rng = SM.makeRng(opts.seed != null ? opts.seed : 902);
    }
    var q = sys.q.map(function (x) { return x / scale; });
    var v0 = sys.V0 / scale;
    var kf = sys.k_fill;
    var epsCtrl = opts.eps_ctrl != null ? opts.eps_ctrl : sys.eps_ctrl;
    var nOut = opts.n_out != null ? opts.n_out : 90;

    var tFirst = 1e-11;
    var outs = [0.0];
    var i;
    for (i = 0; i < nOut; i++) {
      outs.push(tFirst * Math.pow(tEnd / tFirst, i / (nOut - 1.0)));
    }
    var profTimes = opts.profile_times
      || [tEnd * 1e-6, tEnd * 1e-3, tEnd * 3e-2, tEnd];
    var profLeft = profTimes.slice().sort(function (a, b) { return a - b; });
    var profiles = [];

    var times = [0.0];
    var rec = [0.0];
    var t = 0.0;
    var dt = 1e-12;
    var oi = 1;
    var removed = 0.0;
    var guard = 0;
    var stepI = 0;
    while (t < tEnd && guard < 400000) {
      guard += 1;
      if (mode === 'stoch') {
        var capd = leapCap(sys, q, eta, kf);
        if (dt > capd) dt = capd;
      }
      var neu, fill, worst;
      for (;;) {
        if (mode === 'stoch') {
          fill = poisson(kf * eta[0] * dt, rng);
          if (fill > eta[0]) fill = eta[0];
          var work = eta.slice();
          work[0] -= fill;
          neu = leap(sys, q, work, dt, rng);
        } else {
          var sr = sweep(sys, q, eta, dt, stepI % 2 === 0);
          neu = sr[0];
          fill = sr[1];
        }
        worst = 0.0;
        var floor = 1e-6 * v0;
        for (var k2 = 0; k2 < m; k2++) {
          var ref = eta[k2] < q[k2] - eta[k2] ? eta[k2] : q[k2] - eta[k2];
          var d = Math.abs(neu[k2] - eta[k2]) / (ref + floor);
          if (d > worst) worst = d;
        }
        if (worst <= epsCtrl || dt <= 1e-13) break;
        dt *= 0.5;
      }
      eta = neu;
      removed += fill;
      t += dt;
      stepI += 1;
      if (worst < 0.5 * epsCtrl) dt *= 1.3;
      if (t + dt > tEnd) dt = Math.max(tEnd - t, 1e-13);
      while (oi < outs.length && t >= outs[oi]) {
        times.push(t);
        rec.push(removed / v0);
        oi += 1;
      }
      while (profLeft.length && t >= profLeft[0]) {
        profiles.push({ t_s: t, theta: eta.map(function (x, kk) {
          return x / q[kk]; }) });
        profLeft.shift();
      }
    }
    while (profLeft.length) {
      profiles.push({ t_s: t, theta: eta.map(function (x, kk) {
        return x / q[kk]; }) });
      profLeft.shift();
    }
    times.push(t);
    rec.push(removed / v0);
    var sumEta = 0.0;
    for (i = 0; i < m; i++) sumEta += eta[i];
    return { t_s: times, recovery: rec, final_eta: eta,
             final_theta: eta.map(function (x, kk) { return x / q[kk]; }),
             profiles: profiles, steps: guard, mode: mode,
             recovered_fraction: rec[rec.length - 1], V0_column: v0,
             mass_check: (sumEta + removed) / v0 - 1.0,
             completed: t >= 0.999 * tEnd,
             centers_nm: sys.centers,
             depth_nm: sys.centers.map(function (c) {
               return sys.R_nm - c; }) };
  }

  function effectiveBarrier(p, targetRecovery, tExp, opts) {
    opts = opts || {};
    var lo = opts.lo != null ? opts.lo : 0.3;
    var hi = opts.hi != null ? opts.hi : 3.0;
    for (var i = 0; i < 36; i++) {
      var mid = 0.5 * (lo + hi);
      var cg = {};
      if (opts.cg) for (var k in opts.cg) cg[k] = opts.cg[k];
      cg.E_m_bulk_eV = mid;
      var sys = build(p, { vo_total: opts.vo_total, d_um: opts.d_um,
                           kin: opts.kin, cg: cg,
                           dist_fractions: opts.dist_fractions,
                           eps: opts.eps });
      var r = run(sys, tExp, { n_out: 12 });
      if (r.recovered_fraction > targetRecovery) lo = mid;
      else hi = mid;
    }
    return 0.5 * (lo + hi);
  }

  return { build: build, run: run, effectiveBarrier: effectiveBarrier,
           expm2: expm2 };
});
