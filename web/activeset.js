/* Active-set gas/solid Gibbs equilibrium for the browser.
   Line-for-line port of solidgas/activeset.py: gas species in logarithmic
   element-potential coordinates, condensed phases enumerated as active sets,
   excluded phases exactly zero, KKT reduced costs in two bases, one Gibbs
   reference per charge. tests/test_activeset_port.py pins this against the
   Python production solver, which is itself pinned against the 80-digit
   oracle. */

(function (root) {
'use strict';

var TOL_NEWTON = 1e-13;
var TOL_KKT_KJ = 1e-9;
var TOL_DEG_KJ = 1e-6;
var TRACE_MOL = 1e-12;
var EXP_CAP = 700.0;

function Solver(D) {
  this.D = D;
  this.R = D.R_KJ;
  this.RATM = D.R_ATM;
}

/* ------------------------------------------------------------------- mu0 */

Solver.prototype.mu0Gas = function (name, T) {
  var e = this.D.shomate[name];
  var c = T <= e.brk ? e.lo : e.hi;
  var t = T / 1000.0;
  var Ht = e.dHf + (c[0] * t + c[1] * t * t / 2 + c[2] * t * t * t / 3
                    + c[3] * t * t * t * t / 4 - c[4] / t + c[5] - c[7]);
  var St = c[0] * Math.log(t) + c[1] * t + c[2] * t * t / 2
           + c[3] * t * t * t / 3 - c[4] / (2 * t * t) + c[6];
  return Ht - T * St / 1000.0;
};

function cpInt(form, c, T0, T1) {
  var F = form === 'poly'
    ? function (T) { return c[0] * T + c[1] * T * T / 2 + c[2] * T * T * T / 3 - c[3] / T; }
    : function (T) { return c[0] * T + 2 * c[1] * Math.sqrt(T) - c[2] / T - c[3] / (2 * T * T); };
  return F(T1) - F(T0);
}
function cpOverTInt(form, c, T0, T1) {
  var F = form === 'poly'
    ? function (T) { return c[0] * Math.log(T) + c[1] * T + c[2] * T * T / 2 - c[3] / (2 * T * T); }
    : function (T) { return c[0] * Math.log(T) - 2 * c[1] / Math.sqrt(T) - c[2] / (2 * T * T) - c[3] / (3 * T * T * T); };
  return F(T1) - F(T0);
}

Solver.prototype.mu0Solid = function (name, T) {
  var d = this.D.waldner[name];
  var H = d.dHf298, S = d.S298, lo = 298.15;
  var segs = [], i, r;
  for (i = 0; i < d.ranges.length; i++) {
    r = d.ranges[i];
    if (T <= r.T0) break;
    var hi = T < r.T1 ? T : r.T1;
    var from = lo > r.T0 ? lo : r.T0;
    if (hi > from) { segs.push({ r: r, t0: from, t1: hi, Ttr: r.T0 }); lo = hi; }
    if (T <= r.T1) break;
  }
  if (segs.length && segs[segs.length - 1].t1 < T) {
    r = d.ranges[d.ranges.length - 1];
    segs.push({ r: { form: r.form, c: r.c, dHt: 0 },
                t0: segs[segs.length - 1].t1, t1: T, Ttr: r.T0 });
  }
  for (i = 0; i < segs.length; i++) {
    var s = segs[i];
    H += s.r.dHt + cpInt(s.r.form, s.r.c, s.t0, s.t1);
    S += (s.r.dHt ? s.r.dHt / s.Ttr : 0) + cpOverTInt(s.r.form, s.r.c, s.t0, s.t1);
  }
  return (H - T * S) / 1000.0;
};

Solver.prototype.moleWeight = function (solid) {
  var st = this.D.stoich[solid];
  return st[0] * this.D.atomic_weight.Ti + st[1] * this.D.atomic_weight.O;
};

/* -------------------------------------------------------------- the charge */

Solver.prototype.buildCharge = function (feedRel, T_K, P, V, massG, initialSolid, TchargeK) {
  var D = this.D, g, tot = 0;
  for (g in feedRel) {
    if (!D.gas_formula[g]) throw new Error('unknown gas species ' + g);
    if (feedRel[g] < 0) throw new Error('negative feed amount for ' + g);
    tot += Number(feedRel[g]) || 0;
  }
  if (!(tot > 0)) throw new Error('gas feed sums to zero');
  if (!(massG > 0)) throw new Error('solid mass must be positive');
  if (!D.stoich[initialSolid]) throw new Error('unknown initial solid ' + initialSolid);
  var Tc = TchargeK || T_K;
  var ng0 = P * V / (this.RATM * Tc);
  var y0 = {}, n0 = {};
  D.gases.forEach(function (s) {
    y0[s] = (Number(feedRel[s]) || 0) / tot;
    n0[s] = y0[s] * ng0;
  });
  var st = D.stoich[initialSolid];
  var nSolid0 = massG / this.moleWeight(initialSolid);
  var nTi = nSolid0 * st[0];
  var rho = st[1] / st[0];
  var b = { C: 0, H: 0, N: 0, O: 0, Ti: nTi }, bOgas = 0;
  D.gases.forEach(function (s) {
    var f = D.gas_formula[s], e;
    for (e in f) b[e] += f[e] * n0[s];
    bOgas += (f.O || 0) * n0[s];
  });
  b.O += st[1] * nSolid0;
  return { T: T_K, Tcharge: Tc, P: P, V: V, mass: massG,
           initialSolid: initialSolid, nSolid0: nSolid0, ng0: ng0,
           y0: y0, n0: n0, nTi: nTi, rho: rho, b: b, bOgasfeed: bOgas };
};

/* ------------------------------------------------------- linear elimination */

function linSolve(A, rhs) {
  var n = rhs.length, M = [], i, r, c, k;
  for (i = 0; i < n; i++) M.push(A[i].concat([rhs[i]]));
  for (c = 0; c < n; c++) {
    var p = c;
    for (r = c + 1; r < n; r++) if (Math.abs(M[r][c]) > Math.abs(M[p][c])) p = r;
    if (Math.abs(M[p][c]) < 1e-300) return null;
    var tmp = M[c]; M[c] = M[p]; M[p] = tmp;
    var piv = M[c][c];
    for (r = 0; r < n; r++) {
      if (r === c) continue;
      var f = M[r][c] / piv;
      if (f) for (k = c; k <= n; k++) M[r][k] -= f * M[c][k];
    }
  }
  var x = [];
  for (r = 0; r < n; r++) x.push(M[r][n] / M[r][r]);
  return x;
}

/* ---------------------------------------------------------- one active set */

Solver.prototype.solveCandidate = function (active, ch, mu, registry) {
  var D = this.D, self = this;
  var T = ch.T, P = ch.P, b = ch.b, RT = this.R * T, rho = ch.rho;
  var tj = {}, oj = {}, dj = {};
  active.forEach(function (s) {
    tj[s] = D.stoich[s][0]; oj[s] = D.stoich[s][1];
    dj[s] = rho * tj[s] - oj[s];
  });

  var lam = {};
  if (active.length >= 2) {
    var a1 = active[0], a2 = active[1];
    var det0 = tj[a1] * oj[a2] - tj[a2] * oj[a1];
    lam.Ti = (mu[a1] * oj[a2] - mu[a2] * oj[a1]) / det0;
    lam.O = (tj[a1] * mu[a2] - tj[a2] * mu[a1]) / det0;
  }

  var have = { C: b.C > 0, H: b.H > 0, N: b.N > 0 };
  var mFixed = null, availO = null, oFree = false;
  if (active.length === 1) {
    var s0 = active[0];
    mFixed = b.Ti / tj[s0];
    availO = ch.bOgasfeed + dj[s0] * mFixed;
    oFree = availO > 0;
  }

  function elemOk(e) {
    if (e === 'C' || e === 'H' || e === 'N') return have[e];
    return oFree || active.length >= 2;
  }
  var species = D.gases.filter(function (g) {
    var f = D.gas_formula[g], e;
    for (e in f) if (!elemOk(e)) return false;
    return true;
  });
  if (!species.length) {
    return { active: active.slice(), feasible: false, kkt_ok: false,
             reason: 'no gas species with inventory' };
  }

  var elems = ['C', 'H', 'N'].filter(function (e) { return have[e]; });
  if (oFree) elems.push('O');
  var nE = elems.length;

  var gi = {};
  species.forEach(function (s) {
    var v = mu[s] / RT + Math.log(P);
    if (active.length >= 2 && D.gas_formula[s].O) {
      v -= D.gas_formula[s].O * lam.O / RT;
    }
    gi[s] = v;
  });

  var rhs = { C: b.C, H: b.H, N: b.N };
  if (oFree) rhs.O = availO;
  var scale = {};
  elems.forEach(function (e) { scale[e] = Math.max(rhs[e], 1e-300); });

  /* seed: feed-shaped gas, floored; theta by normal equations */
  var yg = {}, tot = 0;
  species.forEach(function (s) { yg[s] = Math.max(ch.y0[s] || 0, 0); tot += yg[s]; });
  species.forEach(function (s) {
    yg[s] = tot > 0 ? yg[s] / tot : 1.0 / species.length;
    yg[s] = Math.max(yg[s], 1e-8);
  });
  tot = 0; species.forEach(function (s) { tot += yg[s]; });
  species.forEach(function (s) { yg[s] /= tot; });
  var AtA = [], Atr = [], i, k;
  for (i = 0; i < nE; i++) { AtA.push(new Array(nE).fill(0)); Atr.push(0); }
  species.forEach(function (s) {
    var r = Math.log(yg[s]) + gi[s], f = D.gas_formula[s];
    for (i = 0; i < nE; i++) {
      var ae = f[elems[i]] || 0;
      if (!ae) continue;
      Atr[i] += ae * r;
      for (k = 0; k < nE; k++) AtA[i][k] += ae * (f[elems[k]] || 0);
    }
  });
  var th = nE ? (linSolve(AtA, Atr) || new Array(nE).fill(0)) : [];
  var x = th.concat([Math.log(ch.ng0)]);

  function gasY(xv) {
    var y = {};
    species.forEach(function (s) {
      var v = -gi[s], f = D.gas_formula[s];
      for (var q = 0; q < nE; q++) v += (f[elems[q]] || 0) * xv[q];
      y[s] = Math.exp(Math.min(v, EXP_CAP));
    });
    return y;
  }
  function residual(xv) {
    var y = gasY(xv);
    var Ng = Math.exp(Math.min(xv[nE], EXP_CAP));
    var res = [];
    for (var q = 0; q < nE; q++) {
      var e = elems[q], t2 = 0;
      species.forEach(function (s) { t2 += (D.gas_formula[s][e] || 0) * y[s]; });
      res.push((t2 * Ng - rhs[e]) / scale[e]);
    }
    var sy = 0;
    species.forEach(function (s) { sy += y[s]; });
    res.push(sy - 1.0);
    return { res: res, y: y, Ng: Ng };
  }
  function maxAbs(v) {
    var m = 0;
    for (var q = 0; q < v.length; q++) if (Math.abs(v[q]) > m) m = Math.abs(v[q]);
    return m;
  }

  var st = residual(x), rn = maxAbs(st.res), it = 0;
  while (rn > TOL_NEWTON && it < 200) {
    var J = [];
    for (i = 0; i <= nE; i++) J.push(new Array(nE + 1).fill(0));
    for (i = 0; i < nE; i++) {
      var e = elems[i];
      for (k = 0; k < nE; k++) {
        var f2 = elems[k], acc = 0;
        species.forEach(function (s) {
          acc += (D.gas_formula[s][e] || 0) * (D.gas_formula[s][f2] || 0) * st.y[s];
        });
        J[i][k] = acc * st.Ng / scale[e];
      }
      var acc2 = 0;
      species.forEach(function (s) { acc2 += (D.gas_formula[s][e] || 0) * st.y[s]; });
      J[i][nE] = acc2 * st.Ng / scale[e];
    }
    for (k = 0; k < nE; k++) {
      var f3 = elems[k], acc3 = 0;
      species.forEach(function (s) { acc3 += (D.gas_formula[s][f3] || 0) * st.y[s]; });
      J[nE][k] = acc3;
    }
    J[nE][nE] = 0;
    var negres = st.res.map(function (v) { return -v; });
    var step = linSolve(J, negres);
    if (!step) {
      return { active: active.slice(), feasible: false, kkt_ok: false,
               reason: 'singular Jacobian' };
    }
    var damp = 1.0, ok = false, xn, stn, rnn;
    for (var h = 0; h < 60; h++) {
      xn = x.map(function (v, q) { return v + damp * step[q]; });
      stn = residual(xn);
      rnn = maxAbs(stn.res);
      if (rnn < rn) { ok = true; break; }
      damp /= 2;
    }
    if (!ok) {
      return { active: active.slice(), feasible: false, kkt_ok: false,
               reason: 'stalled at residual ' + rn.toExponential(2) };
    }
    x = xn; st = stn; rn = rnn; it += 1;
  }
  if (rn > TOL_NEWTON) {
    return { active: active.slice(), feasible: false, kkt_ok: false,
             reason: 'no convergence, residual ' + rn.toExponential(2) };
  }

  var nGas = {};
  D.gases.forEach(function (g) { nGas[g] = 0; });
  species.forEach(function (s) { nGas[s] = st.y[s] * st.Ng; });

  var gasO = 0;
  species.forEach(function (s) { gasO += (D.gas_formula[s].O || 0) * nGas[s]; });
  var m = {};
  if (active.length === 1) {
    m[active[0]] = mFixed;
  } else {
    var b1 = active[0], b2 = active[1];
    var det = tj[b1] * dj[b2] - tj[b2] * dj[b1];
    if (det === 0) {
      return { active: active.slice(), feasible: false, kkt_ok: false,
               reason: 'degenerate stoichiometry pair' };
    }
    var rO = gasO - ch.bOgasfeed;
    m[b1] = (b.Ti * dj[b2] - tj[b2] * rO) / det;
    m[b2] = (tj[b1] * rO - dj[b1] * b.Ti) / det;
  }

  if (oFree) lam.O = RT * x[elems.indexOf('O')];
  if (active.length === 1 && ('O' in lam)) {
    lam.Ti = (mu[active[0]] - oj[active[0]] * lam.O) / tj[active[0]];
  }
  for (i = 0; i < nE; i++) lam[elems[i]] = RT * x[i];

  var rCosts = {}, lamODef = ('O' in lam);
  (registry || D.solids).forEach(function (s) {
    if (active.indexOf(s) >= 0) return;
    var ts = D.stoich[s][0], os = D.stoich[s][1];
    var ds = rho * ts - os;
    if (!lamODef) {
      rCosts[s] = { per_formula: ds > 0 ? -Infinity : null,
                    per_mol_O: ds > 0 ? -Infinity : null, dO: ds };
      return;
    }
    var rf = mu[s] - ts * lam.Ti - os * lam.O;
    rCosts[s] = { per_formula: rf, per_mol_O: ds !== 0 ? rf / ds : null, dO: ds };
  });

  var feasible = active.every(function (s) { return m[s] !== null && m[s] > 0; });
  var reason = feasible ? 'ok'
    : 'negative condensed amount: ' + active.map(function (s) {
        return s + '=' + (m[s] === null ? 'null' : m[s].toExponential(3));
      }).join(', ');

  var kktOk = feasible, worst = null, s2;
  for (s2 in rCosts) {
    var v2 = rCosts[s2].per_formula;
    if (v2 === null) continue;
    if (worst === null || v2 < worst) worst = v2;
    if (v2 < -TOL_KKT_KJ) kktOk = false;
  }

  return { active: active.slice(), feasible: feasible, kkt_ok: kktOk,
           reason: reason, y: st.y, n_gas: nGas, m: m, lam: lam,
           r_costs: rCosts, worst_inactive_r: worst,
           iterations: it, newton_residual: rn };
};

/* --------------------------------------------------------------- full solve */

Solver.prototype.enumerateCandidates = function (solids) {
  var D = this.D;
  var ss = solids.slice().sort(function (a, b) {
    var ra = D.stoich[a][1] / D.stoich[a][0], rb = D.stoich[b][1] / D.stoich[b][0];
    if (ra !== rb) return rb - ra;
    return a < b ? -1 : 1;
  });
  var out = ss.map(function (s) { return [s]; });
  for (var i = 0; i < ss.length; i++) {
    for (var j = i + 1; j < ss.length; j++) out.push([ss[i], ss[j]]);
  }
  return out;
};

Solver.prototype.tripleProbes = function (solids, mu) {
  var D = this.D, hits = [];
  var ss = solids.slice().sort(function (a, b) {
    var ra = D.stoich[a][1] / D.stoich[a][0], rb = D.stoich[b][1] / D.stoich[b][0];
    if (ra !== rb) return rb - ra;
    return a < b ? -1 : 1;
  });
  for (var i = 0; i < ss.length; i++) {
    for (var j = i + 1; j < ss.length; j++) {
      for (var k = j + 1; k < ss.length; k++) {
        var a1 = ss[i], a2 = ss[j], a3 = ss[k];
        var t1 = D.stoich[a1][0], o1 = D.stoich[a1][1];
        var t2 = D.stoich[a2][0], o2 = D.stoich[a2][1];
        var t3 = D.stoich[a3][0], o3 = D.stoich[a3][1];
        var det = t1 * o2 - t2 * o1;
        var lamTi = (mu[a1] * o2 - mu[a2] * o1) / det;
        var lamO = (t1 * mu[a2] - t2 * mu[a1]) / det;
        var res = mu[a3] - t3 * lamTi - o3 * lamO;
        if (Math.abs(res) < TOL_DEG_KJ) {
          hits.push({ phases: [a1, a2, a3], residual_kJ: res });
        }
      }
    }
  }
  return hits;
};

Solver.prototype.gAbsolute = function (nGas, m, mu, T, P) {
  var RT = this.R * T, Ng = 0, G = 0, s;
  for (s in nGas) if (nGas[s] > 0) Ng += nGas[s];
  for (s in nGas) {
    if (nGas[s] > 0) G += nGas[s] * (mu[s] + RT * Math.log(nGas[s] / Ng * P));
  }
  for (s in m) G += m[s] * mu[s];
  return G;
};

Solver.prototype.solve = function (opts) {
  var D = this.D, self = this;
  var T = opts.T_K, P = opts.P_atm || 1.0, V = opts.V_cm3 || 25.0;
  var mass = opts.mass_g !== undefined ? opts.mass_g : 0.100;
  var initialSolid = opts.initial_solid || 'TiO2';
  /* Canonical phase order at the door: every output is then independent of
     how the caller ordered the input list. */
  var solids = (opts.solids || D.solids).slice().sort(function (a, b) {
    var ra = D.stoich[a][1] / D.stoich[a][0], rb = D.stoich[b][1] / D.stoich[b][0];
    if (ra !== rb) return rb - ra;
    return a < b ? -1 : 1;
  });
  var ch = this.buildCharge(opts.feed, T, P, V, mass, initialSolid,
                            opts.T_charge_K);
  var mu = {};
  D.gases.forEach(function (g) { mu[g] = self.mu0Gas(g, T); });
  solids.forEach(function (s) { mu[s] = self.mu0Solid(s, T); });
  if (!(initialSolid in mu)) mu[initialSolid] = this.mu0Solid(initialSolid, T);

  var results = [], passing = [];
  this.enumerateCandidates(solids).forEach(function (cand) {
    var r = self.solveCandidate(cand, ch, mu, solids);
    results.push(r);
    if (r.feasible && r.kkt_ok) passing.push(r);
  });
  if (!passing.length) {
    throw new Error('no candidate active set satisfies KKT');
  }
  passing.sort(function (a, b) { return a.active.length - b.active.length; });
  var win = passing[0];

  var degenerate = [];
  var s;
  for (s in win.r_costs) {
    var v = win.r_costs[s].per_mol_O;
    if (v !== null && isFinite(v) && Math.abs(v) < TOL_DEG_KJ) degenerate.push(s);
  }
  var probes = this.tripleProbes(solids, mu);

  var m = {};
  solids.forEach(function (q) { m[q] = 0.0; });        /* excluded: exactly 0 */
  for (s in win.m) m[s] = win.m[s];
  var nGas = win.n_gas;

  var Gtot = this.gAbsolute(nGas, m, mu, T, P);
  var RT = this.R * T;
  var Gref = ch.nSolid0 * mu[initialSolid];
  D.gases.forEach(function (g) {
    if (ch.n0[g] > 0) Gref += ch.n0[g] * (mu[g] + RT * Math.log(ch.n0[g] / ch.ng0 * P));
  });
  var Grel = Gtot - Gref;
  var gFloor = 1e-12 * Math.max(Math.abs(Gtot), Math.abs(Gref));

  var balance = {};
  ['C', 'H', 'N', 'O', 'Ti'].forEach(function (e) {
    var t2 = 0;
    D.gases.forEach(function (g) { t2 += (D.gas_formula[g][e] || 0) * nGas[g]; });
    for (var q in m) {
      var st2 = D.stoich[q];
      t2 += (e === 'Ti' ? st2[0] : e === 'O' ? st2[1] : 0) * m[q];
    }
    balance[e] = t2 - ch.b[e];
  });

  var solidO = 0, solidTi = 0, ti3 = 0, gasO = 0, Ngtot = 0;
  for (s in m) {
    solidO += D.stoich[s][1] * m[s];
    solidTi += D.stoich[s][0] * m[s];
    ti3 += D.ti3[s] * m[s];
  }
  D.gases.forEach(function (g) {
    gasO += (D.gas_formula[g].O || 0) * nGas[g];
    Ngtot += nGas[g];
  });

  var traces = {};
  D.gases.forEach(function (g) {
    if (nGas[g] > 0 && nGas[g] < TRACE_MOL) traces[g] = Math.log(nGas[g]) / Math.LN10;
  });
  for (s in m) {
    if (m[s] > 0 && m[s] < TRACE_MOL) traces[s] = Math.log(m[s]) / Math.LN10;
  }

  var gasFractions = {}, deltaN = {}, phaseSplit = {};
  D.gases.forEach(function (g) {
    gasFractions[g] = Ngtot > 0 ? nGas[g] / Ngtot : 0;
    deltaN[g] = nGas[g] - ch.n0[g];
  });
  for (s in m) phaseSplit[s] = D.stoich[s][0] * m[s] / ch.nTi * 100.0;

  var rc = {};
  for (s in win.r_costs) {
    rc[s] = { per_formula_kJ: win.r_costs[s].per_formula,
              per_mol_O_kJ: win.r_costs[s].per_mol_O,
              oxygen_removed_per_extent: win.r_costs[s].dO };
  }

  var feedY = {};
  D.gases.forEach(function (g) { if (ch.y0[g] > 0) feedY[g] = ch.y0[g]; });

  return {
    system: 'Ti-O',
    mode: D.mode_note,
    method: 'gibbs_min',
    solver: 'active-set element-potential Newton, double precision',
    T_C: T - 273.15, T_K: T, P_atm: P, V_cm3: V,
    T_charge_K: ch.Tcharge, mass_g: mass,
    initial_solid: initialSolid,
    feed_mole_fractions: feedY,
    n_gas_charge_mol: ch.ng0, n_Ti_mol: ch.nTi,
    active_condensed_phases: win.active,
    species_moles: nGas,
    solid_moles: m,
    gas_fractions: gasFractions,
    phase_split_pct: phaseSplit,
    lambda_kJ_per_mol: win.lam,
    inactive_phase_reduced_costs: rc,
    reduced_cost_basis: 'per_formula: mu0_j - t_j*lambda_Ti - o_j*lambda_O; '
      + 'per_mol_O: divided by def_j = rho_host*t_j - o_j',
    G_total_kJ: Gtot,
    common_reference_G_kJ: Gref,
    G_rel_kJ: Grel,
    G_rel_below_double_resolution: Math.abs(Grel) < gFloor,
    reduced_pct: 100.0 * ti3 / ch.nTi,
    x_in_TiO2x: solidTi > 0 ? 2.0 - solidO / solidTi : null,
    oxygen_to_gas_mol_O: gasO - ch.bOgasfeed,
    conversion_CO2_pct: ch.n0.CO2 > 0
      ? (ch.n0.CO2 - nGas.CO2) / ch.n0.CO2 * 100.0 : null,
    delta_n_gas: deltaN,
    balance_by_element: balance,
    log10_trace_amounts: traces,
    boundary_or_degenerate: degenerate.length > 0 || probes.length > 0,
    degenerate_phases: degenerate,
    triple_probe_hits: probes,
    solver_tolerance: { newton_scaled_residual: TOL_NEWTON,
                        kkt_kJ: TOL_KKT_KJ,
                        degeneracy_kJ_per_mol_O: TOL_DEG_KJ },
    reporting_detection_limit: {
      moles: 'exact zeros for excluded phases; doubles elsewhere',
      G_rel_kJ_floor: gFloor },
    converged: true,
    iterations: win.iterations,
    newton_residual: win.newton_residual,
    candidates: results.map(function (r) {
      return { active: r.active, feasible: r.feasible, kkt_ok: r.kkt_ok,
               reason: r.reason,
               worst_inactive_r_kJ: r.worst_inactive_r !== undefined
                 ? r.worst_inactive_r : null,
               solid_moles: r.m || null };
    }),
  };
};

/* ------------------------------------------------- the reduction line

   Mirror of the same section in solidgas/activeset.py. Titanium never
   enters the gas here, so solid and gas trade nothing but oxygen and the
   whole question is one comparison: does the gas sit above or below the
   potential at which the first reduced phase replaces the host. Both
   sides are closed form, which is why this curve is free where bisecting
   the solver over feed composition would be sixty solves a point.

   Free O2 is carried. With hydrogen in the feed it is worth nothing; take
   the hydrogen away and a CO2/CO mixture at 1650 C runs 0.2 mol% O2 and
   it is worth 5.8 kJ per mol O. It brings the total pressure in with it,
   since a partial pressure is not a property of a mixture on its own. */

function sigma(x) {
  if (x >= 0) return 1 / (1 + Math.exp(-x));
  var e = Math.exp(x);
  return e / (1 + e);
}

Solver.prototype.couplePotential = function (red, ox, T) {
  return this.mu0Gas(ox, T) - this.mu0Gas(red, T);
};

/* Moles of free O2 per mole of everything that is not O2. Null where the
   arithmetic says the gas would be all O2. */
Solver.prototype.o2PerUnit = function (mu, T, P) {
  var k = Math.exp(Math.min((2 * mu - this.mu0Gas('O2', T)) / (this.R * T),
                            EXP_CAP));
  var f = k / P;
  if (f >= 1) return null;
  return f / (1 - f);
};

/* Oxygen potential of a C-H-O gas equilibrated on its own, kJ per mol O.
   The feed before the solid touches it. Null when no couple exists, or
   when the feed is short of enough oxygen to make every carbon a CO. */
Solver.prototype.gasOxygenPotential = function (feedRel, T_K, P, iters) {
  var D = this.D, nC = 0, nH = 0, nO = 0, nI = 0, self = this;
  P = P || 1;
  Object.keys(feedRel).forEach(function (g) {
    var amt = feedRel[g];
    if (!(amt > 0)) return;
    var f = D.gas_formula[g];
    nC += amt * (f.C || 0);
    nH += amt * (f.H || 0);
    nO += amt * (f.O || 0);
    if (!f.C && !f.H && !f.O) nI += amt;
  });
  if (nC <= 0 && nH <= 0) return null;
  if (nO <= nC) return null;

  var rt = this.R * T_K;
  var a = this.couplePotential('CO', 'CO2', T_K);
  var b = this.couplePotential('H2', 'H2O', T_K);
  var base = nC + nH / 2 + nI;
  var oxygenAt = function (mu) {
    var w = self.o2PerUnit(mu, T_K, P);
    return nC * (2 - sigma((a - mu) / rt))
         + (nH / 2) * sigma((mu - b) / rt)
         + (w === null ? Infinity : 2 * base * w);
  };

  var lo = Math.min(a, b) - 300 * rt, hi = Math.max(a, b) + 300 * rt, mid;
  if (oxygenAt(hi) < nO) return null;
  var n = iters || 200, i;
  for (i = 0; i < n; i++) {
    mid = 0.5 * (lo + hi);
    if (oxygenAt(mid) < nO) lo = mid; else hi = mid;
  }
  return 0.5 * (lo + hi);
};

/* The H2/CO2 feed that carries a given oxygen potential, as a mole
   fraction: the same oxygen balance run backwards, so for a feed that
   really is H2 and CO2 it returns the fraction it was mixed at. For
   anything else, the mixture that would be equally reducing. */
Solver.prototype.equivalentCO2Fraction = function (muO, T_K, P) {
  var rt = this.R * T_K;
  P = P || 1;
  var sc = sigma((this.couplePotential('CO', 'CO2', T_K) - muO) / rt);
  var sh = sigma((muO - this.couplePotential('H2', 'H2O', T_K)) / rt);
  var w = this.o2PerUnit(muO, T_K, P);
  if (w === null) return 1;
  var den = sc + sh;
  if (!(den > 0)) return 1;
  return Math.min(Math.max((sh + 2 * w) / den, 0), 1);
};

/* Oxygen potential at which `higher` gives way to `lower`, kJ per mol O. */
Solver.prototype.criticalPotential = function (lower, higher, T) {
  var sh = this.D.stoich[higher], sl = this.D.stoich[lower];
  var nHi = sl[0] / sh[0];
  var dO = nHi * sh[1] - sl[1];
  if (!(dO > 0)) throw new Error(lower + ' is not more reduced than ' + higher);
  return (nHi * this.mu0Solid(higher, T) - this.mu0Solid(lower, T)) / dO;
};

/* The feed at which the host stops surviving, at one temperature. Every
   reduced phase is compared straight to the host, because the question is
   which one appears first out of an untouched host. */
Solver.prototype.reductionBoundary = function (T_K, host, solids, P) {
  var D = this.D, self = this;
  host = host || 'TiO2';
  var reg = solids || D.solids;
  var rhoHost = D.stoich[host][1] / D.stoich[host][0];
  var first = null, muCrit = -Infinity;
  reg.forEach(function (s) {
    if (D.stoich[s][1] / D.stoich[s][0] >= rhoHost) return;
    var c = self.criticalPotential(s, host, T_K);
    if (c > muCrit) { first = s; muCrit = c; }
  });
  if (first === null) return null;
  return { phase: first, host: host, mu_crit_kJ_per_mol_O: muCrit,
           y_CO2: this.equivalentCO2Fraction(muCrit, T_K, P || 1) };
};

/* Where a feed lands on that same axis, with the two ways it can fail to
   land anywhere named rather than swallowed. `literal` marks the feeds
   the axis was built for - H2 and CO2, plus any inert - where the marker
   is the mixing fraction itself and not an equivalent. */
Solver.prototype.feedOnBoundaryAxis = function (feedRel, T_K, P) {
  var D = this.D;
  P = P || 1;
  var mu = this.gasOxygenPotential(feedRel, T_K, P);
  if (mu === null) {
    var anyC = false, anyH = false;
    Object.keys(feedRel).forEach(function (g) {
      if (!(feedRel[g] > 0)) return;
      if (D.gas_formula[g].C) anyC = true;
      if (D.gas_formula[g].H) anyH = true;
    });
    return { y_CO2: null, mu_O_kJ_per_mol: null, literal: false,
             why: (anyC || anyH)
               ? 'this feed is short of the oxygen to make every carbon a CO'
               : 'this feed carries no CO/CO₂ or H₂/H₂O couple' };
  }
  var literal = Object.keys(feedRel).every(function (g) {
    if (!(feedRel[g] > 0)) return true;
    var f = D.gas_formula[g];
    return g === 'CO2' || g === 'H2' || (!f.C && !f.H && !f.O);
  });
  return { y_CO2: this.equivalentCO2Fraction(mu, T_K, P),
           mu_O_kJ_per_mol: mu, literal: literal, why: null };
};


root.Solver = Solver;
})(typeof module !== 'undefined' && module.exports ? module.exports
   : (window.ActiveSet = {}));
