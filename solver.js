/* Thermodynamics for the browser. Same data and same equations as the Python
   package - the coefficients arrive from thermo_data.json rather than being
   retyped, and solver_check.js pins the answers against it. */

(function (root) {
'use strict';

function Engine(D) {
  this.D = D;
  this.R = D.R_KJ;
  this.RATM = D.R_ATM;
}

/* ---------------------------------------------------------------- NIST side */

Engine.prototype.coeffs = function (name, T) {
  var e = this.D.shomate[name];
  return T <= e.brk ? e.lo : e.hi;
};

Engine.prototype.nistMu0 = function (name, T) {
  var e = this.D.shomate[name], c = this.coeffs(name, T), t = T / 1000;
  var H = e.dHf + (c[0]*t + c[1]*t*t/2 + c[2]*t*t*t/3 + c[3]*t*t*t*t/4 - c[4]/t + c[5] - c[7]);
  var S = c[0]*Math.log(t) + c[1]*t + c[2]*t*t/2 + c[3]*t*t*t/3 - c[4]/(2*t*t) + c[6];
  return H - T * S / 1000;
};

/* ------------------------------------------------------------- Waldner side */

function intCp(r, T0, T1) {
  var c = r.c;
  var F = r.form === 'poly'
    ? function (T) { return c[0]*T + c[1]*T*T/2 + c[2]*T*T*T/3 - c[3]/T; }
    : function (T) { return c[0]*T + 2*c[1]*Math.sqrt(T) - c[2]/T - c[3]/(2*T*T); };
  return F(T1) - F(T0);
}
function intCpOverT(r, T0, T1) {
  var c = r.c;
  var F = r.form === 'poly'
    ? function (T) { return c[0]*Math.log(T) + c[1]*T + c[2]*T*T/2 - c[3]/(2*T*T); }
    : function (T) { return c[0]*Math.log(T) - 2*c[1]/Math.sqrt(T) - c[2]/(2*T*T) - c[3]/(3*T*T*T); };
  return F(T1) - F(T0);
}

Engine.prototype.segments = function (name, T) {
  var rs = this.D.waldner[name].ranges, out = [], lo = 298.15, i;
  for (i = 0; i < rs.length; i++) {
    var r = rs[i];
    if (T <= r.lo) break;
    var hi = Math.min(T, r.hi), from = Math.max(lo, r.lo);
    if (hi > from) { out.push({ r: r, t0: from, t1: hi, Ttr: r.lo }); lo = hi; }
    if (T <= r.hi) break;
  }
  if (out.length && out[out.length - 1].t1 < T) {
    var last = rs[rs.length - 1];
    out.push({ r: { form: last.form, c: last.c, dHt: 0 },
               t0: out[out.length - 1].t1, t1: T, Ttr: last.lo });
  }
  return out;
};

Engine.prototype.waldnerMu0 = function (name, T) {
  var e = this.D.waldner[name], segs = this.segments(name, T);
  var H = e.dHf, S = e.S298, i;
  for (i = 0; i < segs.length; i++) {
    var s = segs[i];
    H += s.r.dHt + intCp(s.r, s.t0, s.t1);
    S += (s.r.dHt ? s.r.dHt / s.Ttr : 0) + intCpOverT(s.r, s.t0, s.t1);
  }
  return (H - T * S) / 1000;
};

Engine.prototype.solidMu0 = function (name, T, source) {
  if ((source || 'waldner') === 'waldner' && this.D.waldner[name]) return this.waldnerMu0(name, T);
  return this.nistMu0(name, T);
};

/* ------------------------------------------------------------ oxygen potential */

Engine.prototype.stoich = function (name) {
  var f = this.D.formulas[name];
  return { Ti: f.Ti || 0, O: f.O || 0 };
};

Engine.prototype.muOCriticalPair = function (higher, lower, T, source) {
  var hi = this.stoich(higher), lo = this.stoich(lower);
  var n = lo.Ti / hi.Ti, dO = n * hi.O - lo.O;
  if (dO <= 0) return NaN;
  return (n * this.solidMu0(higher, T, source) - this.solidMu0(lower, T, source)) / dO;
};

Engine.prototype.muOCritical = function (phase, T, source) {
  var hi = this.stoich('TiO2'), lo = this.stoich(phase);
  var n = lo.Ti / hi.Ti, dO = n * hi.O - lo.O;
  if (dO <= 0) return NaN;
  return (n * this.solidMu0('TiO2', T, source) - this.solidMu0(phase, T, source)) / dO;
};

Engine.prototype.muOFromRatio = function (T, ratio) {  /* ratio = H2O/H2 */
  return this.nistMu0('H2O', T) - this.nistMu0('H2', T) + this.R * T * Math.log(ratio);
};

Engine.prototype.muOToPO2 = function (muO, T) {
  return Math.exp((2 * muO - this.nistMu0('O2', T)) / (this.R * T));
};
Engine.prototype.pO2ToMuO = function (pO2, T) {
  return 0.5 * (this.nistMu0('O2', T) + this.R * T * Math.log(pO2));
};

/* Which reduced phase forms first out of untouched rutile, and by what margin. */
Engine.prototype.firstPhase = function (muO, T, phases, source) {
  var best = null, bestCrit = -Infinity, i;
  for (i = 0; i < phases.length; i++) {
    if (phases[i] === 'TiO2') continue;
    var c = this.muOCritical(phases[i], T, source);
    if (c > bestCrit) { bestCrit = c; best = phases[i]; }
  }
  return { phase: best, crit: bestCrit, margin: muO - bestCrit };
};

/* ------------------------------------------------------ C-H-O gas equilibrium

   Five species, three elements, so two independent reactions:
       CO2 + H2  <-> CO + H2O          (no change in moles)
       CO  + 3H2 <-> CH4 + H2O         (loses two)
   Solved as a damped Newton on the two extents, in the variables that stay
   positive. Seeded from the shift reaction alone, which is a quadratic. */

Engine.prototype.gasEquilibrium = function (bC, bH, bO, T, P) {
  /* No carbon means no reaction left to equilibrate: H2 and H2O are fixed by
     the element balance alone. Falling through to the general path would ask
     for a methane extent between zero and zero. */
  if (bC <= 0) {
    var nH2 = bH / 2 - bO, nH2O = bO;
    if (nH2 < 0) { nH2O = bH / 2; nH2 = 0; }   /* more oxygen than hydrogen */
    var Nt = nH2 + nH2O;
    if (!(Nt > 0)) return null;
    var nn = { CO2: 0, H2: nH2, CO: 0, H2O: nH2O, CH4: 0 };
    return { n: nn, ntot: Nt, y: { CO2: 0, H2: nH2 / Nt, CO: 0, H2O: nH2O / Nt, CH4: 0 },
             K1: NaN, K2: NaN };
  }
  var K1 = Math.exp(-(this.nistMu0('CO', T) + this.nistMu0('H2O', T)
                    - this.nistMu0('CO2', T) - this.nistMu0('H2', T)) / (this.R * T));
  var K2 = Math.exp(-(this.nistMu0('CH4', T) + this.nistMu0('H2O', T)
                    - this.nistMu0('CO', T) - 3 * this.nistMu0('H2', T)) / (this.R * T));

  /* Composition from the two extents. c is CO, e is CH4; the rest follow from
     the element balance, which is linear. */
  function comp(c, e) {
    return { CO2: bC - c - e,
             H2:  bH / 2 - bO + 2 * bC - c - 4 * e,
             CO:  c,
             H2O: bO - 2 * bC + c + 2 * e,
             CH4: e };
  }
  /* Range of c that keeps every species positive, at fixed e. */
  function cRange(e) {
    var lo = Math.max(0, 2 * bC - bO - 2 * e);
    var hi = Math.min(bC - e, bH / 2 - bO + 2 * bC - 4 * e);
    return [lo, hi];
  }
  /* Shift equilibrium at fixed e. Monotonic in c, so bisection cannot miss it. */
  function solveShift(e) {
    var r = cRange(e), lo = r[0], hi = r[1], i, mid, n;
    if (!(hi > lo)) return null;
    var pad = (hi - lo) * 1e-14;
    lo += pad; hi -= pad;
    for (i = 0; i < 200; i++) {
      mid = 0.5 * (lo + hi);
      n = comp(mid, e);
      if (n.CO * n.H2O - K1 * n.CO2 * n.H2 > 0) hi = mid; else lo = mid;
    }
    return 0.5 * (lo + hi);
  }
  /* Methanation residual along that curve. Negative with no methane, positive
     once carbon monoxide or hydrogen is exhausted, so it brackets a root. */
  function f2(e) {
    var c = solveShift(e);
    if (c === null) return null;
    var n = comp(c, e), N = n.CO2 + n.H2 + n.CO + n.H2O + n.CH4;
    return { v: n.CH4 * n.H2O * N * N - K2 * n.CO * Math.pow(n.H2, 3) * P * P,
             c: c, n: n, N: N };
  }

  /* The feasible window in e is not always [0, eMax]. A methane-rich charge
     with little oxygen has no solution at e = 0 at all: the element balance
     needs nearly all the carbon sitting in CH4 before CO2, CO and H2O can be
     positive together. Bracket the window first, then bisect inside it. */
  var eMax = Math.min(bC, (bH / 2 - bO + 2 * bC) / 4);
  if (!(eMax > 0)) return null;
  var N = 400, i, k, r, grid = [], any = -1;
  for (k = 0; k <= N; k++) {
    var e = eMax * k / N * (1 - 1e-12);
    var ok = f2(e);
    grid.push(ok ? { e: e, v: ok.v } : null);
    if (ok && any < 0) any = k;
  }
  if (any < 0) return null;
  var last = any;
  for (k = any; k <= N; k++) if (grid[k]) last = k;

  /* narrow the feasible edges so the bracket sits strictly inside */
  function edge(kIn, kOut) {
    var a = grid[kIn].e, b = eMax * kOut / N;
    for (i = 0; i < 80; i++) {
      var m = 0.5 * (a + b);
      if (f2(m)) a = m; else b = m;
    }
    return a;
  }
  var eLo = any === 0 ? grid[0].e : edge(any, any - 1);
  var eHi = last === N ? grid[last].e : edge(last, last + 1);

  var lo = eLo, hi = eHi;
  var rLo = f2(lo), rHi = f2(hi);
  if (!rLo || !rHi) return null;
  if (rLo.v > 0) { hi = lo; }                 /* already past the root */
  else if (rHi.v < 0) { lo = hi; }            /* methane never turns over */
  else {
    for (i = 0; i < 200; i++) {
      var mid = 0.5 * (lo + hi);
      r = f2(mid);
      if (r === null || r.v > 0) hi = mid; else lo = mid;
    }
  }
  r = f2(0.5 * (lo + hi)) || rLo;
  var y = {}, self = this;
  ['CO2', 'H2', 'CO', 'H2O', 'CH4'].forEach(function (k) { y[k] = r.n[k] / r.N; });
  return { n: r.n, ntot: r.N, y: y, K1: K1, K2: K2 };
};

root.Engine = Engine;
})(typeof module !== 'undefined' && module.exports ? module.exports : (window.Thermo = {}));

/* ------------------------------------------------------------ buffered solve

   When a reduced phase is stable, the solid and gas exchange oxygen until the
   gas potential sits exactly on the phase boundary - the two-phase buffer. The
   unknown is how much oxygen moved, and the residual is monotonic in it, so
   bisection is enough and no optimiser is needed. */

(function (ns) {
var E = ns.Engine;

E.prototype.buffered = function (feed, T, P, V, nTi, phases, source) {
  var self = this, g = P * V / (this.RATM * T);

  /* Order the series by oxygen content and keep only the members that lie
     below rutile. Reduction walks down this ladder one pair at a time. */
  var lad = phases.slice().filter(function (p) { return self.D.formulas[p]; })
    .sort(function (a, b) {
      return self.stoich(b).O / self.stoich(b).Ti - self.stoich(a).O / self.stoich(a).Ti;
    });

  var bC = 0, bH = 0, bO = 0;
  if (feed === 'H2') { bH = 2 * g; }
  else if (feed === 'CH4') { bC = g; bH = 4 * g; }
  else if (feed === 'CO2H2') { bC = g / 2; bH = g; bO = g; }

  /* Oxygen the gas has taken, as a function of how far down the ladder we are.
     Converting everything to phase k costs nTi*(2 - O/Ti) in total. */
  function totalRemoved(k) {
    var st = self.stoich(lad[k]);
    return nTi * (2 - st.O / st.Ti);
  }
  function gasMuO(x) {
    var r = self.gasEquilibrium(bC, bH, bO + x, T, P);
    if (!r || !(r.y.H2 > 1e-25) || !(r.y.H2O > 1e-25)) return null;
    return self.muOFromRatio(T, r.y.H2O / r.y.H2);
  }
  /* Solve gasMuO(x) = target on [a, b]; the gas potential rises with x. */
  function crossing(a, bq, target) {
    var i, m, v;
    for (i = 0; i < 160; i++) {
      m = 0.5 * (a + bq);
      v = gasMuO(m);
      if (v === null) { bq = m; continue; }
      if (v < target) a = m; else bq = m;
    }
    return 0.5 * (a + bq);
  }

  var x = 0, hiPhase = lad[0], loPhase = null, crit = NaN, k, capped = false;
  for (k = 1; k < lad.length; k++) {
    crit = self.muOCriticalPair(lad[k - 1], lad[k], T, source);
    var xLo = totalRemoved(k - 1), xHi = totalRemoved(k);
    var mLo = gasMuO(xLo + (xHi - xLo) * 1e-12);
    if (mLo !== null && mLo > crit) {           /* cannot even start this step */
      hiPhase = lad[k - 1]; loPhase = lad[k]; x = xLo;
      return { phase: label(stateAt(x)), assemblage: stateAt(x), crit: crit, x: x,
               ti3: ti3At(x), conv: 0, O_fraction: x / (2 * nTi),
               muO: mLo, margin: mLo - crit, reduces: k > 1,
               solidLimited: false, gas: self.gasEquilibrium(bC, bH, bO + x, T, P) };
    }
    var mHi = gasMuO(xHi * (1 - 1e-12));
    if (mHi !== null && mHi < crit) { x = xHi; hiPhase = lad[k]; continue; }  /* pair used up */
    x = crossing(xLo, xHi, crit);
    hiPhase = lad[k - 1]; loPhase = lad[k];
    break;
  }
  if (k >= lad.length) { capped = true; x = totalRemoved(lad.length - 1); }

  /* Ti(3+) from the oxygen actually removed. Along a two-phase segment the mix
     is linear in x, and every Magneli member carries two Ti(3+) per formula. */
  /* Where x sits on the ladder: which two phases coexist, in what proportion,
     and the Ti(3+) that follows. The mix is linear in x along a segment, and
     every Magneli member carries two Ti(3+) per formula unit. */
  function stateAt(xx) {
    var i, prev = 0, prevPhase = lad[0];
    for (i = 1; i < lad.length; i++) {
      var lim = totalRemoved(i);
      if (xx <= lim + 1e-30) {
        var st = self.stoich(lad[i]), stp = self.stoich(prevPhase);
        var f = (lim - prev) > 0 ? (xx - prev) / (lim - prev) : 1;
        var ti3lo = 2 / st.Ti, ti3hi = prevPhase === 'TiO2' ? 0 : 2 / stp.Ti;
        return { ti3: (ti3hi + f * (ti3lo - ti3hi)) * 100,
                 upper: prevPhase, lower: lad[i], fraction: f };
      }
      prev = lim; prevPhase = lad[i];
    }
    var lastP = lad[lad.length - 1];
    return { ti3: 2 / self.stoich(lastP).Ti * 100, upper: lastP, lower: lastP, fraction: 1 };
  }
  function ti3At(xx) { return stateAt(xx).ti3; }

  /* Name the assemblage the way a phase diagram would: the dominant phase,
     with its partner when the two genuinely coexist. */
  function label(stt) {
    if (stt.fraction >= 1 - 1e-9) return stt.lower;
    if (stt.fraction <= 1e-9) return stt.upper;
    return stt.fraction < 0.5 ? stt.upper + ' + ' + stt.lower
                              : stt.lower + ' + ' + stt.upper;
  }

  var gas = self.gasEquilibrium(bC, bH, bO + x, T, P);
  var stt = stateAt(x);
  return { phase: label(stt), assemblage: stt, crit: crit, x: x,
           ti3: stt.ti3, O_fraction: x / (2 * nTi),
           conv: 0,
           muO: gas && gas.y.H2 > 1e-25 && gas.y.H2O > 1e-25
                ? self.muOFromRatio(T, gas.y.H2O / gas.y.H2) : crit,
           margin: 0, reduces: x > 0, solidLimited: capped, gas: gas };
};

/* Sealed inert headspace: the only place oxygen can go is O2 in the gas. */
E.prototype.inert = function (T, P, V, nTi, pO2Initial, phases, source) {
  var first = this.firstPhase(-Infinity, T, phases, source);
  var pEq = this.muOToPO2(first.crit, T);
  var nO2 = Math.max(pEq - (pO2Initial || 0), 0) * V / (this.RATM * T);
  var st = this.stoich(first.phase), dO = 2 * st.Ti - st.O;
  var nPhase = 2 * nO2 / dO;
  var capped = nPhase * st.Ti > nTi;
  if (capped) nPhase = nTi / st.Ti;
  return { phase: first.phase, crit: first.crit, pO2_eq: pEq, nO2: nO2,
           ti3: Math.min(2 * nPhase / nTi, 1) * 100,
           O_fraction: nPhase * dO / (2 * nTi),
           reduces: (pO2Initial || 0) < pEq, solidLimited: capped };
};
})(typeof module !== 'undefined' && module.exports ? module.exports : window.Thermo);
