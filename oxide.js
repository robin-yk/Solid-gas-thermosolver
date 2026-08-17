/* Gas-solid equilibrium by full-composition Gibbs minimisation.

   The engine is the Python package `solidgas` ported one to one:
   reaction-extent coordinates, the constant part of G removed analytically,
   and an analytic gradient. It is not an oxygen-potential shortcut and it is
   not a bisection on a phase pair - it minimises the total Gibbs energy over
   the whole species set at fixed elemental balance, and the phases fall out.

   What it cannot do is stated on the page next to the results, not here.

   Verified against the Python reference to five significant figures at every
   condition the page ships with; `tests/test_oxide_port.py` holds that. */

const Oxide = (function () {
  'use strict';

  const D = JSON.parse(document.getElementById('oxide-data').textContent);
  const R = D.R_KJ;                       // kJ / mol / K
  const LOG_FLOOR = -110;                 // e^-110 = 1.7e-48
  const LOG_CEIL = Math.log(10);

  /* ---------------------------------------------------------------- thermo */

  /* Body copied from the WGS tools page. One token differs: the coefficients
     arrive as an argument instead of through `thermo[s].c`, because this data
     carries a low and a high range per species where that page carries one.
     The arithmetic is untouched, and the formation term is folded into F by
     export_oxide.py so this still returns the absolute potential. */
  function shomate(c, T) {
    const [A, B, C, D_, E, F, G] = c, t = T / 1000,
      Cp = A + B * t + C * t * t + D_ * t ** 3 + E / (t * t),
      H = A * t + B * t * t / 2 + C * t ** 3 / 3 + D_ * t ** 4 / 4 - E / t + F,
      S = A * Math.log(t) + B * t + C * t * t / 2 + D_ * t ** 3 / 3 - E / (2 * t * t) + G;
    return { Cp, H, S, mu: H * 1000 - T * S };
  }

  /* Sibling, not a generalisation. CALPHAD condensed phases use a different Cp
     functional form and are integrated from 298.15 K with a transition
     enthalpy entered at each range boundary; folding that into shomate() would
     break the gas path, which is the one thing here that is already checked. */
  function cpInt(form, c, T0, T1) {
    const F = form === 'poly'
      ? T => c[0] * T + c[1] * T * T / 2 + c[2] * T ** 3 / 3 - c[3] / T
      : T => c[0] * T + 2 * c[1] * Math.sqrt(T) - c[2] / T - c[3] / (2 * T * T);
    return F(T1) - F(T0);
  }

  function cpIntOverT(form, c, T0, T1) {
    const F = form === 'poly'
      ? T => c[0] * Math.log(T) + c[1] * T + c[2] * T * T / 2 - c[3] / (2 * T * T)
      : T => c[0] * Math.log(T) - 2 * c[1] / Math.sqrt(T) - c[2] / (2 * T * T)
             - c[3] / (3 * T ** 3);
    return F(T1) - F(T0);
  }

  function segments(e, T) {
    const out = []; let lo = 298.15;
    for (const r of e.ranges) {
      if (T <= r.T0) break;
      const hi = Math.min(T, r.T1);
      if (hi > Math.max(lo, r.T0)) { out.push({ t0: Math.max(lo, r.T0), t1: hi, r }); lo = hi; }
      if (T <= r.T1) break;
    }
    if (out.length && out[out.length - 1].t1 < T) {          // extrapolate the last
      const r = e.ranges[e.ranges.length - 1];
      out.push({ t0: out[out.length - 1].t1, t1: T, r: Object.assign({}, r, { dHt: 0 }) });
    }
    return out;
  }

  function calphad(e, T) {
    let H = e.dHf298, S = e.S298;
    for (const s of segments(e, T)) {
      H += s.r.dHt + cpInt(s.r.form, s.r.c, s.t0, s.t1);
      S += (s.r.dHt ? s.r.dHt / s.r.T0 : 0) + cpIntOverT(s.r.form, s.r.c, s.t0, s.t1);
    }
    return { H, S, mu: (H - T * S) / 1000 };                 // kJ/mol
  }

  /* Chemical potential in kJ/mol, whichever table the species lives in. */
  function mu0(name, T) {
    if (D.calphad[name]) return calphad(D.calphad[name], T).mu;
    const e = D.shomate[name];
    return shomate(T <= e.brk ? e.lo : e.hi, T).mu / 1000;
  }

  /* ------------------------------------------------------------ small algebra */

  function inv(M) {
    const n = M.length, A = M.map((r, i) => r.concat(M.map((_, j) => (i === j ? 1 : 0))));
    for (let c = 0; c < n; c++) {
      let p = c;
      for (let r = c + 1; r < n; r++) if (Math.abs(A[r][c]) > Math.abs(A[p][c])) p = r;
      if (Math.abs(A[p][c]) < 1e-12) throw new Error('component set is singular');
      [A[c], A[p]] = [A[p], A[c]];
      const d = A[c][c];
      for (let j = 0; j < 2 * n; j++) A[c][j] /= d;
      for (let r = 0; r < n; r++) {
        if (r === c) continue;
        const f = A[r][c];
        if (f) for (let j = 0; j < 2 * n; j++) A[r][j] -= f * A[c][j];
      }
    }
    return A.map(r => r.slice(n));
  }

  /* ---------------------------------------------------------------- system */

  class GibbsSystem {
    constructor(gases, solids, components) {
      this.gases = gases.slice();
      this.solids = solids.slice();
      this.species = this.gases.concat(this.solids);

      const els = [];
      for (const s of this.species)
        for (const e of Object.keys(D.formulas[s])) if (!els.includes(e)) els.push(e);
      this.elements = els.sort();

      this.E = this.species.map(s => this.elements.map(e => D.formulas[s][e] || 0));
      this.isGas = this.species.map(s => this.gases.includes(s));
      this.comp = components.map(c => this.species.indexOf(c));
      this.free = this.species.map((_, i) => i).filter(i => !this.comp.includes(i));
      if (this.comp.length !== this.elements.length)
        throw new Error('need one component per element');

      // Ek is elements x components; invert it once.
      const Ek = this.elements.map((_, e) => this.comp.map(i => this.E[i][e]));
      this.Ainv = inv(Ek);
      // dK[c][f] : how much of component c one unit of free species f consumes.
      this.dK = this.comp.map((_, ci) => this.free.map((fi) => {
        let v = 0;
        for (let e = 0; e < this.elements.length; e++) v += this.Ainv[ci][e] * this.E[fi][e];
        return -v;
      }));
    }

    /* Element vector for a charge written as species moles. Built by name: the
       element order is whatever sorting gave, and it differs between Ti-O and
       Ce-O because Ce sorts before H. */
    bFrom(moles) {
      const b = this.elements.map(() => 0);
      for (const name of Object.keys(moles))
        for (const [e, k] of Object.entries(D.formulas[name]))
          b[this.elements.indexOf(e)] += k * moles[name];
      return b;
    }

    n0(b) {
      const n = this.species.map(() => 0);
      const mx = Math.max(...b.map(Math.abs), 1e-300);
      this.comp.forEach((si, ci) => {
        let v = 0;
        for (let e = 0; e < this.elements.length; e++) v += this.Ainv[ci][e] * b[e];
        n[si] = Math.abs(v) < 1e-14 * mx ? 0 : v;            // snap the round-off
      });
      return n;
    }

    deltaN(xi) {
      const d = this.species.map(() => 0);
      this.free.forEach((si, f) => { d[si] = xi[f]; });
      this.comp.forEach((si, ci) => {
        let v = 0;
        for (let f = 0; f < this.free.length; f++) v += this.dK[ci][f] * xi[f];
        d[si] = v;
      });
      return d;
    }

    nOf(xi, b) { const z = this.n0(b), d = this.deltaN(xi); return z.map((v, i) => v + d[i]); }

    /* Reaction Gibbs energy for forming each free species from the components.
       This is the whole standard-state dependence of G on the extents, formed
       directly rather than as a difference of two much larger sums - which is
       what lets an extent of 1e-38 be seen at all. */
    dGForm(T) {
      const mu = this.species.map(s => mu0(s, T));
      return this.free.map((si, f) => {
        let v = mu[si];
        this.comp.forEach((ci, k) => { v += this.dK[k][f] * mu[ci]; });
        return v;
      });
    }

    /* Change in the gas mixing term, in units of RT. Every piece goes through
       log1p rather than a difference of two logs, so a small change against a
       large reference survives. */
    dmix(n0, dn, n, P) {
      let s = 0, ng0 = 0, d = 0;
      for (let i = 0; i < this.species.length; i++) {
        if (!this.isGas[i]) continue;
        if (n[i] < 0) return Infinity;
        ng0 += n0[i]; d += dn[i];
        if (n0[i] > 0 && Math.abs(dn[i]) < 0.5 * n0[i]) {
          s += n0[i] * Math.log1p(dn[i] / n0[i]) + dn[i] * Math.log(n[i]);
        } else {
          if (n[i] > 0) s += n[i] * Math.log(n[i]);
          if (n0[i] > 0) s -= n0[i] * Math.log(n0[i]);
        }
      }
      const ng = ng0 + d;
      if (ng <= 0) return Infinity;
      if (ng0 > 0 && Math.abs(d) < 0.5 * ng0) s -= ng0 * Math.log1p(d / ng0) + d * Math.log(ng);
      else s -= ng * Math.log(ng) - (ng0 > 0 ? ng0 * Math.log(ng0) : 0);
      return s + d * Math.log(P);
    }

    /* Gibbs energy relative to the unreacted charge, in kJ. */
    GRel(xi, b, T, P) {
      const z = this.n0(b), dn = this.deltaN(xi), n = z.map((v, i) => v + dn[i]);
      const floor = -1e-12 * Math.max(1, Math.max(...z.map(Math.abs)));
      for (let i = 0; i < n.length; i++) if (n[i] < floor) return 1e10;
      for (let i = 0; i < n.length; i++) if (n[i] < 0) n[i] = 0;
      const m = this.dmix(z, dn, n, P);
      if (!isFinite(m)) return 1e10;
      const g = this.dGForm(T);
      let lin = 0;
      for (let f = 0; f < xi.length; f++) lin += g[f] * xi[f];
      return lin + R * T * m;
    }

    /* dG/dxi, analytic. The d(ng)/dxi terms cancel identically because the mole
       fractions sum to one. Finite differences are useless here - the whole of
       GRel can sit at 1e-19 kJ - so the gradient is written out. */
    gradGRel(xi, b, T, P) {
      const n = this.nOf(xi, b).map(v => (v > 0 ? v : 0));
      let ng = 0;
      for (let i = 0; i < n.length; i++) if (this.isGas[i]) ng += n[i];
      const ly = n.map((v, i) => (this.isGas[i] && ng > 0
        ? Math.log(Math.max(v, 1e-300) / ng * P) : 0));
      const g = this.dGForm(T);
      return this.free.map((si, f) => {
        let v = g[f];
        if (this.isGas[si]) v += R * T * ly[si];
        this.comp.forEach((ci, k) => {
          if (this.isGas[ci]) v += R * T * this.dK[k][f] * ly[ci];
        });
        return v;
      });
    }

    gasFractions(n) {
      let tot = 0;
      for (let i = 0; i < n.length; i++) if (this.isGas[i]) tot += n[i];
      const y = {};
      this.species.forEach((s, i) => { if (this.isGas[i]) y[s] = tot > 0 ? n[i] / tot : 0; });
      return y;
    }

    reducedPercent(n, nMetal) {
      if (nMetal <= 0) return 0;
      let red = 0;
      this.solids.forEach(p => {
        const m = D.metal_phases[p];
        if (m) red += n[this.species.indexOf(p)] * m[1];
      });
      return Math.min(red / nMetal, 1) * 100;
    }

    phaseSplit(n, nMetal) {
      const out = {};
      if (nMetal <= 0) return out;
      this.solids.forEach(p => {
        const m = D.metal_phases[p];
        if (m) out[p] = n[this.species.indexOf(p)] * m[0] / nMetal * 100;
      });
      return out;
    }
  }

  /* -------------------------------------------------------------- optimiser */

  const clamp = (v, a, b) => (v < a ? a : v > b ? b : v);

  /* Projected gradient with a Barzilai-Borwein step and backtracking, over the
     log of the extents. SciPy's SLSQP does not exist in a browser and this is
     what replaces it. Two things make the substitution safe: the change of
     variables has already removed every equality constraint, so only the
     component non-negativities remain, and GRel returns 1e10 on an infeasible
     composition - so a line search that only accepts a decrease can never step
     outside the feasible set.

     Reports its own iteration count and gradient norm; the page shows them. */
  function descend(sys, u0, b, T, P, maxIter) {
    const nf = u0.length;
    let u = u0.slice(), f = sys.GRel(u.map(Math.exp), b, T, P);
    if (!isFinite(f) || f >= 1e9) return null;
    let g = sys.gradGRel(u.map(Math.exp), b, T, P).map((v, i) => v * Math.exp(u[i]));
    let step = 1e-2, pu = null, pg = null, it = 0, g0 = null;
    for (; it < maxIter; it++) {
      if (pu) {
        let sy = 0, ss = 0;
        for (let i = 0; i < nf; i++) {
          const ds = u[i] - pu[i], dg = g[i] - pg[i];
          sy += ds * dg; ss += ds * ds;
        }
        if (sy > 1e-300 && isFinite(ss / sy)) step = clamp(ss / sy, 1e-12, 1e6);
      }
      let t = step, moved = false;
      for (let k = 0; k < 80; k++) {
        const un = u.map((v, i) => clamp(v - t * g[i], LOG_FLOOR, LOG_CEIL));
        const fn = sys.GRel(un.map(Math.exp), b, T, P);
        if (isFinite(fn) && fn < f) { pu = u; pg = g; u = un; f = fn; moved = true; break; }
        t *= 0.5;
      }
      if (!moved) break;
      g = sys.gradGRel(u.map(Math.exp), b, T, P).map((v, i) => v * Math.exp(u[i]));
      let gn = 0; for (let i = 0; i < nf; i++) gn += g[i] * g[i];
      gn = Math.sqrt(gn);
      if (g0 === null) g0 = gn;
      if (gn < 1e-12 * Math.max(g0, 1e-300)) break;
    }
    let gf = 0; for (let i = 0; i < nf; i++) gf += g[i] * g[i];
    return { u, f, it, gnorm: Math.sqrt(gf) };
  }

  /* Minimise from every seed and keep the lowest. The basis is scaled to order
     one first: G is homogeneous of degree one in n, so that is exact, and it
     keeps the tolerances meaningful against a charge of a millimole. */
  function minimise(sys, b, T, seeds, P, maxIter) {
    const scale = b.reduce((a, v) => a + Math.abs(v), 0);
    if (!(scale > 0)) throw new Error('empty charge');
    const bs = b.map(v => v / scale);
    let best = null, tried = 0, iters = 0;

    // Coarse pass: every seed, few iterations. This is only choosing a basin -
    // running any of them to convergence would be waste, since all but one gets
    // thrown away.
    for (const seed of seeds) {
      const u0 = seed.map(v => clamp(Math.log(Math.max(v / scale, Math.exp(LOG_FLOOR))),
                                     LOG_FLOOR, LOG_CEIL));
      const r = descend(sys, u0, bs, T, P, 150);
      tried++;
      if (!r) continue;
      iters += r.it;
      if (best === null || r.f < best.f) best = r;
    }
    if (!best) return null;

    // Polish: the winning basin, run out properly. Projected gradient gets the
    // last digits more slowly than SLSQP does, and the port has to match the
    // Python answer to five figures, so it is given the room to.
    //
    // Restarted rather than run once. The Barzilai-Borwein step is built from
    // the previous two iterates, and when it stalls it stalls with a step that
    // no amount of backtracking rescues; beginning again from the same point
    // throws that history away and the descent resumes. It stops when a whole
    // restart buys nothing.
    for (let round = 0; round < 6; round++) {
      const pol = descend(sys, best.u, bs, T, P, maxIter || 3000);
      if (!pol) break;
      iters += pol.it;
      const gained = best.f - pol.f;
      if (pol.f < best.f) best = pol;
      if (!(gained > 1e-15 * Math.abs(best.f))) break;
    }
    const xi = best.u.map(v => Math.exp(clamp(v, LOG_FLOOR, LOG_CEIL)) * scale);
    const n = sys.nOf(xi, b).map(v => (v > 0 ? v : 0));
    let resid = 0;
    for (let e = 0; e < sys.elements.length; e++) {
      let v = 0;
      for (let i = 0; i < sys.species.length; i++) v += sys.E[i][e] * n[i];
      resid = Math.max(resid, Math.abs(v - b[e]));
    }
    return {
      n, xi, G_rel_kJ: best.f * scale, residual: resid,
      method: 'gibbs_min', formulation: 'reaction_extent',
      seeds: tried, iterations: iters, gnorm: best.gnorm,
    };
  }

  /* Extent starts: a ladder of gas conversion crossed with each phase basin.
     A start within a couple of decades of the answer is what the optimiser
     needs, and there is no way to know in advance which decade. */
  function buildSeeds(sys, ng, nMetal) {
    const ref = sys.free.map(i => (sys.isGas[i] ? ng : nMetal));
    const gasPos = sys.free.map((i, j) => (sys.isGas[i] ? j : -1)).filter(j => j >= 0);
    const solPos = sys.free.map((i, j) => (sys.isGas[i] ? -1 : j)).filter(j => j >= 0);
    const out = [];
    for (const gf of [1e-6, 0.05, 0.2, 0.4, 0.7]) {
      const base = sys.free.map(() => 1e-48);
      gasPos.forEach(j => { base[j] = ref[j] * gf; });
      out.push(base.slice());
      for (const j of solPos)
        for (const m of [1e-1, 1e-4, 1e-8, 1e-14, 1e-22, 1e-32, 1e-42]) {
          const v = base.slice(); v[j] = ref[j] * m; out.push(v);
        }
    }
    return out;
  }

  /* ------------------------------------------------------- verification only

     A second, independent route to the same answer, for the verification tab
     and nowhere else. It is not a minimisation: the metal does not enter the
     gas below about 2000 K, so the solid exchanges nothing with it but oxygen,
     and each condensed phase's critical oxygen potential becomes a plain
     function of temperature with no optimiser in it.

     That shortcut is legitimate here and it is fast, and it is still not what
     the page reports. Two phases in contact fix the gas oxygen potential
     exactly, so where the minimisation puts two phases together its gas has to
     reproduce a number computed from the phase data alone - which is a test of
     the minimiser using nothing the minimiser produced. Reporting the shortcut
     instead would turn that test into a number compared against itself. */
  const Verify = {
    /* Oxygen potential of a mixture, per mole of O, from the redox couples.
       Molecular O2 is not used: in a C-H-O gas its equilibrium amount is around
       1e-40 and sits at the optimiser's floor, so its value is arbitrary. */
    muOFromGas(y, T) {
      const vals = [];
      for (const [red, ox] of [['H2', 'H2O'], ['CO', 'CO2']]) {
        const a = y[red] || 0, b = y[ox] || 0;
        if (a > 1e-280 && b > 1e-280)
          vals.push(mu0(ox, T) - mu0(red, T) + R * T * Math.log(b / a));
      }
      if (!vals.length) return { muO: null, spread: null };
      return {
        muO: vals.reduce((a, v) => a + v, 0) / vals.length,
        spread: Math.max(...vals) - Math.min(...vals),
      };
    },

    counts(name, metal) {
      const f = D.formulas[name];
      return [f[metal] || 0, f.O || 0];
    },

    /* Critical potential for higher -> lower, per mole of O removed. Both are
       line compounds, so this is a pure function of temperature. */
    muOCritical(lower, higher, T, metal) {
      const [th, oh] = this.counts(higher, metal), [tl, ol] = this.counts(lower, metal);
      const nHi = tl / th, dO = nHi * oh - ol;
      if (dO <= 0) throw new Error(lower + ' is not more reduced than ' + higher);
      return (nHi * mu0(higher, T) - mu0(lower, T)) / dO;
    },

    ladder(solids, T, metal) {
      const ord = solids.slice().sort((a, b) => {
        const [ta, oa] = this.counts(a, metal), [tb, ob] = this.counts(b, metal);
        return ob / tb - oa / ta;
      });
      const out = [];
      for (let i = 0; i < ord.length - 1; i++)
        out.push({ higher: ord[i], lower: ord[i + 1],
                   crit: this.muOCritical(ord[i + 1], ord[i], T, metal) });
      return out;
    },

    /* pO2, in atm, carrying a given oxygen potential. The mu0(O2) term is the
       one that is easy to drop, and dropping it moves this by fourteen orders
       of magnitude - though never a margin, which is a difference and cancels. */
    muOToPO2(muO, T) { return Math.exp((2 * muO - mu0('O2', T)) / (R * T)); },
  };

  return { D, R, mu0, shomate, calphad, GibbsSystem, minimise, buildSeeds,
           Verify, LOG_FLOOR };
})();
