/* Canonical defect statistical mechanics for oxygen vacancies in rutile.
   Line-for-line port of solidgas/statmech.py: same loops, same operation
   order, same sfc32 RNG, so a same-seed run is bit-reproducible against the
   Python reference (verified by tests/test_statmech_port.py). Runs on the
   page, inside a Worker, and under node. */

(function (root, factory) {
  'use strict';
  var api = factory();
  if (typeof module !== 'undefined' && module.exports) module.exports = api;
  root.StatMech = api;
})(typeof self !== 'undefined' ? self : this, function () {
  'use strict';

  var KB_EV = 8.617333262e-5;
  var N_AVO = 6.02214076e23;

  /* ------------------------------------------------------------- RNG */

  function makeRng(seed) {
    /* sfc32 seeded by splitmix32. State lives in a Uint32Array: V8 then
       keeps the words as raw int32 instead of heap-boxed closure slots,
       which is ~4x faster at an identical output stream. */
    var st = new Uint32Array(4);
    var state = seed >>> 0;
    for (var k = 0; k < 4; k++) {
      state = (state + 0x9E3779B9) >>> 0;
      var z = state;
      z = Math.imul(z ^ (z >>> 16), 0x21F0AAAD) >>> 0;
      z = Math.imul(z ^ (z >>> 15), 0x735A2D97) >>> 0;
      st[k] = (z ^ (z >>> 15)) >>> 0;
    }
    var rng = function () {
      var a = st[0], b = st[1], c = st[2], d = st[3];
      var t = (a + b) >>> 0;
      a = (b ^ (b >>> 9)) >>> 0;
      b = (c + ((c << 3) >>> 0)) >>> 0;
      c = (((c << 21) >>> 0) | (c >>> 11)) >>> 0;
      d = (d + 1) >>> 0;
      t = (t + d) >>> 0;
      c = (c + t) >>> 0;
      st[0] = a; st[1] = b; st[2] = c; st[3] = d;
      return t / 4294967296;
    };
    for (var i = 0; i < 12; i++) rng();
    return rng;
  }

  /* -------------------------------------------------------- geometry */

  function siteDensities(p) {
    var a = p.lattice_constants_A.a, c = p.lattice_constants_A.c;
    var cellA2 = c * a * Math.sqrt(2.0);
    return [1e16 / cellA2, 4e16 / cellA2];
  }

  function densityGCm3(p) {
    var a = p.lattice_constants_A.a, c = p.lattice_constants_A.c;
    var vcell = a * a * c * 1e-24;
    return 2.0 * p.molar_mass_g_mol / (N_AVO * vcell);
  }

  function geometry(p, opts) {
    opts = opts || {};
    var ssLayers = opts.ss_layers != null ? opts.ss_layers
      : p.defaults.subsurface_layers;
    var area;
    if (opts.bet_m2_g != null) {
      area = opts.bet_m2_g * 1e4;
    } else {
      var d = opts.d_um != null ? opts.d_um : p.defaults.particle_diameter_um;
      area = 6.0 / (densityGCm3(p) * d * 1e-4);
    }
    var sig = siteDensities(p);
    var nO = 2.0 / p.molar_mass_g_mol * 1e6;
    var nS = area * sig[0] / N_AVO * 1e6;
    var nSS = area * sig[1] * ssLayers / N_AVO * 1e6;
    return { area_m2_g: area * 1e-4, N_s: nS, N_ss: nSS,
             N_b: nO - nS - nSS, N_O_total: nO, ss_layers: ssLayers };
  }

  /* --------------------------------------------------------- lattice */

  function buildLattice(rows, rowSites, ssLayers, bulkLayers, ordering) {
    var along = ordering.pair_eV.along_row;
    var cross = ordering.pair_eV.cross_row;
    if (along.length && rowSites < 2 * along.length + 1) {
      throw new Error('row too short for the along-row interaction range');
    }
    var hasCross = false;
    for (var q = 0; q < cross.length; q++) if (cross[q] !== 0.0) hasCross = true;
    if (hasCross && rows < 3) {
      throw new Error('cross-row interactions need at least 3 rows');
    }
    var nSurf = rows * rowSites;
    var nSS = ssLayers * rows * rowSites;
    var nBulk = bulkLayers * rows * rowSites;
    var n = nSurf + nSS + nBulk;
    var cls = new Uint8Array(n);
    var i;
    for (i = 0; i < n; i++) {
      cls[i] = i < nSurf ? 0 : (i < nSurf + nSS ? 1 : 2);
    }
    var nbrs = new Array(nSurf);
    function sid(r, k) {
      return (((r % rows) + rows) % rows) * rowSites
        + (((k % rowSites) + rowSites) % rowSites);
    }
    for (var r = 0; r < rows; r++) {
      for (var k = 0; k < rowSites; k++) {
        var lst = [];
        for (var dz = 1; dz <= along.length; dz++) {
          var v = along[dz - 1];
          if (v !== 0.0) {
            lst.push([sid(r, k + dz), v]);
            lst.push([sid(r, k - dz), v]);
          }
        }
        for (var off = 0; off < cross.length; off++) {
          var w = cross[off];
          if (w === 0.0) continue;
          if (off === 0) {
            lst.push([sid(r + 1, k), w]);
            lst.push([sid(r - 1, k), w]);
          } else {
            lst.push([sid(r + 1, k + off), w]);
            lst.push([sid(r + 1, k - off), w]);
            lst.push([sid(r - 1, k + off), w]);
            lst.push([sid(r - 1, k - off), w]);
          }
        }
        nbrs[sid(r, k)] = lst;
      }
    }
    /* flat neighbour arrays: same iteration order, no per-pair objects in
       the sampling loop (and the same layout the Python reference uses) */
    var total = 0;
    for (i = 0; i < nSurf; i++) total += nbrs[i].length;
    var nbrOff = new Int32Array(nSurf + 1);
    var nbrIdx = new Int32Array(total);
    var nbrV = new Float64Array(total);
    var pos = 0;
    for (i = 0; i < nSurf; i++) {
      var lst2 = nbrs[i];
      for (var t2 = 0; t2 < lst2.length; t2++) {
        nbrIdx[pos] = lst2[t2][0];
        nbrV[pos] = lst2[t2][1];
        pos++;
      }
      nbrOff[i + 1] = pos;
    }
    return { n: n, n_surf: nSurf, n_ss: nSS, n_bulk: nBulk,
             rows: rows, row_sites: rowSites, cls: cls,
             nbr_off: nbrOff, nbr_idx: nbrIdx, nbr_v: nbrV };
  }

  function inter(i, occ, lat) {
    if (i >= lat.n_surf) return 0.0;
    var s = 0.0;
    var off = lat.nbr_off, idx = lat.nbr_idx, vv = lat.nbr_v;
    for (var t = off[i]; t < off[i + 1]; t++) {
      if (occ[idx[t]]) s += vv[t];
    }
    return s;
  }

  function totalEnergy(lat, occ, epsCls) {
    var e = 0.0, i, t;
    for (i = 0; i < lat.n; i++) {
      if (occ[i]) e += epsCls[lat.cls[i]];
    }
    var p = 0.0;
    var off = lat.nbr_off, idx = lat.nbr_idx, vv = lat.nbr_v;
    for (i = 0; i < lat.n_surf; i++) {
      if (occ[i]) {
        for (t = off[i]; t < off[i + 1]; t++) {
          if (occ[idx[t]]) p += vv[t];
        }
      }
    }
    return e + 0.5 * p;
  }

  /* -------------------------------------------------------------- MC */

  function gapHist(lat, occ, hist) {
    var rows = lat.rows, ls = lat.row_sites;
    for (var r = 0; r < rows; r++) {
      var base = r * ls;
      var ks = [];
      for (var k = 0; k < ls; k++) if (occ[base + k]) ks.push(k);
      var m = ks.length;
      if (m < 2) continue;
      for (var t = 0; t < m; t++) {
        var gap = (t + 1 < m) ? ks[t + 1] - ks[t] : ks[0] + ls - ks[t];
        hist[gap] += 1;
      }
    }
  }

  function adjustFilling(lat, occ, vac, target, rng) {
    var n = lat.n;
    while (vac.length > target) {
      var iv = Math.floor(rng() * vac.length);
      occ[vac[iv]] = 0;
      vac.splice(iv, 1);
    }
    while (vac.length < target) {
      var s = Math.floor(rng() * n);
      if (!occ[s]) {
        occ[s] = 1;
        vac.push(s);
      }
    }
  }

  function mcRun(lat, epsCls, kt, occ, vac, sweepsEq, sweepsSample, rng, opts) {
    opts = opts || {};
    var widomEvery = opts.widom_every != null ? opts.widom_every : 10;
    var histEvery = opts.hist_every != null ? opts.hist_every : 25;
    var blocks = opts.blocks != null ? opts.blocks : 20;
    var n = lat.n;
    var cls = lat.cls;
    var nSurf = lat.n_surf;
    var off = lat.nbr_off, idx = lat.nbr_idx, vv = lat.nbr_v;
    var nv = vac.length;
    var beta = 1.0 / kt;
    var cnt = new Int32Array(3);
    var i;
    for (i = 0; i < nv; i++) cnt[cls[vac[i]]] += 1;
    /* single-element typed arrays keep the mutated accumulators unboxed */
    var es = new Float64Array(1);
    es[0] = totalEnergy(lat, occ, epsCls);
    var eps0 = new Float64Array(3);
    eps0[0] = epsCls[0]; eps0[1] = epsCls[1]; eps0[2] = epsCls[2];

    function attempt() {
      var iv = (rng() * nv) | 0;
      var v = vac[iv];
      var s = (rng() * n) | 0;
      if (occ[s]) return;
      occ[v] = 0;
      var iS = 0.0, iV = 0.0, t;
      if (s < nSurf) {
        for (t = off[s]; t < off[s + 1]; t++) {
          if (occ[idx[t]]) iS += vv[t];
        }
      }
      if (v < nSurf) {
        for (t = off[v]; t < off[v + 1]; t++) {
          if (occ[idx[t]]) iV += vv[t];
        }
      }
      var de = (eps0[cls[s]] - eps0[cls[v]]) + (iS - iV);
      if (de <= 0.0 || rng() < Math.exp(-de * beta)) {
        occ[s] = 1;
        vac[iv] = s;
        cnt[cls[v]] -= 1;
        cnt[cls[s]] += 1;
        es[0] += de;
      } else {
        occ[v] = 1;
      }
    }

    var moves = sweepsEq * n;
    for (i = 0; i < moves; i++) attempt();

    var bsize = Math.max(1, Math.floor(sweepsSample / blocks));
    var thBlocks = [];
    var widSum = 0.0, widN = 0;
    var hist = new Array(lat.row_sites + 1);
    for (i = 0; i <= lat.row_sites; i++) hist[i] = 0;
    var eSum = 0.0;
    for (var b = 0; b < blocks; b++) {
      var acc0 = 0, acc1 = 0, acc2 = 0;
      for (var sw = 0; sw < bsize; sw++) {
        for (i = 0; i < n; i++) attempt();
        acc0 += cnt[0];
        acc1 += cnt[1];
        acc2 += cnt[2];
        eSum += es[0];
        var gsw = b * bsize + sw;
        if (widomEvery && gsw % widomEvery === 0) {
          var w = 0.0;
          for (i = 0; i < n; i++) {
            if (!occ[i]) {
              w += Math.exp(-beta * (epsCls[cls[i]] + inter(i, occ, lat)));
            }
          }
          widSum += w;
          widN += 1;
        }
        if (histEvery && gsw % histEvery === 0) gapHist(lat, occ, hist);
      }
      thBlocks.push([acc0 / (bsize * lat.n_surf),
                     lat.n_ss ? acc1 / (bsize * lat.n_ss) : 0.0,
                     lat.n_bulk ? acc2 / (bsize * lat.n_bulk) : 0.0]);
    }

    var theta = [0.0, 0.0, 0.0];
    var k2;
    for (i = 0; i < thBlocks.length; i++) {
      for (k2 = 0; k2 < 3; k2++) theta[k2] += thBlocks[i][k2];
    }
    for (k2 = 0; k2 < 3; k2++) theta[k2] /= blocks;
    var err = [0.0, 0.0, 0.0];
    if (blocks > 1) {
      for (k2 = 0; k2 < 3; k2++) {
        var s2 = 0.0;
        for (i = 0; i < thBlocks.length; i++) {
          var dd = thBlocks[i][k2] - theta[k2];
          s2 += dd * dd;
        }
        err[k2] = Math.sqrt(s2 / (blocks - 1) / blocks);
      }
    }
    var muW = widN ? -kt * Math.log((widSum / widN) / (nv + 1)) : null;
    var eCheck = totalEnergy(lat, occ, epsCls);
    var vacSorted = vac.slice().sort(function (x, y) { return x - y; });
    return { theta_s: theta[0], theta_ss: theta[1], theta_b: theta[2],
             err_s: err[0], err_ss: err[1], err_b: err[2],
             mu_eV: muW, E_mean_eV: eSum / (blocks * bsize),
             E_drift_eV: Math.abs(es[0] - eCheck),
             hist: hist, n_vac: nv, vac_sorted: vacSorted };
  }

  function isothermScan(p, tC, opts) {
    opts = opts || {};
    var d = {};
    var key;
    for (key in p.defaults.mc) d[key] = p.defaults.mc[key];
    if (opts.mc) for (key in opts.mc) d[key] = opts.mc[key];
    var ordering = opts.ordering != null ? opts.ordering : p.surface_ordering;
    if (opts.zero_pairs) {
      ordering = { pair_eV: { along_row: [], cross_row: [] } };
    }
    var lat = buildLattice(d.rows, d.row_sites, d.ss_layers, d.bulk_layers,
                           ordering);
    var epsCls = opts.eps != null ? opts.eps
      : energetics(p, opts.preset).eps;
    var kt = KB_EV * (tC + 273.15);
    var rng = makeRng(opts.seed != null ? opts.seed : d.seed);
    var occ = new Uint8Array(lat.n);
    var vac = [];
    var quality = opts.quality != null ? opts.quality : 1.0;
    var swEq = Math.max(1, Math.round(d.equil_sweeps * quality));
    var swWarm = Math.max(1, Math.round(d.warm_equil_sweeps * quality));
    var swSm = Math.max(d.blocks, Math.round(d.sample_sweeps * quality));
    var pts = [];
    for (var idx = 0; idx < d.fillings.length; idx++) {
      adjustFilling(lat, occ, vac, d.fillings[idx], rng);
      var res = mcRun(lat, epsCls, kt, occ, vac,
                      idx === 0 ? swEq : swWarm, swSm, rng,
                      { widom_every: d.widom_every, hist_every: d.hist_every,
                        blocks: d.blocks });
      res.filling = d.fillings[idx];
      delete res.vac_sorted;
      pts.push(res);
      if (opts.progress) opts.progress(idx + 1, d.fillings.length);
    }
    return pts;
  }

  /* -------------------------------------------------------- matching */

  function fd(mu, epsI, kt) {
    var x = (epsI - mu) / kt;
    if (x > 700.0) return 0.0;
    if (x < -700.0) return 1.0;
    return 1.0 / (1.0 + Math.exp(x));
  }

  function energetics(p, preset) {
    var ve = p.vacancy_energetics;
    var key = preset || ve.default_preset;
    var q = ve.presets[key];
    return { preset: key, eps: [q.surface_eV, q.subsurface_eV, q.bulk_eV] };
  }

  function thetaSInterp(points, epsS, kt) {
    var pts = [];
    var i;
    for (i = 0; i < points.length; i++) {
      if (points[i].mu_eV != null) pts.push(points[i]);
    }
    pts.sort(function (a, b) { return a.mu_eV - b.mu_eV; });
    var mus = [], ths = [];
    var run = 0.0;
    for (i = 0; i < pts.length; i++) {
      if (pts[i].theta_s > run) run = pts[i].theta_s;
      mus.push(pts[i].mu_eV);
      ths.push(run);
    }
    return function (mu) {
      if (!mus.length) return fd(mu, epsS, kt);
      if (mu <= mus[0]) {
        var t = fd(mu, epsS, kt);
        return t < ths[0] ? t : ths[0];
      }
      if (mu >= mus[mus.length - 1]) return ths[ths.length - 1];
      var lo = 0, hi = mus.length - 1;
      while (hi - lo > 1) {
        var mid = (lo + hi) >> 1;
        if (mus[mid] <= mu) lo = mid; else hi = mid;
      }
      var w = (mu - mus[lo]) / (mus[hi] - mus[lo]);
      return ths[lo] + w * (ths[hi] - ths[lo]);
    };
  }

  function distribute(p, tC, voTotal, geom, opts) {
    opts = opts || {};
    var kt = KB_EV * (tC + 273.15);
    var eps = opts.eps != null ? opts.eps : energetics(p, opts.preset).eps;
    var cap = p.saturation.surface_cap_coverage;
    var thetaSFn = opts.theta_s_fn || function (mu) { return fd(mu, eps[0], kt); };

    function parts(mu) {
      var ts = thetaSFn(mu);
      if (ts > cap) ts = cap;
      return [ts, fd(mu, eps[1], kt), fd(mu, eps[2], kt)];
    }
    function total(mu) {
      var t = parts(mu);
      return geom.N_s * t[0] + geom.N_ss * t[1] + geom.N_b * t[2];
    }

    var lo = -3.0, hi = 3.0;
    var it;
    for (it = 0; it < 40; it++) {
      if (total(lo) <= voTotal) break;
      lo *= 2.0;
    }
    for (it = 0; it < 40; it++) {
      if (total(hi) >= voTotal) break;
      hi *= 2.0;
    }
    for (it = 0; it < 200; it++) {
      var mid = 0.5 * (lo + hi);
      if (total(mid) < voTotal) lo = mid; else hi = mid;
    }
    var mu = 0.5 * (lo + hi);
    var t3 = parts(mu);
    var ts = t3[0], tss = t3[1], tb = t3[2];
    var us = geom.N_s * ts, uss = geom.N_ss * tss, ub = geom.N_b * tb;
    var tot = us + uss + ub;
    var molTi = 1e6 / p.molar_mass_g_mol;
    var x = voTotal / molTi;
    var sat = p.saturation;
    var warnings = [];
    if (ts >= cap - 1e-12) {
      warnings.push({ kind: 'surface_cap',
        text: 'Surface pinned at the ' + (cap * 100).toFixed(0)
          + '% reconstruction threshold (' + us.toFixed(2) + ' of '
          + geom.N_s.toFixed(2) + ' umol-O/g capacity); coverage beyond it '
          + 'is not modelled.' });
    }
    if (tss > sat.subsurface_dilute_warn) {
      warnings.push({ kind: 'subsurface_dense',
        text: 'Subsurface layer occupancy ' + (tss * 100).toFixed(0)
          + '% is far beyond the dilute-defect regime; read this as a '
          + 'heavily reduced near-surface shell (extended defects), not '
          + 'isolated vacancies.' });
    }
    if (x >= sat.x_max_shear) {
      warnings.push({ kind: 'x_max',
        text: 'x = ' + x.toFixed(4) + ' in TiO2-x exceeds the rutile '
          + 'point-defect range (~' + sat.x_max_shear.toFixed(3)
          + '); crystallographic shear planes are expected.' });
    } else if (x >= 0.8 * sat.x_max_shear) {
      warnings.push({ kind: 'x_near',
        text: 'x = ' + x.toFixed(4) + ' in TiO2-x is near the rutile '
          + 'point-defect range (~' + sat.x_max_shear.toFixed(3) + ').' });
    }
    if (voTotal > 0 && Math.abs(tot - voTotal) > 1e-6 * voTotal) {
      warnings.push({ kind: 'unresolved',
        text: 'Inventory not matched within tolerance; check the isotherm '
          + 'range.' });
    }
    return { method: 'canonical_statmech',
             T_C: tC, VO_total_umol_g: voTotal,
             mu_V_eV: mu,
             theta: { surface: ts, subsurface: tss, bulk: tb },
             umol_g: { surface: us, subsurface: uss, bulk: ub },
             fractions: { surface: tot ? us / tot : 0.0,
                          subsurface: tot ? uss / tot : 0.0,
                          bulk: tot ? ub / tot : 0.0 },
             matched_umol_g: tot,
             x_TiO2mx: x, Ti3_frac: 2.0 * x,
             surface_cap_bound: ts >= cap - 1e-12,
             warnings: warnings };
  }

  function spacingSummary(points, thetaS) {
    var best = null;
    var i;
    for (i = 0; i < points.length; i++) {
      var q = points[i];
      if (q.mu_eV == null) continue;
      var d = Math.abs(q.theta_s - thetaS);
      if (best === null || d < best[0]) best = [d, q];
    }
    if (best === null) return null;
    var pt = best[1];
    var tot = 0, mode = 0, modeN = -1, near = 0;
    for (var gap = 1; gap < pt.hist.length; gap++) {
      var c = pt.hist[gap];
      tot += c;
      if (c > modeN) { modeN = c; mode = gap; }
      if (gap <= 2) near += c;
    }
    if (tot === 0) return null;
    return { at_theta_s: pt.theta_s, modal_gap_sites: mode,
             P_gap_le_2: near / tot, pairs_sampled: tot };
  }

  /* ------------------------------------------------------- CO2 layer */

  function co2KineticsParams(p, kin) {
    var base = p.co2_kinetics;
    var d = { Ea_eV: { surface: base.Ea_eV.surface,
                       subsurface: base.Ea_eV.subsurface,
                       bulk: base.Ea_eV.bulk },
              A_eff_per_s_atm: { surface: base.A_eff_per_s_atm.surface,
                                 subsurface: base.A_eff_per_s_atm.subsurface,
                                 bulk: base.A_eff_per_s_atm.bulk },
              p_order_m: base.p_order_m,
              exposure_s: base.defaults.exposure_s,
              p_CO2_atm: base.defaults.p_CO2_atm,
              T_reox_C: base.defaults.T_reox_C };
    if (kin) {
      for (var k in kin) {
        if (k === 'Ea_eV' || k === 'A_eff_per_s_atm') {
          for (var c in kin[k]) d[k][c] = kin[k][c];
        } else {
          d[k] = kin[k];
        }
      }
    }
    return d;
  }

  function accessibility(p, umolByClass, kin) {
    var d = co2KineticsParams(p, kin);
    var kt = KB_EV * (d.T_reox_C + 273.15);
    var pf = Math.pow(d.p_CO2_atm, d.p_order_m);
    var probs = {};
    var acc = 0.0, tot = 0.0;
    var order = ['surface', 'subsurface', 'bulk'];
    for (var i = 0; i < order.length; i++) {
      var c = order[i];
      var k = d.A_eff_per_s_atm[c] * pf * Math.exp(-d.Ea_eV[c] / kt);
      var x = k * d.exposure_s;
      var prob = x > 700.0 ? 1.0 : 1.0 - Math.exp(-x);
      probs[c] = prob;
      acc += umolByClass[c] * prob;
      tot += umolByClass[c];
    }
    return { P: probs, accessible_umol_g: acc,
             f_recoverable: tot ? acc / tot : 0.0, params: d };
  }

  function dielectric(p, voTotal) {
    var c = p.dielectric_correlation;
    if (c.a_per_umol_g == null || c.b == null) return null;
    var e2 = c.a_per_umol_g * voTotal + c.b;
    var out = { eps2: e2 };
    if (c.eps_prime) out.tan_delta = e2 / c.eps_prime;
    return out;
  }

  function sweep(p, tC, geom, voList, opts) {
    opts = opts || {};
    var rows = [];
    for (var i = 0; i < voList.length; i++) {
      var d = distribute(p, tC, voList[i], geom,
                         { theta_s_fn: opts.theta_s_fn, preset: opts.preset,
                           eps: opts.eps });
      var a = accessibility(p, d.umol_g, opts.kin);
      rows.push({ VO_total_umol_g: voList[i],
                  surface: d.umol_g.surface,
                  subsurface: d.umol_g.subsurface,
                  bulk: d.umol_g.bulk,
                  accessible: a.accessible_umol_g,
                  f_recoverable: a.f_recoverable,
                  mu_V_eV: d.mu_V_eV });
    }
    return rows;
  }

  function runFull(p, opts) {
    opts = opts || {};
    var tC = opts.T_C != null ? opts.T_C : p.defaults.T_C;
    var vo = opts.VO_total != null ? opts.VO_total
      : p.defaults.VO_total_umol_g;
    var geom = geometry(p, { d_um: opts.d_um, bet_m2_g: opts.bet_m2_g,
                             ss_layers: opts.ss_layers });
    var kt = KB_EV * (tC + 273.15);
    var eps = opts.eps != null ? opts.eps : energetics(p, opts.preset).eps;
    var pts = isothermScan(p, tC, { seed: opts.seed, quality: opts.quality,
                                    eps: eps, zero_pairs: opts.zero_pairs,
                                    ordering: opts.ordering, mc: opts.mc,
                                    progress: opts.progress });
    var fn = thetaSInterp(pts, eps[0], kt);
    var out = distribute(p, tC, vo, geom, { theta_s_fn: fn, eps: eps });
    out.geometry = geom;
    out.isotherm = pts;
    out.spacing = spacingSummary(pts, out.theta.surface);
    out.accessibility = accessibility(p, out.umol_g, opts.kin);
    return out;
  }

  return { KB_EV: KB_EV, N_AVO: N_AVO,
           makeRng: makeRng, siteDensities: siteDensities,
           densityGCm3: densityGCm3, geometry: geometry,
           buildLattice: buildLattice, inter: inter,
           totalEnergy: totalEnergy, adjustFilling: adjustFilling,
           mcRun: mcRun, isothermScan: isothermScan, fd: fd,
           energetics: energetics, thetaSInterp: thetaSInterp,
           distribute: distribute, spacingSummary: spacingSummary,
           co2KineticsParams: co2KineticsParams,
           accessibility: accessibility, dielectric: dielectric,
           sweep: sweep, runFull: runFull };
});
