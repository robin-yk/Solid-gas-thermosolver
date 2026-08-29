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
    var nRes = opts.resolved_layers != null ? opts.resolved_layers
      : (p.defaults.resolved_layers || 1);
    nRes = Math.max(1, Math.round(nRes));
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
    /* One resolved layer IS the declared shell, so the legacy capacity
       is reproduced exactly. Past one, resolution is the point and each
       resolved layer is a single d110 oxygen slab. */
    var slabs = [], i;
    if (nRes === 1) slabs.push(ssLayers);
    else for (i = 0; i < nRes; i++) slabs.push(1);
    var nLayers = [], nSS = 0.0, slabSum = 0;
    for (i = 0; i < slabs.length; i++) {
      nLayers.push(area * sig[1] * slabs[i] / N_AVO * 1e6);
      nSS += nLayers[i];
      slabSum += slabs[i];
    }
    var d110nm = p.lattice_constants_A.a / Math.sqrt(2) / 10;
    return { area_m2_g: area * 1e-4, N_s: nS, N_ss: nSS,
             N_b: nO - nS - nSS, N_O_total: nO, ss_layers: ssLayers,
             layer_nm: d110nm, resolved_layers: nRes,
             slabs_per_layer: slabs, N_layers: nLayers,
             shell_nm: (1 + 4 * slabSum) / 4 * d110nm };
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
    var surfVac = [];
    for (i = 0; i < vacSorted.length; i++) {
      if (vacSorted[i] < lat.n_surf) surfVac.push(vacSorted[i]);
    }
    return { theta_s: theta[0], theta_ss: theta[1], theta_b: theta[2],
             err_s: err[0], err_ss: err[1], err_b: err[2],
             mu_eV: muW, E_mean_eV: eSum / (blocks * bsize),
             E_drift_eV: Math.abs(es[0] - eCheck),
             hist: hist, n_vac: nv, vac_sorted: vacSorted,
             surf_vac: surfVac };
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
      /* the final surface configuration travels with the point: the
         particle view draws the sampled arrangement, not synthetic noise */
      res.rows = d.rows;
      res.row_sites = d.row_sites;
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

  function sigmoid(z) {
    if (z >= 0.0) return z < 700.0 ? 1.0 / (1.0 + Math.exp(-z)) : 1.0;
    var e = z > -700.0 ? Math.exp(z) : 0.0;
    return e / (1.0 + e);
  }

  function thetaOfMu(mu, epsI, omegaI, kt) {
    /* mu = eps + omega*theta + kT ln(theta/(1-theta)), solved in the
       logit where the sigmoid bounds the root between (mu-eps-omega)/kT
       and (mu-eps)/kT with no search for a bracket. At omega = 0 the
       bracket collapses and fd is returned directly, so the legacy
       partition reproduces bit for bit rather than to a tolerance.
       Line-for-line port of solidgas/statmech.py::theta_of_mu. */
    if (omegaI === 0.0) return fd(mu, epsI, kt);
    var hi = (mu - epsI) / kt;
    var lo = hi - omegaI / kt;
    var t;
    if (lo > hi) { t = lo; lo = hi; hi = t; }
    var it, mid, th;
    for (it = 0; it < 200; it++) {
      mid = 0.5 * (lo + hi);
      th = sigmoid(mid);
      if (epsI + omegaI * th + kt * mid < mu) lo = mid; else hi = mid;
    }
    return sigmoid(0.5 * (lo + hi));
  }

  function energetics(p, preset, resolvedLayers, coveragePenalty) {
    var ve = p.vacancy_energetics;
    var key = preset || ve.default_preset;
    var q = ve.presets[key];
    var eps = [q.surface_eV, q.subsurface_eV, q.bulk_eV];
    var d = p.defaults || {};
    var nRes = resolvedLayers != null ? resolvedLayers
      : (d.resolved_layers || 1);
    nRes = Math.max(1, Math.round(nRes));
    var prof = ve.layer_profile || {};
    var span = prof.span_fraction || [0.0];
    if (nRes > span.length) {
      throw new Error('layer_profile resolves at most ' + span.length
                      + ' layers');
    }
    var epsLayers = [], i;
    for (i = 0; i < nRes; i++) {
      epsLayers.push(eps[1] + span[i] * (eps[2] - eps[1]));
    }
    var cp = coveragePenalty != null ? coveragePenalty
      : (d.coverage_penalty || 'off');
    var on = !(cp === null || cp === false || cp === 'off');
    var om = (ve.coverage_penalty || {}).omega_eV || {};
    var omegaS = 0.0, omegaB = 0.0, omegaLayers = [];
    for (i = 0; i < nRes; i++) {
      omegaLayers.push(on ? (om['layer' + (i + 1)] || 0.0) : 0.0);
    }
    if (on) { omegaS = om.surface || 0.0; omegaB = om.bulk || 0.0; }
    var all = [omegaS].concat(omegaLayers, [omegaB]);
    for (i = 0; i < all.length; i++) {
      if (all[i] < 0.0) {
        throw new Error('only repulsive coverage penalties are admitted');
      }
    }
    return { preset: key, eps: eps, resolved_layers: nRes,
             eps_layers: epsLayers, omega_surface: omegaS,
             omega_layers: omegaLayers, omega_bulk: omegaB,
             coverage_penalty: on ? 'on' : 'off' };
  }

  function classesOf(p, geom, opts) {
    /* The class table the deep sum runs over. eps as a bare triple is the
       legacy call and still names surface / one subsurface / bulk. */
    var en = energetics(p, opts.preset, geom.resolved_layers || 1,
                        opts.omega);
    if (opts.eps == null) return en;
    var out = {}, k;
    for (k in en) if (en.hasOwnProperty(k)) out[k] = en[k];
    out.eps = opts.eps.slice();
    if (en.eps_layers.length === 1) {
      out.eps_layers = [opts.eps[1]];
    } else {
      var span = p.vacancy_energetics.layer_profile.span_fraction;
      out.eps_layers = [];
      for (var i = 0; i < en.eps_layers.length; i++) {
        out.eps_layers.push(opts.eps[1] + span[i]
                            * (opts.eps[2] - opts.eps[1]));
      }
    }
    return out;
  }

  function deepOf(en, geom, mu, kt) {
    /* At one resolved layer with no interaction this is the legacy
       expression term for term, including the order of the additions. */
    var caps = geom.N_layers || [geom.N_ss];
    var total = 0.0, i;
    for (i = 0; i < caps.length; i++) {
      total += caps[i] * thetaOfMu(mu, en.eps_layers[i],
                                   en.omega_layers[i], kt);
    }
    return total + geom.N_b * thetaOfMu(mu, en.eps[2], en.omega_bulk, kt);
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
      if (mu >= mus[mus.length - 1]) {
        var n = mus.length;
        if (n < 2 || mus[n - 1] <= mus[n - 2]) return ths[n - 1];
        var tiny = 1e-12;
        var t0 = Math.min(1 - tiny, Math.max(tiny, ths[n - 2]));
        var t1 = Math.min(1 - tiny, Math.max(tiny, ths[n - 1]));
        var l0 = Math.log(t0 / (1 - t0));
        var l1 = Math.log(t1 / (1 - t1));
        var slope = Math.max(0, (l1 - l0) / (mus[n - 1] - mus[n - 2]));
        var z = l1 + slope * (mu - mus[n - 1]);
        return z > 700 ? 1 : 1 / (1 + Math.exp(-z));
      }
      var lo = 0, hi = mus.length - 1;
      while (hi - lo > 1) {
        var mid = (lo + hi) >> 1;
        if (mus[mid] <= mu) lo = mid; else hi = mid;
      }
      var w = (mu - mus[lo]) / (mus[hi] - mus[lo]);
      return ths[lo] + w * (ths[hi] - ths[lo]);
    };
  }

  function phaseBoundaries(p, tC, geom, opts) {
    /* Closed-form mu and inventory boundaries of the phase construction:
       (1x2) onset, reconstruction complete, CS-precipitation ceiling. */
    opts = opts || {};
    var kt = KB_EV * (tC + 273.15);
    var en = classesOf(p, geom, opts);
    var eps = en.eps;
    var sp = p.surface_phases;
    var thT = sp.theta_transition;
    var thR = sp.theta_reconstructed_eff;
    var thSol = p.saturation.x_max_shear / 2.0;
    var muT = eps[0] + en.omega_surface * thT
      + kt * Math.log(thT / (1.0 - thT));
    /* the CS coexistence is a condition on the BULK class, so the bulk
       interaction belongs in it */
    var muCs = eps[2] + en.omega_bulk * thSol
      + kt * Math.log(thSol / (1.0 - thSol));
    var deepT = deepOf(en, geom, muT, kt);
    var deepCs = deepOf(en, geom, muCs, kt);
    return { mu_t_eV: muT, mu_cs_eV: muCs,
             theta_transition: thT, theta_reconstructed_eff: thR,
             theta_sol: thSol,
             VO_onset_umol_g: geom.N_s * thT + deepT,
             VO_recon_complete_umol_g: geom.N_s * thR + deepT,
             VO_cs_umol_g: geom.N_s * thR + deepCs };
  }

  function distribute(p, tC, voTotal, geom, opts) {
    /* Phase-aware matching: (1x1) lattice gas, first-order (1x2)
       coexistence riser, reconstructed line phase, mu pinned at the CS
       coexistence estimate with lever-rule precipitation. Line-for-line
       port of solidgas/statmech.py::distribute. */
    opts = opts || {};
    var kt = KB_EV * (tC + 273.15);
    var en = classesOf(p, geom, opts);
    var eps = en.eps;
    var b = phaseBoundaries(p, tC, geom, { eps: eps, omega: opts.omega });
    var nS = geom.N_s;
    var thT = b.theta_transition;
    var thR = b.theta_reconstructed_eff;

    function deep(mu) {
      return deepOf(en, geom, mu, kt);
    }
    function surf(mu) {
      return thetaOfMu(mu, eps[0], en.omega_surface, kt);
    }

    var uExt = 0.0;
    var regime, phase, mu, ts, us, phi;
    var it, mid;
    if (voTotal <= b.VO_onset_umol_g) {
      regime = 'dilute'; phase = '1x1';
      var lo = -3.0, hi = b.mu_t_eV;
      for (it = 0; it < 60; it++) {
        if (nS * surf(lo) + deep(lo) <= voTotal) break;
        lo = 2.0 * lo - hi;
      }
      for (it = 0; it < 200; it++) {
        mid = 0.5 * (lo + hi);
        if (nS * surf(mid) + deep(mid) < voTotal) lo = mid;
        else hi = mid;
      }
      mu = 0.5 * (lo + hi);
      ts = surf(mu);
      us = nS * ts;
      phi = 0.0;
    } else if (voTotal < b.VO_recon_complete_umol_g) {
      regime = 'surface_coexistence'; phase = 'coexistence';
      mu = b.mu_t_eV;
      us = voTotal - deep(mu);
      ts = us / nS;
      phi = (ts - thT) / (thR - thT);
    } else {
      var muLo = b.mu_t_eV, muHi = b.mu_cs_eV;
      if (voTotal <= b.VO_cs_umol_g) {
        regime = 'reconstructed'; phase = '1x2';
        var target = voTotal - nS * thR;
        for (it = 0; it < 200; it++) {
          mid = 0.5 * (muLo + muHi);
          if (deep(mid) < target) muLo = mid; else muHi = mid;
        }
        mu = 0.5 * (muLo + muHi);
      } else {
        regime = 'cs_coexistence'; phase = '1x2';
        mu = b.mu_cs_eV;
        uExt = voTotal - nS * thR - deep(mu);
      }
      ts = thR;
      us = nS * thR;
      phi = 1.0;
    }
    var caps = geom.N_layers || [geom.N_ss];
    var thLayers = [], uLayers = [], uss = 0.0, li;
    for (li = 0; li < caps.length; li++) {
      thLayers.push(thetaOfMu(mu, en.eps_layers[li], en.omega_layers[li],
                              kt));
      uLayers.push(caps[li] * thLayers[li]);
      uss += uLayers[li];
    }
    /* one resolved layer IS the shell, so its occupancy is reported
       directly rather than divided back out */
    var tss = caps.length === 1 ? thLayers[0]
      : (geom.N_ss ? uss / geom.N_ss : 0.0);
    var tb = thetaOfMu(mu, eps[2], en.omega_bulk, kt);
    var ub = geom.N_b * tb;
    var tot = us + uss + ub + uExt;
    var molTi = 1e6 / p.molar_mass_g_mol;
    var x = voTotal / molTi;
    var xDiss = (us + uss + ub) / molTi;
    var sat = p.saturation;
    var warnings = [];
    if (phase === 'coexistence') {
      warnings.push({ kind: 'surface_phase',
        text: 'Surface is in (1x1)/(1x2) two-phase coexistence at mu_V = '
          + mu.toFixed(3) + ' eV: ' + (phi * 100).toFixed(0)
          + '% of the area is reconstructed to the added-row Ti2O3 phase '
          + '(lever rule between theta = ' + thT.toFixed(2)
          + ' and the line-phase deficiency 0.5). The (1x2) branch is a '
          + 'zero-width line compound; cross-linked variants are not '
          + 'distinguished.' });
    } else if (phase === '1x2') {
      warnings.push({ kind: 'surface_phase',
        text: 'Surface is fully reconstructed to the (1x2) added-row Ti2O3 '
          + 'phase: areal O deficiency 0.5 bridging-ML = ' + us.toFixed(2)
          + ' umol-O/g at this geometry (Onishi-Iwasawa stoichiometry). '
          + 'Surface phases beyond theta_eff = 0.5 are not modelled; the '
          + '(1x2) branch is a zero-width line compound.' });
    }
    if (tss > sat.subsurface_dilute_warn) {
      warnings.push({ kind: 'subsurface_dense',
        text: 'Subsurface layer occupancy ' + (tss * 100).toFixed(0)
          + '%: the dilute-defect approximation fails for this class; read '
          + 'it as a heavily reduced near-surface shell (extended defects), '
          + 'not isolated vacancies.' });
    }
    if (regime === 'cs_coexistence') {
      warnings.push({ kind: 'cs_precipitation',
        text: 'Total inventory exceeds the estimated point-defect '
          + 'solubility (x_sol ~ ' + sat.x_max_shear.toFixed(3)
          + ', approximate): mu_V pins at the rutile / CS-phase coexistence '
          + 'estimate and ' + uExt.toFixed(1) + ' umol-O/g precipitates as '
          + 'extended defects (crystallographic shear planes).' });
    } else if (x >= 0.8 * sat.x_max_shear) {
      warnings.push({ kind: 'x_near',
        text: 'x = ' + x.toFixed(4) + ' in TiO2-x is near the rutile '
          + 'point-defect range (~' + sat.x_max_shear.toFixed(3)
          + '); the CS ceiling sits at ' + b.VO_cs_umol_g.toFixed(1)
          + ' umol-O/g.' });
    }
    if (voTotal > 0 && Math.abs(tot - voTotal) > 1e-6 * voTotal) {
      warnings.push({ kind: 'unresolved',
        text: 'Inventory not matched within tolerance; check the solver '
          + 'bracket.' });
    }
    return { method: 'canonical_statmech',
             T_C: tC, VO_total_umol_g: voTotal,
             mu_V_eV: mu,
             theta: { surface: ts, subsurface: tss, bulk: tb },
             layers: (function () {
               var out = [], i;
               for (i = 0; i < caps.length; i++) {
                 out.push({ index: i + 1, N_umol_g: caps[i],
                            eps_eV: en.eps_layers[i],
                            omega_eV: en.omega_layers[i],
                            theta: thLayers[i], umol_g: uLayers[i] });
               }
               return out;
             }()),
             energetics: { preset: en.preset,
                           resolved_layers: en.resolved_layers,
                           coverage_penalty: en.coverage_penalty,
                           omega_surface: en.omega_surface,
                           omega_bulk: en.omega_bulk },
             surface_phase: { phase: phase, phi_reconstructed: phi },
             regime: regime,
             umol_g: { surface: us, subsurface: uss, bulk: ub },
             extended_defects_umol_g: uExt,
             fractions: { surface: tot ? us / tot : 0.0,
                          subsurface: tot ? uss / tot : 0.0,
                          bulk: tot ? ub / tot : 0.0,
                          extended: tot ? uExt / tot : 0.0 },
             matched_umol_g: tot,
             x_TiO2mx: x, x_dissolved: xDiss, Ti3_frac: 2.0 * x,
             phase_boundaries: b,
             surface_reconstruction_regime: phase !== '1x1',
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

  /* ---------------------------------------------- thermodynamic bridge */

  function phaseValidity(active, reducedCosts, isReduced) {
    var redActive = [];
    var i;
    for (i = 0; i < active.length; i++) {
      if (isReduced[active[i]]) redActive.push(active[i]);
    }
    if (redActive.length) {
      var tier = active.indexOf('TiO2') >= 0 ? 'boundary' : 'invalid';
      return { tier: tier, reduced_active: redActive,
               nearest_phase: redActive[0],
               margin_per_mol_O_kJ: tier === 'boundary' ? 0.0 : null };
    }
    var bestS = null, bestR = null;
    var names = [];
    for (var s in reducedCosts) names.push(s);
    names.sort();
    for (i = 0; i < names.length; i++) {
      s = names[i];
      if (!isReduced[s]) continue;
      var rc = reducedCosts[s];
      var r = (rc && typeof rc === 'object') ? rc.per_mol_O_kJ : rc;
      if (r == null || !isFinite(r)) continue;
      if (bestR === null || r < bestR) { bestR = r; bestS = s; }
    }
    return { tier: 'stable', reduced_active: [],
             nearest_phase: bestS, margin_per_mol_O_kJ: bestR };
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

  function accessibility(p, umolByClass, kin, totalUmol) {
    /* totalUmol, when given, is the denominator of f_recoverable - pass
       the full inventory so precipitated extended defects count against
       the recoverable fraction. */
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
    if (totalUmol != null) tot = totalUmol;
    return { P: probs, accessible_umol_g: acc,
             f_recoverable: tot ? acc / tot : 0.0, params: d };
  }

  function dielectric(p, voTotal) {
    /* power law from the log10-log10 regression of this work's data */
    var c = p.dielectric_correlation;
    if (c.log10_intercept == null || c.exponent == null || !(voTotal > 0)) {
      return null;
    }
    var e2 = Math.pow(10.0, c.log10_intercept
                      + c.exponent * Math.log10(voTotal));
    var rng = c.fit_range_umol_g;
    var out = { eps2: e2,
                extrapolated: !!(rng && (voTotal < rng[0]
                                         || voTotal > rng[1])) };
    if (c.eps_prime) out.tan_delta = e2 / c.eps_prime;
    return out;
  }

  function sweep(p, tC, geom, voList, opts) {
    opts = opts || {};
    var rows = [];
    for (var i = 0; i < voList.length; i++) {
      var d = distribute(p, tC, voList[i], geom,
                         { preset: opts.preset, eps: opts.eps });
      var a = accessibility(p, d.umol_g, opts.kin, voList[i]);
      rows.push({ VO_total_umol_g: voList[i],
                  surface: d.umol_g.surface,
                  subsurface: d.umol_g.subsurface,
                  bulk: d.umol_g.bulk,
                  extended: d.extended_defects_umol_g,
                  accessible: a.accessible_umol_g,
                  f_recoverable: a.f_recoverable,
                  mu_V_eV: d.mu_V_eV,
                  regime: d.regime });
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
    var eps = opts.eps != null ? opts.eps : energetics(p, opts.preset).eps;
    var pts = isothermScan(p, tC, { seed: opts.seed, quality: opts.quality,
                                    eps: eps, zero_pairs: opts.zero_pairs,
                                    ordering: opts.ordering, mc: opts.mc,
                                    progress: opts.progress });
    /* decoupled by design: the layer partition is the analytic
       phase-aware solve on the formation energies alone; the sampled
       isotherm feeds only the ordering exhibit, taken at the coverage of
       the surviving (1x1) patches once coexistence is reached */
    var out = distribute(p, tC, vo, geom, { eps: eps });
    out.geometry = geom;
    out.isotherm = pts;
    out.spacing = spacingSummary(pts,
      Math.min(out.theta.surface, p.surface_phases.theta_transition));
    out.accessibility = accessibility(p, out.umol_g, opts.kin, vo);
    return out;
  }

  return { KB_EV: KB_EV, N_AVO: N_AVO,
           makeRng: makeRng, siteDensities: siteDensities,
           thetaOfMu: thetaOfMu,
           densityGCm3: densityGCm3, geometry: geometry,
           buildLattice: buildLattice, inter: inter,
           totalEnergy: totalEnergy, adjustFilling: adjustFilling,
           mcRun: mcRun, isothermScan: isothermScan, fd: fd,
           energetics: energetics, thetaSInterp: thetaSInterp,
           phaseBoundaries: phaseBoundaries,
           distribute: distribute, spacingSummary: spacingSummary,
           co2KineticsParams: co2KineticsParams,
           accessibility: accessibility, dielectric: dielectric,
           sweep: sweep, phaseValidity: phaseValidity, runFull: runFull };
});
