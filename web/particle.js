/* Vacancy distribution in one oxide particle, from a titrated inventory.

   Line-for-line port of solidgas/particle/: same branches, same bracket
   algebra, same layer stack summed in the same order, so the browser and
   the Python engine agree to the last bits (verified by
   tests/test_particle_port.py against the same fifty-digit reference).

   The physics is in the Python module headers and in
   docs/particle-model.md; the short version is that an oxygen vacancy
   leaves two electrons behind, they localise as two Ti(3+) on a cation
   sublattice half the size of the oxygen one, and counting that
   sublattice is what caps reduction at Ti2O3 and turns the dilute
   isotherm from p(O2)^(-1/2) into p(O2)^(-1/6). */

(function (root, factory) {
  'use strict';
  var api = factory();
  if (typeof module !== 'undefined' && module.exports) module.exports = api;
  root.Particle = api;
})(typeof self !== 'undefined' ? self : this, function () {
  'use strict';

  var KB_EV = 8.617333262e-5;
  var N_AVO = 6.02214076e23;

  var POLARON_RATIO = 4.0;
  var THETA_CAP = 1.0 / POLARON_RATIO;
  var ITERS = 80;

  /* --------------------------------------------------------- free energy */

  function cap(ratio) { return 1.0 / (ratio === undefined ? POLARON_RATIO : ratio); }

  function phi(theta, ratio) {
    ratio = ratio === undefined ? POLARON_RATIO : ratio;
    if (theta <= 0.0) return -Infinity;
    if (theta >= 1.0 / ratio) return Infinity;
    var p = ratio * theta;
    return Math.log(theta / (1.0 - theta)) + 2.0 * Math.log(p / (1.0 - p));
  }

  function phiOxygenOnly(theta) {
    if (theta <= 0.0) return -Infinity;
    if (theta >= 1.0) return Infinity;
    return Math.log(theta / (1.0 - theta));
  }

  function brackets(ratio) {
    var half = 0.5 / ratio;
    return { phiHalf: Math.log(half / (1.0 - half)),
             logRatio2: 2.0 * Math.log(ratio),
             slack: -Math.log(1.0 - half) - 2.0 * Math.log(0.5) };
  }

  function thetaOfPhi(x, compensated, ratio) {
    compensated = compensated === undefined ? true : compensated;
    ratio = ratio === undefined ? POLARON_RATIO : ratio;
    if (!compensated) {
      if (x > 700.0) return 1.0;
      if (x < -700.0) return 0.0;
      return 1.0 / (1.0 + Math.exp(-x));
    }
    var b = brackets(ratio), top = 1.0 / ratio, lo, hi, mid, i, t;

    if (x <= b.phiHalf) {
      hi = (x - b.logRatio2) / 3.0;
      if (hi < -745.0) return 0.0;
      lo = Math.max((x - b.logRatio2 - b.slack) / 3.0, -745.0);
      for (i = 0; i < ITERS; i++) {
        mid = 0.5 * (lo + hi);
        if (phi(Math.exp(mid), ratio) < x) lo = mid; else hi = mid;
      }
      return Math.exp(0.5 * (lo + hi));
    }

    hi = Math.log(0.5 * top);
    lo = hi - 1.0;
    while (lo > -745.0 && phi(top - Math.exp(lo), ratio) < x) {
      lo = hi - 2.0 * (hi - lo);
    }
    if (lo <= -745.0) {
      if (phi(top - Math.exp(-745.0), ratio) < x) return top;
      lo = -745.0;
    }
    for (i = 0; i < ITERS; i++) {
      mid = 0.5 * (lo + hi);
      t = top - Math.exp(mid);
      if (t <= 0.0 || phi(t, ratio) > x) lo = mid; else hi = mid;
    }
    return top - Math.exp(0.5 * (lo + hi));
  }

  function muOfTheta(theta, eForm, kt, compensated, ratio) {
    var f = (compensated === undefined || compensated)
      ? phi(theta, ratio) : phiOxygenOnly(theta);
    return eForm + kt * f;
  }

  function thetaOfMu(mu, eForm, kt, compensated, ratio) {
    return thetaOfPhi((mu - eForm) / kt, compensated, ratio);
  }

  function diluteTheta(mu, eForm, kt, ratio) {
    ratio = ratio === undefined ? POLARON_RATIO : ratio;
    return Math.exp(((mu - eForm) / kt - 2.0 * Math.log(ratio)) / 3.0);
  }

  function freeEnergyPerSite(theta, eForm, kt, compensated, ratio) {
    ratio = ratio === undefined ? POLARON_RATIO : ratio;
    if (theta <= 0.0) return 0.0;
    var ideal = theta * Math.log(theta) + (1.0 - theta) * Math.log(1.0 - theta);
    var f = eForm * theta + kt * ideal;
    if (compensated === false) return f;
    var p = ratio * theta;
    if (p >= 1.0) return Infinity;
    var catTerm = p * Math.log(p) + (1.0 - p) * Math.log(1.0 - p);
    return f + kt * catTerm * (2.0 / ratio);
  }

  /* ------------------------------------------------------------ material */

  function segregation(p, preset) {
    var ve = p.vacancy_energetics;
    var key = preset || ve.default_preset;
    var q = ve.presets[key];
    var eS = q.surface_eV, eSs = q.subsurface_eV, eB = q.bulk_eV;
    var d110 = p.lattice_constants_A.a / Math.SQRT2 / 10.0;
    var de = eB - eS, rem = eB - eSs;
    if (!(de > 0.0 && rem > 0.0 && rem < de)) {
      throw new Error('the ' + key + ' triple does not decay monotonically '
                      + 'towards the bulk, so E(z) = E_bulk - dE exp(-z/xi) '
                      + 'cannot pass through all three');
    }
    var xi = -d110 / Math.log(rem / de);
    return { preset: key, label: q.label || key,
             e_surface_eV: eS, e_subsurface_eV: eSs, e_bulk_eV: eB,
             dE_seg_eV: de, xi_nm: xi, d110_nm: d110,
             xi_in_layers: xi / d110 };
  }

  function formationEnergy(zNm, seg) {
    return seg.e_bulk_eV - seg.dE_seg_eV * Math.exp(-Math.max(0.0, zNm) / seg.xi_nm);
  }

  function densityGCm3(p) {
    var a = p.lattice_constants_A.a, c = p.lattice_constants_A.c;
    return 2.0 * p.molar_mass_g_mol / (N_AVO * a * a * c * 1e-24);
  }

  function bridgingDensityCm2(p) {
    var a = p.lattice_constants_A.a, c = p.lattice_constants_A.c;
    return 1e16 / (c * a * Math.SQRT2);
  }

  function particle(p, opts) {
    opts = opts || {};
    var areaCm2G;
    if (opts.bet_m2_g != null) {
      areaCm2G = opts.bet_m2_g * 1e4;
    } else {
      var dUm = opts.d_um != null ? opts.d_um : p.defaults.particle_diameter_um;
      areaCm2G = 6.0 / (densityGCm3(p) * dUm * 1e-4);
    }
    var nO = 2.0 / p.molar_mass_g_mol * 1e6;
    var nSurface = areaCm2G * bridgingDensityCm2(p) / N_AVO * 1e6;
    var rNm = 3.0 / (densityGCm3(p) * areaCm2G) * 1e7;
    var xShear = (p.saturation || {}).x_max_shear;
    var ph = p.surface_phases || {};
    return { area_m2_g: areaCm2G * 1e-4, radius_nm: rNm,
             diameter_um: 2.0 * rNm * 1e-3,
             N_surface_umol_g: nSurface,
             N_interior_umol_g: nO - nSurface,
             N_oxygen_umol_g: nO,
             segregation: segregation(p, opts.preset),
             x_max_shear: xShear == null ? null : xShear,
             theta_solubility: xShear == null ? null : 0.5 * xShear,
             surface_phases: {
               theta_transition: ph.theta_transition,
               theta_reconstructed_eff: ph.theta_reconstructed_eff,
               transition_interval: ph.transition_interval } };
  }

  function thetaFromLoading(vo, part) { return vo / part.N_oxygen_umol_g; }
  function loadingFromTheta(th, part) { return th * part.N_oxygen_umol_g; }
  function ktEV(tC) { return KB_EV * (tC + 273.15); }

  /* ------------------------------------------------------------- profile */

  function layerGeometry(part) {
    var d = part.segregation.d110_nm;
    var rIn = part.radius_nm - 0.5 * d;
    return { d: d, rIn: rIn, vol: rIn * rIn * rIn / 3.0 };
  }

  function interiorAverage(mu, kt, part, compensated, ratio, relTol) {
    relTol = relTol === undefined ? 1e-15 : relTol;
    var seg = part.segregation, g = layerGeometry(part);
    var thBulk = thetaOfMu(mu, seg.e_bulk_eV, kt, compensated, ratio);
    var excess = 0.0, k = 1, nMax = Math.floor(g.rIn / g.d), ro, ri, v, e, term;
    while (k <= nMax) {
      ro = g.rIn - (k - 1) * g.d;
      ri = Math.max(0.0, g.rIn - k * g.d);
      v = (ro * ro * ro - ri * ri * ri) / 3.0;
      e = formationEnergy(k * g.d, seg);
      term = v * (thetaOfMu(mu, e, kt, compensated, ratio) - thBulk);
      excess += term;
      if (Math.abs(term) <= relTol * Math.max(Math.abs(excess), 1e-300) && k > 4) break;
      k++;
    }
    return thBulk + excess / g.vol;
  }

  function surfaceState(mu, kt, part, compensated, ratio, reconstruction) {
    var seg = part.segregation, eS = seg.e_surface_eV;
    var th11 = thetaOfMu(mu, eS, kt, compensated, ratio);
    var ph = part.surface_phases || {};
    if (reconstruction === false || ph.theta_transition == null
        || ph.theta_reconstructed_eff == null) {
      return { theta: th11, phase: '1x1', reconstructed_fraction: 0.0,
               mu_transition_eV: null };
    }
    var muT = muOfTheta(ph.theta_transition, eS, kt, compensated, ratio);
    if (mu < muT) {
      return { theta: th11, phase: '1x1', reconstructed_fraction: 0.0,
               mu_transition_eV: muT };
    }
    if (mu > muT) {
      return { theta: ph.theta_reconstructed_eff, phase: '1x2',
               reconstructed_fraction: 1.0, mu_transition_eV: muT };
    }
    return { theta: ph.theta_transition, phase: 'coexistence',
             reconstructed_fraction: 0.0, mu_transition_eV: muT };
  }

  function inventory(mu, kt, part, compensated, ratio, reconstruction, surfFraction) {
    var s = surfaceState(mu, kt, part, compensated, ratio, reconstruction);
    var thS = surfFraction == null ? s.theta : surfFraction;
    return part.N_surface_umol_g * thS
      + part.N_interior_umol_g * interiorAverage(mu, kt, part, compensated, ratio);
  }

  function bracketMu(target, kt, part, compensated, ratio, reconstruction) {
    var seg = part.segregation;
    var th = Math.min(Math.max(target / part.N_oxygen_umol_g, 1e-300), 0.5 / ratio);
    var guess = seg.e_bulk_eV + kt * (3.0 * Math.log(th) + 2.0 * Math.log(ratio));
    var lo = guess - 1.0, hi = guess + 1.0, i;
    for (i = 0; i < 60; i++) {
      if (inventory(lo, kt, part, compensated, ratio, reconstruction) <= target) break;
      lo -= (hi - lo);
    }
    for (i = 0; i < 60; i++) {
      if (inventory(hi, kt, part, compensated, ratio, reconstruction) >= target) break;
      hi += (hi - lo);
    }
    return [lo, hi];
  }

  function equilibrate(voUmolG, tC, part, opts) {
    opts = opts || {};
    var compensated = opts.compensated === undefined ? true : opts.compensated;
    var ratio = opts.ratio === undefined ? POLARON_RATIO : opts.ratio;
    var reconstruction = opts.reconstruction === undefined ? true : opts.reconstruction;
    var shearLimit = opts.shear_limit === undefined ? true : opts.shear_limit;

    var kt = ktEV(tC), seg = part.segregation;
    var capTotal = part.N_oxygen_umol_g / ratio;
    if (voUmolG > capTotal) {
      throw new Error('a loading of ' + voUmolG + ' umol-O/g is past the cation '
        + 'sublattice itself (' + capTotal.toFixed(1) + ' umol-O/g is theta = 1/'
        + ratio + ' everywhere, which is Ti2O3, not reduced rutile)');
    }

    var muCs = null;
    if (shearLimit && part.theta_solubility) {
      muCs = muOfTheta(part.theta_solubility, seg.e_bulk_eV, kt, compensated, ratio);
    }
    var ph = part.surface_phases || {};
    var thT = reconstruction ? ph.theta_transition : null;
    var thR = reconstruction ? ph.theta_reconstructed_eff : null;
    var muT = thT == null ? null
      : muOfTheta(thT, seg.e_surface_eV, kt, compensated, ratio);

    var mu = null, limitedBy = null, extended = 0.0, thS = null, frac = 0.0;

    if (muT != null && (muCs == null || muT <= muCs)) {
      var loInv = inventory(muT, kt, part, compensated, ratio, reconstruction, thT);
      var hiInv = inventory(muT, kt, part, compensated, ratio, reconstruction, thR);
      if (loInv <= voUmolG && voUmolG <= hiInv) {
        mu = muT;
        frac = hiInv > loInv ? (voUmolG - loInv) / (hiInv - loInv) : 0.0;
        thS = thT + frac * (thR - thT);
        limitedBy = 'surface reconstruction';
      }
    }

    if (mu === null) {
      var b = bracketMu(voUmolG, kt, part, compensated, ratio, reconstruction);
      var lo = b[0], hi = b[1], mid;
      for (var i = 0; i < 200; i++) {
        mid = 0.5 * (lo + hi);
        if (inventory(mid, kt, part, compensated, ratio, reconstruction) < voUmolG) {
          lo = mid;
        } else {
          hi = mid;
        }
      }
      mu = 0.5 * (lo + hi);
    }

    if (muCs != null && mu > muCs) {
      mu = muCs; thS = null; limitedBy = 'crystallographic shear';
      extended = voUmolG - inventory(mu, kt, part, compensated, ratio, reconstruction);
    }

    var surf = surfaceState(mu, kt, part, compensated, ratio, reconstruction);
    if (thS != null) {
      surf = { theta: thS, phase: 'coexistence', reconstructed_fraction: frac,
               mu_transition_eV: surf.mu_transition_eV };
    }
    var thI = interiorAverage(mu, kt, part, compensated, ratio);
    var heldS = part.N_surface_umol_g * surf.theta;
    var heldI = part.N_interior_umol_g * thI;

    return {
      mu_eV: mu, kt_eV: kt, T_C: tC,
      compensated: compensated, polaron_ratio: ratio,
      loading_umol_g: voUmolG,
      theta_particle: voUmolG / part.N_oxygen_umol_g,
      surface: surf,
      theta_interior_mean: thI,
      theta_bulk_plateau: thetaOfMu(mu, seg.e_bulk_eV, kt, compensated, ratio),
      theta_first_layer: thetaOfMu(mu, formationEnergy(seg.d110_nm, seg), kt,
                                   compensated, ratio),
      umol_surface: heldS, umol_interior: heldI,
      umol_extended: extended,
      held_umol_g: heldS + heldI,
      limited_by: limitedBy,
      mu_shear_eV: muCs, mu_transition_eV: muT,
      point_defect_ceiling_umol_g: muCs == null ? null
        : inventory(muCs, kt, part, compensated, ratio, reconstruction)
    };
  }

  function depthProfile(res, part, n, zMaxNm) {
    n = n || 160;
    var seg = part.segregation, kt = res.kt_eV;
    var zmax = zMaxNm != null ? zMaxNm
      : Math.min(part.radius_nm, 40.0 * seg.xi_nm);
    var out = [{ z_nm: 0.0, theta: res.surface.theta, pool: 'bridging row' }];
    for (var i = 0; i < n; i++) {
      var z = 0.5 * seg.d110_nm
        + (zmax - 0.5 * seg.d110_nm) * Math.pow(i / (n - 1.0), 3);
      out.push({ z_nm: z,
                 theta: thetaOfMu(res.mu_eV, formationEnergy(z, seg), kt,
                                  res.compensated, res.polaron_ratio),
                 pool: 'interior' });
    }
    return out;
  }

  function shellSplit(res, part, depthsNm) {
    depthsNm = depthsNm || [1.0, 2.0, 5.0, 10.0];
    var seg = part.segregation, kt = res.kt_eV, mu = res.mu_eV;
    var g = layerGeometry(part), total = res.held_umol_g, rows = [];
    for (var j = 0; j < depthsNm.length; j++) {
      var dep = depthsNm[j];
      if (dep >= g.rIn) continue;
      var acc = 0.0, k = 1, nMax = Math.floor(g.rIn / g.d);
      while (k <= nMax) {
        var zOut = (k - 1) * g.d;
        if (zOut >= dep) break;
        var ro = g.rIn - zOut;
        var ri = Math.max(g.rIn - Math.min(k * g.d, dep), 0.0);
        acc += (ro * ro * ro - ri * ri * ri) / 3.0
          * thetaOfMu(mu, formationEnergy(k * g.d, seg), kt,
                      res.compensated, res.polaron_ratio);
        k++;
      }
      var umol = res.umol_surface + part.N_interior_umol_g * acc / g.vol;
      rows.push({ depth_nm: dep,
                  volume_fraction: 1.0 - Math.pow((g.rIn - dep) / g.rIn, 3),
                  umol_g: umol,
                  inventory_fraction: total ? umol / total : 0.0 });
    }
    return rows;
  }

  /* ----------------------------------------------------------- transport */

  function diffusivity(tC, eM, nu0, aNm) {
    return aNm * aNm * nu0 * Math.exp(-eM / (KB_EV * (tC + 273.15)));
  }

  function uptake(tau) {
    if (tau <= 0.0) return 0.0;
    var s = 0.0;
    for (var n = 1; n <= 200; n++) {
      var x = n * n * Math.PI * Math.PI * tau;
      if (x > 700.0) break;
      s += Math.exp(-x) / (n * n);
    }
    return Math.max(0.0, Math.min(1.0, 1.0 - 6.0 / (Math.PI * Math.PI) * s));
  }

  function tauForUptake(f) {
    if (!(f > 0.0 && f < 1.0)) {
      throw new Error('uptake fraction must be strictly between 0 and 1');
    }
    var lo = 1e-12, hi = 10.0, mid;
    for (var i = 0; i < 200; i++) {
      mid = Math.sqrt(lo * hi);
      if (uptake(mid) < f) lo = mid; else hi = mid;
    }
    return Math.sqrt(lo * hi);
  }

  function timescales(tC, part, p, opts) {
    opts = opts || {};
    var kin = (p || {}).kinetics_cgmc || {};
    var eM = opts.e_m_eV != null ? opts.e_m_eV : kin.E_m_bulk_eV;
    var nu0 = opts.nu0_per_s != null ? opts.nu0_per_s : kin.nu0_per_s;
    var aNm = opts.a_nm != null ? opts.a_nm : p.lattice_constants_A.c / 10.0;
    var fractions = opts.fractions || [0.5, 0.9, 0.99];
    var d = diffusivity(tC, eM, nu0, aNm);
    var scale = part.radius_nm * part.radius_nm / d;
    var times = {};
    for (var i = 0; i < fractions.length; i++) {
      times[fractions[i]] = tauForUptake(fractions[i]) * scale;
    }
    return { T_C: tC, E_m_eV: eM, nu0_per_s: nu0, a_nm: aNm,
             radius_nm: part.radius_nm, D_nm2_per_s: d, D_cm2_per_s: d * 1e-14,
             tau_unit_s: scale, times_s: times };
  }

  function barrierForTime(tC, part, tTargetS, nu0, aNm, fraction) {
    fraction = fraction === undefined ? 0.99 : fraction;
    var kt = KB_EV * (tC + 273.15);
    var need = tauForUptake(fraction) * part.radius_nm * part.radius_nm / tTargetS;
    return kt * Math.log(aNm * aNm * nu0 / need);
  }

  function verdict(tC, part, exposureS, p, opts) {
    opts = opts || {};
    var fraction = opts.fraction === undefined ? 0.99 : opts.fraction;
    var ts = timescales(tC, part, p,
                        { e_m_eV: opts.e_m_eV, nu0_per_s: opts.nu0_per_s,
                          a_nm: opts.a_nm, fractions: [fraction] });
    var tHom = ts.times_s[fraction];
    var eBreak = barrierForTime(tC, part, exposureS, ts.nu0_per_s, ts.a_nm, fraction);
    ts.fraction = fraction;
    ts.exposure_s = exposureS;
    ts.t_homogenise_s = tHom;
    ts.ratio = exposureS ? tHom / exposureS : Infinity;
    ts.equilibrium_reachable = tHom < exposureS;
    ts.E_m_that_would_break_it_eV = eBreak;
    ts.barrier_margin_eV = eBreak - ts.E_m_eV;
    return ts;
  }

  /* --------------------------------------------------------- space charge */

  var E_CHARGE = 1.602176634e-19;
  var EPS0_F_PER_CM = 8.8541878128e-14;
  var EPS_R_RUTILE = 100.0;

  function oxygenDensityCm3(p) {
    var a = p.lattice_constants_A.a, c = p.lattice_constants_A.c;
    return 4.0 / (a * a * c * 1e-24);
  }

  function debyeLengthNm(thetaV, p, tC, epsR, ratio) {
    epsR = epsR === undefined ? EPS_R_RUTILE : epsR;
    ratio = ratio === undefined ? POLARON_RATIO : ratio;
    var z2n = (4.0 * thetaV + 0.5 * ratio * thetaV) * oxygenDensityCm3(p);
    if (z2n <= 0.0) return Infinity;
    var ktJ = KB_EV * (tC + 273.15) * E_CHARGE;
    return Math.sqrt(epsR * EPS0_F_PER_CM * ktJ / (E_CHARGE * E_CHARGE * z2n)) * 1e7;
  }

  function assess(res, part, p, epsR) {
    epsR = epsR === undefined ? EPS_R_RUTILE : epsR;
    var xi = part.segregation.xi_nm;
    var lam = debyeLengthNm(res.theta_bulk_plateau, p, res.T_C, epsR,
                            res.polaron_ratio);
    var r = lam / xi, verdictText;
    if (r < 0.2) {
      verdictText = 'screening is much shorter than the segregation depth, '
        + 'so charge follows chemistry and pointwise electroneutrality is safe';
    } else if (r < 2.0) {
      verdictText = 'screening and segregation happen over the same distance, '
        + 'so they cannot be separated; the outer few nanometres of this '
        + 'profile carry an error that this model does not quantify, and the '
        + 'ceria study finds that error runs towards MORE surface enrichment';
    } else {
      verdictText = 'screening is longer than the segregation depth, so a '
        + 'genuine space-charge layer extends past the chemically perturbed '
        + 'shell and pointwise electroneutrality is not a safe assumption here';
    }
    return { eps_r: epsR, theta_bulk_plateau: res.theta_bulk_plateau,
             debye_nm_at_bulk: lam,
             debye_nm_at_first_layer: debyeLengthNm(
               Math.max(res.theta_first_layer, 1e-30), p, res.T_C, epsR,
               res.polaron_ratio),
             xi_nm: xi, ratio_debye_over_xi: r,
             radius_nm: part.radius_nm,
             fraction_of_radius: lam / part.radius_nm,
             verdict: verdictText };
  }

  return { KB_EV: KB_EV, N_AVO: N_AVO,
           POLARON_RATIO: POLARON_RATIO, THETA_CAP: THETA_CAP,
           cap: cap, phi: phi, phiOxygenOnly: phiOxygenOnly,
           thetaOfPhi: thetaOfPhi, muOfTheta: muOfTheta,
           thetaOfMu: thetaOfMu, diluteTheta: diluteTheta,
           freeEnergyPerSite: freeEnergyPerSite,
           segregation: segregation, formationEnergy: formationEnergy,
           densityGCm3: densityGCm3, bridgingDensityCm2: bridgingDensityCm2,
           particle: particle, thetaFromLoading: thetaFromLoading,
           loadingFromTheta: loadingFromTheta, ktEV: ktEV,
           interiorAverage: interiorAverage, surfaceState: surfaceState,
           inventory: inventory, equilibrate: equilibrate,
           depthProfile: depthProfile, shellSplit: shellSplit,
           diffusivity: diffusivity, uptake: uptake,
           tauForUptake: tauForUptake, timescales: timescales,
           barrierForTime: barrierForTime, verdict: verdict,
           debyeLengthNm: debyeLengthNm, assess: assess,
           EPS_R_RUTILE: EPS_R_RUTILE };
});
