/* Steady-state redox of an oxygen carrier, in the browser.

   The whole model is one crossing. Reduction consumes lattice oxygen so
   its rate falls as the solid empties; oxidation consumes vacancies so its
   rate rises. Two opposed monotone curves meet once, and where they meet
   is the steady state. At unit site orders the crossing has a closed form
   and the only material quantity in it is the RATIO of the two rate
   constants:

       theta(face) = K / (1 + K),        K = k(reduction) / k(oxidation)

   with no temperature, no pressure and no prefactor left in it. That is
   why K is the axis this page is drawn on: every other parameterisation -
   fitted rate constants, a vacancy formation energy fed through
   Bronsted-Evans-Polanyi slopes - is a way of predicting K, and each adds
   an assumption that the closed form above does not need.

   Two things then decide whether a carrier can be RUN at that crossing:

       the solid must stay one phase        theta(particle) <= theta_sol
       the face must stay below its ceiling theta(face)     <= theta_max

   and theta(face) = E * theta(particle), where E is the surface
   enrichment a layer-resolved partition supplies and the uniform model
   silently sets to 1.

   This file mirrors the closed forms of solidgas/redox.py and nothing
   else: the Python module's bisection, its descriptor axis and its BEP
   overlay stay there. The forms here are certified against that solver by
   tests/test_redox_port.py, which runs the real bisection - with the
   enrichment entered as the average -> surface mapping - and compares. */

(function (root, factory) {
  'use strict';
  var api = factory();
  if (typeof module !== 'undefined' && module.exports) module.exports = api;
  root.Redox = api;
}(typeof self !== 'undefined' ? self : this, function () {
  'use strict';

  /* theta* = K/(1+K). Written as 1/(1+1/K) for large K so the ratio, not
     the difference, sets the precision of the answer near saturation. */
  function thetaOfRatio(K) {
    if (!(K >= 0.0)) return 0.0;
    if (!isFinite(K)) return 1.0;
    return K <= 1.0 ? K / (1.0 + K) : 1.0 / (1.0 + 1.0 / K);
  }

  /* the inverse, K = theta/(1-theta): what ratio a given occupancy needs */
  function ratioOfTheta(th) {
    if (!(th > 0.0)) return 0.0;
    return th < 1.0 ? th / (1.0 - th) : Infinity;
  }

  function usableWindow(enrichment, thetaSol, thetaMax) {
    var thPhase = Math.min(1.0, enrichment * thetaSol);
    var kPhase = ratioOfTheta(thPhase);
    var kFace = ratioOfTheta(thetaMax);
    return { enrichment: enrichment,
             K_max: Math.min(kPhase, kFace),
             K_second_phase: kPhase,
             K_no_steady_state: kFace,
             binding: kPhase < kFace ? 'second_phase' : 'no_steady_state' };
  }

  /* The enrichment at which the two limits swap places. Below it the solid
     changes phase first; above it the face saturates first, and widening
     the single-phase range stops buying anything. */
  function crossoverEnrichment(thetaSol, thetaMax) {
    return thetaMax / thetaSol;
  }

  function steady(K, opts) {
    opts = opts || {};
    var thSol = opts.theta_sol, thMax = opts.theta_surface_max;
    var E = opts.enrichment != null ? opts.enrichment : 1.0;
    var band = opts.bands || { theta_lo: 0.25, theta_hi: 0.75 };
    var thFace = thetaOfRatio(K);
    var thPart = E > 0.0 ? thFace / E : thFace;
    var w = usableWindow(E, thSol, thMax);
    /* Past the face ceiling the oxidation rate stops rising with the
       inventory - the face is pinned at theta_max by its own surface
       phase - so the imbalance keeps its sign and the particle reduces
       without settling. This is the CEILING speaking, not the solver: the
       mass-action crossing above still has a root, and what the flag says
       is that the face cannot go there to supply it. */
    var noSteadyState = K > w.K_no_steady_state;
    /* the bands are a statement about the solid, so they read the
       particle average - on an enriched particle the face can be well
       into the balanced band while the solid is barely reduced */
    var regime = thPart < band.theta_lo ? 'reduction_limited'
      : thPart > band.theta_hi ? 'oxidation_limited' : 'balanced';
    return { K: K, enrichment: E,
             theta_surface: thFace, theta: thPart,
             /* both half-rates equal theta(face) at the crossing, in units
                of k(oxidation); only the ratio was given, so only a ratio
                of rates can come back */
             rate_over_k_ox: thFace,
             regime: regime,
             theta_solubility: thSol,
             x_equiv: thPart * 2.0,
             beyond_solubility: thPart > thSol,
             no_steady_state: noSteadyState,
             usable: !noSteadyState && thPart <= thSol,
             window: w };
  }

  /* ------------------------------------------------------------ figures */

  /* The two curves the crossing is read off, in units of k(oxidation).
     They are the rate laws themselves at unit site orders - nothing here
     is fitted or interpolated. */
  function crossingData(K, n) {
    n = n || 201;
    var rows = [], i, th;
    for (i = 0; i < n; i++) {
      th = i / (n - 1);
      rows.push({ theta: th, r_red: K * (1.0 - th), r_ox: th });
    }
    return { K: K, theta_star: thetaOfRatio(K),
             rate_over_k_ox: thetaOfRatio(K), rows: rows };
  }

  /* The second-phase boundary over the enrichment axis, stopped where it
     meets the face ceiling: above that enrichment the curve bounds
     nothing and drawing it further would claim a limit that is not there. */
  function phaseData(opts) {
    opts = opts || {};
    var thSol = opts.theta_sol, thMax = opts.theta_surface_max;
    var kLo = opts.K_lo || 1e-4, kHi = opts.K_hi || 1e4;
    var eLo = opts.E_lo || 1.0, eHi = opts.E_hi || 1e4;
    var n = opts.n || 241;
    var crossE = crossoverEnrichment(thSol, thMax);
    var band = [], i, e, t;
    for (i = 0; i < n; i++) {
      e = Math.pow(10.0, Math.log10(eHi) * i / (n - 1));
      if (e > crossE) break;
      t = e * thSol;
      band.push({ E: e, K: t / (1.0 - t) });
    }
    band.push({ E: crossE, K: ratioOfTheta(thMax) });
    return { theta_sol: thSol, theta_surface_max: thMax,
             K_lo: kLo, K_hi: kHi, E_lo: eLo, E_hi: eHi,
             K_no_steady_state: ratioOfTheta(thMax),
             crossover_E: crossE,
             second_phase: band,
             points: opts.points || [],
             current: opts.current || null };
  }

  return { thetaOfRatio: thetaOfRatio, ratioOfTheta: ratioOfTheta,
           usableWindow: usableWindow,
           crossoverEnrichment: crossoverEnrichment,
           steady: steady,
           crossingData: crossingData, phaseData: phaseData };
}));
