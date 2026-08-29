/* The two figures the steady-state redox model exists to produce.

   crossing   the two half-reaction rates against the vacancy fraction,
              meeting once. The whole model is this picture.
   phase      where a carrier can actually be run, over the ratio of rate
              constants and the surface enrichment. The second axis is the
              one the uniform assumption sets to 1 without saying so.

   Both are drawn from rows computed by solidgas/redox.py and written out
   by scripts/export_redox_figures.py, so nothing here evaluates a rate
   law or a boundary condition - a figure cannot disagree with the engine
   it never calls. Style, colour roles and geometry are the kit's. */

(function (root, factory) {
  'use strict';
  var kit = (typeof module !== 'undefined' && module.exports)
    ? require('./figkit.js') : root.FigKit;
  if (!kit) throw new Error('figures_redox.js needs figkit.js');
  var api = factory(kit);
  if (typeof module !== 'undefined' && module.exports) module.exports = api;
  root.RedoxFigures = api;
})(typeof self !== 'undefined' ? self : this, function (K) {
  'use strict';

  var C = K.C, LW = K.LW, TICK = K.TICK, T = K.FIGTYPE;
  var lin = K.lin, lg = K.lg, frame = K.frame;
  var axisX = K.axisX, axisY = K.axisY, series = K.series, tint = K.tint;

  /* reduction consumes lattice oxygen and oxidation consumes vacancies, so
     the gas role carries the reducing half and the surface role the
     oxidising one - the same two things they carry everywhere else */
  var RED = C.gas, OX = C.surface;

  function crossing(D) {
    var d = D.crossing;
    var f = K.square();
    var p = f.pane;
    var X = lin(0, 1, p.x0, p.x1);
    var top = 1.0;
    var Y = lin(0, top, p.y1, p.y0);
    frame(f, X, Y);
    axisX(f, X, p.y1, [0, 0.25, 0.5, 0.75, 1],
          'vacancy fraction on the reacting face, θ',
          function (v) { return v.toFixed(2); });
    axisY(f, Y, p.x0, [0, 0.25, 0.5, 0.75, 1], 'rate, in units of k(oxidation)',
          function (v) { return v.toFixed(2); });

    var red = d.rows.map(function (r) { return [r.theta, r.r_red]; });
    var ox = d.rows.map(function (r) { return [r.theta, r.r_ox]; });
    series(f, red, X, Y, RED, LW.curve, null);
    series(f, ox, X, Y, OX, LW.curve, '5 2.5');

    var cx = X(d.theta_star), cy = Y(d.rate_over_k_ox);
    f.line(cx, p.y1, cx, cy, tint(C.ink, 0.35), LW.hair, '2 2');
    f.line(p.x0, cy, cx, cy, tint(C.ink, 0.35), LW.hair, '2 2');
    f.dot(cx, cy, 2.5, C.ink);
    f.text(cx + 6, cy - 5, 'steady state', { size: T.small });
    f.text(cx + 6, cy + 5, 'θ* = ' + d.theta_star.toFixed(3),
           { size: T.small, fill: C.structure });

    /* direct labels: a legend would repeat what the curves already say */
    var iRed = Math.round(0.18 * (d.rows.length - 1));
    var iOx = Math.round(0.80 * (d.rows.length - 1));
    f.text(X(d.rows[iRed].theta) + 4, Y(d.rows[iRed].r_red) - 5,
           'reduction ∝ K(1 − θ)', { size: T.small, fill: RED });
    f.text(X(d.rows[iOx].theta) - 4, Y(d.rows[iOx].r_ox) + 11,
           'oxidation ∝ θ', { size: T.small, fill: OX, anchor: 'end' });
    f.text(p.x0 + 6, p.y0 + 11, 'K = k(red)/k(ox) = ' + d.K,
           { size: T.small });
    return f.done();
  }

  function phase(D) {
    var d = D.phase;
    var f = K.square();
    var p = f.pane;
    var X = lg(d.K_lo, d.K_hi, p.x0, p.x1);
    var Y = lg(d.E_lo, d.E_hi, p.y1, p.y0);

    /* Two regions are shaded because they are the ones a carrier cannot
       be run in; what is left is usable and needs no wash of colour.

       second phase   right of the boundary curve and left of the face
                      ceiling, and only below the enrichment where those
                      two meet - above it the curve bounds nothing.
       no steady      right of the face ceiling, at every enrichment. */
    var i, q;
    var poly = '';
    for (i = 0; i < d.second_phase.length; i++) {
      q = d.second_phase[i];
      poly += (i ? ' L' : 'M') + X(q.K) + ',' + Y(q.E);
    }
    poly += ' L' + X(d.K_no_steady_state) + ',' + Y(d.crossover_E)
      + ' L' + X(d.K_no_steady_state) + ',' + Y(d.E_lo) + ' Z';
    f.path(poly, 'none', 0, null, tint(C.bulk, 0.13));
    f.rect(X(d.K_no_steady_state), p.y0,
           p.x1 - X(d.K_no_steady_state), p.y1 - p.y0,
           tint(C.subsurface, 0.13));

    frame(f, X, Y);
    axisX(f, X, p.y1, K.decades(d.K_lo, d.K_hi), 'K = k(reduction) / k(oxidation)',
          K.expLabel);
    axisY(f, Y, p.x0, K.decades(d.E_lo, d.E_hi),
          'surface enrichment, θ(face) / θ(particle)', K.expLabel);

    var curve = [];
    for (i = 0; i < d.second_phase.length; i++) {
      if (d.second_phase[i].E > d.E_hi) break;
      curve.push([d.second_phase[i].K, d.second_phase[i].E]);
    }
    series(f, curve, X, Y, C.bulk, LW.curve, '5 2.5');
    f.line(X(d.K_no_steady_state), p.y0,
           X(d.K_no_steady_state), p.y1, C.subsurface, LW.curve);

    f.text(X(1.2e-3), Y(2200), 'usable', { size: T.small, weight: 'bold' });
    f.text(X(0.14), Y(9), 'second phase',
           { size: T.small, anchor: 'middle', fill: C.bulk });
    f.text(X(0.14), Y(5.2), 'vacancies order away',
           { size: T.small, anchor: 'middle', fill: C.structure });
    f.text(X(12), Y(60), 'no steady state',
           { size: T.small, fill: C.subsurface });

    for (i = 0; i < d.points.length; i++) {
      q = d.points[i];
      var px = X(Math.max(q.K, d.K_lo)), py = Y(q.E);
      var below = q.E < 3.0;              /* keep clear of the bottom spine */
      f.dot(px, py, 2.5, C.ink);
      f.text(px + 6, py + (below ? -6 : 3), q.label + ', K ≤ '
             + (q.K < 0.01 ? q.K.toFixed(4) : q.K.toFixed(2)),
             { size: T.small });
    }
    return f.done();
  }

  return { crossing: crossing, phase: phase,
           FIGURES: { crossing: crossing, phase: phase } };
});
