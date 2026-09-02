/* The two figures the population model exists to produce.

   The first is the claim: the measured inventory rises eightfold across
   the series while the isolated population passes through a maximum. Both
   are drawn on one axis so the divergence is a shape and not a sentence.
   The second is the fit itself, five measured rates against five fitted
   ones, because a three-parameter fit to five numbers should be shown
   rather than summarised. */

(function (root, factory) {
  'use strict';
  var kit = (typeof module !== 'undefined' && module.exports)
    ? require('./figkit.js') : root.FigKit;
  if (!kit) throw new Error('figures_population.js needs figkit.js');
  var api = factory(kit);
  if (typeof module !== 'undefined' && module.exports) module.exports = api;
  root.PopulationFigures = api;
})(typeof self !== 'undefined' ? self : this, function (K) {
  'use strict';

  var C = K.C, LW = K.LW, T = K.FIGTYPE;
  var lin = K.lin, lg = K.lg, frame = K.frame, axisX = K.axisX, axisY = K.axisY;
  var tint = K.tint;

  /* One shape and one hue per assigned population, so the three series
     stay apart where they cross and still read printed in greyscale. */
  var SERIES = [
    { key: 'iso', shape: 'circle', col: C.surface,
      text: 'Isolated accessible' },
    { key: 'assoc', shape: 'square', col: C.subsurface,
      text: 'Associated surface-region' },
    { key: 'below', shape: 'triangle', col: C.extended,
      text: 'Below active region' }
  ];
  var ISO = C.surface;

  function partitionData(rows) {
    return {
      points: rows.map(function (q) {
        return { sample: q.sample, total: q.nV_total_umol_g,
                 iso: q.n_iso_umol_g, assoc: q.n_assoc_umol_g,
                 below: q.n_below_umol_g };
      })
    };
  }

  function ratesData(rows) {
    return {
      points: rows.map(function (q) {
        return { sample: q.sample, measured: q.r_CO_measured,
                 fitted: q.r_CO_fitted, iso: q.n_iso_umol_g };
      })
    };
  }

  /* ---------------------------------------------------------- drawing */

  /* The claim, on one pair of axes: the measured total runs eightfold
     across the series while the isolated population turns over. Plotting
     the assignment against the total rather than against a sample name
     puts the turnover on the axis - the blue guide rises, peaks near
     100 µmol g⁻¹ and falls, while the grey series climbs the whole way.
     Both axes are logarithmic because the assignments span five decades
     and a stacked bar would hide the smallest of them entirely. */
  function population(D) {
    var f = K.square();
    var p = { x0: f.pane.x0, y0: f.pane.y0, x1: f.pane.x1, y1: f.pane.y1 };
    if (!D || !D.points.length) return f.done();

    var xs = D.points.map(function (q) { return q.total; })
      .filter(function (v) { return isFinite(v) && v > 0; });
    var ys = [];
    D.points.forEach(function (q) {
      SERIES.forEach(function (S) {
        if (isFinite(q[S.key]) && q[S.key] > 0) ys.push(q[S.key]);
      });
    });
    if (!xs.length || !ys.length) return f.done();

    var X = lg(K.decadeFloor(Math.min.apply(null, xs)),
               K.decadeCeil(Math.max.apply(null, xs)), p.x0, p.x1);
    var Y = lg(K.decadeFloor(Math.min.apply(null, ys)),
               K.decadeCeil(Math.max.apply(null, ys)), p.y1, p.y0);

    frame(f, [p.x0, p.x1], [p.y0, p.y1]);

    /* the guide runs through the isolated series only, and it is a guide:
       it interpolates the five computed points and asserts nothing between
       them, which is why it is dashed and carries no markers of its own */
    var guide = D.points.filter(function (q) { return q.iso > 0; })
      .sort(function (a, b) { return a.total - b.total; })
      .map(function (q) { return [X(q.total), Y(q.iso)]; });
    if (guide.length > 1) {
      f.path(K.smooth(guide), ISO, LW.curve, '8,5');
    }

    SERIES.forEach(function (S) {
      D.points.forEach(function (q) {
        var v = q[S.key];
        if (!isFinite(v) || v <= 0) return;
        K.marker(f, S.shape, X(q.total), Y(v), K.MARK, S.col);
      });
    });

    axisX(f, X, p.y1, K.decades(X.d0, X.d1), 'Total Vₒ (µmol g⁻¹)',
          K.powLabel);
    axisY(f, Y, p.x0, K.decades(Y.d0, Y.d1), 'Assigned Vₒ (µmol g⁻¹)',
          K.powLabel);
    K.legend(f, p.x0 + 10, p.y0 + 17, SERIES.map(function (S) {
      return { col: S.col, marker: S.shape, text: S.text };
    }));
    return f.done();
  }

  function rates(D) {
    var f = K.square();
    var p = { x0: f.pane.x0, y0: f.pane.y0, x1: f.pane.x1, y1: f.pane.y1 };
    if (!D || !D.points.length) return f.done();

    /* the rates span 86-fold and the collapse at the last sample is the
       fact the model exists to explain, so the axis has to be logarithmic:
       on a linear one that bar is a line along the floor */
    var vals = [];
    D.points.forEach(function (q) {
      if (q.measured > 0) vals.push(q.measured);
      if (q.fitted > 0) vals.push(q.fitted);
    });
    if (!vals.length) return f.done();
    var Y = lg(K.decadeFloor(Math.min.apply(null, vals)),
               K.decadeCeil(Math.max.apply(null, vals)), p.y1, p.y0);
    var n = D.points.length;
    var X = lin(-0.5, n - 0.5, p.x0, p.x1);
    var w = (p.x1 - p.x0) / n * 0.30;

    frame(f, [p.x0, p.x1], [p.y0, p.y1]);
    D.points.forEach(function (q, i) {
      var cx = X(i);
      var y0 = Y(Y.d0);
      f.rect(cx - w - 1.5, Y(q.measured), w, y0 - Y(q.measured),
             tint(C.ink, 0.18), C.ink, LW.axis);
      f.rect(cx + 1.5, Y(q.fitted), w, y0 - Y(q.fitted), ISO, ISO, LW.axis);
    });
    axisX(f, X, p.y1, D.points.map(function (_, i) { return i; }),
          'reduction treatment', function (i) { return D.points[i].sample; });
    axisY(f, Y, p.x0, K.decades(Y.d0, Y.d1),
          'initial CO rate (µmol g⁻¹ s⁻¹)', K.powLabel);
    K.legend(f, p.x0 + 10, p.y0 + 17, [
      { col: tint(C.ink, 0.18), edge: C.ink, swatch: true, text: 'measured' },
      { col: ISO, swatch: true, text: 'fitted, ∝ n(isolated)' }
    ]);
    return f.done();
  }

  return { population: population, rates: rates,
           partitionData: partitionData, ratesData: ratesData,
           FIGURES: { population: population, rates: rates } };
});
