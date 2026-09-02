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

  var ISO = C.surface, ASSOC = C.subsurface, BELOW = C.extended;

  function partitionData(rows) {
    return {
      bars: rows.map(function (q) {
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

  function population(D) {
    var f = K.square();
    var p = { x0: f.pane.x0 + 4, y0: 22, x1: f.pane.x1, y1: f.pane.y1 };
    if (!D || !D.bars.length) return f.done();

    var hi = Math.max.apply(null, D.bars.map(function (b) { return b.total; }));
    /* four orders of magnitude between the smallest split and the largest
       total, so the count axis is logarithmic and the bars are stacked on it */
    var Y = lg(0.05, hi * 1.6, p.y1, p.y0);
    var n = D.bars.length;
    var X = lin(-0.5, n - 0.5, p.x0, p.x1);
    var w = (p.x1 - p.x0) / n * 0.24;

    frame(f, [p.x0, p.x1], [p.y0, p.y1]);

    D.bars.forEach(function (b, i) {
      var cx = X(i);
      var trio = [[b.iso, ISO], [b.assoc, ASSOC], [b.below, BELOW]];
      trio.forEach(function (t, j) {
        var v = Math.max(t[0], 0.05);
        var y = Y(v), y0 = Y(0.05);
        var x = cx + (j - 1) * (w + 1.6) - w / 2;
        f.rect(x, y, w, y0 - y, t[0] > 0.05 ? t[1] : tint(t[1], 0.3),
               t[1], LW.axis);
      });
      f.dot(cx, Y(b.total), K.MARK.r, C.paper, C.ink);
    });

    axisX(f, X, p.y1, D.bars.map(function (_, i) { return i; }), 'reduction treatment',
          function (i) { return D.bars[i].sample; });
    axisY(f, Y, p.x0, K.decades(0.05, hi * 1.6), 'µmol-O g⁻¹', K.expLabel);
    K.legend(f, p.x0 + 7, p.y0 + 11, [
      { col: ISO, swatch: true, text: 'isolated, active' },
      { col: ASSOC, swatch: true, text: 'associated' },
      { col: BELOW, swatch: true, text: 'below the active region' },
      { col: C.ink, text: 'measured total' }
    ]);
    return f.done();
  }

  function rates(D) {
    var f = K.square();
    var p = { x0: f.pane.x0 + 4, y0: 22, x1: f.pane.x1, y1: f.pane.y1 };
    if (!D || !D.points.length) return f.done();

    var vals = [];
    D.points.forEach(function (q) { vals.push(q.measured, q.fitted); });
    var hi = Math.max.apply(null, vals);
    var Y = lin(0, hi * 1.15, p.y1, p.y0);
    var n = D.points.length;
    var X = lin(-0.5, n - 0.5, p.x0, p.x1);
    var w = (p.x1 - p.x0) / n * 0.30;

    frame(f, [p.x0, p.x1], [p.y0, p.y1]);
    D.points.forEach(function (q, i) {
      var cx = X(i);
      f.rect(cx - w - 1, Y(q.measured), w, Y(0) - Y(q.measured),
             tint(C.ink, 0.18), C.ink, LW.axis);
      f.rect(cx + 1, Y(q.fitted), w, Y(0) - Y(q.fitted), ISO, ISO, LW.axis);
    });
    axisX(f, X, p.y1, D.points.map(function (_, i) { return i; }),
          'reduction treatment', function (i) { return D.points[i].sample; });
    axisY(f, Y, p.x0, K.niceTicks(0, hi * 1.15, 5),
          'initial CO rate (µmol g⁻¹ s⁻¹)');
    K.legend(f, p.x0 + 7, p.y0 + 11, [
      { col: C.ink, swatch: true, text: 'measured' },
      { col: ISO, swatch: true, text: 'fitted, ∝ n(isolated)' }
    ]);
    return f.done();
  }

  return { population: population, rates: rates,
           partitionData: partitionData, ratesData: ratesData,
           FIGURES: { population: population, rates: rates } };
});
