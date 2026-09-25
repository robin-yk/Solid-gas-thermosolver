/* The figure the one-model calculation exists to produce: the apparent CO
   turnover frequency of each sample, as the range it spans over every
   parameter point with a countable site denominator, against the single
   value the supplement's fixed denominator gives. The scenario chosen on
   the page is one filled marker per sample, so the reader sees where one
   set of assumptions falls inside the range. */

(function (root, factory) {
  'use strict';
  var kit = (typeof module !== 'undefined' && module.exports)
    ? require('./figkit.js') : root.FigKit;
  if (!kit) throw new Error('figures_distribution.js needs figkit.js');
  var api = factory(kit);
  if (typeof module !== 'undefined' && module.exports) module.exports = api;
  root.DistributionFigures = api;
})(typeof self !== 'undefined' ? self : this, function (K) {
  'use strict';

  var C = K.C, LW = K.LW, T = K.FIGTYPE;

  /* D.points: [{sample, min, median, max, fixed, nBelow, nOk, cur, curBelow}] */
  function tofRange(D) {
    var f = K.square({ wide: true });
    var p = { x0: 70, y0: 16, x1: f.W - 16, y1: f.pane.y1 };
    if (!D || !D.points.length) return f.done();
    var vals = [];
    D.points.forEach(function (q) {
      [q.min, q.max, q.fixed, q.cur].forEach(function (v) {
        if (isFinite(v) && v > 0) vals.push(v);
      });
    });
    var Y = K.lg(K.decadeFloor(Math.min.apply(null, vals)),
                 K.decadeCeil(Math.max.apply(null, vals)), p.y1, p.y0);
    var n = D.points.length;
    var X = K.lin(-0.5, n - 0.5, p.x0, p.x1);
    var cap = 7, gap = 12;

    K.frame(f, [p.x0, p.x1], [p.y0, p.y1]);
    var inside = function (v) { return Math.min(Math.max(Y(v), p.y0 + 2), p.y1 - 2); };
    D.points.forEach(function (q, i) {
      var x = X(i) - gap / 2;
      /* the range over the countable cases, with its median */
      f.line(x, Y(q.min), x, Y(q.max), C.ink, LW.curve);
      f.line(x - cap, Y(q.min), x + cap, Y(q.min), C.ink, LW.axis);
      f.line(x - cap, Y(q.max), x + cap, Y(q.max), C.ink, LW.axis);
      f.line(x - 1.6 * cap, Y(q.median), x + 1.6 * cap, Y(q.median), C.ink, LW.guide);
      /* the fixed-denominator value and the chosen scenario, side by side */
      K.marker(f, 'square', X(i) + gap, Y(q.fixed), K.MARK, C.extended);
      if (isFinite(q.cur) && q.cur > 0) {
        K.marker(f, 'circle', x, inside(q.cur), K.MARK, C.surface, !q.curBelow);
      }
      if (q.nBelow) {
        f.text(X(i), p.y1 - 8, q.nBelow + ' below threshold',
               { size: T.note, anchor: 'middle', fill: C.extended });
      }
    });
    K.axisX(f, X, p.y1, D.points.map(function (_, i) { return i; }),
            'reduction treatment', function (i) { return D.points[i].sample; });
    K.axisY(f, Y, p.x0, K.decades(Y.d0, Y.d1),
            'apparent TOF (s⁻¹)', K.powLabel);
    K.legend(f, p.x1 - 250, p.y0 + 17, [
      { col: C.ink, text: 'range over ' + D.nPoints + ' parameter points' },
      { col: C.surface, marker: 'circle', filled: true, text: 'scenario chosen on the page' },
      { col: C.extended, marker: 'square', text: 'fixed ' + D.fixed + ' µmol g⁻¹ (SI 2a)' }
    ]);
    return f.done();
  }

  return { tofRange: tofRange, FIGURES: { tofRange: tofRange } };
});
