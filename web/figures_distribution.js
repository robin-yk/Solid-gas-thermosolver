/* Apparent CO TOF for five explicit vacancy-distribution assumptions. */

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

  /* D.points: [{sample, scenarios: [{id, tof, below}]}] */
  function tofRange(D) {
    var f = K.square({ wide: true });
    var p = { x0: 70, y0: 36, x1: f.W - 16, y1: f.pane.y1 };
    if (!D || !D.points.length) return f.done();
    var vals = [];
    D.points.forEach(function (q) {
      q.scenarios.forEach(function (s) {
        var v = s.tof;
        if (isFinite(v) && v > 0) vals.push(v);
      });
    });
    var Y = K.lg(K.decadeFloor(Math.min.apply(null, vals)),
                 K.decadeCeil(Math.max.apply(null, vals)), p.y1, p.y0);
    var n = D.points.length;
    var X = K.lin(-0.5, n - 0.5, p.x0, p.x1);
    var colors = ['#0072B2', '#E69F00', '#009E73', '#CC79A7', '#D55E00'];
    var shapes = ['circle', 'square', 'triangle', 'circle', 'square'];

    K.frame(f, [p.x0, p.x1], [p.y0, p.y1]);
    D.points.forEach(function (q, i) {
      q.scenarios.forEach(function (s, j) {
        if (!isFinite(s.tof) || s.tof <= 0) return;
        f.raw('<g data-sample="' + K.esc(q.sample) + '" data-scenario="' + s.id + '" data-tof="' + s.tof + '">');
        f.el('title', {}, K.esc(q.sample + ', ' + s.id + ': ' + s.tof + ' s⁻¹'));
        K.marker(f, shapes[j], X(i) + (j - 2) * 12, Y(s.tof), K.MARK, colors[j], !s.below);
        f.raw('</g>');
      });
    });
    K.axisX(f, X, p.y1, D.points.map(function (_, i) { return i; }),
            'sample', function (i) { return D.points[i].sample; });
    K.axisY(f, Y, p.x0, K.decades(Y.d0, Y.d1),
            'TOF (s⁻¹)', K.powLabel);
    D.points[0].scenarios.forEach(function (s, j) {
      K.legend(f, p.x0 + j * 80, 19, [{col: colors[j], marker: shapes[j], filled: true, text: s.id}]);
    });
    return f.done();
  }

  return { tofRange: tofRange, FIGURES: { tofRange: tofRange } };
});
