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

  function regionFractions(pools, cap) {
    var capacities = [cap.bridging + cap.basal_L1,
      cap.L1_subbridging + cap.subsurface_L2_4];
    capacities.push(2e6 / 79.866 - capacities[0] - capacities[1]);
    var amounts = [pools.bridging + pools.reconstructed_row + pools.basal_L1,
      pools.L1_subbridging + pools.subsurface_L2_4, pools.bulk];
    return ['Surface', 'Subsurface', 'Bulk'].map(function (region, i) {
      return {region: region, amount: amounts[i], capacity: capacities[i], percent: 100 * amounts[i] / capacities[i]};
    });
  }

  function vacancyProfile(D) {
    var f = K.square({wide: true});
    var p = {x0: 70, y0: 36, x1: f.W - 16, y1: f.pane.y1};
    var colors = ['#0072B2', '#E69F00', '#009E73', '#CC79A7', '#D55E00'];
    var shapes = ['circle', 'square', 'triangle', 'circle', 'square'];
    var vals = D.scenarios.reduce(function (a, s) { return a.concat(s.points.map(function (r) {return r.percent;})); }, []).filter(function(v) {return v > 0;});
    var Y = K.lg(K.decadeFloor(Math.min.apply(null, vals) / 1.2), 100, p.y1, p.y0);
    var mid = (p.x0 + p.x1) / 2;
    var axes = [K.lin(0, 1.3, p.x0, mid - 24), K.lg(1.3, 450, mid + 24, p.x1)];
    K.frame(f, [p.x0, mid - 24], [p.y0, p.y1]);
    K.frame(f, [mid + 24, p.x1], [p.y0, p.y1]);
    D.scenarios.forEach(function (s, j) {
      K.legend(f, p.x0 + j * 80, 19, [{col: colors[j], marker: shapes[j], filled: true, text: s.id}]);
      axes.forEach(function(X, panel) {
        var pts = s.points.filter(function(r) {return r.percent > 0 && (panel ? r.depth_nm >= 1.3 : r.depth_nm < 1.3);});
        f.el('polyline', {points:pts.map(function(r) {return X(r.depth_nm) + ',' + Y(r.percent);}).join(' '), fill:'none', stroke:colors[j], 'stroke-width':1.3});
        pts.forEach(function (r) {
          f.raw('<g data-scenario="' + s.id + '" data-depth="' + r.depth_nm + '" data-percent="' + r.percent + '">');
          f.el('title', {}, K.esc(s.id + ', ' + r.depth_nm.toPrecision(4) + ' nm: ' + r.percent.toPrecision(4) + '%'));
          if (!panel) K.marker(f, shapes[j], X(r.depth_nm), Y(r.percent), K.MARK * .65, colors[j], true);
          f.raw('</g>');
        });
      });
    });
    K.axisX(f, axes[0], p.y1, [0, .4, .8, 1.2], 'Depth (nm)', function(v) {return String(v);});
    K.axisX(f, axes[1], p.y1, [2, 10, 100, 450], 'Depth (nm)', function(v) {return String(v);});
    K.axisY(f, Y, p.x0, K.decades(Y.d0, Y.d1), 'Vacancy fraction (%)', K.powLabel);
    return f.done();
  }
  return { tofRange: tofRange, regionFractions: regionFractions, vacancyProfile: vacancyProfile, FIGURES: { tofRange: tofRange } };
});
