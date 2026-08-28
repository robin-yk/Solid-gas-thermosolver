/* The defect workspace figures.

   Five figures, each one publication output rather than a screen chart:
   four square 3.5 in panels and one landscape schematic. Every one is
   drawn through web/figkit.js, so the style, the colour roles and the
   coordinate system are the kit's and nothing is set here. What the page
   previews is the file that downloads.

   Each builder takes the same model object, so the download menu can be
   wired once:

     M = { T_C, vo, geom, out, sweep: {rows, xmax}, bounds,
           cg: {r, sys, tEnd, target_pct}, acc, surf }

   No builder computes physics. They read solver output and draw it. */

(function (root, factory) {
  'use strict';
  var kit = (typeof module !== 'undefined' && module.exports)
    ? require('./figkit.js') : root.FigKit;
  if (!kit) throw new Error('figures_defect.js needs figkit.js');
  var api = factory(kit);
  if (typeof module !== 'undefined' && module.exports) module.exports = api;
  root.DefectFigures = api;
})(typeof self !== 'undefined' ? self : this, function (K) {
  'use strict';

  var C = K.C, LW = K.LW, MARK = K.MARK, T = K.FIGTYPE;
  var lin = K.lin, lg = K.lg, frame = K.frame;
  var axisX = K.axisX, axisY = K.axisY, series = K.series;
  var legend = K.legend, expLabel = K.expLabel, tint = K.tint;

  /* the classes, in depth order, with the line style that carries them
     into greyscale */
  var CLASSES = [
    { key: 'surface', label: 'surface', col: C.surface, dash: null },
    { key: 'subsurface', label: 'subsurface', col: C.subsurface,
      dash: '5 2.5' },
    { key: 'bulk', label: 'bulk', col: C.bulk, dash: '1.6 2' },
    { key: 'extended', label: 'extended', col: C.extended, dash: '7 2 1.6 2' }
  ];

  function pane(f) { return f.pane; }

  /* Round the axis top up to a four-tick number, so the labels read
     0/25/50/75/100 rather than 0/17.4/34.8. */
  function niceMax(v) {
    if (!(v > 0)) return 1;
    var mag = Math.pow(10, Math.floor(Math.log10(v / 4)));
    var mult = [1, 1.25, 1.5, 2, 2.5, 3, 4, 5, 6, 7.5, 10];
    for (var i = 0; i < mult.length; i++) {
      if (4 * mult[i] * mag >= v) return 4 * mult[i] * mag;
    }
    return 40 * mag;
  }

  function ticks4(max) {
    return [0, max / 4, max / 2, 3 * max / 4, max];
  }

  function fmt(v) {
    if (v === 0) return '0';
    if (Math.abs(v) >= 100) return String(Math.round(v));
    if (Math.abs(v) >= 10) return v.toFixed(0);
    if (Math.abs(v) >= 1) return v.toFixed(1);
    return v.toFixed(2);
  }

  /* ------------------------------------------------- vacancy distribution */

  /* A square bar plot over the four classes. The numbers live in the table
     beside the figure, not on the bars: a value printed twice is a value
     that can disagree with itself. */
  function distribution(M) {
    var f = K.square();
    var p = pane(f);
    var rows = CLASSES.map(function (c) {
      return { c: c, v: c.key === 'extended'
        ? (M.out.extended_defects_umol_g || 0) : M.out.umol_g[c.key] };
    });
    var vmax = 0;
    rows.forEach(function (r) { if (r.v > vmax) vmax = r.v; });
    var ymax = niceMax(vmax * 1.08);
    var Y = lin(0, ymax, p.y1, p.y0);
    var n = rows.length;
    var slot = (p.x1 - p.x0) / n;
    var bw = slot * 0.46;

    frame(f, [p.x0, p.x1], Y);
    rows.forEach(function (r, i) {
      var xc = p.x0 + (i + 0.5) * slot;
      f.rect(xc - bw / 2, Y(r.v), bw, p.y1 - Y(r.v),
        tint(r.c.col, 0.34), r.c.col, LW.axis);
      f.text(xc, p.y1 + T.tick + 3, r.c.label,
        { size: T.tick, anchor: 'middle' });
    });
    axisY(f, Y, p.x0, ticks4(ymax), 'vacancy inventory (µmol-O g⁻¹)', fmt);
    f.text((p.x0 + p.x1) / 2, p.y1 + 2 * T.tick + 10.5, 'site class',
      { size: T.body, anchor: 'middle' });
    return f.done();
  }

  /* --------------------------------------------------- accessible inventory */

  /* Total against accessible along the inventory axis. Phase boundaries are
     thin light rules and only the two that carry the mechanism are named. */
  function accessible(M) {
    var f = K.square();
    var p = pane(f);
    var rows = M.sweep.rows;
    var xmax = M.sweep.xmax;
    var vmax = 0;
    rows.forEach(function (r) {
      if (r.VO_total_umol_g <= xmax && r.VO_total_umol_g > vmax) {
        vmax = r.VO_total_umol_g;
      }
    });
    var ymax = niceMax(vmax * 1.04);
    var X = lin(0, xmax, p.x0, p.x1), Y = lin(0, ymax, p.y1, p.y0);

    if (M.bounds) {
      [[M.bounds.VO_recon_complete_umol_g, '(1×2) complete'],
       [M.bounds.VO_cs_umol_g, 'CS ceiling']].forEach(function (b) {
        if (b[0] > 0 && b[0] < xmax) {
          f.line(X(b[0]), p.y0, X(b[0]), p.y1, C.structure, LW.guide,
            '3 2.2');
          f.text(X(b[0]) + 3, p.y0 + T.small + 2, b[1],
            { size: T.small, fill: C.structure });
        }
      });
    }
    var inb = rows.filter(function (r) { return r.VO_total_umol_g <= xmax; });
    series(f, inb.map(function (r) {
      return [r.VO_total_umol_g, r.VO_total_umol_g]; }), X, Y, C.ink, 2.3);
    series(f, inb.map(function (r) {
      return [r.VO_total_umol_g, r.accessible]; }), X, Y, C.surface, 2.3,
      '5 2.5');
    if (M.vo != null && M.acc) {
      f.dot(X(M.vo), Y(M.acc.accessible_umol_g), MARK, C.surface);
    }
    frame(f, X, Y);
    axisX(f, X, p.y1, ticks4(xmax), 'total vacancy content (µmol-O g⁻¹)',
      fmt);
    axisY(f, Y, p.x0, ticks4(ymax), 'inventory (µmol-O g⁻¹)', fmt);
    legend(f, p.x0 + 12, p.y0 + 26, [
      { col: C.ink, text: 'total' },
      { col: C.surface, dash: '5 2.5', text: 'CO₂-accessible' }
    ]);
    return f.done();
  }

  /* ---------------------------------------------------------- CGMC recovery */

  /* Recovered fraction against exposure on a log time axis. Direct labels
     rather than a legend: two curves, each with somewhere to put its name. */
  function recovery(M) {
    var f = K.square();
    var p = pane(f);
    var r = M.cg.r;
    var lo = -6, hi = Math.max(lo + 1, Math.ceil(Math.log10(M.cg.tEnd)));
    var X = lg(Math.pow(10, lo), Math.pow(10, hi), p.x0, p.x1);
    var Y = lin(0, 100, p.y1, p.y0);
    var tmin = Math.pow(10, lo);
    var pts = [];
    for (var i = 0; i < r.t_s.length; i++) {
      /* the run starts below the drawn decade; clipping here keeps the
         curve inside its box instead of running out of the panel */
      if (r.t_s[i] >= tmin) pts.push([r.t_s[i], r.recovery[i] * 100]);
    }
    frame(f, X, Y);
    series(f, pts, X, Y, C.subsurface, 2.3);
    /* name the curve where it is rising, and the measured point under
       the marker, so the two annotations cannot meet */
    var mid = null;
    for (var j = 0; j < pts.length; j++) {
      if (pts[j][1] >= 50) { mid = pts[j]; break; }
    }
    mid = mid || pts[pts.length - 1];
    if (mid) {
      f.text(X(mid[0]) + 7, Y(mid[1]) + 4,
        'E\u2098 = ' + M.cg.E_m_eV.toFixed(2) + ' eV',
        { size: T.small, fill: C.subsurface });
    }
    if (M.cg.target_pct != null && M.cg.target_s != null) {
      f.dot(X(M.cg.target_s), Y(M.cg.target_pct), MARK, C.ink);
      f.text(X(M.cg.target_s) - 5, Y(M.cg.target_pct) + T.small + 5,
        'measured', { size: T.small, anchor: 'end' });
    }
    axisX(f, X, p.y1, [1e-6, 1e-4, 1e-2, 1, 100], 'exposure time (s)',
      expLabel);
    axisY(f, Y, p.x0, [0, 25, 50, 75, 100], 'recovered (%)');
    return f.done();
  }

  /* ------------------------------------------------------------ radial profile */

  /* Occupancy against depth from the surface, at a few representative
     times. Four curves at most, so the panel stays readable. */
  function radial(M) {
    var f = K.square();
    var p = pane(f);
    var r = M.cg.r;
    var d = r.depth_nm || r.centers_nm;
    var all = r.profiles || [];
    /* four snapshots at most, spread over the run rather than the first
       four, which would all sit before anything has moved */
    var profs = [];
    var take = Math.min(4, all.length);
    for (var j = 0; j < take; j++) {
      profs.push(all[Math.round(j * (all.length - 1) / (take - 1 || 1))]);
    }
    /* floor the occupancy axis just under the smallest value drawn, so a
       drained particle is not a flat line along an arbitrary decade */
    /* six decades below the largest occupancy drawn. A drained particle
       runs to arbitrarily small numbers and an axis that chases them
       turns the panel into noise. */
    var hiTh = 0;
    profs.forEach(function (pr) {
      pr.theta.forEach(function (th) { if (th > hiTh) hiTh = th; });
    });
    var top = Math.pow(10, Math.ceil(Math.log10(hiTh || 1)));
    if (!(top > 0)) top = 1;
    var floor = top * 1e-6;
    var X = lg(Math.max(d[0] * 0.7, 0.05), d[d.length - 1] * 1.15,
      p.x0, p.x1);
    var Y = lg(floor, top, p.y1, p.y0);
    frame(f, X, Y);
    profs.forEach(function (pr, i) {
      var dash = [null, '5 2.5', '1.6 2', '7 2 1.6 2'][i];
      series(f, pr.theta.map(function (th, k) {
        return [d[k], Math.max(th, floor)]; }), X, Y, C.bulk, 2.3, dash);
      f.text(p.x1 - 5, p.y0 + 12 + i * (T.small + 3), tfmt(pr.t_s),
        { size: T.small, anchor: 'end', fill: C.bulk });
    });
    var dec = [];
    for (var e = Math.round(Math.log10(floor));
         e <= Math.round(Math.log10(top)); e += 2) {
      dec.push(Math.pow(10, e));
    }
    axisX(f, X, p.y1, [0.1, 1, 10, 100, 1000], 'depth from surface (nm)',
      function (v) { return v < 1 ? String(v) : String(Math.round(v)); });
    axisY(f, Y, p.x0, dec, 'vacancy occupancy', expLabel);
    return f.done();
  }

  function tfmt(t) {
    if (t >= 1) return Math.round(t) + ' s';
    if (t >= 1e-3) return (t * 1e3).toFixed(0) + ' ms';
    return t.toExponential(0) + ' s';
  }

  return {
    CLASSES: CLASSES,
    distribution: distribution, accessible: accessible,
    recovery: recovery, radial: radial,
    FIGURES: { distribution: distribution, accessible: accessible,
               recovery: recovery, radial: radial }
  };
});
