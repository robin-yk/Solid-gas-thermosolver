/* The two figures the surface-capacity bound exists to produce.

   Both draw the same comparison from opposite ends. The first puts the
   measured oxygen removal next to the oxygen the (110) surface could
   hold, on one axis, so the gap is a length rather than a ratio in a
   sentence. The second turns that gap into the bound itself: the largest
   fraction of the measured inventory that could be surface oxygen, and
   by complement the smallest fraction that must lie below it.

   Nothing here is a partition. The bar's upper band is "not surface",
   not "subsurface" and not "bulk" - this repository does not distinguish
   them and the figure must not imply that it does. */

(function (root, factory) {
  'use strict';
  var kit = (typeof module !== 'undefined' && module.exports)
    ? require('./figkit.js') : root.FigKit;
  if (!kit) throw new Error('figures_surface.js needs figkit.js');
  var api = factory(kit);
  if (typeof module !== 'undefined' && module.exports) module.exports = api;
  root.SurfaceFigures = api;
})(typeof self !== 'undefined' ? self : this, function (K) {
  'use strict';

  var C = K.C, LW = K.LW, T = K.FIGTYPE;
  var lin = K.lin, frame = K.frame, axisX = K.axisX, axisY = K.axisY;
  var series = K.series, tint = K.tint;

  var SURF = C.surface, DEEP = C.extended;

  /* --------------------------------------------------------- payloads */

  function capacityData(rows, d_um, markInventory) {
    var sel = rows.filter(function (r) {
      return Math.abs(r.diameter_um - d_um) < 1e-9;
    }).sort(function (a, b) {
      return a.measured_umol_O_g - b.measured_umol_O_g;
    });
    if (!sel.length) return null;
    return {
      d_um: d_um,
      area_m2_g: sel[0].area_m2_g,
      ceiling: sel[0].reconstruction_deficiency_umol_O_g,
      full_row: sel[0].bridging_capacity_umol_O_g,
      bars: sel.map(function (r) {
        return { measured: r.measured_umol_O_g,
                 surface_max: r.reconstruction_deficiency_umol_O_g,
                 below_min: r.measured_umol_O_g
                   - r.reconstruction_deficiency_umol_O_g };
      }),
      mark: markInventory == null ? null : Number(markInventory)
    };
  }

  function boundData(rows, d_um, markInventory) {
    var sel = rows.filter(function (r) {
      return Math.abs(r.diameter_um - d_um) < 1e-9;
    }).sort(function (a, b) {
      return a.measured_umol_O_g - b.measured_umol_O_g;
    });
    if (!sel.length) return null;
    var ceiling = sel[0].reconstruction_deficiency_umol_O_g;
    var full = sel[0].bridging_capacity_umol_O_g;
    var lo = sel[0].measured_umol_O_g, hi = sel[sel.length - 1].measured_umol_O_g;
    var recon = [], row = [];
    for (var i = 0; i <= 120; i++) {
      var v = lo + (hi - lo) * i / 120;
      recon.push([v, Math.min(ceiling / v, 1) * 100]);
      row.push([v, Math.min(full / v, 1) * 100]);
    }
    var at = null;
    if (markInventory != null && markInventory > 0) {
      at = { measured: Number(markInventory),
             pct: Math.min(ceiling / markInventory, 1) * 100,
             pct_full: Math.min(full / markInventory, 1) * 100 };
    }
    return { d_um: d_um, ceiling: ceiling, full_row: full,
             lo: lo, hi: hi, recon: recon, row: row, at: at };
  }

  /* ---------------------------------------------------------- drawing */

  function capacity(D) {
    var f = K.square();
    var p = f.pane;
    if (!D) return f.done();

    var top = Math.max.apply(null, D.bars.map(function (b) {
      return b.measured;
    }));
    var ys = lin(0, top * 1.08, p.y1, p.y0);
    var n = D.bars.length;
    var xs = lin(-0.5, n - 0.5, p.x0, p.x1);
    var wid = (p.x1 - p.x0) / n * 0.62;

    frame(f, [p.x0, p.x1], [p.y0, p.y1]);

    D.bars.forEach(function (b, i) {
      var cx = xs(i);
      var yTop = ys(b.measured), yCut = ys(b.surface_max), y0 = ys(0);
      f.rect(cx - wid / 2, yTop, wid, yCut - yTop,
             tint(DEEP, 0.22), DEEP, LW.axis);
      f.rect(cx - wid / 2, yCut, wid, y0 - yCut, SURF, SURF, LW.axis);
      if (D.mark != null && Math.abs(b.measured - D.mark) < 1e-9) {
        f.line(cx, yTop - 9, cx, yTop - 2.5, C.ink, LW.curve);
      }
    });

    /* the whole bridging row, which nothing has been observed to strip */
    f.line(p.x0, ys(D.full_row), p.x1, ys(D.full_row), C.ink, LW.axis, '3 2');
    f.text(p.x1 - 3, ys(D.full_row) - 4, 'whole bridging row',
           { size: T.small, anchor: 'end' });

    axisX(f, xs, p.y1, D.bars.map(function (_, i) { return i; }),
          'oxygen removed by titration (µmol-O g⁻¹)',
          function (i) { return String(Math.round(D.bars[i].measured)); });
    axisY(f, ys, p.x0, K.niceTicks(0, top * 1.08, 5),
          'µmol-O g⁻¹');
    K.legend(f, p.x0 + 8, p.y0 + 12, [
      { col: SURF, text: 'could be surface (at most)' },
      { col: DEEP, swatch: true, text: 'must be below the surface' }
    ]);
    return f.done();
  }

  function bound(D) {
    var f = K.square();
    var p = f.pane;
    if (!D) return f.done();

    var top = Math.min(100, Math.max.apply(null, D.row.map(function (q) {
      return q[1];
    })) * 1.15);
    var xs = lin(D.lo, D.hi, p.x0, p.x1);
    var ys = lin(0, top, p.y1, p.y0);
    frame(f, [p.x0, p.x1], [p.y0, p.y1]);

    series(f, D.row, xs, ys, C.ink, LW.axis, '3 2');
    series(f, D.recon, xs, ys, SURF, LW.curve);

    if (D.at) {
      f.line(xs(D.at.measured), ys(0), xs(D.at.measured), ys(D.at.pct),
             C.ink, LW.axis, '1 2');
      f.dot(xs(D.at.measured), ys(D.at.pct), K.MARK.r, SURF);
      f.text(xs(D.at.measured) + 5, ys(D.at.pct) - 5,
             D.at.pct.toFixed(1) + '% at most',
             { size: T.small });
    }

    axisX(f, xs, p.y1, K.niceTicks(D.lo, D.hi, 5),
          'oxygen removed by titration (µmol-O g⁻¹)');
    axisY(f, ys, p.x0, K.niceTicks(0, top, 5),
          'largest possible surface share (%)');
    K.legend(f, p.x0 + 8, p.y0 + 12, [
      { col: SURF, text: '(1×2) reconstructed surface' },
      { col: C.ink, dash: '3 2', text: 'whole bridging row' }
    ]);
    return f.done();
  }

  return { capacity: capacity, bound: bound,
           capacityData: capacityData, boundData: boundData,
           FIGURES: { capacity: capacity, bound: bound } };
});
