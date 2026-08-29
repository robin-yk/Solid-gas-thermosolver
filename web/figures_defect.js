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
  var TICK = K.TICK;
  var lin = K.lin, lg = K.lg, frame = K.frame, n = K.n2;
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

  /* ------------------------------------------------ particle schematic */

  /* Placement, not physics. The subsurface class has one occupancy and no
     structure inside it, so which cells carry the drawn dots is a drawing
     choice; it is made by a fixed generator here rather than by the
     engine's RNG so the figure redraws identically from its own inputs. */
  function placer(seed) {
    var s = seed >>> 0;
    return function () {
      s = (s * 1664525 + 1013904223) >>> 0;
      return s / 4294967296;
    };
  }

  function pickCells(n, k, rnd) {
    var idx = [], t;
    for (t = 0; t < n; t++) idx.push(t);
    for (t = 0; t < k && t < n; t++) {
      var sw = t + Math.floor(rnd() * (n - t));
      var tmp = idx[t]; idx[t] = idx[sw]; idx[sw] = tmp;
    }
    return idx.slice(0, Math.min(k, n));
  }

  /* A vacancy on one oxygen site. In accessibility mode the fill strength
     is the refill probability and a locked site is drawn hollow, so the
     same dot carries both what is there and whether CO2 can reach it. */
  function vacDot(f, x, y, col, prob, accMode, r) {
    if (!accMode) return f.dot(x, y, r, col);
    if (prob < 0.02) return f.dot(x, y, r, C.paper, col);
    return f.dot(x, y, r, tint(col, Math.max(0.18, prob)));
  }

  /* The cross-section: a true-scale disc, the shell magnified until one
     drawn dot is one oxygen site, and the same three occupancies against
     depth. Landscape, because the disc and the magnified window have to
     sit side by side for the magnification to read as a magnification. */
  function particle(M) {
    var f = K.square({ wide: true });
    var p = M.pv;
    var cell = p.cell_nm, t1 = p.layer_nm, nLay = p.n_layers;
    var acc = p.accMode;
    var core = Math.min(0.30, 0.05
      + 0.25 * Math.sqrt(Math.max(0, p.theta.bulk)));

    /* ---- a: the particle at true scale, where the shell cannot be seen */
    var cx = 90, cy = 100, R = 66;
    f.dot(cx, cy, R, tint(C.structure, 0.30));
    f.dot(cx, cy, R - 1.4, tint(C.bulk, core));
    f.el('circle', { cx: cx, cy: cy, r: R, fill: 'none',
      stroke: C.surface, 'stroke-width': 1.1 });
    f.el('circle', { cx: cx, cy: cy, r: R - 1.3, fill: 'none',
      stroke: C.subsurface, 'stroke-width': 0.9 });
    f.panel(4, 12, 'a');
    f.text(cx, cy - 1, 'bulk occupation',
      { size: T.small, anchor: 'middle' });
    f.text(cx, cy + 9, 'θb = ' + p.theta_bulk_txt,
      { size: T.small, anchor: 'middle', fill: C.structure });
    f.text(cx, cy - R - 7, 'shell ' + p.shell_nm.toFixed(2)
      + ' nm, not resolvable here',
      { size: T.small, anchor: 'middle', fill: C.structure });
    var dy = cy + R + 12;
    f.line(cx - R, dy, cx + R, dy, C.ink, LW.hair);
    f.line(cx - R, dy - 3, cx - R, dy + 3, C.ink, LW.hair);
    f.line(cx + R, dy - 3, cx + R, dy + 3, C.ink, LW.hair);
    f.text(cx, dy + 11, 'd = ' + p.d_um.toFixed(2) + ' μm',
      { size: T.tick, anchor: 'middle' });

    /* ---- b: the shell, magnified until one dot is one oxygen site */
    var x0 = 214, y0 = 20, fw = 284, fh = 96;
    var nCols = 14;
    var sitePx = fw / nCols;
    var pxNm = sitePx / cell;               /* isotropic: same both ways */
    var hPx = t1 * pxNm / 2;                /* O sub-row spacing */
    var gasH = 15;
    f.rect(x0, y0, fw, fh, C.paper);
    f.rect(x0, y0, fw, gasH, tint(C.gas, 0.10));
    f.text(x0 + fw - 4, y0 + 10, 'gas',
      { size: T.small, anchor: 'end', fill: C.structure });
    f.panel(x0, 12, 'b');

    var surfY = y0 + gasH + 10;
    var subY = function (i) { return surfY + i * hPx; };
    var rows = 2 * nLay;
    var nSub = Math.floor((y0 + fh - 4 - surfY) / hPx);
    f.rect(x0 + 0.6, surfY - hPx / 2, fw - 1.2, (rows + 1) * hPx,
      tint(C.subsurface, 0.09));

    /* the oxygen sublattice, one dashed rule per sub-row */
    function latticeRow(y, period, off) {
      f.el('line', { x1: x0 + 0.6, x2: x0 + fw - 0.6, y1: y, y2: y,
        stroke: C.structure, 'stroke-width': 1.1, opacity: 0.45,
        'stroke-dasharray': '1.1 ' + (period - 1.1).toFixed(2),
        'stroke-dashoffset': (-off).toFixed(2) });
    }
    latticeRow(surfY, sitePx, 0.5 * sitePx - 0.55);
    for (var i = 1; i <= nSub; i++) {
      latticeRow(subY(i), sitePx / 2, 0.25 * sitePx - 0.55);
    }

    var siteX = function (j) { return x0 + (j + 0.5) * sitePx; };
    var siteX2 = function (j, s) {
      return x0 + (j + 0.25 + 0.5 * s) * sitePx;
    };
    /* the bridging row carries the sampled arrangement */
    for (var j = 0; j < nCols; j++) {
      if (p.surf_occ[j % p.surf_occ.length]) {
        vacDot(f, siteX(j), surfY, C.surface, p.P.surface, acc, 2.0);
      }
    }
    /* the subsurface class carries its expected count on seeded sites */
    var rnd = placer(20240517);
    var nSites = rows * nCols * 2;
    pickCells(nSites, Math.round(nSites * p.theta.subsurface), rnd)
      .forEach(function (ix) {
        var s = ix % 2, rest = (ix - s) / 2;
        var jj = rest % nCols, rr = (rest - jj) / nCols;
        vacDot(f, siteX2(jj, s), subY(1 + rr), C.subsurface,
          p.P.subsurface, acc, 1.7);
      });
    /* the bulk class spans the whole interior and the window shows a
       sliver of it, so it stays a class-average tint: placing individual
       deep vacancies would invent depth structure the partition does not
       carry */
    var coreTop = subY(rows) + hPx / 2;
    f.rect(x0 + 0.6, coreTop, fw - 1.2,
      Math.max(0, y0 + fh - 1.1 - coreTop), tint(C.bulk, core));

    var brX = x0 + 4;
    f.path('M' + n(brX + 3) + ',' + n(surfY - hPx / 2) + ' L' + n(brX) + ','
      + n(surfY - hPx / 2) + ' L' + n(brX) + ',' + n(coreTop) + ' L'
      + n(brX + 3) + ',' + n(coreTop), C.subsurface, LW.guide);
    f.text(brX + 6, coreTop + 8, 'shell ' + p.shell_nm.toFixed(2)
      + ' nm · 1 + ' + (4 * nLay) + ' O per cell',
      { size: T.small, fill: C.subsurface });
    var sbY = y0 + fh - 8;
    f.line(x0 + 8, sbY, x0 + 8 + pxNm, sbY, C.ink, 1.4);
    /* the frame goes on last so the fills cannot creep over it */
    f.line(x0 + 8, sbY - 2.5, x0 + 8, sbY + 2.5, C.ink, 1.4);
    f.line(x0 + 8 + pxNm, sbY - 2.5, x0 + 8 + pxNm, sbY + 2.5, C.ink, 1.4);
    f.text(x0 + 12 + pxNm, sbY + 3, '1 nm', { size: T.small });
    f.rect(x0, y0, fw, fh, 'none', C.ink, LW.axis);

    /* ---- c: the same three occupancies against depth */
    var px0 = 250, pw = 248, papTop = 148, papBot = 210;
    var linW = 0.46 * pw, logW = 0.44 * pw;
    var lgMax = Math.log10(p.R_nm);
    var XD = function (d) {
      if (d <= 1) return px0 + d * linW;
      return px0 + linW + 10 + Math.log10(d) / lgMax * logW;
    };
    /* six decades: a thick declared shell drives the bulk class into the
       1e-5 range, which four decades would flatten onto the floor */
    var YT = function (th) {
      return papBot - (Math.log10(Math.max(th, 1e-6)) + 6) / 6
        * (papBot - papTop);
    };
    f.panel(218, 140, 'c');
    var yTicks = [[1, '100%'], [0.01, '1%'], [0.0001, '0.01%'],
     [1e-6, '10⁻⁴%']];
    var bx = px0 + linW + 5;
    f.path('M' + n(bx - 2) + ',' + n(papBot + 3) + ' L' + n(bx + 1) + ','
      + n(papBot - 3), C.structure, LW.hair);
    f.path('M' + n(bx + 2) + ',' + n(papBot + 3) + ' L' + n(bx + 5) + ','
      + n(papBot - 3), C.structure, LW.hair);
    var tBr = t1 / 4;
    [[0, tBr, p.theta.surface, C.surface],
     [tBr, p.shell_nm, p.theta.subsurface, C.subsurface],
     [p.shell_nm, p.R_nm, p.theta.bulk, C.bulk]].forEach(function (s, si, a) {
      var y = YT(s[2]);
      if (si === 2) {
        f.line(XD(s[0]), y, px0 + linW, y, s[3], 2.3);
        f.line(px0 + linW + 10, y, XD(s[1]), y, s[3], 2.3);
      } else {
        f.line(XD(s[0]), y, XD(s[1]), y, s[3], 2.3);
      }
      if (si) {
        f.line(XD(s[0]), YT(a[si - 1][2]), XD(s[0]), y, C.structure,
          LW.hair);
      }
    });
    f.rect(px0, papTop, pw, papBot - papTop, 'none', C.ink, LW.axis);
    yTicks.forEach(function (g) {
      /* the outer two land on the frame itself, so they are labelled but
         not ticked - a tick drawn on the spine is not a tick */
      var y = YT(g[0]);
      if (y > papTop + 0.5 && y < papBot - 0.5) {
        f.line(px0, y, px0 + TICK.major, y, C.ink, LW.axis);
        f.line(px0 + pw, y, px0 + pw - TICK.major, y, C.ink, LW.axis);
      }
      f.text(px0 - 3, y + 2.8, g[1], { size: T.tick, anchor: 'end' });
    });
    [[0, '0'], [p.shell_nm, p.shell_nm.toFixed(2)], [1, '1'],
     [10, '10'], [100, '100'], [p.R_nm, Math.round(p.R_nm) + '']]
      .forEach(function (g) {
        if (g[0] > p.R_nm) return;
        var x = XD(g[0]);
        if (x > px0 + 0.5 && x < px0 + pw - 0.5) {
          f.line(x, papBot, x, papBot - TICK.major, C.ink, LW.axis);
          f.line(x, papTop, x, papTop + TICK.major, C.ink, LW.axis);
        }
        f.text(x, papBot + 9, g[1], { size: T.tick, anchor: 'middle' });
      });
    f.text(px0 + pw / 2, papBot + 21, 'depth from surface (nm)',
      { size: T.body, anchor: 'middle' });
    f.text(XD(p.shell_nm) + 5, papTop + 20, p.shell_pct.toFixed(1)
      + '% of the inventory within ' + p.shell_nm.toFixed(2) + ' nm',
      { size: T.small });

    K.legend(f, 6, 206, [
      { swatch: true, col: C.surface, text: 'surface, 1 O per cell' },
      { swatch: true, col: C.subsurface,
        text: 'subsurface, ' + (4 * nLay) + ' O per cell' },
      { swatch: true, col: C.bulk, text: 'bulk, drawn as a tint' }
    ], 11);
    f.text(6, 243, acc ? 'fill = P(CO₂ refill), hollow = locked'
      : 'one dot = one vacancy on one O site',
      { size: T.small, fill: C.structure });
    return f.done();
  }

  return {
    CLASSES: CLASSES,
    distribution: distribution, accessible: accessible,
    recovery: recovery, radial: radial, particle: particle,
    FIGURES: { distribution: distribution, accessible: accessible,
               recovery: recovery, radial: radial, particle: particle }
  };
});
