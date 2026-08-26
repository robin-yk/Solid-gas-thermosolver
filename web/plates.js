/* Manuscript figure plates, drawn to print size from frozen solver output.

   Every plate is a string of SVG in points, sized in millimetres, so what
   the browser shows is what the page will print: Nature column widths, a
   5 pt type floor, flat fills, no gradients, no transparency, no external
   reference. The kit runs unchanged under node, which is how the gate in
   tests/test_plates.py measures type size and width without a browser.

   Colour carries physical role and nothing else - structure is grey. */

(function (root, factory) {
  'use strict';
  var api = factory();
  if (typeof module !== 'undefined' && module.exports) module.exports = api;
  root.Plates = api;
})(typeof self !== 'undefined' ? self : this, function () {
  'use strict';

  var PT_PER_MM = 72 / 25.4;
  var WIDTH_MM = { single: 89, onehalf: 120, double: 183 };
  var MAX_H_MM = 247;
  var FONT = 'Arial, Helvetica, sans-serif';
  var TYPE = { panel: 8, body: 7, tick: 6.5, small: 5.5 };
  var C = {
    surface: '#0072B2', subsurface: '#D55E00', bulk: '#CC79A7',
    extended: '#E69F00', gas: '#009E73', structure: '#6E6E6E',
    ink: '#000000', paper: '#FFFFFF'
  };
  /* a flat 14% tint of a role hue, for fills under a line of the same role */
  function tint(hex, f) {
    f = f == null ? 0.14 : f;
    var n = parseInt(hex.slice(1), 16);
    var r = (n >> 16) & 255, g = (n >> 8) & 255, b = n & 255;
    var m = function (v) { return Math.round(255 - (255 - v) * f); };
    return '#' + [m(r), m(g), m(b)].map(function (v) {
      return ('0' + v.toString(16)).slice(-2);
    }).join('');
  }

  /* ------------------------------------------------------------ build */

  function esc(s) {
    return String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;')
      .replace(/>/g, '&gt;');
  }
  function n2(v) { return (Math.round(v * 100) / 100).toString(); }
  var SUB = '\u2080\u2081\u2082\u2083\u2084\u2085\u2086\u2087'
    + '\u2088\u2089';
  /* Ti10O19 -> Ti(10 subscript)O(19 subscript); digits that do not follow an
     element letter, as in a 1:1 feed ratio, are left alone */
  function chem(s) {
    return String(s).replace(/([A-Za-z])(\d+)/g, function (_, a, d) {
      return a + d.replace(/\d/g, function (c) { return SUB[Number(c)]; });
    });
  }
  function attrs(o) {
    var s = '';
    for (var k in o) {
      if (o[k] == null || o[k] === '') continue;
      s += ' ' + k + '="' + (typeof o[k] === 'number' ? n2(o[k]) : o[k]) + '"';
    }
    return s;
  }

  function plate(widthKey, heightMm) {
    var wmm = WIDTH_MM[widthKey];
    if (!wmm) throw new Error('unknown plate width: ' + widthKey);
    if (heightMm > MAX_H_MM) throw new Error('plate taller than the page');
    var parts = [];
    var api = {
      widthKey: widthKey, wmm: wmm, hmm: heightMm,
      W: wmm * PT_PER_MM, H: heightMm * PT_PER_MM,
      raw: function (s) { parts.push(s); return api; },
      el: function (name, o, inner) {
        parts.push('<' + name + attrs(o)
          + (inner == null ? '/>' : '>' + inner + '</' + name + '>'));
        return api;
      },
      line: function (x1, y1, x2, y2, col, w, dash) {
        return api.el('line', { x1: x1, y1: y1, x2: x2, y2: y2,
          stroke: col || C.ink, 'stroke-width': w || 0.6,
          'stroke-dasharray': dash });
      },
      rect: function (x, y, w, h, fill, stroke, sw) {
        return api.el('rect', { x: x, y: y, width: Math.max(0, w),
          height: Math.max(0, h), fill: fill || 'none',
          stroke: stroke, 'stroke-width': sw || (stroke ? 0.6 : null) });
      },
      path: function (d, col, w, dash, fill) {
        return api.el('path', { d: d, fill: fill || 'none',
          stroke: col || C.ink, 'stroke-width': w || 0.9,
          'stroke-dasharray': dash, 'stroke-linejoin': 'round' });
      },
      dot: function (x, y, r, fill, stroke) {
        return api.el('circle', { cx: x, cy: y, r: r, fill: fill || C.ink,
          stroke: stroke, 'stroke-width': stroke ? 0.6 : null });
      },
      text: function (x, y, s, o) {
        o = o || {};
        var size = o.size || TYPE.body;
        if (size < 5) throw new Error('type below the 5 pt floor: ' + size);
        return api.el('text', { x: x, y: y, 'font-family': FONT,
          'font-size': size, 'font-weight': o.weight,
          'text-anchor': o.anchor, fill: o.fill || C.ink,
          transform: o.rotate == null ? null
            : 'rotate(' + o.rotate + ' ' + n2(x) + ' ' + n2(y) + ')'
        }, esc(s));
      },
      panel: function (x, y, ch) {
        return api.text(x, y, ch, { size: TYPE.panel, weight: 'bold' });
      },
      done: function () {
        return '<svg xmlns="http://www.w3.org/2000/svg" '
          + 'width="' + wmm + 'mm" height="' + heightMm + 'mm" '
          + 'viewBox="0 0 ' + n2(api.W) + ' ' + n2(api.H) + '" '
          + 'data-width-mm="' + wmm + '" data-height-mm="' + heightMm + '">'
          + '<rect x="0" y="0" width="' + n2(api.W) + '" height="'
          + n2(api.H) + '" fill="' + C.paper + '"/>'
          + parts.join('') + '</svg>';
      }
    };
    return api;
  }

  /* ----------------------------------------------------------- scales */

  function lin(d0, d1, r0, r1) {
    var f = function (v) { return r0 + (v - d0) / (d1 - d0) * (r1 - r0); };
    f.d0 = d0; f.d1 = d1; f.r0 = r0; f.r1 = r1;
    f.clamp = function (v) { return f(Math.max(d0, Math.min(d1, v))); };
    return f;
  }
  function log10(v) { return Math.log(v) / Math.LN10; }
  function lg(d0, d1, r0, r1) {
    var l0 = log10(d0), l1 = log10(d1);
    var f = function (v) {
      return r0 + (log10(Math.max(v, d0 * 1e-12)) - l0) / (l1 - l0)
        * (r1 - r0);
    };
    f.d0 = d0; f.d1 = d1; f.r0 = r0; f.r1 = r1;
    f.clamp = function (v) { return f(Math.max(d0, Math.min(d1, v))); };
    return f;
  }
  function niceTicks(d0, d1, want) {
    var span = d1 - d0, raw = span / (want || 4);
    var mag = Math.pow(10, Math.floor(Math.log(raw) / Math.LN10));
    var step = [1, 2, 2.5, 5, 10].map(function (m) { return m * mag; })
      .filter(function (v) { return v >= raw; })[0] || 10 * mag;
    var out = [], v = Math.ceil(d0 / step) * step;
    for (; v <= d1 + 1e-9; v += step) {
      out.push(Math.round(v / step) * step);
    }
    return out;
  }
  function decades(d0, d1) {
    var out = [];
    for (var e = Math.ceil(log10(d0)); e <= Math.floor(log10(d1)); e++) {
      out.push(Math.pow(10, e));
    }
    return out;
  }
  function expLabel(v) {
    var e = Math.round(log10(v));
    if (e >= 0 && e <= 3) return String(Math.round(v));
    return '10' + String(e).replace(/-/g, '−')
      .split('').map(function (ch) {
        return '⁰¹²³⁴⁵⁶⁷⁸⁹'
          [Number(ch)] || (ch === '−' ? '⁻' : ch);
      }).join('');
  }

  /* ------------------------------------------------------------- axes */

  function axisX(p, sc, y, ticks, label, fmt) {
    p.line(sc.r0, y, sc.r1, y, C.ink, 0.7);
    ticks.forEach(function (t) {
      var x = sc(t);
      p.line(x, y, x, y + 2.4, C.ink, 0.7);
      p.text(x, y + 9.5, fmt ? fmt(t) : String(t),
        { size: TYPE.tick, anchor: 'middle' });
    });
    if (label) {
      p.text((sc.r0 + sc.r1) / 2, y + 21, label,
        { size: TYPE.body, anchor: 'middle' });
    }
    return p;
  }
  function axisY(p, sc, x, ticks, label, fmt) {
    p.line(x, sc.r0, x, sc.r1, C.ink, 0.7);
    ticks.forEach(function (t) {
      var y = sc(t);
      p.line(x - 2.4, y, x, y, C.ink, 0.7);
      p.text(x - 4, y + 2.3, fmt ? fmt(t) : String(t),
        { size: TYPE.tick, anchor: 'end' });
    });
    if (label) {
      p.text(x - 26, (sc.r0 + sc.r1) / 2, label,
        { size: TYPE.body, anchor: 'middle', rotate: -90 });
    }
    return p;
  }
  function series(p, pts, xs, ys, col, w, dash) {
    var d = '', started = false;
    pts.forEach(function (q) {
      if (q[1] == null || !isFinite(q[1])) { started = false; return; }
      d += (started ? 'L' : 'M') + n2(xs(q[0])) + ',' + n2(ys(q[1]));
      started = true;
    });
    if (d) p.path(d, col, w || 1.0, dash);
    return p;
  }
  function legend(p, x, y, items, gap) {
    gap = gap || 9.5;
    items.forEach(function (it, i) {
      var yy = y + i * gap;
      if (it.dash) {
        p.line(x, yy - 2.2, x + 9, yy - 2.2, it.col, 1.0, it.dash);
      } else if (it.swatch) {
        p.rect(x, yy - 5.2, 7, 6, it.col);
      } else {
        p.line(x, yy - 2.2, x + 9, yy - 2.2, it.col, 1.4);
      }
      p.text(x + 12, yy, it.text, { size: TYPE.small });
    });
    return p;
  }

  /* ----------------------------------------------------------- plates */

  var PLATES = {};

  /* d1 - site capacity and the phase ladder */
  PLATES.d1 = function (D) {
    var g = D.plot.defect.geometry, b = D.plot.defect.boundaries;
    var sw = D.plot.defect.sweep, fl = D.plot.defect.flagship;
    var p = plate('double', 58);

    /* a: capacity of the three classes, log axis, with the measured
       inventory drawn across them */
    var ax0 = 34, ax1 = 196, ay = 116;
    var xs = lg(1, 1e5, ax0, ax1);
    p.panel(8, 16, 'a');
    p.text(ax0, 16, 'Oxygen-site capacity per class',
      { size: TYPE.body, weight: 'bold' });
    var caps = [
      { k: 'surface', label: 'bridging O', v: g.N_s },
      { k: 'subsurface', label: 'subsurface shell', v: g.N_ss },
      { k: 'bulk', label: 'bulk O', v: g.N_b }
    ];
    var xv = xs(fl.VO_total_umol_g);
    caps.forEach(function (c, i) {
      var y = 32 + i * 27;
      var xr = xs(c.v);
      p.rect(ax0, y, xr - ax0, 12, tint(C[c.k], 0.30), C[c.k], 0.7);
      p.text(ax0 + 2, y - 3, c.label, { size: TYPE.small, fill: C[c.k] });
      var v = c.v >= 1000 ? String(Math.round(c.v)) : c.v.toFixed(2);
      /* the value sits inside its own bar: outside it would collide with
         the measured-inventory line */
      p.text(xr - 3, y + 8.8, v, { size: TYPE.tick, anchor: 'end' });
    });
    p.line(xv, 26, xv, ay, C.ink, 0.8, '2.5 1.8');
    p.text(xv, 22, 'measured ' + fl.VO_total_umol_g + ' µmol-O g⁻¹',
      { size: TYPE.small, anchor: 'middle' });
    axisX(p, xs, ay, decades(1, 1e5),
      'capacity (µmol-O g⁻¹)', expLabel);

    /* b: class inventory against total */
    var bx0 = 268, bx1 = 508, by0 = 116, by1 = 30;
    var vmax = sw[sw.length - 1].vo;
    var ymax = 0;
    sw.forEach(function (r) {
      ymax = Math.max(ymax, r.surface, r.subsurface, r.bulk, r.extended);
    });
    ymax = Math.ceil(ymax / 20) * 20;
    /* log inventory axis: the ladder spans two decades and a linear axis
       hides the reconstruction riser entirely */
    var X = lg(1, vmax, bx0, bx1), Y = lin(0, ymax, by0, by1);
    p.panel(238, 16, 'b');
    p.text(bx0, 16, 'Class inventory against total vacancy content',
      { size: TYPE.body, weight: 'bold' });

    /* regime bands, named once at the top */
    var edges = [1, b.VO_onset_umol_g, b.VO_recon_complete_umol_g,
                 b.VO_cs_umol_g, vmax];
    var bands = [
      { t: '(1×1)', f: '#FFFFFF' },
      { t: 'coexist.', f: tint(C.surface, 0.10) },
      { t: 'reconstructed', f: tint(C.subsurface, 0.08) },
      { t: 'CS precip.', f: tint(C.extended, 0.14) }
    ];
    bands.forEach(function (bd, i) {
      p.rect(X(edges[i]), by1, X(edges[i + 1]) - X(edges[i]), by0 - by1,
        bd.f);
    });
    bands.forEach(function (bd, i) {
      var xm = (X(edges[i]) + X(edges[i + 1])) / 2;
      if (X(edges[i + 1]) - X(edges[i]) > 26) {
        p.text(xm, by1 - 3, bd.t, { size: TYPE.small, anchor: 'middle',
          fill: C.structure });
      }
    });
    [b.VO_onset_umol_g, b.VO_recon_complete_umol_g,
     b.VO_cs_umol_g].forEach(function (v) {
      p.line(X(v), by1, X(v), by0, C.structure, 0.6, '2 1.6');
    });
    /* clip to the drawn decade range: a log axis would otherwise carry the
       low-inventory tail out of the panel */
    var inb = sw.filter(function (r) { return r.vo >= X.d0; });
    [['surface', C.surface], ['subsurface', C.subsurface],
     ['bulk', C.bulk], ['extended', C.extended]].forEach(function (kv) {
      series(p, inb.map(function (r) { return [r.vo, r[kv[0]]]; }),
        X, Y, kv[1], 1.1);
    });
    p.line(X(fl.VO_total_umol_g), by1, X(fl.VO_total_umol_g), by0,
      C.ink, 0.8, '2.5 1.8');
    p.dot(X(fl.VO_total_umol_g), Y(fl.umol_g.subsurface), 1.9, C.ink);
    axisX(p, X, by0, [1, 3, 10, 30, 100, 220],
      'total vacancy content (µmol-O g⁻¹)',
      function (v) { return String(v); });
    axisY(p, Y, bx0, [0, 20, 40, 60, 80, 100].filter(function (v) {
      return v <= ymax; }), 'class inventory (µmol-O g⁻¹)');
    p.rect(bx0 + 8, by1 + 3, 62, 40, C.paper);
    legend(p, bx0 + 11, by1 + 12, [
      { col: C.surface, text: 'surface' },
      { col: C.subsurface, text: 'subsurface' },
      { col: C.bulk, text: 'bulk' },
      { col: C.extended, text: 'extended (CS)' }
    ]);
    return p.done();
  };

  /* d3 - grand-potential construction of the surface phase set */
  PLATES.d3 = function (D) {
    var o = D.plot.defect.omega;
    var p = plate('onehalf', 54);
    var x0 = 40, x1 = 326, y0 = 110, y1 = 26;
    var mu0 = o.curve[0].mu, mu1 = o.curve[o.curve.length - 1].mu;
    var wmin = 0, wmax = 0;
    o.curve.forEach(function (q) {
      wmin = Math.min(wmin, q.gas, q.line);
      wmax = Math.max(wmax, q.gas, q.line);
    });
    var X = lin(mu0, mu1, x0, x1), Y = lin(wmin, wmax, y0, y1);
    p.text(x0, 16, 'Surface grand potential per bridging site',
      { size: TYPE.body, weight: 'bold' });

    /* the region past the re-crossing, where the excluded branch would win */
    if (o.mu_recross < mu1) {
      p.rect(X(o.mu_recross), y1, X(mu1) - X(o.mu_recross), y0 - y1,
        tint(C.structure, 0.10));
      p.text((X(o.mu_recross) + X(mu1)) / 2, y1 + 10,
        'where the excluded branch would win',
        { size: TYPE.small, anchor: 'middle', fill: C.structure });
    }
    /* the ideal gas branch: solid up to the transition, dashed beyond */
    series(p, o.curve.filter(function (q) { return q.mu <= o.mu_t; })
      .map(function (q) { return [q.mu, q.gas]; }), X, Y, C.surface, 1.2);
    series(p, o.curve.filter(function (q) { return q.mu >= o.mu_t; })
      .map(function (q) { return [q.mu, q.gas]; }), X, Y, C.surface, 1.0,
      '3 2');
    series(p, o.curve.map(function (q) { return [q.mu, q.line]; }),
      X, Y, C.subsurface, 1.2);

    [o.mu_t, o.mu_recross].forEach(function (v) {
      if (v > mu1) return;
      p.line(X(v), y0, X(v), y1, C.structure, 0.6, '2 1.6');
    });
    p.dot(X(o.mu_t), Y(o.omega_t), 2.0, C.ink);
    p.text(X(o.mu_t) - 4, Y(o.omega_t) + 13, 'μₜ, first-order transition',
      { size: TYPE.small, anchor: 'end' });
    p.text(X(o.mu_recross) + 4, y0 - 8, 'recrossing',
      { size: TYPE.small, anchor: 'start', fill: C.structure });
    axisX(p, X, y0, [-0.5, -0.25, 0, 0.25, 0.5, 0.75],
      'μᵥ (eV)');
    axisY(p, Y, x0, [0, -0.2, -0.4, -0.6, -0.8].filter(function (v) {
      return v <= wmax && v >= wmin;
    }), 'ωₛ (eV per site)', function (v) { return v.toFixed(1); });
    p.rect(x0 + 6, y0 - 36, 104, 32, C.paper);
    legend(p, x0 + 9, y0 - 26, [
      { col: C.surface, text: '(1×1) lattice gas' },
      { col: C.subsurface, text: chem('added-row Ti2O3 line phase') },
      { col: C.surface, dash: '3 2', text: 'excluded high-coverage branch' }
    ]);
    return p.done();
  };

  /* e1 - rutile stability under three gas charges */
  PLATES.e1 = function (D) {
    var m = D.plot.equilibrium.margin, feeds = D.plot.equilibrium.feeds;
    var p = plate('double', 56);

    var ax0 = 52, ax1 = 236, ay0 = 112, ay1 = 30;
    var Ts = m.rows.map(function (r) { return r.T_C; });
    var tmin = Math.min.apply(null, Ts), tmax = Math.max.apply(null, Ts);
    var rmax = Math.ceil(Math.max.apply(null, m.rows.map(function (r) {
      return r.r_kJ; })) / 20) * 20;
    var X = lin(tmin, tmax, ax0, ax1), Y = lin(0, rmax, ay0, ay1);
    p.panel(6, 14, 'a');
    p.text(ax0, 16, 'Stability margin of ' + chem(m.phase)
      + ', ' + m.label, { size: TYPE.body, weight: 'bold' });
    p.rect(ax0, Y(rmax), ax1 - ax0, ay0 - Y(rmax), tint(C.surface, 0.07));
    p.text(ax1 - 3, ay1 + 9, 'rutile stable',
      { size: TYPE.small, anchor: 'end', fill: C.structure });
    series(p, m.rows.map(function (r) { return [r.T_C, r.r_kJ]; }),
      X, Y, C.surface, 1.2);
    m.rows.forEach(function (r) {
      p.dot(X(r.T_C), Y(r.r_kJ), 1.8, C.surface);
    });
    p.line(ax0, ay0, ax1, ay0, C.ink, 0.7);
    axisX(p, X, ay0, [500, 700, 900, 1100, 1300, 1500],
      'temperature (°C)');
    axisY(p, Y, ax0, [0, 20, 40, 60, 80],
      'KKT reduced cost (kJ per mol O)');

    var bx0 = 308, bx1 = 500;
    var Yb = lg(1e-36, 1e1, ay0, ay1);
    var Xb = lin(tmin, tmax, bx0, bx1);
    p.panel(266, 14, 'b');
    p.text(bx0, 16, 'Reduced titanium fraction',
      { size: TYPE.body, weight: 'bold' });
    var cols = { rwgs_1_1: C.gas, pure_h2: C.subsurface,
                 pure_n2: C.structure };
    feeds.forEach(function (fd) {
      var pts = fd.points.filter(function (q) { return q.reduced_pct > 0; })
        .map(function (q) { return [q.T_C, q.reduced_pct]; });
      if (!pts.length) return;
      series(p, pts, Xb, Yb, cols[fd.case], 1.2);
      pts.forEach(function (q) {
        p.dot(Xb(q[0]), Yb(q[1]), 1.6, cols[fd.case]);
      });
    });
    /* the feed that returns exact zero cannot sit on a log axis: say so */
    feeds.forEach(function (fd) {
      if (fd.points.some(function (q) { return q.reduced_pct > 0; })) return;
      p.line(bx0, ay0 - 3, bx1, ay0 - 3, cols[fd.case], 1.2, '3 2');
      p.text(bx0 + 3, ay0 - 6,
        chem(fd.label) + ': exact zero at every temperature',
        { size: TYPE.small, fill: cols[fd.case] });
    });
    axisX(p, Xb, ay0, [500, 900, 1300], 'temperature (°C)');
    axisY(p, Yb, bx0, [1e-36, 1e-30, 1e-24, 1e-18, 1e-12, 1e-6, 1],
      'reduced Ti (%)', expLabel);
    /* the exact-zero feed already carries its own annotation, so the key
       names only the two curves that can sit on a log axis */
    legend(p, bx1 - 59, ay0 - 44, feeds
      .filter(function (fd) {
        return fd.points.some(function (q) { return q.reduced_pct > 0; });
      })
      .map(function (fd) {
        return { col: cols[fd.case], text: chem(fd.label) };
      }));
    return p.done();
  };


  /* d2 - where the measured inventory sits */
  PLATES.d2 = function (D) {
    var g = D.plot.defect.geometry, fl = D.plot.defect.flagship;
    var ord = D.plot.defect.ordering;
    var p = plate('double', 64);
    var t1 = g.layer_nm, nLay = g.ss_layers;

    /* a: slab cross-section, one dot per vacancy, isotropic scale */
    var x0 = 30, y0 = 26, fw = 176, fh = 80;
    var nCols = 20, sitePx = fw / nCols;
    var pxNm = sitePx / 0.29587;                 /* c along [001], nm */
    var hPx = t1 * pxNm / 2;                     /* O sub-row spacing */
    p.panel(8, 14, 'a');
    p.text(x0, 16, 'Outer shell, to scale',
      { size: TYPE.body, weight: 'bold' });
    p.rect(x0, y0, fw, fh, C.paper, C.structure, 0.6);
    p.rect(x0, y0, fw, 9, tint(C.gas, 0.18));
    p.text(x0 + fw - 3, y0 + 6.6, 'gas',
      { size: TYPE.small, anchor: 'end', fill: C.structure });
    var surfY = y0 + 16;
    var nSub = Math.floor((y0 + fh - 5 - surfY) / hPx);
    var shellRows = 2 * nLay;
    p.rect(x0 + 1, surfY - hPx / 2, fw - 2, (shellRows + 1) * hPx,
      tint(C.subsurface, 0.09));
    function row(y, period, off) {
      p.el('line', { x1: x0 + 1, y1: y, x2: x0 + fw - 1, y2: y,
        stroke: C.structure, 'stroke-width': 1.3, 'stroke-opacity': 0.35,
        'stroke-dasharray': '1.2 ' + (period - 1.2).toFixed(2),
        'stroke-dashoffset': (-off).toFixed(2) });
    }
    row(surfY, sitePx, 0.5 * sitePx - 0.6);
    for (var i = 1; i <= nSub; i++) {
      row(surfY + i * hPx, sitePx / 2, 0.25 * sitePx - 0.6);
    }
    var rng = 1;
    function nextRand() {                     /* deterministic placement */
      rng = (rng * 1103515245 + 12345) % 2147483648;
      return rng / 2147483648;
    }
    /* the surface is on the added-row line phase here, which removes one
       bridging oxygen per two sites: the arrangement is periodic, not the
       sampled (1x1) configuration of the ordering exhibit */
    var recon = fl.surface_phase && fl.surface_phase.phase === '1x2';
    for (var j = 0; j < nCols; j++) {
      var vac = recon ? (j % 2 === 0)
        : ord.surf_vac.indexOf(j) >= 0;
      if (vac) p.dot(x0 + (j + 0.5) * sitePx, surfY, 2.1, C.surface);
    }
    var nSites = shellRows * nCols * 2;
    var want = Math.round(nSites * fl.theta.subsurface);
    var order = [];
    for (var k = 0; k < nSites; k++) order.push(k);
    for (var t = 0; t < want; t++) {
      var sw = t + Math.floor(nextRand() * (nSites - t));
      var tmp = order[t]; order[t] = order[sw]; order[sw] = tmp;
    }
    order.slice(0, want).forEach(function (ix) {
      var s2 = ix % 2, rest = (ix - s2) / 2;
      var jj = rest % nCols, rr = (rest - jj) / nCols;
      p.dot(x0 + (jj + 0.25 + 0.5 * s2) * sitePx, surfY + (1 + rr) * hPx,
        1.8, C.subsurface);
    });
    var coreTop = surfY + (shellRows + 0.5) * hPx;
    p.rect(x0 + 1, coreTop, fw - 2, y0 + fh - 1 - coreTop,
      tint(C.bulk, 0.06));
    var brX = x0 + 4;
    p.path('M' + (brX + 3) + ',' + n2(surfY - hPx / 2) + ' L' + brX + ','
      + n2(surfY - hPx / 2) + ' L' + brX + ',' + n2(coreTop) + ' L'
      + (brX + 3) + ',' + n2(coreTop), C.subsurface, 1.0);
    p.text(brX + 6, coreTop + 9, 'shell ' + g.shell_nm.toFixed(2)
      + ' nm, 1 + ' + (4 * nLay) + ' O per cell',
      { size: TYPE.small, fill: C.subsurface });
    p.line(x0 + 6, y0 + fh - 6 - pxNm, x0 + 6, y0 + fh - 6, C.ink, 1.4);
    p.text(x0 + 10, y0 + fh - 7, '1 nm', { size: TYPE.small });

    /* b: class occupancy against depth */
    var bx0 = 250, bx1 = 372, by0 = 118, by1 = 30;
    var dmax = 1000;
    var XD = lg(0.02, dmax, bx0, bx1), YT = lg(1e-4, 1, by0, by1);
    p.panel(224, 14, 'b');
    p.text(bx0, 16, 'Occupancy against depth',
      { size: TYPE.body, weight: 'bold' });
    var segs = [[0.02, t1 / 4, fl.theta.surface, C.surface],
                [t1 / 4, g.shell_nm, fl.theta.subsurface, C.subsurface],
                [g.shell_nm, dmax, fl.theta.bulk, C.bulk]];
    segs.forEach(function (q, si) {
      p.line(XD(q[0]), YT(q[2]), XD(q[1]), YT(q[2]), q[3], 1.6);
      if (si) {
        p.line(XD(q[0]), YT(segs[si - 1][2]), XD(q[0]), YT(q[2]),
          C.structure, 0.6);
      }
    });
    axisX(p, XD, by0, [0.1, 1, 10, 100, 1000], 'depth (nm)',
      function (v) { return v < 1 ? '0.1' : String(v); });
    p.line(XD(g.shell_nm), by0, XD(g.shell_nm), by1, C.structure, 0.6,
      '2 1.6');
    p.text(XD(g.shell_nm) + 3, by1 + 8, 'shell edge',
      { size: TYPE.small, fill: C.structure });
    axisY(p, YT, bx0, [1e-4, 1e-2, 1], 'vacancy occupancy', expLabel);

    /* c: the partition */
    var cx0 = 418, cx1 = 508;
    p.panel(392, 14, 'c');
    p.text(cx0, 16, 'Partition', { size: TYPE.body, weight: 'bold' });
    var vmax = Math.max(fl.umol_g.surface, fl.umol_g.subsurface,
                        fl.umol_g.bulk);
    var XV = lin(0, vmax, cx0, cx1);
    [['surface', 'surface'], ['subsurface', 'subsurface'],
     ['bulk', 'bulk']].forEach(function (kv, i) {
      var y = 34 + i * 26;
      p.rect(cx0, y, XV(fl.umol_g[kv[0]]) - cx0, 11,
        tint(C[kv[0]], 0.32), C[kv[0]], 0.7);
      p.text(cx0, y - 3, kv[1], { size: TYPE.small, fill: C[kv[0]] });
      p.text(cx0 + 2, y + 8.2, fl.umol_g[kv[0]].toFixed(1),
        { size: TYPE.tick });
    });
    axisX(p, XV, 112, [0, 25, 50], 'µmol-O g⁻¹');
    p.text(cx0, 146, 'μᵥ = ' + fl.mu_V_eV.toFixed(3) + ' eV, common to all '
      + 'classes', { size: TYPE.small });
    p.text(cx0, 155, 'shell holds ' + fl.shell_pct.toFixed(1)
      + '% of the inventory', { size: TYPE.small });
    return p.done();
  };

  /* d4 - coarse-grained reoxidation transport */
  PLATES.d4 = function (D) {
    var c = D.plot.defect.cgmc;
    var p = plate('double', 62);
    var ax0 = 40, ax1 = 244, ay0 = 116, ay1 = 30;
    var tmin = 1e-6, tmax = c.target_s;
    var X = lg(tmin, tmax, ax0, ax1), Y = lin(0, 100, ay0, ay1);
    p.panel(8, 14, 'a');
    p.text(ax0, 16, 'Recovered fraction against exposure',
      { size: TYPE.body, weight: 'bold' });
    [['neutral', C.gas, c.E_m_eV], ['fitted', C.subsurface, c.E_m_fit_eV]]
      .forEach(function (kv) {
        var rows = c[kv[0]].rows.filter(function (r) {
          return r.t_s >= tmin; });
        series(p, rows.map(function (r) { return [r.t_s, r.rec * 100]; }),
          X, Y, kv[1], 1.2);
      });
    p.dot(X(c.target_s), Y(c.target_pct), 2.2, C.ink);
    p.text(X(c.target_s) - 4, Y(c.target_pct) - 5,
      'measured ' + c.target_pct + '% at ' + c.target_s + ' s',
      { size: TYPE.small, anchor: 'end' });
    axisX(p, X, ay0, [1e-6, 1e-4, 1e-2, 1, 100], 'exposure time (s)',
      expLabel);
    axisY(p, Y, ax0, [0, 25, 50, 75, 100], 'recovered (%)');
    legend(p, ax0 + 92, ay0 - 44, [
      { col: C.gas, text: 'Eₘ = ' + c.E_m_eV.toFixed(2)
        + ' eV, neutral vacancy' },
      { col: C.subsurface, text: 'Eₘ = ' + c.E_m_fit_eV.toFixed(2)
        + ' eV, fitted' }
    ]);

    /* b: what a measured recovery implies for the barrier */
    var bx0 = 322, bx1 = 500;
    var cal = c.calibration;
    var emin = cal[0].E_m_eV, emax = cal[cal.length - 1].E_m_eV;
    var XB = lin(emin, emax, bx0, bx1), YB = lin(0, 100, ay0, ay1);
    p.panel(286, 14, 'b');
    p.text(bx0, 16, 'Recovery at ' + c.target_s + ' s against barrier',
      { size: TYPE.body, weight: 'bold' });
    series(p, cal.map(function (q) { return [q.E_m_eV, q.rec * 100]; }),
      XB, YB, C.subsurface, 1.2);
    cal.forEach(function (q) {
      p.dot(XB(q.E_m_eV), YB(q.rec * 100), 1.5, C.subsurface);
    });
    p.line(bx0, YB(c.target_pct), XB(c.E_m_fit_eV), YB(c.target_pct),
      C.ink, 0.8, '2.5 1.8');
    p.line(XB(c.E_m_fit_eV), YB(c.target_pct), XB(c.E_m_fit_eV), ay0,
      C.ink, 0.8, '2.5 1.8');
    p.dot(XB(c.E_m_fit_eV), YB(c.target_pct), 2.2, C.ink);
    p.text(XB(c.E_m_fit_eV) + 4, YB(c.target_pct) - 5,
      c.E_m_fit_eV.toFixed(2) + ' eV', { size: TYPE.small });
    p.line(XB(c.E_m_eV), ay0, XB(c.E_m_eV), ay1, C.gas, 0.9, '3 2');
    p.text(XB(c.E_m_eV) + 3, ay1 + 8, 'neutral vacancy',
      { size: TYPE.small, fill: C.gas });
    axisX(p, XB, ay0, [1, 1.5, 2, 2.5], 'migration barrier Eₘ (eV)',
      function (v) { return v.toFixed(1); });
    axisY(p, YB, bx0, [0, 25, 50, 75, 100], 'recovered (%)');
    return p.done();
  };

  /* e2 - active-set selection */
  PLATES.e2 = function (D) {
    var m = D.plot.equilibrium.method;
    var p = plate('onehalf', 62);
    var x0 = 60, x1 = 316, y0 = 104, y1 = 34;
    var costs = m.costs.slice().sort(function (a, b) {
      return a.r_kJ - b.r_kJ; });
    var vmax = Math.max.apply(null, costs.map(function (q) {
      return q.r_kJ; }));
    var X = lin(0, vmax * 1.12, x0, x1);
    var bh = (y0 - y1) / costs.length * 0.62;
    p.text(x0 - 22, 16, 'Optimality of the selected assemblage at '
      + m.T_C + ' °C', { size: TYPE.body, weight: 'bold' });
    p.text(x0 - 22, 26, chem(m.winner.join(' + ')) + ' survives; every '
      + 'excluded phase costs more', { size: TYPE.small,
        fill: C.structure });
    costs.forEach(function (q, i) {
      var y = y1 + (i + 0.5) * (y0 - y1) / costs.length - bh / 2;
      p.rect(x0, y, X(q.r_kJ) - x0, bh, tint(C.surface, 0.30), C.surface,
        0.6);
      p.text(x0 - 4, y + bh - 1, chem(q.phase),
        { size: TYPE.small, anchor: 'end' });
      p.text(X(q.r_kJ) + 3, y + bh - 1, q.r_kJ.toFixed(1),
        { size: TYPE.small });
    });
    p.line(x0, y1, x0, y0, C.ink, 0.9);
    axisX(p, X, y0, niceTicks(0, vmax * 1.12, 4),
      'KKT reduced cost (kJ per mol O)',
      function (v) { return String(Math.round(v)); });
    var nOk = m.candidates.filter(function (c) {
      return c.feasible && c.kkt_ok; }).length;
    var nFeas = m.candidates.filter(function (c) {
      return c.feasible; }).length;
    p.text(x0 - 22, y0 + 30, m.candidates.length + ' candidate sets '
      + 'enumerated, ' + nFeas + ' feasible, ' + nOk + ' optimal',
      { size: TYPE.small });
    p.text(x0 - 22, y0 + 39, 'excluded phases are returned as exact zero',
      { size: TYPE.small, fill: C.structure });
    return p.done();
  };

  /* s1 - declared shell thickness */
  PLATES.s1 = function (D) {
    var L = D.plot.defect.shell_ladder;
    var p = plate('onehalf', 60);
    var x0 = 44, x1 = 246, y0 = 110, y1 = 28;
    var tot = D.plot.defect.flagship.VO_total_umol_g;
    var Y = lin(0, tot, y0, y1);
    var bw = (x1 - x0) / L.length * 0.5;
    p.text(x0, 16, 'Shell thickness against the shell/bulk split',
      { size: TYPE.body, weight: 'bold' });
    L.forEach(function (r, i) {
      var xc = x0 + (i + 0.5) * (x1 - x0) / L.length;
      var acc = 0;
      [['surface', C.surface], ['subsurface', C.subsurface],
       ['bulk', C.bulk]].forEach(function (kv) {
        var v = r[kv[0]];
        p.rect(xc - bw / 2, Y(acc + v), bw, Y(acc) - Y(acc + v),
          tint(kv[1], 0.34), kv[1], 0.6);
        acc += v;
      });
      p.text(xc, Y(acc) - 4, r.shell_pct.toFixed(0) + '%',
        { size: TYPE.small, anchor: 'middle' });
      p.text(xc, y0 + 9.5, String(r.n_lay),
        { size: TYPE.tick, anchor: 'middle' });
      p.text(xc, y0 + 18, r.shell_nm.toFixed(2),
        { size: TYPE.small, anchor: 'middle', fill: C.structure });
    });
    p.line(x0, y0, x1, y0, C.ink, 0.7);
    p.text((x0 + x1) / 2, y0 + 29, 'subsurface layers / shell depth (nm)',
      { size: TYPE.body, anchor: 'middle' });
    axisY(p, Y, x0, [0, 25, 50, 75, 95].filter(function (v) {
      return v <= tot; }), 'class inventory (µmol-O g⁻¹)');
    legend(p, x1 + 12, y1 + 10, [
      { col: C.surface, swatch: true, text: 'surface' },
      { col: C.subsurface, swatch: true, text: 'subsurface' },
      { col: C.bulk, swatch: true, text: 'bulk' }
    ]);
    return p.done();
  };

  /* s2 - surface vacancy arrangement */
  PLATES.s2 = function (D) {
    var o = D.plot.defect.ordering;
    var p = plate('onehalf', 60);
    var x0 = 40, y0 = 26;
    var cell = 5.2, nR = o.rows, nC = o.row_sites;
    p.panel(8, 14, 'a');
    p.text(x0, 16, 'Sampled bridging row',
      { size: TYPE.body, weight: 'bold' });
    var occ = {};
    o.surf_vac.forEach(function (v) { occ[v] = 1; });
    for (var r = 0; r < nR; r++) {
      for (var k = 0; k < nC; k++) {
        var cx = x0 + k * cell, cy = y0 + r * cell;
        if (occ[r * nC + k]) p.dot(cx, cy, 1.8, C.surface);
        else p.dot(cx, cy, 0.7, C.structure);
      }
    }
    p.text(x0, y0 + nR * cell + 8,
      'θₛ = ' + (o.theta_s * 100).toFixed(1) + '%, '
      + o.rows + ' × ' + o.row_sites + ' sites',
      { size: TYPE.small });

    var bx0 = 214, bx1 = 330, by0 = 104, by1 = 30;
    var h = o.hist || [];
    var n = Math.min(h.length, 10);
    var hmax = 0;
    for (var i = 1; i <= n; i++) hmax = Math.max(hmax, h[i - 1] || 0);
    var Y = lin(0, hmax, by0, by1);
    var bw2 = (bx1 - bx0) / n * 0.66;
    p.panel(190, 14, 'b');
    p.text(bx0, 16, 'Separation along the row',
      { size: TYPE.body, weight: 'bold' });
    for (var i2 = 1; i2 <= n; i2++) {
      var xc2 = bx0 + (i2 - 0.5) * (bx1 - bx0) / n;
      var v2 = h[i2 - 1] || 0;
      p.rect(xc2 - bw2 / 2, Y(v2), bw2, by0 - Y(v2),
        tint(C.surface, 0.34), C.surface, 0.6);
      if (i2 % 2 === 1) {
        p.text(xc2, by0 + 9.5, String(i2),
          { size: TYPE.tick, anchor: 'middle' });
      }
    }
    p.line(bx0, by0, bx1, by0, C.ink, 0.7);
    axisY(p, Y, bx0, [0, hmax], 'pairs sampled',
      function (v) { return String(Math.round(v)); });
    p.text((bx0 + bx1) / 2, by0 + 21, 'gap (bridging sites)',
      { size: TYPE.body, anchor: 'middle' });
    p.text(bx0, by1 - 4, 'modal gap ' + o.spacing.modal_gap_sites
      + ' sites, P(gap ≤ 2) = ' + o.spacing.P_gap_le_2.toFixed(2),
      { size: TYPE.small });
    return p.done();
  };

  /* s3 - geometry sensitivity */
  PLATES.s3 = function (D) {
    var b = D.plot.defect.boundaries, bb = D.plot.defect.bet;
    var p = plate('single', 52);
    var x0 = 36, x1 = 236, y0 = 92, y1 = 30;
    var X = lg(1, 400, x0, x1);
    p.text(x0, 15, 'Boundary ladder, two geometries',
      { size: TYPE.body, weight: 'bold' });
    var keys = [['VO_onset_umol_g', 'onset'],
                ['VO_recon_complete_umol_g', 'complete'],
                ['VO_cs_umol_g', 'CS ceiling']];
    keys.forEach(function (kv, i) {
      var y = 30 + i * 19;
      p.line(X(b[kv[0]]), y, X(bb.boundaries[kv[0]]), y, C.structure, 0.8);
      p.dot(X(b[kv[0]]), y, 2.2, C.surface);
      p.dot(X(bb.boundaries[kv[0]]), y, 2.2, C.paper, C.subsurface);
      p.text(x0, y - 5, kv[1], { size: TYPE.small, fill: C.structure });
    });
    axisX(p, X, y0, [1, 10, 100, 400], 'µmol-O g⁻¹',
      function (v) { return String(v); });
    legend(p, x0 + 4, y1 + 30, [
      { col: C.surface, swatch: true,
        text: 'sphere, ' + D.slots['defect.geometry.d_um'] + ' µm' },
      { col: C.subsurface, swatch: true,
        text: 'BET, ' + bb.area_m2_g.toFixed(2) + ' m² g⁻¹' }
    ]);
    p.text(x0, y0 + 22, 'agreement to '
      + bb.max_shift_pct.toFixed(1) + '% on every boundary',
      { size: TYPE.small });
    return p.done();
  };

  var SHORT_CASE = { rwgs_1_1: 'CO2/H2 = 1:1', h2_rich: 'CO2/H2 = 1:3',
    product_mix: 'product mix', pure_h2: 'pure H2', pure_n2: 'sealed N2' };

  /* s4 - validity of the partition */
  PLATES.s4 = function (D) {
    var v = D.plot.equilibrium.validity;
    var p = plate('onehalf', 58);
    var x0 = 108, y0 = 30, cw = 34, ch = 15;
    var temps = v.temps;
    p.text(20, 16, 'Where the partition is valid',
      { size: TYPE.body, weight: 'bold' });
    var fills = { stable: tint(C.surface, 0.30),
                  boundary: tint(C.extended, 0.36),
                  outside: tint(C.structure, 0.30) };
    temps.forEach(function (t, j) {
      p.text(x0 + (j + 0.5) * cw, y0 - 4, String(t),
        { size: TYPE.small, anchor: 'middle' });
    });
    p.text(x0 + temps.length * cw / 2, y0 - 14, 'temperature (°C)',
      { size: TYPE.small, anchor: 'middle', fill: C.structure });
    v.cases.forEach(function (c, i) {
      var y = y0 + i * ch;
      p.text(x0 - 4, y + 10.5, chem(SHORT_CASE[c.case] || c.case),
        { size: TYPE.small, anchor: 'end' });
      c.tiers.forEach(function (q) {
        var j = temps.indexOf(q.T_C);
        p.rect(x0 + j * cw + 1, y + 2, cw - 2, ch - 4, fills[q.tier],
          C.paper, 0.8);
      });
    });
    legend(p, x0, y0 + v.cases.length * ch + 16, [
      { col: fills.stable, swatch: true, text: 'rutile alone, valid' },
      { col: fills.boundary, swatch: true,
        text: 'rutile with a reduced phase, boundary' },
      { col: fills.outside, swatch: true,
        text: 'no rutile, outside the model' }
    ]);
    return p.done();
  };

  /* s5 - numerical verification of the defect stack */
  PLATES.s5 = function (D) {
    var V = D.plot.verification;
    var p = plate('double', 62);
    var ax0 = 52, ax1 = 236, ay0 = 112, ay1 = 30;
    var X = lin(350, 1050, ax0, ax1), Y = lg(1e-16, 1e-8, ay0, ay1);
    p.panel(8, 14, 'a');
    p.text(ax0, 16, 'Partition against the reference',
      { size: TYPE.body, weight: 'bold' });
    p.line(ax0, Y(1e-9), ax1, Y(1e-9), C.extended, 1.0, '3 2');
    p.text(ax1 - 3, Y(1e-9) - 4, 'enforced tolerance',
      { size: TYPE.small, anchor: 'end', fill: C.extended });
    V.grid_dev.forEach(function (q) {
      p.dot(X(q.T_C), Y(Math.max(q.dev_eV, 1e-16)), 1.8, C.surface);
    });
    axisX(p, X, ay0, [400, 600, 800, 1000], 'temperature (°C)');
    axisY(p, Y, ax0, [1e-16, 1e-12, 1e-8], 'deviation in μᵥ (eV)',
      expLabel);
    p.text(ax0 + 4, ay0 - 5, 'points at the floor of double precision',
      { size: TYPE.small, fill: C.structure });

    var bx0 = 330, bx1 = 470;
    var XB = lg(1e-16, 1e-1, bx0, bx1);
    p.panel(276, 14, 'b');
    p.text(bx0, 16, 'Deviation ledger',
      { size: TYPE.body, weight: 'bold' });
    V.layers.forEach(function (L, i) {
      var y = 34 + i * 22;
      p.text(bx0, y - 3, L.layer, { size: TYPE.small });
      p.rect(bx0, y, XB(Math.max(L.measured, 1e-16)) - bx0, 6,
        tint(C.surface, 0.34), C.surface, 0.6);
      p.line(XB(L.tolerance), y - 1.5, XB(L.tolerance), y + 7.5,
        C.extended, 1.1);
    });
    axisX(p, XB, ay0, [1e-15, 1e-10, 1e-5], 'measured against tolerance',
      expLabel);
    p.text(bx0, ay0 + 32, V.n_tests + ' tests, none skipped. '
      + 'Vertical marks are the enforced tolerances.',
      { size: TYPE.small });
    return p.done();
  };

  /* s6 - equilibrium verification */
  PLATES.s6 = function (D) {
    var V = D.plot.equilibrium.verification;
    var p = plate('double', 62);
    var ax0 = 48, ax1 = 240, ay0 = 112, ay1 = 30;
    var Ts = V.traces.map(function (q) { return q.T_C; });
    var tmin = Math.min.apply(null, Ts), tmax = Math.max.apply(null, Ts);
    var lo = Math.floor(Math.min.apply(null,
      V.traces.map(function (q) { return q.log10_mol; })) / 10) * 10;
    var X = lin(tmin, tmax, ax0, ax1), Y = lin(lo, 0, ay0, ay1);
    p.panel(8, 14, 'a');
    p.text(ax0, 16, 'Trace amounts, sealed inert charge',
      { size: TYPE.body, weight: 'bold' });
    var byS = {};
    V.traces.forEach(function (q) {
      (byS[q.species] = byS[q.species] || []).push(q); });
    var cols = [C.surface, C.subsurface, C.gas, C.bulk];
    Object.keys(byS).sort().forEach(function (sp, i) {
      var col = cols[i % cols.length];
      var pts = byS[sp].sort(function (a, b) { return a.T_C - b.T_C; });
      series(p, pts.map(function (q) { return [q.T_C, q.log10_mol]; }),
        X, Y, col, 1.1);
      pts.forEach(function (q) { p.dot(X(q.T_C), Y(q.log10_mol), 1.6, col); });
      /* the traces converge, so the keys are stacked rather than placed
         on the curve ends */
      p.text(ax1 - 2, ay1 + 8 + i * 9, chem(sp),
        { size: TYPE.small, anchor: 'end', fill: col });
    });
    axisX(p, X, ay0, [500, 900, 1300], 'temperature (°C)');
    axisY(p, Y, ax0, [lo, lo / 2, 0], 'log₁₀ (amount / mol)',
      function (v) { return String(Math.round(v)); });

    var bx0 = 330, bx1 = 500;
    var XB = lin(tmin, tmax, bx0, bx1);
    var rmin = -100;
    var YB = lg(1e-100, 1e-60, ay0, ay1);
    p.panel(276, 14, 'b');
    p.text(bx0, 16, 'Elemental closure',
      { size: TYPE.body, weight: 'bold' });
    V.residuals.forEach(function (q) {
      p.dot(XB(q.T_C), YB(Math.max(q.worst_abs_mol, 1e-100)), 1.9, C.gas);
    });
    axisX(p, XB, ay0, [500, 900, 1300], 'temperature (°C)');
    axisY(p, YB, bx0, [1e-100, 1e-80, 1e-60],
      'worst residual (mol)', expLabel);
    p.text(bx0, ay1 - 4, 'zero residuals are drawn at the axis floor',
      { size: TYPE.small, fill: C.structure });
    return p.done();
  };

  function render(id, data) {
    if (!PLATES[id]) throw new Error('no plate ' + id);
    return PLATES[id](data);
  }

  return { plate: plate, render: render, PLATES: PLATES, C: C, tint: tint,
           chem: chem, niceTicks: niceTicks,
           TYPE: TYPE, WIDTH_MM: WIDTH_MM, PT_PER_MM: PT_PER_MM,
           lin: lin, lg: lg, decades: decades, axisX: axisX, axisY: axisY,
           series: series, legend: legend, expLabel: expLabel };
});
