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

  function render(id, data) {
    if (!PLATES[id]) throw new Error('no plate ' + id);
    return PLATES[id](data);
  }

  return { plate: plate, render: render, PLATES: PLATES, C: C, tint: tint,
           chem: chem,
           TYPE: TYPE, WIDTH_MM: WIDTH_MM, PT_PER_MM: PT_PER_MM,
           lin: lin, lg: lg, decades: decades, axisX: axisX, axisY: axisY,
           series: series, legend: legend, expLabel: expLabel };
});
