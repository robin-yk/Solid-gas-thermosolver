/* The figure kit: one drawing style, two products.

   Everything quantitative in this repository is drawn through here - the
   workspace figures on the page, at
   a square 3.5 in. The coordinate system is typographic points in both,
   so a stroke width written here is the stroke width that reaches the
   page, and the browser preview and the print export share one geometry.

   House style, in one place and nowhere else: a closed rectangular axis
   box on all four sides, ticks pointing inward, tick labels on the bottom
   and left only, no grid, white ground, one stroke weight for data and
   another for spines, ticks and guides. Colour carries physical role and
   nothing else.

   Runs unchanged under node, which is how the gates measure type size and
   width without a browser. */

(function (root, factory) {
  'use strict';
  var api = factory();
  if (typeof module !== 'undefined' && module.exports) module.exports = api;
  root.FigKit = api;
})(typeof self !== 'undefined' ? self : this, function () {
  'use strict';

  var PT_PER_MM = 72 / 25.4;
  var WIDTH_MM = { single: 89, onehalf: 120, double: 183 };
  var MAX_H_MM = 247;
  var FONT = 'Arial, Helvetica, sans-serif';
  /* Two type scales, both in points at print size. A plate is a column
     of a printed page and sits at the journal floor; a workspace figure
     is a 3.5 in square read on its own and carries the larger set the
     specification asks for. Offsets below are derived from these, so a
     figure never needs a hand-placed label. */
  var TYPE = { panel: 8, body: 7, tick: 6.5, small: 5.5 };
  var FIGTYPE = { panel: 9, body: 9, tick: 8, small: 8 };
  /* Print style, in points - the plate coordinate system is points, so
     these are the numbers that reach the page. Spines and ticks carry one
     weight, data curves another, guides a third; nothing sets a stroke
     width by hand. */
  var LW = { axis: 0.9, curve: 2.2, guide: 0.8, hair: 0.6 };
  var TICK = { major: 4, minor: 2 };
  var MARK = 2.5;                 /* radius, so a 5 pt marker */
  /* One hue per physical role. Extended defects are the one neutral
     role, so structure - axes, boundaries, anything held fixed - is a
     lighter grey than they are and never competes with them. */
  var C = {
    surface: '#0072B2', subsurface: '#D55E00', bulk: '#6A51A3',
    extended: '#777777', gas: '#009E73', structure: '#AAAAAA',
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

  /* A figure in points. The SVG declares its print size in real units -
     millimetres for a plate, inches for a workspace figure - and the
     viewBox is the same number in points, so one user unit is one point
     and nothing rescales between the screen and the page. */
  function figure(W, H, o) {
    o = o || {};
    var parts = [];
    var api = {
      kind: o.kind || 'figure',
      W: W, H: H,
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
      path: function (d, col, w, dash, fill, cls) {
        return api.el('path', { 'class': cls, d: d, fill: fill || 'none',
          stroke: col || C.ink, 'stroke-width': w || LW.curve,
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
          + 'width="' + o.widthAttr + '" height="' + o.heightAttr + '" '
          + 'viewBox="0 0 ' + n2(api.W) + ' ' + n2(api.H) + '" '
          + (o.dataAttrs ? o.dataAttrs + ' ' : '')
          + 'data-kind="' + api.kind + '">'
          + '<rect x="0" y="0" width="' + n2(api.W) + '" height="'
          + n2(api.H) + '" fill="' + C.paper + '"/>'
          + parts.join('') + '</svg>';
      }
    };
    return api;
  }

  /* The manuscript plate: a Nature column width in millimetres. */
  function plate(widthKey, heightMm) {
    var wmm = WIDTH_MM[widthKey];
    if (!wmm) throw new Error('unknown plate width: ' + widthKey);
    if (heightMm > MAX_H_MM) throw new Error('plate taller than the page');
    var api = figure(wmm * PT_PER_MM, heightMm * PT_PER_MM, {
      widthAttr: wmm + 'mm', heightAttr: heightMm + 'mm',
      dataAttrs: 'data-width-mm="' + wmm + '" data-height-mm="'
        + heightMm + '"',
      kind: 'plate'
    });
    api.widthKey = widthKey;
    api.wmm = wmm;
    api.hmm = heightMm;
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
  /* Every decade inside the range, endpoints included. log10 of an exact
     decade lands a few ulp short of the integer, so a bound written as
     1e3 used to lose its own tick; the tolerance and the literal
     construction keep an axis from dropping the label at its own edge. */
  function decades(d0, d1) {
    var out = [];
    var e1 = Math.floor(log10(d1) + 1e-9);
    for (var e = Math.ceil(log10(d0) - 1e-9); e <= e1; e++) {
      out.push(Number('1e' + e));
    }
    return out;
  }

  /* The enclosing decades of a value, so a log axis brackets its data
     without a hand-set range: a figure that is recomputed from the user's
     own numbers must not be able to draw a point outside its own frame. */
  function decadeFloor(v) { return Number('1e' + Math.floor(log10(v))); }
  function decadeCeil(v) { return Number('1e' + Math.ceil(log10(v))); }

  function expLabel(v) {
    var e = Math.round(log10(v));
    if (e >= 0 && e <= 3) return String(Math.round(v));
    return '10' + String(e).replace(/-/g, '−')
      .split('').map(function (ch) {
        return '⁰¹²³⁴⁵⁶⁷⁸⁹'
          [Number(ch)] || (ch === '−' ? '⁻' : ch);
      }).join('');
  }

  /* Always the exponent form. An axis spanning many decades reads as one
     ladder only if every rung is written the same way; expLabel's plain
     integers are for the two- or three-decade axes. */
  function powLabel(v) {
    var e = Math.round(log10(v));
    return '10' + String(e).replace(/-/g, '−').split('')
      .map(function (ch) {
        return '⁰¹²³⁴⁵⁶⁷⁸⁹'
          [Number(ch)] || (ch === '−' ? '⁻' : ch);
      }).join('');
  }

  /* ------------------------------------------------------------- axes */

  /* Full rectangular axis box: four spines, no grid. Call it after the
     panel fills and before the axes, so nothing paints over a spine. The
     box is stashed on the plate so axisX and axisY mirror their ticks onto
     the opposite spine without being told where it is. Either argument may
     be a scale or an explicit [lo, hi] pixel pair. */
  function frame(p, X, Y) {
    var xr = X.r0 == null ? X : [X.r0, X.r1];
    var yr = Y.r0 == null ? Y : [Y.r0, Y.r1];
    var x0 = Math.min(xr[0], xr[1]), x1 = Math.max(xr[0], xr[1]);
    var y0 = Math.min(yr[0], yr[1]), y1 = Math.max(yr[0], yr[1]);
    p._box = { x0: x0, y0: y0, x1: x1, y1: y1 };
    p.rect(x0, y0, x1 - x0, y1 - y0, 'none', C.ink, LW.axis);
    return p;
  }

  /* Minor ticks at the pixel midpoint of each major interval: the
     arithmetic midpoint on a linear axis, the geometric one on a log
     axis, which is where a minor tick belongs on either. */
  function minorsOf(px, lo, hi) {
    var s = px.slice().sort(function (a, b) { return a - b; });
    var out = [];
    for (var i = 1; i < s.length; i++) {
      var m = 0.5 * (s[i - 1] + s[i]);
      if (m > lo + 0.5 && m < hi - 0.5) out.push(m);
    }
    return out;
  }

  /* A tick outside the scale's own domain is drawn off the axis, which
     is always a bug rather than an intention. */
  function inDomain(sc, t) {
    if (sc.d0 == null) return true;
    var lo = Math.min(sc.d0, sc.d1), hi = Math.max(sc.d0, sc.d1);
    var pad = 1e-9 * (Math.abs(hi) + Math.abs(lo) + 1);
    return t >= lo - pad && t <= hi + pad;
  }

  function axisX(p, sc, y, ticks, label, fmt) {
    var T = p.type || TYPE;
    var bx = p._box;
    ticks = ticks.filter(function (t) { return inDomain(sc, t); });
    var top = bx ? bx.y0 : null;
    if (!bx) p.line(sc.r0, y, sc.r1, y, C.ink, LW.axis);
    var px = [];
    ticks.forEach(function (t) {
      var x = sc(t);
      px.push(x);
      p.line(x, y, x, y - TICK.major, C.ink, LW.axis);
      if (top != null) p.line(x, top, x, top + TICK.major, C.ink, LW.axis);
      /* the outermost label is centred on a tick that sits on the frame,
         so half of it hangs past the figure edge and gets clipped. Nudge
         only those back inside; every interior label is untouched. */
      var s = fmt ? fmt(t) : String(t);
      var half = 0.28 * T.tick * s.length;
      var tx = Math.min(Math.max(x, half + 1), p.W - half - 1);
      p.text(tx, y + T.tick + 3, s, { size: T.tick, anchor: 'middle' });
    });
    minorsOf(px, bx ? bx.x0 : sc.r0, bx ? bx.x1 : sc.r1)
      .forEach(function (x) {
        p.line(x, y, x, y - TICK.minor, C.ink, LW.axis);
        if (top != null) {
          p.line(x, top, x, top + TICK.minor, C.ink, LW.axis);
        }
      });
    if (label) {
      p.text((sc.r0 + sc.r1) / 2, y + T.tick + T.body + 7.5, label,
        { size: T.body, anchor: 'middle' });
    }
    return p;
  }
  function axisY(p, sc, x, ticks, label, fmt) {
    var T = p.type || TYPE;
    var bx = p._box;
    ticks = ticks.filter(function (t) { return inDomain(sc, t); });
    var right = bx ? bx.x1 : null;
    if (!bx) p.line(x, sc.r0, x, sc.r1, C.ink, LW.axis);
    var px = [];
    ticks.forEach(function (t) {
      var y = sc(t);
      px.push(y);
      p.line(x, y, x + TICK.major, y, C.ink, LW.axis);
      if (right != null) {
        p.line(right - TICK.major, y, right, y, C.ink, LW.axis);
      }
      p.text(x - 4, y + T.tick / 2 - 0.95, fmt ? fmt(t) : String(t),
        { size: T.tick, anchor: 'end' });
    });
    minorsOf(px, bx ? bx.y0 : sc.r1, bx ? bx.y1 : sc.r0)
      .forEach(function (y) {
        p.line(x, y, x + TICK.minor, y, C.ink, LW.axis);
        if (right != null) {
          p.line(right - TICK.minor, y, right, y, C.ink, LW.axis);
        }
      });
    if (label) {
      p.text(x - (2 * T.tick + T.body + 6), (sc.r0 + sc.r1) / 2, label,
        { size: T.body, anchor: 'middle', rotate: -90 });
    }
    return p;
  }
  /* Open markers, three shapes. A scatter that carries three populations
     has to survive greyscale and a printed 3.5 in., so the series are told
     apart by shape first and colour second. Every shape is drawn hollow on
     the paper ground and sized to the same visual area as a circle of the
     same r, so no series looks heavier than another. */
  function marker(p, shape, x, y, r, col, filled) {
    var fill = filled ? col : C.paper;
    if (shape === 'square') {
      var s = 1.772 * r;                       /* equal area to the circle */
      return p.rect(x - s / 2, y - s / 2, s, s, fill, col, LW.axis);
    }
    if (shape === 'triangle') {
      var a = 2.694 * r, h = 0.866 * a;        /* equal area, centroid at y */
      return p.path('M' + n2(x) + ',' + n2(y - 2 * h / 3)
        + 'L' + n2(x + a / 2) + ',' + n2(y + h / 3)
        + 'L' + n2(x - a / 2) + ',' + n2(y + h / 3) + 'Z',
        col, LW.axis, null, fill, 'marker');
    }
    return p.el('circle', { cx: x, cy: y, r: r, fill: fill,
      stroke: col, 'stroke-width': LW.axis });
  }

  /* A guide through computed points: Catmull-Rom in pixel space, emitted
     as cubic Beziers. It interpolates every point exactly, so the guide
     cannot claim a value the model did not produce - it only smooths the
     path between them, which is what a dashed guide is for. */
  function smooth(pts) {
    if (pts.length < 3) {
      return pts.map(function (q, i) {
        return (i ? 'L' : 'M') + n2(q[0]) + ',' + n2(q[1]);
      }).join('');
    }
    var d = 'M' + n2(pts[0][0]) + ',' + n2(pts[0][1]);
    for (var i = 0; i < pts.length - 1; i++) {
      var p0 = pts[i === 0 ? 0 : i - 1], p1 = pts[i], p2 = pts[i + 1];
      var p3 = pts[Math.min(i + 2, pts.length - 1)];
      var c1x = p1[0] + (p2[0] - p0[0]) / 6, c1y = p1[1] + (p2[1] - p0[1]) / 6;
      var c2x = p2[0] - (p3[0] - p1[0]) / 6, c2y = p2[1] - (p3[1] - p1[1]) / 6;
      d += 'C' + n2(c1x) + ',' + n2(c1y) + ' ' + n2(c2x) + ',' + n2(c2y)
        + ' ' + n2(p2[0]) + ',' + n2(p2[1]);
    }
    return d;
  }

  function series(p, pts, xs, ys, col, w, dash) {
    var d = '', started = false;
    pts.forEach(function (q) {
      if (q[1] == null || !isFinite(q[1])) { started = false; return; }
      d += (started ? 'L' : 'M') + n2(xs(q[0])) + ',' + n2(ys(q[1]));
      started = true;
    });
    if (d) p.path(d, col, w || LW.curve, dash, null, 'curve');
    return p;
  }
  function legend(p, x, y, items, gap) {
    var T = p.type || TYPE;
    gap = gap || T.small + 4;
    items.forEach(function (it, i) {
      var yy = y + i * gap;
      if (it.marker) {
        marker(p, it.marker, x + 4.5, yy - 2.4, MARK + 0.5, it.col, it.filled);
      } else if (it.dash) {
        p.line(x, yy - 2.2, x + 9, yy - 2.2, it.col, 1.0, it.dash);
      } else if (it.swatch) {
        /* edge is for a filled shape whose outline carries meaning: a
           tinted bar with an ink stroke reads as black in the legend
           unless the key is drawn the same way as the bar. */
        p.rect(x, yy - 5.2, 7, 6, it.col, it.edge, it.edge ? LW.axis : null);
      } else {
        p.line(x, yy - 2.2, x + 9, yy - 2.2, it.col, 1.4);
      }
      p.text(x + 12, yy, it.text, { size: T.small });
    });
    return p;
  }

  /* --------------------------------------------------------- products */

  /* The workspace figure: exactly 3.5 in square, one user unit to the
     point, so the same drawing code serves the preview and the 600 dpi
     export without a change of coordinates. */
  var SQUARE_IN = 3.5;
  var SQUARE_PT = SQUARE_IN * 72;              /* 252 */

  function square(o) {
    o = o || {};
    var wIn = o.wide ? 2 * SQUARE_IN : SQUARE_IN;
    var f = figure(wIn * 72, SQUARE_PT, {
      widthAttr: wIn + 'in', heightAttr: SQUARE_IN + 'in',
      kind: o.wide ? 'wide' : 'square'
    });
    f.type = FIGTYPE;
    /* the plotting panel is square and the margins are what is left */
    f.pane = { x0: 44, y0: 10, x1: 44 + 200, y1: 210 };
    return f;
  }

  /* ---------------------------------------------------------- export */

  /* A 600 dpi raster of a square figure is 2100 px on a side. The browser
     draws the SVG into a canvas at that size; toBlob gives a PNG with no
     resolution chunk, so pngWithDPI splices a pHYs in - 600 dpi is 23622
     pixels per metre - and the file then opens at print size rather than
     at 2100 px of nothing. */
  var DPI = 600;
  var PX_PER_M = Math.round(DPI / 0.0254);     /* 23622 */

  var CRC = (function () {
    var t = [], c, n, k;
    for (n = 0; n < 256; n++) {
      c = n;
      for (k = 0; k < 8; k++) c = (c & 1) ? (0xEDB88320 ^ (c >>> 1)) : (c >>> 1);
      t[n] = c >>> 0;
    }
    return function (bytes, from, to) {
      var c2 = 0xFFFFFFFF;
      for (var i = from; i < to; i++) {
        c2 = t[(c2 ^ bytes[i]) & 255] ^ (c2 >>> 8);
      }
      return (c2 ^ 0xFFFFFFFF) >>> 0;
    };
  }());

  function u32(v) {
    return [(v >>> 24) & 255, (v >>> 16) & 255, (v >>> 8) & 255, v & 255];
  }

  /* Insert a pHYs chunk straight after IHDR. Anything else in the file is
     copied through untouched, so this cannot disturb the image data. */
  function pngWithDPI(buf, pxPerM) {
    var b = new Uint8Array(buf);
    var at = 8 + 4 + 4 + 25 - 12 + 12;         /* signature + full IHDR */
    at = 8 + (8 + 13 + 4);
    var body = [0, 0, 0, 9, 0x70, 0x48, 0x59, 0x73]
      .concat(u32(pxPerM)).concat(u32(pxPerM)).concat([1]);
    var chunk = new Uint8Array(body.length + 4);
    chunk.set(body, 0);
    chunk.set(u32(CRC(chunk, 4, body.length)), body.length);
    var out = new Uint8Array(b.length + chunk.length);
    out.set(b.subarray(0, at), 0);
    out.set(chunk, at);
    out.set(b.subarray(at), at + chunk.length);
    return out;
  }

  /* Rasterise an SVG string at a fixed pixel size. Browser only. */
  function toPNG(svg, pxW, pxH) {
    return new Promise(function (resolve, reject) {
      var img = new Image();
      var url = 'data:image/svg+xml;charset=utf-8,' + encodeURIComponent(svg);
      img.onload = function () {
        var cv = document.createElement('canvas');
        cv.width = pxW;
        cv.height = pxH;
        var cx = cv.getContext('2d');
        cx.fillStyle = C.paper;
        cx.fillRect(0, 0, pxW, pxH);
        cx.drawImage(img, 0, 0, pxW, pxH);
        cv.toBlob(function (blob) {
          if (!blob) { reject(new Error('canvas produced no PNG')); return; }
          blob.arrayBuffer().then(function (ab) {
            resolve(new Blob([pngWithDPI(ab, PX_PER_M)],
              { type: 'image/png' }));
          }).catch(reject);
        }, 'image/png');
      };
      img.onerror = function () { reject(new Error('SVG did not rasterise')); };
      img.src = url;
    });
  }

  /* Pixel size of a figure at 600 dpi, from the inches it declares. */
  function pxAt600(svg) {
    var w = /width="([\d.]+)in"/.exec(svg);
    var h = /height="([\d.]+)in"/.exec(svg);
    if (!w || !h) throw new Error('figure does not declare its print size');
    return [Math.round(parseFloat(w[1]) * DPI),
            Math.round(parseFloat(h[1]) * DPI)];
  }

  function saveBlob(blob, name) {
    var a = document.createElement('a');
    var url = URL.createObjectURL(blob);
    a.href = url;
    a.download = name;
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
    setTimeout(function () { URL.revokeObjectURL(url); }, 4000);
  }

  function downloadSVG(svg, name) {
    saveBlob(new Blob([svg], { type: 'image/svg+xml;charset=utf-8' }),
      name + '.svg');
  }

  function downloadPNG(svg, name) {
    var px = pxAt600(svg);
    return toPNG(svg, px[0], px[1]).then(function (blob) {
      saveBlob(blob, name + '.png');
    });
  }

  return {
    PT_PER_MM: PT_PER_MM, WIDTH_MM: WIDTH_MM, MAX_H_MM: MAX_H_MM,
    FONT: FONT, TYPE: TYPE, FIGTYPE: FIGTYPE, LW: LW, TICK: TICK,
    MARK: MARK, C: C,
    DPI: DPI, PX_PER_M: PX_PER_M, SQUARE_IN: SQUARE_IN, SQUARE_PT: SQUARE_PT,
    tint: tint, esc: esc, n2: n2, chem: chem, expLabel: expLabel,
    powLabel: powLabel, decadeFloor: decadeFloor, decadeCeil: decadeCeil,
    figure: figure, plate: plate, square: square,
    lin: lin, lg: lg, niceTicks: niceTicks, decades: decades,
    frame: frame, axisX: axisX, axisY: axisY, series: series,
    marker: marker, smooth: smooth,
    inDomain: inDomain,
    legend: legend, minorsOf: minorsOf,
    toPNG: toPNG, pxAt600: pxAt600, pngWithDPI: pngWithDPI,
    downloadSVG: downloadSVG, downloadPNG: downloadPNG, saveBlob: saveBlob
  };
});
