/* The two figures the equilibrium workspace exists to produce.

   optimality  why this assemblage and not one of the others, at the
               condition on the panel: one bar per excluded phase, the
               KKT reduced cost it would have to give up to appear.
   margin      how far that answer is from a different one across
               temperature, with the assemblage ladder underneath it.

   Both draw from the record the Python package produced, and are drawn
   on the same kit at the workspace's square size. Neither evaluates a
   free energy: solidgas/activeset.py and its port do that, and the two
   data builders below only reshape what the solver already returned, so
   a figure cannot disagree with the engine it never calls.

   One deliberate difference from the plate. e2 plots the reduced cost
   per mol of oxygen, which is the comparable quantity at the single
   frozen condition it reports. Per mol O divides by the oxygen
   deficiency relative to the host, def_j = rho_host*t_j - o_j, and that
   denominator changes sign as soon as a candidate is less deficient
   than the host - over the conditions this page can be set to, 41% of
   reduced costs come back negative for that reason alone. The KKT test
   is stated per formula unit, where r_j > 0 means the phase is absent,
   and that is what a live figure has to draw. The per-mol-O column is
   still in the table on the KKT tab, next to the deficiency it is
   divided by. */

(function (root, factory) {
  'use strict';
  var kit = (typeof module !== 'undefined' && module.exports)
    ? require('./figkit.js') : root.FigKit;
  if (!kit) throw new Error('figures_thermo.js needs figkit.js');
  var api = factory(kit);
  if (typeof module !== 'undefined' && module.exports) module.exports = api;
  root.ThermoFigures = api;
})(typeof self !== 'undefined' ? self : this, function (K) {
  'use strict';

  var C = K.C, LW = K.LW, T = K.FIGTYPE;
  var lin = K.lin, frame = K.frame, axisX = K.axisX, axisY = K.axisY;
  var series = K.series, tint = K.tint, chem = K.chem;

  /* the phase that would appear first is the one worth a second colour;
     everything else is the same blue e2 draws every excluded phase in */
  var FAR = C.surface, NEAR = C.subsurface;

  /* --------------------------------------------------------- payloads */

  /* Reshape one solver result into the bar figure's payload. Phases whose
     reduced cost is not finite are dropped and counted: an oxygen-free
     gas leaves the oxygen potential unbounded below, so the formation
     direction has no cost to report and a bar would be a lie about it. */
  function optimalityData(R) {
    var rc = R.inactive_phase_reduced_costs || {};
    var tol = (R.solver_tolerance || {}).degeneracy_kJ_per_mol_O || 0;
    var costs = [], unbounded = [];
    Object.keys(rc).forEach(function (ph) {
      var pf = rc[ph].per_formula_kJ;
      if (pf === null || !isFinite(pf)) { unbounded.push(ph); return; }
      /* null per mol O means the candidate and the host differ by no
         oxygen at all, so there is no deficiency to divide by - not that
         the phase is sitting on a boundary. Math.abs(null) is 0, so this
         has to be tested before the tolerance, not after. */
      var po = rc[ph].per_mol_O_kJ;
      costs.push({ phase: ph, r_kJ: pf, r_per_mol_O_kJ: po,
                   degenerate: po !== null && isFinite(po)
                     && Math.abs(po) < tol });
    });
    costs.sort(function (a, b) { return a.r_kJ - b.r_kJ; });
    var cand = R.candidates || [];
    return {
      T_C: R.T_C,
      winner: (R.active_condensed_phases || []).slice(),
      costs: costs,
      unbounded: unbounded,
      candidates: {
        total: cand.length,
        feasible: cand.filter(function (c) { return c.feasible; }).length,
        optimal: cand.filter(function (c) {
          return c.feasible && c.kkt_ok; }).length
      }
    };
  }

  /* Reshape a temperature sweep. Each row carries the smallest reduced
     cost at that temperature and the phase holding it; the segments are
     the runs of one assemblage, and their edges are where the margin
     goes to zero - the curve's own dips are the same crossings, resolved
     only as far as the grid. */
  function marginData(rows, currentT_C) {
    var pts = rows.map(function (r) {
      var rc = r.inactive_phase_reduced_costs || {};
      var best = null, who = null;
      Object.keys(rc).forEach(function (ph) {
        var pf = rc[ph].per_formula_kJ;
        if (pf === null || !isFinite(pf)) return;
        if (best === null || pf < best) { best = pf; who = ph; }
      });
      return { T_C: r.T_C, margin_kJ: best, nearest: who,
               winner: (r.active_condensed_phases || []).slice() };
    });
    var segs = [], last = null;
    pts.forEach(function (q) {
      var key = q.winner.join('+');
      if (key !== last) {
        if (segs.length) segs[segs.length - 1].to_C = q.T_C;
        segs.push({ from_C: q.T_C, to_C: q.T_C, phases: q.winner.slice() });
        last = key;
      } else {
        segs[segs.length - 1].to_C = q.T_C;
      }
    });
    /* the point the panel is actually set to, so the reader can find it */
    var cur = null;
    if (currentT_C != null && pts.length) {
      cur = pts.reduce(function (a, b) {
        return Math.abs(b.T_C - currentT_C) < Math.abs(a.T_C - currentT_C)
          ? b : a;
      });
    }
    return {
      T_lo: pts.length ? pts[0].T_C : 0,
      T_hi: pts.length ? pts[pts.length - 1].T_C : 0,
      rows: pts, segments: segs,
      current: cur && { T_C: currentT_C, margin_kJ: cur.margin_kJ,
                        winner: cur.winner.slice(), nearest: cur.nearest }
    };
  }

  /* ---------------------------------------------------------- figures */

  function optimality(D) {
    var d = D.optimality;
    var f = K.square();
    /* The left margin holds a formula, not a number, so it is wider than
       the kit's default and there is no y axis to carry a label. Bars
       fill the box top to bottom, so the condition goes in a header band
       above the frame. */
    var p = { x0: f.pane.x0 + 34, y0: f.head, x1: f.pane.x1, y1: f.pane.y1 };
    var n = d.costs.length;
    if (!n) {
      f.text(f.W / 2, f.H / 2, 'no excluded phase has a finite reduced cost',
             { size: T.small, anchor: 'middle', fill: C.structure });
      return f.done();
    }
    var vmax = d.costs[n - 1].r_kJ;
    var X = lin(0, vmax * 1.15, p.x0, p.x1);
    var row = (p.y1 - p.y0) / n;
    var bh = row * 0.58;

    d.costs.forEach(function (q, i) {
      var near = i === 0;
      var col = near ? NEAR : FAR;
      var y = p.y0 + (i + 0.5) * row - bh / 2;
      f.rect(p.x0, y, X(q.r_kJ) - p.x0, bh, tint(col, 0.30), col, 0.6);
      f.text(p.x0 - 6, y + bh - 2, chem(q.phase),
             { size: T.small, anchor: 'end',
               weight: near ? 'bold' : null });
      f.text(X(q.r_kJ) + 5, y + bh - 2,
             q.r_kJ < 10 ? q.r_kJ.toFixed(2) : q.r_kJ.toFixed(1),
             { size: T.small, fill: col });
      if (near) {
        f.text(X(q.r_kJ) + 5, y + bh + T.note + 2,
               q.degenerate ? 'on the boundary' : 'nearest boundary',
               { size: T.note, fill: col });
      }
    });

    frame(f, X, [p.y0, p.y1]);
    axisX(f, X, p.y1, K.niceTicks(0, vmax * 1.15, 4),
          'KKT reduced cost (kJ)',
          function (v) { return vmax < 4 ? v.toFixed(1) : String(Math.round(v)); });

    var note = d.costs.length
      ? 'smallest margin: ' + chem(d.costs[0].phase) + ', '
        + (d.costs[0].r_kJ < 10 ? d.costs[0].r_kJ.toFixed(2)
                                : d.costs[0].r_kJ.toFixed(1)) + ' kJ'
      : 'no finite margin';
    if (d.unbounded.length) {
      note += ' \u00b7 ' + d.unbounded.length + ' unbounded (no gas oxygen)';
    }
    f.header(chem(d.winner.join(' + ')) + ' stable at ' + Math.round(d.T_C)
             + ' °C', note);
    return f.done();
  }

  function margin(D) {
    var d = D.margin;
    var f = K.square();
    /* the same header band the bar figure carries, so the pair reads as
       one answer given two ways */
    var p = { x0: f.pane.x0, y0: f.head, x1: f.pane.x1, y1: f.pane.y1 };
    var STRIP = 22;                       /* the assemblage ladder */
    var yStrip = p.y1 - STRIP;
    var vals = d.rows.map(function (q) { return q.margin_kJ; })
      .filter(function (v) { return v != null && isFinite(v); });
    var vmax = vals.length ? Math.max.apply(null, vals) : 1;
    var X = lin(d.T_lo, d.T_hi, p.x0, p.x1);
    var Y = lin(0, vmax * 1.08, yStrip, p.y0);

    /* Every run of one assemblage gets a cell, and the two-phase runs -
       the tie lines the ladder rests on - are shaded darker than the
       one-phase fields between them, so the alternation is readable
       without a label on every segment; most are a few points wide. The
       curve dips to zero at each of these edges, which is the same
       crossing seen from above. */
    d.segments.forEach(function (g) {
      var x0 = X(g.from_C), x1 = X(g.to_C);
      f.rect(x0, yStrip, x1 - x0, STRIP, g.phases.length > 1
             ? tint(C.surface, 0.34) : tint(C.structure, 0.16));
      if (g.from_C > d.T_lo) {
        f.line(x0, yStrip, x0, p.y1, C.paper, LW.axis);
      }
    });
    f.line(p.x0, yStrip, p.x1, yStrip, C.structure, LW.hair);
    /* a cell is named only if the name fits inside it; the ticks come up
       4 pt from the bottom spine, so the text sits clear above them */
    d.segments.forEach(function (g) {
      var x0 = X(g.from_C), x1 = X(g.to_C);
      var label = chem(g.phases.join('+'));
      if (x1 - x0 < label.length * 0.55 * T.note + 6) return;
      f.text((x0 + x1) / 2, yStrip + STRIP - 7, label,
             { size: T.note, anchor: 'middle' });
    });

    series(f, d.rows.map(function (q) { return [q.T_C, q.margin_kJ]; }),
           X, Y, C.surface, LW.curve);

    frame(f, X, [p.y0, p.y1]);
    axisX(f, X, p.y1, K.niceTicks(d.T_lo, d.T_hi, 4), 'temperature (°C)',
          function (v) { return String(Math.round(v)); });
    axisY(f, Y, p.x0, K.niceTicks(0, vmax * 1.08, 4),
          'reduction margin (kJ)',
          function (v) { return vmax < 4 ? v.toFixed(1) : String(Math.round(v)); });

    var c = d.current;
    if (c && c.margin_kJ != null) {
      var qx = X(Math.min(Math.max(c.T_C, d.T_lo), d.T_hi));
      var qy = Y(Math.min(c.margin_kJ, vmax * 1.08));
      /* The curve runs through the marker, so a box directly above or
         below it is never clear on a steep limb; the room is diagonal.
         Each corner is tried in turn and the first free one takes the
         word - and on a sawtooth, where a spike can fill all four, the
         ring goes down alone, because the header already names the
         condition and a word laid over the curve costs more than it
         says. Free is measured on the curve's segments, not its sample
         points: on a steep limb two neighbours are a point and a half
         apart in x and thirty apart in y, and a box tested only at the
         points is jumped clean over. */
      var pad = (d.T_hi - d.T_lo) * 30 / (p.x1 - p.x0);
      var nearPts = d.rows.filter(function (q) {
        return q.margin_kJ != null && isFinite(q.margin_kJ)
          && Math.abs(q.T_C - c.T_C) <= pad;
      }).map(function (q) {
        return [X(q.T_C), Y(Math.min(q.margin_kJ, vmax * 1.08))];
      });
      var clear = function (bx0, bx1, top, bot) {
        if (top < p.y0 + 1 || bot > yStrip - 1) return false;
        if (bx0 < p.x0 + 1 || bx1 > p.x1 - 1) return false;
        for (var i = 1; i < nearPts.length; i++) {
          var a = nearPts[i - 1], b = nearPts[i];
          if (Math.max(a[0], b[0]) < bx0 || Math.min(a[0], b[0]) > bx1) {
            continue;
          }
          if (Math.min(a[1], b[1]) < bot && Math.max(a[1], b[1]) > top) {
            return false;
          }
        }
        return true;
      };
      f.dot(qx, qy, 6, C.paper, C.gas);
      f.dot(qx, qy, 2.2, C.gas);
      var W = 9 * 0.6 * T.note, corner = null;
      [[9, -16, 'start'], [-9, -16, 'end'],
       [9, 20, 'start'], [-9, 20, 'end']].forEach(function (o) {
        if (corner) return;
        var bx0 = o[0] > 0 ? qx + o[0] : qx + o[0] - W;
        if (clear(bx0, bx0 + W, qy + o[1] - T.note, qy + o[1] + 4)) corner = o;
      });
      if (corner) {
        f.text(qx + corner[0], qy + corner[1], 'current T',
               { size: T.note, anchor: corner[2], fill: C.gas,
                 weight: 'bold' });
      }
      f.header(chem(c.winner.join(' + ')) + ' stable at ' + Math.round(c.T_C)
               + ' °C', 'nearest reduced phase ' + chem(c.nearest || '—')
               + ', margin ' + (c.margin_kJ < 10 ? c.margin_kJ.toFixed(2)
                                                  : c.margin_kJ.toFixed(1))
               + ' kJ');
    }
    return f.done();
  }

  /* --------------------------------------------- conversion over T */

  /* How much of the CO2 is gone, at every temperature the panel admits.
     The rows are the sweep the margin figure already ran, so this costs
     nothing extra and cannot report a temperature the curve does not.

     One thing the number hides and the band shows. Conversion is counted
     against the feed, so oxygen that left the gas is counted whether it
     went into CO or into the solid. Where the assemblage is no longer the
     host alone, part of the answer is the solid being eaten, and the
     shaded run says exactly where that starts. */
  function conversionData(rows, currentT_C, host, exactPct) {
    var pts = rows.map(function (r) {
      return { T_C: r.T_C, conv_pct: r.conversion_CO2_pct,
               winner: (r.active_condensed_phases || []).slice() };
    }).filter(function (q) {
      return q.conv_pct != null && isFinite(q.conv_pct);
    });
    var segs = [], last = null;
    pts.forEach(function (q) {
      var key = q.winner.join('+');
      if (key !== last) {
        if (segs.length) segs[segs.length - 1].to_C = q.T_C;
        segs.push({ from_C: q.T_C, to_C: q.T_C, phases: q.winner.slice() });
        last = key;
      } else {
        segs[segs.length - 1].to_C = q.T_C;
      }
    });
    /* The curve is a grid; the panel is not. Reading the marker off the
       nearest sample would put a number in the header the workspace never
       reported - eleven degrees away is a fifth of a percent here - so the
       caller passes the solve it already has and the grid is only drawn. */
    var cur = null;
    if (currentT_C != null && pts.length) {
      cur = pts.reduce(function (a, b) {
        return Math.abs(b.T_C - currentT_C) < Math.abs(a.T_C - currentT_C)
          ? b : a;
      });
    }
    var host2 = host || 'TiO2';
    return {
      T_lo: pts.length ? pts[0].T_C : 0,
      T_hi: pts.length ? pts[pts.length - 1].T_C : 0,
      rows: pts, segments: segs, host: host2,
      any_reduction: segs.some(function (g) {
        return g.phases.length !== 1 || g.phases[0] !== host2;
      }),
      current: cur && {
        T_C: currentT_C,
        conv_pct: exactPct != null && isFinite(exactPct) ? exactPct
                                                         : cur.conv_pct,
        winner: cur.winner.slice() }
    };
  }

  function conversion(D) {
    var d = D.conversion;
    var f = K.square();
    var p = { x0: f.pane.x0, y0: f.head, x1: f.pane.x1, y1: f.pane.y1 };
    var vals = d.rows.map(function (q) { return q.conv_pct; });
    var vmax = vals.length ? Math.max.apply(null, vals) : 1;
    var vmin = vals.length ? Math.min.apply(null, vals) : 0;
    var top = Math.max(vmax * 1.08, 1e-9);
    var bot = Math.min(vmin, 0);
    var X = lin(d.T_lo, d.T_hi, p.x0, p.x1);
    var Y = lin(bot, top, p.y1, p.y0);

    /* the runs where something other than the host alone is standing */
    d.segments.forEach(function (g) {
      if (g.phases.length === 1 && g.phases[0] === d.host) return;
      var x0 = X(g.from_C), x1 = X(g.to_C);
      f.rect(x0, p.y0, x1 - x0, p.y1 - p.y0, tint(C.subsurface, 0.16));
      if (x1 - x0 > 14 * 0.6 * T.note + 6) {
        f.text((x0 + x1) / 2, p.y0 + T.note + 5, 'solid reducing',
               { size: T.note, anchor: 'middle', fill: C.subsurface });
      }
    });
    if (bot < 0) f.line(p.x0, Y(0), p.x1, Y(0), C.structure, LW.hair);

    series(f, d.rows.map(function (q) { return [q.T_C, q.conv_pct]; }),
           X, Y, C.surface, LW.curve);

    frame(f, X, [p.y0, p.y1]);
    axisX(f, X, p.y1, K.niceTicks(d.T_lo, d.T_hi, 4), 'temperature (°C)',
          function (v) { return String(Math.round(v)); });
    axisY(f, Y, p.x0, K.niceTicks(bot, top, 4),
          'CO₂ converted (%)',
          function (v) { return top < 4 ? v.toFixed(1) : String(Math.round(v)); });

    var c = d.current;
    if (c && c.conv_pct != null) {
      var qx = X(Math.min(Math.max(c.T_C, d.T_lo), d.T_hi));
      var qy = Y(Math.min(Math.max(c.conv_pct, bot), top));
      f.line(qx, Y(bot), qx, qy, C.gas, LW.hair, '3 3');
      f.dot(qx, qy, 6, C.paper, C.gas);
      f.dot(qx, qy, 2.2, C.gas);
      f.header('CO₂ conversion ' + c.conv_pct.toFixed(1) + '% at '
               + Math.round(c.T_C) + ' °C', d.any_reduction
               ? 'shaded: includes oxygen transfer from the solid'
               : 'no oxygen transfer from the solid');
    }
    return f.done();
  }

  /* ------------------------------------------- where reduction starts */

  /* The line the feed has to stay above. Below it the host gives way to
     the phase named in the header; above it the host survives.

     The y axis is a feed, not a potential, because a feed is what anyone
     reading this can set on a mass-flow controller. Every gas can be put
     on it: a mixture of H2 and CO2 lands at the fraction it was mixed at,
     and anything else lands at the H2/CO2 mixture that would be equally
     reducing. The marker is the panel's own feed before the solid touches
     it, which is the quantity the question is about - once a closed charge
     has equilibrated against a reducible solid the gas is buffered onto
     this very line and no longer says how reducing it started out. */
  function boundaryData(curve, feed, currentT_C, exactAt) {
    var pts = (curve || []).filter(function (q) {
      return q && q.y_CO2 > 0 && isFinite(q.y_CO2);
    }).map(function (q) {
      return { T_C: q.T_C, pct: q.y_CO2 * 100.0, phase: q.phase };
    });
    /* same rule as the conversion figure: the header quotes the closed
       form at the panel's own temperature, not the nearest grid point -
       on a curve this steep those differ by eight percent */
    var at = null;
    if (exactAt && exactAt.y_CO2 > 0 && isFinite(exactAt.y_CO2)) {
      at = { T_C: currentT_C, pct: exactAt.y_CO2 * 100.0,
             phase: exactAt.phase };
    } else if (currentT_C != null && pts.length) {
      at = pts.reduce(function (a, b) {
        return Math.abs(b.T_C - currentT_C) < Math.abs(a.T_C - currentT_C)
          ? b : a;
      });
    }
    var fd = null;
    if (feed && feed.y_CO2 > 0 && isFinite(feed.y_CO2)) {
      fd = { T_C: feed.T_C, pct: feed.y_CO2 * 100.0,
             literal: !!feed.literal,
             ratio: at ? feed.y_CO2 * 100.0 / at.pct : null };
    }
    return {
      T_lo: pts.length ? pts[0].T_C : 0,
      T_hi: pts.length ? pts[pts.length - 1].T_C : 0,
      rows: pts, host: (feed && feed.host) || 'TiO2',
      at: at, feed: fd, why_no_feed: (feed && feed.why) || null
    };
  }

  function boundary(D) {
    var d = D.boundary;
    var f = K.square();
    var p = { x0: f.pane.x0 + 8, y0: f.head, x1: f.pane.x1, y1: f.pane.y1 };
    var vals = d.rows.map(function (q) { return q.pct; });
    if (d.feed) vals = vals.concat([d.feed.pct]);
    var lo = Math.pow(10, Math.floor(log10(Math.min.apply(null, vals))) - 0.3);
    var hi = Math.pow(10, Math.ceil(log10(Math.max.apply(null, vals))) + 0.3);
    var X = lin(d.T_lo, d.T_hi, p.x0, p.x1);
    var Y = K.lg(lo, hi, p.y1, p.y0);

    /* everything above the line is the host surviving; fill it so the
       reader never has to work out which side is which */
    var above = ['M', p.x0, p.y0].join(' ');
    d.rows.forEach(function (q) { above += ' L ' + X(q.T_C) + ' ' + Y(q.pct); });
    above += ' L ' + p.x1 + ' ' + p.y0 + ' Z';
    f.path(above, 'none', 0, null, tint(C.surface, 0.13));

    series(f, d.rows.map(function (q) { return [q.T_C, q.pct]; }),
           X, Y, C.surface, LW.curve);

    /* Where the marker is, and how far its stem runs, worked out before
       anything is written so the words can get out of its way. */
    var q = d.feed, qx = null, qy = null, stem = null;
    if (q) {
      qx = X(Math.min(Math.max(q.T_C, d.T_lo), d.T_hi));
      qy = Y(Math.min(Math.max(q.pct, lo), hi));
      if (d.at) {
        var qb = Y(Math.min(Math.max(d.at.pct, lo), hi));
        stem = [Math.min(qy, qb), Math.max(qy, qb)];
      }
    }

    /* One word in each field. Anywhere along the line a horizontal word
       set close to it gets crossed, because over its own width the line
       climbs a decade or more. The right-hand edge is the one place that
       cannot happen: the curve leaves the panel at a fixed height, and
       the two words sit above and below that exit, each in its own
       field, reading against the spine rather than against the slope.
       A word is dropped if the marker's stem would run through it - the
       header and the marker have already said which side the feed is
       on, and a struck-through label says it worse. */
    var last = d.rows[d.rows.length - 1];
    if (last) {
      var ye = Y(last.pct);
      var edge = function (y, txt, col) {
        var w = txt.length * 0.55 * T.note + 6;
        var hitX = qx !== null && qx > p.x1 - 12 - w && qx < p.x1 - 1;
        var hitY = stem && stem[0] < y + 4 && stem[1] > y - T.note;
        if (hitX && (hitY || (qy !== null && qy > y - T.note && qy < y + 4))) {
          return;
        }
        f.text(p.x1 - 12, y, txt, { size: T.note, anchor: 'end', fill: col });
      };
      edge(Math.max(ye - 10, p.y0 + T.note + 6), chem(d.host) + ' stable', C.surface);
      edge(Math.min(ye + 22, p.y1 - 8), 'reduction', C.structure);
    }

    frame(f, X, [p.y0, p.y1]);
    axisX(f, X, p.y1, K.niceTicks(d.T_lo, d.T_hi, 4), 'temperature (°C)',
          function (v) { return String(Math.round(v)); });
    axisY(f, Y, p.x0, thinDecades(lo, hi, 6), 'CO₂ in the feed (mol%)',
          K.expLabel);

    if (q) {
      if (stem) f.line(qx, stem[0], qx, stem[1], C.subsurface, LW.hair, '3 3');
      f.rect(qx - 5, qy - 5, 10, 10, C.subsurface, C.paper, LW.hair);
      var right = qx > p.x0 + (p.x1 - p.x0) * 0.62;
      f.text(qx + (right ? -10 : 10), qy + 4, 'feed',
             { size: T.note, anchor: right ? 'end' : 'start',
               fill: C.subsurface });
    }
    if (d.at) {
      f.header(chem(d.host) + ' reduces to ' + chem(d.at.phase) + ' below '
               + fmtPct(d.at.pct) + ' CO₂', null);
      /* with no feed marked there is nothing to say about one, and a
         warning about a feed the caller never supplied reads as an error */
      var why = q
        ? (q.ratio >= 1
            ? 'feed: ' + fmtTimes(q.ratio) + ' the boundary'
            : 'feed: ' + fmtTimes(1 / q.ratio)
              + ' below the boundary, the solid reduces')
        : d.why_no_feed;
      if (why) {
        f.text(8, 8 + T.body + 4 + T.note * 0.85, why,
               { size: T.note, fill: q ? C.structure : C.subsurface });
      }
    }
    return f.done();
  }

  function log10(v) { return Math.log(v) / Math.LN10; }

  /* Decade ticks, thinned to every second or fifth one when the span is
     tall enough that labelling all of them would collide. */
  function thinDecades(lo, hi, want) {
    var all = K.decades(lo, hi);
    var step = 1;
    while (Math.ceil(all.length / step) > want) step += 1;
    return all.filter(function (v, i) {
      return (all.length - 1 - i) % step === 0;
    });
  }

  function fmtPct(v) {
    if (v >= 1) return v.toFixed(1) + '%';
    if (v >= 0.01) return v.toFixed(3) + '%';
    return v.toPrecision(2) + '%';
  }

  function fmtTimes(r) {
    if (r >= 1e4) return r.toPrecision(2).replace(/e\+?(-?\d+)/,
      '×10^$1') + '×';
    if (r >= 100) return String(Math.round(r)) + '×';
    if (r >= 10) return r.toFixed(0) + '×';
    return r.toFixed(1) + '×';
  }

  return { optimality: optimality, margin: margin,
           conversion: conversion, boundary: boundary,
           optimalityData: optimalityData, marginData: marginData,
           conversionData: conversionData, boundaryData: boundaryData,
           FIGURES: { optimality: optimality, margin: margin,
                      conversion: conversion, boundary: boundary } };
});
