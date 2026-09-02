/* The four figures on the equilibrium workspace.

   None of them solves anything of its own. The bar and the reduction
   line come off the result already on screen; the two curves ask
   window.ThermoBridge to re-run the charge on the panel at each
   temperature - the same solver, the same inputs, one argument changed -
   so a figure here cannot report a condition the workspace was never
   set to.

   Two cost nothing and draw on every run: the bar, and the reduction
   line, whose curve is closed form. Two need the temperature sweep and
   sit behind its button, because that is 121 solves, about half a
   second. Once the sweep has been run they follow every later solve on
   their own. */

(function () {
  'use strict';

  var TB = window.ThermoBridge, FIG = window.ThermoFigures;
  var KIT = window.FigKit;
  if (!TB || !FIG || !KIT) return;

  var $ = function (id) { return document.getElementById(id); };
  var btn = $('sweepBtn'), note = $('sweepNote');
  if (!$('figOptimality') || !$('figMargin') || !btn) return;
  var hasNew = $('figBoundary') && $('figConversion');

  var figState = {};
  function mountFigure(hostId, name, svg) {
    var host = $(hostId);
    if (!host) return;
    figState[hostId] = { svg: svg, name: name };
    host.querySelector('.figbox').innerHTML = svg;
    var dl = host.querySelector('.figdl');
    if (dl && !dl.dataset.wired) {
      dl.dataset.wired = '1';
      dl.addEventListener('click', function (ev) {
        var fmt = ev.target && ev.target.getAttribute('data-fmt');
        if (!fmt) return;
        var st = figState[hostId];
        if (!st) return;
        if (fmt === 'svg') {
          KIT.downloadSVG(st.svg, st.name);
        } else {
          ev.target.disabled = true;
          KIT.downloadPNG(st.svg, st.name).catch(function (e) {
            window.alert('PNG export failed: ' + e.message);
          }).then(function () { ev.target.disabled = false; });
        }
        dl.open = false;
      });
    }
  }
  function figFail(hostId, why) {
    var host = $(hostId);
    if (host) host.querySelector('.figbox').textContent = why;
  }

  /* -------------------------------------------------------------- bar */

  function drawOptimality(R) {
    if (!R) { figFail('figOptimality', 'no result to draw'); return; }
    try {
      var d = FIG.optimalityData(R);
      mountFigure('figOptimality', 'equilibrium-optimality',
                  FIG.optimality({ optimality: d }));
      var cap = $('capOptimality');
      if (cap) {
        var near = d.costs[0];
        cap.innerHTML = !near
          ? 'No reduced Ti&ndash;O phase carries a finite margin at this condition.'
          : (near.r_kJ > 0
             ? 'All reduction margins are positive, so no reduced Ti&ndash;O phase is stable. '
               + KIT.chem(near.phase) + ' has the smallest margin and defines the nearest phase boundary.'
             : KIT.chem(near.phase) + ' is at zero margin and is part of the stable assemblage.');
      }
    } catch (e) {
      figFail('figOptimality', 'figure failed: ' + e.message);
    }
  }

  /* ------------------------------------------ the line, drawn for free */

  var BOUND_N = 181;

  /* The same statement the margin curve makes, in the units of a
     mass-flow controller instead of kilojoules: not "54.9 kJ from the
     nearest phase" but "7918 times more CO2 than it would take". The
     curve is two exponentials a point, so unlike the margin it does not
     need the sweep behind it. */
  function drawBoundary() {
    var r = range(), lo = r[0], hi = r[1];
    var curve = [], i, T_C, b;
    for (i = 0; i < BOUND_N; i++) {
      T_C = lo + (hi - lo) * i / (BOUND_N - 1);
      b = TB.boundaryAt(T_C);
      if (b) curve.push({ T_C: T_C, y_CO2: b.y_CO2, phase: b.phase });
    }
    if (!curve.length) {
      figFail('figBoundary', 'nothing in the registry is more reduced than '
              + TB.host());
      return;
    }
    var here = TB.T_C();
    var d = FIG.boundaryData(curve, TB.feedPoint(here), here, TB.boundaryAt(here));
    mountFigure('figBoundary', 'equilibrium-boundary', FIG.boundary({ boundary: d }));
    var cap = $('capBoundary');
    if (cap) {
      var pct = function (v) { return v >= 1 ? v.toFixed(1) : v >= 0.01 ? v.toFixed(3) : v.toPrecision(2); };
      var times = function (r) { return r >= 100 ? String(Math.round(r)) : r >= 10 ? r.toFixed(0) : r.toFixed(1); };
      var txt = d.at
        ? 'At ' + Math.round(here) + ' °C, ' + KIT.chem(d.host) + ' begins to reduce below '
          + pct(d.at.pct) + ' mol% CO₂.'
        : '';
      if (d.at && d.feed && d.feed.ratio != null) {
        txt += d.feed.ratio >= 1
          ? ' The current feed contains ' + pct(d.feed.pct) + ' mol% CO₂, '
            + times(d.feed.ratio) + ' times this boundary.'
          : ' The current feed contains ' + pct(d.feed.pct) + ' mol% CO₂, '
            + times(1 / d.feed.ratio) + ' times below this boundary, so the solid reduces.';
      } else if (d.why_no_feed) {
        txt += ' ' + d.why_no_feed;
      }
      cap.innerHTML = txt;
    }
  }

  /* ------------------------------------------------------------ sweep */

  var SWEEP_N = 121;
  var swept = false;

  /* the range the panel itself admits, so the curve cannot run past a
     temperature the workspace would refuse */
  function range() {
    var el = $('T');
    var lo = parseFloat(el && el.getAttribute('min'));
    var hi = parseFloat(el && el.getAttribute('max'));
    return [isFinite(lo) ? lo : 300, isFinite(hi) ? hi : 1650];
  }

  function drawSwept() {
    var r = range(), lo = r[0], hi = r[1];
    var rows = [], i, T_C;
    for (i = 0; i < SWEEP_N; i++) {
      T_C = lo + (hi - lo) * i / (SWEEP_N - 1);
      var one = TB.solveAt(T_C);
      if (!one) { figFail('figMargin', 'gas feed sums to zero'); return 0; }
      rows.push(one);
    }
    var here = TB.T_C();
    mountFigure('figMargin', 'equilibrium-margin',
                FIG.margin({ margin: FIG.marginData(rows, here) }));
    if (hasNew) {
      /* the marker quotes the solve the workspace is showing, not the
         nearest point of the grid - eleven degrees apart is a fifth of a
         percent of conversion, and the header would not match the panel */
      var cur = TB.last();
      mountFigure('figConversion', 'equilibrium-conversion',
                  FIG.conversion({ conversion: FIG.conversionData(
                    rows, here, TB.host(),
                    cur && cur.conversion_CO2_pct) }));
    }
    return rows.length;
  }

  function sweep() {
    btn.disabled = true;
    var was = btn.textContent;
    btn.textContent = 'sweeping…';
    /* let the label paint before the solver takes the thread */
    window.setTimeout(function () {
      var t0 = (window.performance || Date).now();
      var n = 0;
      try {
        n = drawSwept();
      } catch (e) {
        figFail('figMargin', 'figure failed: ' + e.message);
      }
      var dt = (window.performance || Date).now() - t0;
      btn.disabled = false;
      btn.textContent = was;
      if (n) {
        swept = true;
        var r = range();
        note.textContent = n + ' solves from ' + r[0] + ' to ' + r[1]
          + ' °C in ' + dt.toFixed(0) + ' ms · redrawn with every run';
      }
    }, 0);
  }

  btn.addEventListener('click', sweep);
  TB.onSolve(function (R) {
    drawOptimality(R);
    if (hasNew) {
      try { drawBoundary(); }
      catch (e) { figFail('figBoundary', 'figure failed: ' + e.message); }
    }
    if (swept) { try { drawSwept(); } catch (e) { swept = false; } }
  });
  drawOptimality(TB.last());
  if (hasNew) {
    try { drawBoundary(); }
    catch (e) { figFail('figBoundary', 'figure failed: ' + e.message); }
  }
}());
