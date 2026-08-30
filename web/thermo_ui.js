/* The two figures on the equilibrium workspace.

   Neither solves anything. The bar comes off the result already on
   screen, and the curve asks window.ThermoBridge to re-run the charge on
   the panel at each temperature - the same solver, the same inputs, one
   argument changed - so a figure here cannot report a condition the
   workspace was never set to.

   The sweep is behind a button because it is 121 solves, about half a
   second, and the bar answers the question most visitors are asking. Once
   it has been drawn it follows every later solve on its own. */

(function () {
  'use strict';

  var TB = window.ThermoBridge, FIG = window.ThermoFigures;
  var KIT = window.FigKit;
  if (!TB || !FIG || !KIT) return;

  var $ = function (id) { return document.getElementById(id); };
  var btn = $('sweepBtn'), note = $('sweepNote');
  if (!$('figOptimality') || !$('figMargin') || !btn) return;

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
      mountFigure('figOptimality', 'equilibrium-optimality',
                  FIG.optimality({ optimality: FIG.optimalityData(R) }));
    } catch (e) {
      figFail('figOptimality', 'figure failed: ' + e.message);
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

  function drawMargin() {
    var r = range(), lo = r[0], hi = r[1];
    var rows = [], i, T_C;
    for (i = 0; i < SWEEP_N; i++) {
      T_C = lo + (hi - lo) * i / (SWEEP_N - 1);
      var one = TB.solveAt(T_C);
      if (!one) { figFail('figMargin', 'gas feed sums to zero'); return 0; }
      rows.push(one);
    }
    mountFigure('figMargin', 'equilibrium-margin',
                FIG.margin({ margin: FIG.marginData(rows, TB.T_C()) }));
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
        n = drawMargin();
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
    if (swept) { try { drawMargin(); } catch (e) { swept = false; } }
  });
  drawOptimality(TB.last());
}());
