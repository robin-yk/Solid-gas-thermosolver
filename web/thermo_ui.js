/* The seven figures on the thermodynamic-equilibrium workspace.

   None of them solves anything of its own. The bar and the reduction
   line come off the result already on screen; the two curves ask
   window.ThermoBridge to re-run the charge on the panel at each
   temperature - the same solver, the same inputs, one argument changed -
   so a figure here cannot report a condition the workspace was never
   set to.

   Two cost nothing and draw on every run: the bar and the reduction line.
   Two share the optional 121-point temperature sweep. The operating atlas
   adds a browser-computed gas-only map, a feed-ratio trade-off, and a
   49-point finite gas-solid comparison. */

(function () {
  'use strict';

  var TB = window.ThermoBridge, FIG = window.ThermoFigures;
  var KIT = window.FigKit;
  if (!TB || !FIG || !KIT) return;

  var $ = function (id) { return document.getElementById(id); };
  var btn = $('sweepBtn'), note = $('sweepNote');
  if (!$('figOptimality') || !$('figMargin') || !btn) return;
  var hasNew = $('figBoundary') && $('figConversion');
  var atlasBtn = $('atlasRun'), atlasNote = $('atlasNote');
  var hasAtlas = atlasBtn && $('figOperatingMap') && $('figRatioTradeoff')
    && $('figSolidCoupling');

  var figState = {};
  function mountFigure(hostId, name, svg, csv) {
    var host = $(hostId);
    if (!host) return;
    figState[hostId] = { svg: svg, name: name, csv: csv || null };
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
        } else if (fmt === 'png') {
          ev.target.disabled = true;
          KIT.downloadPNG(st.svg, st.name).catch(function (e) {
            window.alert('PNG export failed: ' + e.message);
          }).then(function () { ev.target.disabled = false; });
        } else if (fmt === 'csv' && st.csv) {
          KIT.saveBlob(new Blob([st.csv], { type: 'text/csv;charset=utf-8' }),
                       st.name + '.csv');
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

  /* ------------------------------------------------ operating atlas */

  var atlasBuilt = false;
  function csvOf(rows, columns) {
    var quote = function (v) {
      var s = v == null ? '' : String(v);
      return /[",\n]/.test(s) ? '"' + s.replace(/"/g, '""') + '"' : s;
    };
    var keys = columns.map(function (c) {
      return typeof c === 'string' ? c : c[0];
    });
    var labels = columns.map(function (c) {
      return typeof c === 'string' ? c : c[1];
    });
    return labels.join(',') + '\n' + rows.map(function (q) {
      return keys.map(function (k) { return quote(q[k]); }).join(',');
    }).join('\n') + '\n';
  }

  function logPoint(lo, hi, i, n, centred) {
    var f = centred ? (i + 0.5) / n : i / (n - 1);
    return Math.exp(Math.log(lo) + (Math.log(hi) - Math.log(lo)) * f);
  }

  function currentRatio(T_C) {
    var pt = TB.feedPoint(T_C);
    if (!pt || !(pt.y_CO2 > 0) || !(pt.y_CO2 < 1)) {
      return { ratio: 1, literal: true, fallback: true };
    }
    return { ratio: (1 - pt.y_CO2) / pt.y_CO2,
             literal: pt.literal, fallback: false };
  }

  function drawAtlas() {
    var t0 = (window.performance || Date).now();
    var T_LO = 400, T_HI = 1500, R_LO = 0.1, R_HI = 1e7;
    var NT = 44, NR = 36, i, j, T_C, ratio, one;
    var cells = [], boundary = [];
    for (i = 0; i < NT; i++) {
      T_C = T_LO + (T_HI - T_LO) * (i + 0.5) / NT;
      var cellBoundary = TB.boundaryAt(T_C);
      var cellBoundaryRatio = cellBoundary && cellBoundary.y_CO2 > 0
        && cellBoundary.y_CO2 < 1
        ? (1 - cellBoundary.y_CO2) / cellBoundary.y_CO2 : null;
      for (j = 0; j < NR; j++) {
        ratio = logPoint(R_LO, R_HI, j, NR, true);
        one = TB.gasOnlyAt(T_C, ratio);
        cells.push({ T_C: T_C, ratio: ratio,
                     conv_pct: one.CO2_conversion_pct,
                     reduction_boundary_ratio: cellBoundaryRatio,
                     TiO2_state: cellBoundaryRatio !== null
                       && ratio > cellBoundaryRatio
                       ? 'reduction favoured' : 'rutile stable' });
      }
    }
    for (i = 0; i < 181; i++) {
      T_C = T_LO + (T_HI - T_LO) * i / 180;
      var b = TB.boundaryAt(T_C);
      if (b && b.y_CO2 > 0 && b.y_CO2 < 1) {
        boundary.push({ T_C: T_C, ratio: (1 - b.y_CO2) / b.y_CO2,
                        phase: b.phase });
      }
    }
    var hereT = TB.T_C(), picked = currentRatio(hereT);
    var exact = TB.gasOnlyAt(hereT, picked.ratio);
    var map = {
      T_lo: T_LO, T_hi: T_HI, ratio_lo: R_LO, ratio_hi: R_HI,
      nT: NT, nR: NR, cells: cells, boundary: boundary, host: TB.host(),
      current: { T_C: hereT, ratio: picked.ratio,
                 conv_pct: exact.CO2_conversion_pct,
                 literal: picked.literal }
    };
    mountFigure('figOperatingMap', 'rwgs-operating-window',
                FIG.operatingMap({ operatingMap: map }),
                csvOf(cells, [
                  ['T_C', 'temperature_C'],
                  ['ratio', 'H2_CO2_ratio'],
                  ['conv_pct', 'gas_only_CO2_conversion_pct'],
                  ['reduction_boundary_ratio',
                   'TiO2_reduction_boundary_H2_CO2_ratio'],
                  'TiO2_state'
                ]));

    var trade = [], NTRADE = 121;
    for (i = 0; i < NTRADE; i++) {
      ratio = logPoint(0.1, 10, i, NTRADE, false);
      one = TB.gasOnlyAt(hereT, ratio);
      trade.push({ ratio: ratio, co2: one.CO2_conversion_pct,
                   h2: one.H2_utilization_pct,
                   co: one.CO_yield_per_mol_feed_pct });
    }
    var tradeCurrent = picked.ratio >= 0.1 && picked.ratio <= 10
      ? { ratio: picked.ratio, co2: exact.CO2_conversion_pct } : null;
    mountFigure('figRatioTradeoff', 'rwgs-feed-ratio-tradeoff',
      FIG.ratioTradeoff({ ratioTradeoff: {
        T_C: hereT, ratio_lo: 0.1, ratio_hi: 10,
        rows: trade, current: tradeCurrent
      } }), csvOf(trade, [
        ['ratio', 'H2_CO2_ratio'], ['co2', 'CO2_conversion_pct'],
        ['h2', 'H2_utilization_pct'],
        ['co', 'CO_per_total_fresh_gas_pct']
      ]));

    var coupled = [], NCOUPLE = 49;
    for (i = 0; i < NCOUPLE; i++) {
      ratio = logPoint(R_LO, R_HI, i, NCOUPLE, false);
      var gas = TB.gasOnlyAt(hereT, ratio);
      var both = TB.solveFeedAt(hereT, { CO2: 1, H2: ratio });
      coupled.push({ ratio: ratio, gas_only: gas.CO2_conversion_pct,
                     gas_solid: both.conversion_CO2_pct,
                     reduced_pct: both.reduced_pct,
                     phases: both.active_condensed_phases.join('+') });
    }
    var at = TB.boundaryAt(hereT);
    var br = at && at.y_CO2 > 0 && at.y_CO2 < 1
      ? (1 - at.y_CO2) / at.y_CO2 : R_HI;
    mountFigure('figSolidCoupling', 'rwgs-gas-solid-coupling',
      FIG.solidCoupling({ solidCoupling: {
        T_C: hereT, ratio_lo: R_LO, ratio_hi: R_HI, rows: coupled,
        boundary_ratio: br, host: TB.host()
      } }), csvOf(coupled, [
        ['ratio', 'H2_CO2_ratio'],
        ['gas_only', 'gas_only_CO2_conversion_pct'],
        ['gas_solid', 'gas_solid_CO2_conversion_pct'],
        ['reduced_pct', 'Ti3_fraction_pct'], 'phases'
      ]));

    var dt = (window.performance || Date).now() - t0;
    atlasBuilt = true;
    atlasNote.textContent = (NT * NR).toLocaleString()
      + ' gas-phase states + ' + NCOUPLE + ' gas–solid equilibria · '
      + dt.toFixed(0) + ' ms · temperature ' + Math.round(hereT) + ' °C'
      + (picked.fallback ? ' · 1:1 feed used because this gas has no RWGS ratio' : '');
  }

  function buildAtlas() {
    if (!hasAtlas) return;
    atlasBtn.disabled = true;
    atlasBtn.textContent = 'calculating…';
    window.setTimeout(function () {
      try { drawAtlas(); }
      catch (e) {
        figFail('figOperatingMap', 'figure failed: ' + e.message);
        atlasNote.textContent = 'Atlas calculation failed: ' + e.message;
      }
      atlasBtn.disabled = false;
      atlasBtn.textContent = 'Recalculate atlas';
    }, 0);
  }

  if (hasAtlas) atlasBtn.addEventListener('click', buildAtlas);

  btn.addEventListener('click', sweep);
  TB.onSolve(function (R) {
    drawOptimality(R);
    if (hasNew) {
      try { drawBoundary(); }
      catch (e) { figFail('figBoundary', 'figure failed: ' + e.message); }
    }
    if (swept) { try { drawSwept(); } catch (e) { swept = false; } }
    if (hasAtlas && atlasBuilt) {
      atlasNote.textContent = 'Conditions changed · recalculate the atlas to update the gas–solid comparison.';
      atlasBtn.textContent = 'Recalculate atlas';
    }
  });
  drawOptimality(TB.last());
  if (hasNew) {
    try { drawBoundary(); }
    catch (e) { figFail('figBoundary', 'figure failed: ' + e.message); }
  }
  if (hasAtlas) window.setTimeout(buildAtlas, 80);
}());
