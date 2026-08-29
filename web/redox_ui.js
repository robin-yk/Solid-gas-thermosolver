/* UI for the steady-state redox workspace.

   Three inputs and two pictures. The engine is web/redox.js, which carries
   the closed forms of solidgas/redox.py and nothing else, and the figures
   are the same two web/figures_redox.js draws for the manuscript - the
   page hands them a payload and gets back the SVG it downloads, so what
   is on screen and what is in the file cannot differ.

   The cycling card shares the workspace because it asks the second half
   of the same question: the steady state says where one exposure lands,
   and the cycle model says what a sequence of them does to a solid that
   remembers. */

(function () {
  'use strict';

  if (!window.Redox || !window.RedoxFigures) return;
  var R = window.Redox, CY = window.Cycling;
  var FIG = window.RedoxFigures, KIT = window.FigKit;
  var P = JSON.parse(document.getElementById('redox-data').textContent);

  var $ = function (id) { return document.getElementById(id); };
  var el = {};
  ['rxK', 'rxKs', 'rxE', 'rxEnote', 'rxSol', 'rxSolNote', 'rxMax',
   'rxEnv', 'rxWindow', 'rxSummary', 'rxKpis', 'rxPhaseNote',
   'cySummary', 'cyKpis', 'cyIn', 'cyOut', 'cyNotes', 'cyAdd', 'cyReset',
   'cyTotal', 'cySplit', 'cyPredKpis', 'cyPredNote'].forEach(function (id) {
     el[id] = $(id);
   });

  var E_RESOLVED = 1844.0;   /* 0.9 um rutile at 600 C, layer-resolved */

  /* -------------------------------------------------------- figures */

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

  /* ----------------------------------------------------- formatting */

  function num(node, fallback) {
    var v = parseFloat(node.value);
    return isFinite(v) ? v : fallback;
  }
  function sig(v, n) {
    if (!isFinite(v)) return '∞';
    if (v === 0) return '0';
    var a = Math.abs(v);
    if (a >= 1e4 || a < 1e-3) return v.toExponential((n || 3) - 1);
    return v.toPrecision(n || 3).replace(/\.?0+$/, '');
  }
  function thf(th) {
    if (th >= 0.001) return (th * 100).toFixed(2) + '%';
    if (th > 0) return th.toExponential(2);
    return '0';
  }
  function kpi(label, value, unit) {
    return '<div class="kpi"><div class="k">' + label + '</div>'
      + '<div class="v">' + value + '</div>'
      + '<div class="u">' + unit + '</div></div>';
  }

  /* ------------------------------------------------- steady state */

  function inputs() {
    var x = num(el.rxSol, 0.008);
    return { K: Math.max(0, num(el.rxK, 0.5)),
             E: Math.max(1, num(el.rxE, 1)),
             theta_sol: x / P.structure.oxygen_per_formula,
             theta_max: Math.min(0.99, Math.max(0.01, num(el.rxMax, 0.5))) };
  }

  var REGIME = {
    reduction_limited: 'Reduction-limited: the crossing is set by how fast '
      + 'the reductant can strip oxygen, and the solid stays oxidised.',
    balanced: 'Balanced: both half-reactions matter at the crossing. This '
      + 'is the band a carrier is worth running in.',
    oxidation_limited: 'Oxidation-limited: the crossing is set by how fast '
      + 'the oxidant can refill, and the solid stays reduced.'
  };

  function solve() {
    var q = inputs();
    var s = R.steady(q.K, { theta_sol: q.theta_sol,
                            theta_surface_max: q.theta_max,
                            enrichment: q.E,
                            bands: P.regime_thresholds });
    var w = s.window;

    el.rxSolNote.innerHTML = 'point-defect x in MO<sub>2−x</sub>, so '
      + 'θ<sub>sol</sub> = ' + sig(q.theta_sol, 3) + ' · '
      + P.structure.label;
    el.rxEnote.innerHTML = q.E === 1
      ? 'The uniform particle. Nothing says this is wrong; it says nothing '
        + 'at all, which is the problem.'
      : 'The face is ' + sig(q.E, 4) + '× as vacant as the particle '
        + 'average, so the solid holds ' + sig(q.E, 4) + '× less inventory '
        + 'to keep it there.';

    /* the verdict is about whether the carrier can be RUN here, which is
       a different question from whether the algebra has a root */
    var ok = s.usable;
    el.rxSummary.className = 'verdict ' + (ok ? 'none' : 'red');
    var head, why;
    if (s.no_steady_state) {
      head = 'No steady state at K = ' + sig(q.K, 3) + '.';
      why = 'The crossing wants θ(face) = ' + thf(s.theta_surface)
        + ', above the ' + thf(q.theta_max) + ' its surface phase pins it '
        + 'at. Past the ceiling the oxidation rate stops rising with the '
        + 'inventory, so the imbalance keeps its sign and the particle '
        + 'reduces until something outside this pair of rate laws stops '
        + 'it. That is an answer about the carrier, not a failure to '
        + 'converge.';
    } else if (s.beyond_solubility) {
      head = 'The crossing is outside the single-phase solid.';
      why = 'It sits at θ(particle) = ' + thf(s.theta)
        + ' (x = ' + sig(s.x_equiv, 3) + ' in MO<sub>2−x</sub>), above the '
        + 'declared solubility θ = ' + sig(q.theta_sol, 3) + '. The '
        + 'vacancies order into a second structure before the carrier gets '
        + 'there, so the number below is where an unshearing solid <i>would'
        + '</i> sit and the real one changes phase instead.';
    } else {
      head = 'Usable steady state at K = ' + sig(q.K, 3) + '.';
      why = 'θ(particle) = ' + thf(s.theta) + ' stays inside the '
        + 'single-phase range and θ(face) = ' + thf(s.theta_surface)
        + ' stays under the ceiling.';
    }
    el.rxSummary.innerHTML = '<b>' + head + '</b> ' + why
      + '<span class="more">' + REGIME[s.regime] + '</span>';

    el.rxKpis.innerHTML =
      kpi('θ on the reacting face', thf(s.theta_surface),
          'K/(1+K), exactly')
      + kpi('θ in the particle', thf(s.theta),
            'x = ' + sig(s.x_equiv, 3) + ' in MO<sub>2−x</sub>')
      + kpi('Turnover', sig(s.rate_over_k_ox, 3),
            'per site, in units of k(oxidation)')
      + kpi('Largest usable K', sig(w.K_max, 3),
            w.binding === 'second_phase'
              ? 'the solid changes phase first'
              : 'the face saturates first');

    var cross = R.crossoverEnrichment(q.theta_sol, q.theta_max);
    el.rxWindow.innerHTML =
      'At this enrichment the two limits are <b>K ≤ '
      + sig(w.K_second_phase, 3) + '</b> before the solid leaves its '
      + 'single-phase range and <b>K ≤ ' + sig(w.K_no_steady_state, 3)
      + '</b> before the face saturates, so the window is K ≤ '
      + sig(w.K_max, 3) + '. Which one binds swaps at E = '
      + sig(cross, 4) + ': below it a wider single-phase range buys room, '
      + 'above it only a higher face ceiling does. The uniform particle '
      + 'is held to K ≤ ' + sig(R.usableWindow(1, q.theta_sol,
                                               q.theta_max).K_max, 3)
      + ' and the layer-resolved one to K ≤ '
      + sig(R.usableWindow(E_RESOLVED, q.theta_sol,
                           q.theta_max).K_max, 3) + '.';

    el.rxEnv.innerHTML = 'At ' + P.conditions.T_C + ' °C, '
      + P.conditions.p_red_atm + ' atm reductant and '
      + P.conditions.p_ox_atm + ' atm oxidant &mdash; none of which enter '
      + 'the crossing. <b>K is the only material quantity in it.</b>'
      + '<span class="more">Temperature, pressure and both prefactors '
      + 'cancel out of θ* = K/(1+K) at unit site orders. They set how fast '
      + 'the carrier reaches the steady state, not where it is.</span>';

    try {
      mountFigure('figCrossing', 'redox-crossing',
        FIG.crossing({ crossing: R.crossingData(q.K, 201) }));
    } catch (e) { figFail('figCrossing', e.message); }

    try {
      var pd = R.phaseData({
        theta_sol: q.theta_sol, theta_surface_max: q.theta_max,
        points: [{ E: 1.0, label: 'uniform',
                   K: R.usableWindow(1.0, q.theta_sol, q.theta_max).K_max },
                 { E: E_RESOLVED, label: 'layer-resolved',
                   K: R.usableWindow(E_RESOLVED, q.theta_sol,
                                     q.theta_max).K_max }],
        current: { K: q.K, E: q.E, usable: ok } });
      mountFigure('figPhase', 'redox-phase', FIG.phase({ phase: pd }));
    } catch (e2) { figFail('figPhase', e2.message); }

    var uni = R.usableWindow(1.0, q.theta_sol, q.theta_max).K_max;
    var res = R.usableWindow(E_RESOLVED, q.theta_sol, q.theta_max).K_max;
    el.rxPhaseNote.innerHTML = 'Reading the diagram vertically is the '
      + 'point. The uniform model is the bottom edge and answers K ≤ '
      + sig(uni, 3) + '; the layer-resolved particle at E = '
      + E_RESOLVED + ' answers K ≤ ' + sig(res, 3) + ', a factor of '
      + sig(res / uni, 3) + '. The turnover at a given K is identical in '
      + 'both — any strictly increasing average→face mapping only '
      + 'reparameterises the same crossing. What the correction moves is '
      + 'the inventory the solid must hold to keep the face there, and so '
      + 'whether it survives being run at its own optimum.';
  }

  /* ------------------------------------------------------- cycling */

  var MEASURED = [{ created: 94.8, recovered: 81.5 },
                  { created: 61.0, recovered: 51.6 }];
  var rows = MEASURED.map(function (r) { return { c: r.created,
                                                  v: r.recovered }; });

  function drawInputs() {
    var h = '<tr><th>Cycle</th><th>V<sub>O</sub> made, μmol-O/g</th>'
      + '<th>O recovered, μmol-O/g</th><th>f = out / in</th><th></th></tr>';
    rows.forEach(function (r, i) {
      var f = r.c ? r.v / r.c : 0;
      h += '<tr><td>' + (i + 1) + '</td>'
        + '<td><input type="number" step="any" min="0" data-row="' + i
        + '" data-key="c" value="' + r.c + '"></td>'
        + '<td><input type="number" step="any" min="0" data-row="' + i
        + '" data-key="v" value="' + r.v + '"></td>'
        + '<td>' + (r.c ? f.toFixed(3) : '—') + '</td>'
        + '<td>' + (rows.length > 2
          ? '<button class="mini" data-drop="' + i + '">remove</button>'
          : '') + '</td></tr>';
    });
    el.cyIn.innerHTML = h;
  }

  function fitCycles() {
    var obs = rows.map(function (r) {
      return { created: r.c, recovered: r.v }; });
    var fit;
    try {
      fit = CY.fit(obs);
    } catch (e) {
      el.cySummary.className = 'verdict red';
      el.cySummary.textContent = e.message;
      el.cyKpis.innerHTML = '';
      el.cyOut.innerHTML = '';
      el.cyNotes.innerHTML = '';
      el.cyPredKpis.innerHTML = '';
      el.cyPredNote.textContent = '';
      return;
    }

    var made = obs.reduce(function (a, o) { return a + o.created; }, 0);
    var got = obs.reduce(function (a, o) { return a + o.recovered; }, 0);

    el.cySummary.className = 'verdict ' + (fit.lock_at_bound ? 'red' : 'none');
    el.cySummary.innerHTML = '<b>f = ' + fit.f_rec.toFixed(3)
      + ', lock = ' + fit.lock.toFixed(3) + '.</b> '
      + 'One exposure recovers ' + (fit.f_rec * 100).toFixed(1) + '% of '
      + 'what it can still reach, and ' + (fit.lock * 100).toFixed(1)
      + '% of what it misses stops responding to later ones. Residual '
      + 'after ' + obs.length + ' cycles: ' + (made - got).toFixed(1)
      + ' μmol-O/g of the ' + made.toFixed(1) + ' made.'
      + '<span class="more">lock = 0 is the memoryless model exactly, so '
      + 'the memory is a switch and not a new assumption. Two cycles '
      + 'constrain two parameters and no more: an unlocking rate is the '
      + 'obvious third and is left out rather than fitted to data that '
      + 'cannot see it.</span>';

    el.cyKpis.innerHTML =
      kpi('f, per exposure', fit.f_rec.toFixed(3),
          'of the reachable pool')
      + kpi('lock', fit.lock.toFixed(3),
            'of the miss, per exposure')
      + kpi('Residual', (made - got).toFixed(1),
            'μmol-O/g still in the solid')
      + kpi('Fit RMS', fit.rms_umol_g.toFixed(2),
            'μmol-O/g over ' + obs.length + ' cycles');

    var h = '<tr><th>Cycle</th><th>reachable before</th><th>locked before</th>'
      + '<th>recovered, model</th><th>recovered, measured</th>'
      + '<th>residual after</th></tr>';
    fit.predicted.forEach(function (r, i) {
      h += '<tr><td>' + r.cycle + '</td><td>'
        + r.pre_accessible_umol_g.toFixed(1) + '</td><td>'
        + r.pre_locked_umol_g.toFixed(1) + '</td><td>'
        + r.recovered_umol_g.toFixed(1) + '</td><td>'
        + obs[i].recovered.toFixed(1) + '</td><td>'
        + r.residual_umol_g.toFixed(1) + '</td></tr>';
    });
    el.cyOut.innerHTML = h;

    var n = '';
    fit.notes.forEach(function (t) {
      n += '<div class="note">' + t + '</div>';
    });
    if (fit.lock_at_bound && !fit.notes.length) {
      n += '<div class="note">Locking is at its bound: the model cannot '
        + 'predict less recovery than it already does.</div>';
    }
    el.cyNotes.innerHTML = n;

    if (!el.cyTotal.value) el.cyTotal.value = made.toFixed(1);
    predict(fit);
  }

  function predict(fit) {
    var totalVo = num(el.cyTotal, 100);
    var split = Math.max(1, Math.min(12, Math.round(num(el.cySplit, 2))));
    var t = CY.twoSamplePrediction(totalVo, fit.f_rec, fit.lock, split);
    el.cyPredKpis.innerHTML =
      kpi('One long reduction', t.fresh.recovered_umol_g.toFixed(1),
          'μmol-O/g recovered')
      + kpi('Built in ' + split + ' cycles',
            t.cycled.recovered_umol_g.toFixed(1), 'μmol-O/g recovered')
      + kpi('Gap', t.recovery_gap_umol_g.toFixed(1),
            'μmol-O/g, ' + ((1 - t.recovery_ratio) * 100).toFixed(1)
            + '% less')
      + kpi('Locked in the cycled one', t.cycled.locked.toFixed(1),
            'μmol-O/g out of reach');
    el.cyPredNote.innerHTML = 'A memoryless inventory predicts a gap of '
      + '<b>zero</b> at any f, because the total is its entire state. So '
      + 'the measurement is a yes/no about the model, not a fit to it — '
      + 'and if the gap is there, its size <i>is</i> the lock parameter. '
      + 'What the gap cannot say is <i>why</i>: migration beyond reach, '
      + 'aggregation into an ordered defect and a local relaxation all '
      + 'produce the same number, and cycle integrals alone do not '
      + 'separate them.';
  }

  /* -------------------------------------------------------- wiring */

  /* K spans eight decades, so the slider is in log10 and the box is the
     number itself; each writes the other and neither is the master. */
  function setK(v) {
    el.rxK.value = parseFloat(v.toPrecision(4));
    el.rxKs.value = Math.log10(Math.max(1e-4, Math.min(1e4, v)));
  }
  el.rxKs.addEventListener('input', function () {
    setK(Math.pow(10, parseFloat(el.rxKs.value)));
    solve();
  });
  el.rxK.addEventListener('input', function () {
    var v = num(el.rxK, 0.5);
    if (v > 0) el.rxKs.value = Math.log10(Math.max(1e-4, Math.min(1e4, v)));
    solve();
  });
  [el.rxE, el.rxSol, el.rxMax].forEach(function (n) {
    n.addEventListener('input', solve);
  });
  Array.prototype.forEach.call(
    document.querySelectorAll('#ws-redox .mini[data-e]'), function (b) {
      b.addEventListener('click', function () {
        el.rxE.value = b.dataset.e;
        solve();
      });
    });

  el.cyIn.addEventListener('input', function (ev) {
    var t = ev.target;
    if (t.dataset.row == null) return;
    var v = parseFloat(t.value);
    rows[+t.dataset.row][t.dataset.key] = isFinite(v) ? v : 0;
    fitCycles();
    var f = rows[+t.dataset.row];
    var cell = t.closest('tr').children[3];
    if (cell) cell.textContent = f.c ? (f.v / f.c).toFixed(3) : '—';
  });
  el.cyIn.addEventListener('click', function (ev) {
    var i = ev.target && ev.target.getAttribute('data-drop');
    if (i == null) return;
    rows.splice(+i, 1);
    drawInputs();
    fitCycles();
  });
  el.cyAdd.addEventListener('click', function () {
    var last = rows[rows.length - 1];
    rows.push({ c: +(last.c * 0.75).toFixed(1),
                v: +(last.v * 0.75).toFixed(1) });
    drawInputs();
    fitCycles();
  });
  el.cyReset.addEventListener('click', function () {
    rows = MEASURED.map(function (r) {
      return { c: r.created, v: r.recovered }; });
    el.cyTotal.value = '';
    drawInputs();
    fitCycles();
  });
  [el.cyTotal, el.cySplit].forEach(function (n) {
    n.addEventListener('input', fitCycles);
  });

  var rxTabs = document.querySelectorAll('#rxTabs button');
  Array.prototype.forEach.call(rxTabs, function (btn) {
    btn.addEventListener('click', function () {
      Array.prototype.forEach.call(rxTabs, function (b) {
        var on = b === btn;
        b.setAttribute('aria-selected', on ? 'true' : 'false');
        var panel = $(b.dataset.panel);
        if (panel) panel.hidden = !on;
      });
    });
  });

  /* ------------------------------------------- workspace switching */

  var started = false;
  var tabs = document.querySelectorAll('#wstabs .wstab');
  Array.prototype.forEach.call(tabs, function (btn) {
    btn.addEventListener('click', function () {
      if (btn.dataset.ws !== 'ws-redox' || started) return;
      started = true;
      setK(num(el.rxK, 0.5));
      solve();
      drawInputs();
      fitCycles();
    });
  });
}());
