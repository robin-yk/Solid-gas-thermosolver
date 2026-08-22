/* UI for the Defect stat mech workspace.

   The Monte Carlo isotherm runs in a Worker assembled from the inlined
   engine source (no request leaves the page); everything downstream of the
   isotherm - inventory, particle size, geometry mode - re-solves instantly
   on the main thread because the isotherm depends only on temperature and
   energetics. Before the first run finishes, results render from the exact
   ideal-surface isotherm and are labelled as a preview. */

(function () {
  'use strict';

  if (!window.StatMech) return;
  var SM = window.StatMech;
  var DATA = JSON.parse(document.getElementById('statmech-data').textContent);

  var $ = function (id) { return document.getElementById(id); };
  var el = {};
  ['smT', 'smVO', 'smGeoMode', 'smGeoVal', 'smGeoUnit', 'smDerived',
   'smPreset', 'smEpsSS', 'smEpsB', 'smV1', 'smV2', 'smV3', 'smQuality',
   'smSeed', 'smRun', 'smStatus', 'smSummary', 'smKpis', 'smChart',
   'smTable', 'smSpacing', 'smNotes',
   'smEaS', 'smEaSS', 'smEaB', 'smAS', 'smASS', 'smAB', 'smKinT', 'smExp',
   'smPCO2', 'smAccKpis', 'smSweep', 'smTip', 'smSweepLegend', 'smAccNote',
   'smCoupledNote', 'smEnv', 'smParticle', 'smPvAcc', 'smPvDl',
   'smPvLegend', 'smPvNote'].forEach(function (id) { el[id] = $(id); });

  /* ------------------------------------------- workspace switching */

  var started = false;
  var tabs = document.querySelectorAll('#wstabs .wstab');
  Array.prototype.forEach.call(tabs, function (btn) {
    btn.addEventListener('click', function () {
      Array.prototype.forEach.call(tabs, function (b) {
        b.classList.toggle('active', b === btn);
      });
      ['ws-thermo', 'ws-defect'].forEach(function (ws) {
        $(ws).style.display = ws === btn.dataset.ws ? '' : 'none';
      });
      if (btn.dataset.ws === 'ws-defect') {
        if (!started) {
          started = true;
          solve();
          launch();
        }
        updateEnv();      // the charge may have been edited over there
      }
    });
  });

  /* ------------------------------------------------ input reading */

  function num(node, fallback) {
    var v = parseFloat(node.value);
    return isFinite(v) ? v : fallback;
  }
  function currentT() { return num(el.smT, DATA.defaults.T_C); }
  function currentEps() {
    return [0.0, num(el.smEpsSS, 0.59), num(el.smEpsB, 1.31)];
  }
  function currentOrdering() {
    return { pair_eV: {
      along_row: [num(el.smV1, 0.40), num(el.smV2, 0.20), num(el.smV3, 0.05)],
      cross_row: DATA.surface_ordering.pair_eV.cross_row } };
  }
  function currentGeom() {
    if (el.smGeoMode.value === 'bet') {
      return SM.geometry(DATA, { bet_m2_g: num(el.smGeoVal, 1.57) });
    }
    return SM.geometry(DATA, {
      d_um: num(el.smGeoVal, DATA.defaults.particle_diameter_um) });
  }
  function currentKin() {
    return { Ea_eV: { surface: num(el.smEaS, 0.45),
                      subsurface: num(el.smEaSS, 1.25),
                      bulk: num(el.smEaB, 2.10) },
             A_eff_per_s_atm: { surface: num(el.smAS, 1e13),
                                subsurface: num(el.smASS, 3e5),
                                bulk: num(el.smAB, 1e5) },
             T_reox_C: num(el.smKinT, 600),
             exposure_s: num(el.smExp, 300),
             p_CO2_atm: num(el.smPCO2, 0.2) };
  }
  function isoSignature() {
    return JSON.stringify([currentT(), currentEps(), currentOrdering(),
                           el.smQuality.value, num(el.smSeed, 12345)]);
  }

  /* --------------------------------------------------- MC run flow */

  var iso = null, isoKey = '', isoInfo = '';
  var worker = null, runToken = 0, running = false, pending = null;

  function setStatus(t) { el.smStatus.textContent = t; }

  function buildWorker() {
    var src = document.getElementById('sm-engine').textContent + '\n'
      + document.getElementById('sm-worker-src').textContent;
    var blob = new Blob([src], { type: 'application/javascript' });
    return new Worker(URL.createObjectURL(blob));
  }

  function stopRun(msg) {
    if (worker) { worker.terminate(); worker = null; }
    running = false;
    el.smRun.textContent = 'Run model';
    if (msg) setStatus(msg);
  }

  function finishRun(sig, pts, ms) {
    running = false;
    el.smRun.textContent = 'Run model';
    iso = pts;
    isoKey = sig;
    isoInfo = 'isotherm: ' + pts.length + ' points, '
      + (ms / 1000).toFixed(1) + ' s';
    setStatus(isoInfo);
    solve();
  }

  function launch() {
    if (running) stopRun(null);
    var sig = isoSignature();
    var token = ++runToken;
    var t0 = Date.now();
    var msg = { cmd: 'run', data: DATA, T_C: currentT(),
                seed: num(el.smSeed, 12345),
                quality: parseFloat(el.smQuality.value) || 1.0,
                eps: currentEps(), ordering: currentOrdering(),
                mc: null, token: token };
    running = true;
    el.smRun.textContent = 'Cancel run';
    setStatus('sampling the surface isotherm…');
    try {
      worker = buildWorker();
      worker.onmessage = function (ev) {
        var m = ev.data;
        if (m.token !== token) return;
        if (m.kind === 'progress') {
          setStatus('sampling the surface isotherm… point '
            + m.done + ' of ' + m.total);
        } else if (m.kind === 'done') {
          finishRun(sig, m.isotherm, Date.now() - t0);
        } else if (m.kind === 'error') {
          stopRun('run failed: ' + m.message);
        }
      };
      worker.postMessage(msg);
    } catch (e) {
      setStatus('sampling on the main thread (workers unavailable)…');
      setTimeout(function () {
        var pts = SM.isothermScan(DATA, msg.T_C, {
          seed: msg.seed, quality: msg.quality, eps: msg.eps,
          ordering: msg.ordering });
        finishRun(sig, pts, Date.now() - t0);
      }, 30);
    }
  }

  function scheduleRun() {
    clearTimeout(pending);
    pending = setTimeout(function () {
      if (!started) return;
      if (isoKey !== isoSignature()) launch();
    }, 700);
  }

  el.smRun.addEventListener('click', function () {
    if (running) stopRun('run cancelled'); else launch();
  });

  /* ------------------------------------------- thermodynamic bridge */

  var env = null, envKey = '', envTimer = null;

  function scheduleEnv() {
    clearTimeout(envTimer);
    envTimer = setTimeout(updateEnv, 400);
  }

  function updateEnv() {
    if (!el.smEnv) return;
    var TB = window.ThermoBridge;
    var tC = currentT();
    if (!TB) {
      env = null;
      el.smEnv.innerHTML = '<div class="readout">The equilibrium workspace '
        + 'is unavailable on this page.</div>';
      return;
    }
    var desc = TB.describe();
    if (!desc) {
      env = null;
      el.smEnv.innerHTML = '<div class="readout">The equilibrium '
        + "workspace's gas charge sums to zero.</div>";
      return;
    }
    var key = tC + '|' + desc;
    if (key === envKey && env) { renderEnv(); return; }
    envKey = key;
    if (tC < 300 || tC > 1650) {
      env = { outOfRange: true, T_C: tC, desc: desc, thermoT: TB.T_C() };
      renderEnv();
      solve();
      return;
    }
    try {
      var R = TB.solveAt(tC);
      if (!R) throw new Error('no solution returned');
      var v = SM.phaseValidity(R.active_condensed_phases,
                               R.inactive_phase_reduced_costs,
                               TB.isReduced());
      env = { T_C: tC, desc: desc, mu_O: R.lambda_kJ_per_mol.O,
              active: R.active_condensed_phases, v: v, thermoT: TB.T_C() };
    } catch (e) {
      env = { error: String(e && e.message ? e.message : e), T_C: tC,
              desc: desc, thermoT: TB.T_C() };
    }
    renderEnv();
    solve();                       // the validity note lives in the results
  }

  function subHtml(s) {
    return String(s).replace(/(\d+)/g, '<sub>$1</sub>');
  }

  function renderEnv() {
    var box = el.smEnv;
    if (!env) return;
    var h = '';
    if (env.outOfRange) {
      h = '<div class="note calm">' + env.T_C + ' °C is outside the '
        + 'equilibrium data range (300–1650 °C); no verdict at this '
        + 'temperature.</div>';
    } else if (env.error) {
      h = '<div class="note">Equilibrium solve failed: ' + esc(env.error)
        + '</div>';
    } else {
      var v = env.v;
      h = '<div style="font-size:15px;line-height:1.8">'
        + '<div><span class="k2">Charge</span> ' + env.desc + '</div>'
        + '<div><span class="k2">Oxygen potential</span> μ<sub>O</sub> = '
        + env.mu_O.toFixed(2) + ' kJ/mol-O at ' + env.T_C + ' °C</div>'
        + '<div><span class="k2">Equilibrium phases</span> '
        + env.active.map(subHtml).join(' + ') + '</div>';
      if (v.tier === 'stable') {
        h += '<div><span class="k2">Nearest reduced</span> '
          + subHtml(v.nearest_phase) + ' at r = +'
          + v.margin_per_mol_O_kJ.toFixed(1) + ' kJ/mol-O</div>';
      }
      h += '</div>';
      if (v.tier === 'stable') {
        h += '<div class="note calm" style="margin:10px 0 4px"><b>Rutile '
          + 'stable</b> under this charge at ' + env.T_C + ' °C, with a +'
          + v.margin_per_mol_O_kJ.toFixed(1) + ' kJ/mol-O margin to '
          + subHtml(v.nearest_phase) + '.</div>';
      } else if (v.tier === 'boundary') {
        h += '<div class="note" style="margin:10px 0 4px"><b>On the '
          + 'reduction boundary.</b> TiO<sub>2</sub> coexists with '
          + v.reduced_active.map(subHtml).join(' + ')
          + ' at equilibrium; reduction begins at this condition.</div>';
      } else {
        h += '<div class="note" style="margin:10px 0 4px;border-left-color:'
          + 'var(--red)"><b>Single-phase rutile model invalid at '
          + 'equilibrium.</b> Thermodynamics predicts '
          + v.reduced_active.map(subHtml).join(' + ')
          + '. Rutile is absent under this charge at ' + env.T_C + ' °C.</div>';
      }
    }
    h += '<div style="display:flex;gap:8px;margin-top:8px;flex-wrap:wrap">'
      + '<button type="button" class="mini" id="smEnvOpen">Edit charge in '
      + 'equilibrium workspace</button>';
    if (env.thermoT != null && Math.abs(env.thermoT - env.T_C) > 1e-9) {
      h += '<button type="button" class="mini" id="smEnvMatch">Use '
        + 'equilibrium T (' + env.thermoT + ' °C)</button>';
    }
    h += '</div>';
    box.innerHTML = h;
    var open = $('smEnvOpen');
    if (open) {
      open.onclick = function () {
        var b = document.querySelector('#wstabs .wstab[data-ws="ws-thermo"]');
        if (b) b.click();
      };
    }
    var match = $('smEnvMatch');
    if (match) {
      match.onclick = function () {
        el.smT.value = env.thermoT;
        solve();
        scheduleRun();
        scheduleEnv();
      };
    }
  }

  /* ------------------------------------------------ solve + render */

  function solve() {
    var geom = currentGeom();
    var tC = currentT();
    var vo = num(el.smVO, DATA.defaults.VO_total_umol_g);
    var eps = currentEps();
    var kt = SM.KB_EV * (tC + 273.15);
    var fresh = !!iso && isoKey === isoSignature();
    var fn = fresh ? SM.thetaSInterp(iso, eps[0], kt) : undefined;
    var out = SM.distribute(DATA, tC, vo, geom,
                            { eps: eps, theta_s_fn: fn });
    out.geometry = geom;
    out.spacing = fresh ? SM.spacingSummary(iso, out.theta.surface) : null;
    render(out, fresh, { tC: tC, eps: eps, fn: fn, geom: geom,
                         kin: currentKin() });
  }

  function f2(v) {
    if (v >= 100) return v.toFixed(0);
    if (v >= 10) return v.toFixed(1);
    return v.toFixed(2);
  }
  function pctf(f) {
    var p = f * 100;
    if (p >= 10) return p.toFixed(1) + '%';
    if (p >= 0.1) return p.toFixed(2) + '%';
    if (p > 0) return p.toExponential(1) + '%';
    return '0%';
  }
  function thf(th) {
    if (th >= 0.001) return (th * 100).toFixed(1) + '%';
    if (th > 0) return th.toExponential(2);
    return '0';
  }
  function esc(s) {
    return String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;');
  }

  var CLASSES = [
    { key: 'surface', label: 'Surface', varname: '--cSurf' },
    { key: 'subsurface', label: 'Subsurface', varname: '--cSub' },
    { key: 'bulk', label: 'Bulk', varname: '--cBulk' }];

  function render(out, fresh, ctx) {
    var geom = out.geometry;
    var acc = SM.accessibility(DATA, out.umol_g, ctx.kin);
    el.smDerived.textContent = 'x in TiO2-x = ' + out.x_TiO2mx.toFixed(4)
      + ' · Ti³⁺/Ti = ' + (out.Ti3_frac * 100).toFixed(2)
      + '% · area ' + geom.area_m2_g.toFixed(2)
      + ' m²/g · surface capacity ' + geom.N_s.toFixed(1)
      + ' μmol-O/g';

    var bad = out.warnings.some(function (w) {
      return w.kind === 'x_max' || w.kind === 'unresolved';
    });
    var parts = CLASSES.map(function (c) {
      return c.label.toLowerCase() + ' ' + f2(out.umol_g[c.key]) + ' ('
        + pctf(out.fractions[c.key]) + ')';
    });
    el.smSummary.className = 'verdict ' + (bad ? 'red' : 'none');
    el.smSummary.innerHTML = '<b>' + f2(out.VO_total_umol_g)
      + ' μmol-O/g at ' + out.T_C + ' °C:</b> ' + parts.join(', ')
      + ' μmol-O/g. Vacancy chemical potential μ<sub>V</sub> = '
      + out.mu_V_eV.toFixed(3) + ' eV.'
      + (fresh ? '' : ' <span style="color:var(--muted)">Ideal-surface '
         + 'preview. The interacting isotherm is being sampled.</span>');

    var k = '';
    CLASSES.forEach(function (c) {
      k += '<div class="kpi"><div class="k"><span class="swatch" style="'
        + 'background:var(' + c.varname + ')"></span>' + c.label
        + ' V<sub>O</sub></div><div class="v">' + f2(out.umol_g[c.key])
        + '</div><div class="u">μmol-O/g · '
        + pctf(out.fractions[c.key]) + ' of inventory</div></div>';
    });
    k += '<div class="kpi"><div class="k">CO<sub>2</sub>-recoverable</div>'
      + '<div class="v">' + pctf(acc.f_recoverable) + '</div>'
      + '<div class="u">' + f2(acc.accessible_umol_g)
      + ' μmol-O/g accessible</div></div>';
    k += '<div class="kpi"><div class="k">μ<sub>V</sub> (vacancy)</div>'
      + '<div class="v">' + out.mu_V_eV.toFixed(3) + '</div>'
      + '<div class="u">eV, common to all classes</div></div>';
    el.smKpis.innerHTML = k;

    drawChart(out);

    var h = '<tr><th>Site class</th><th>θ occupancy</th>'
      + '<th>μmol-O/g</th><th>share of inventory</th></tr>';
    CLASSES.forEach(function (c) {
      h += '<tr><td><span class="swatch" style="background:var(' + c.varname
        + ')"></span>' + c.label + '</td><td>' + thf(out.theta[c.key])
        + '</td><td>' + f2(out.umol_g[c.key]) + '</td><td>'
        + pctf(out.fractions[c.key]) + '</td></tr>';
    });
    h += '<tr><td><b>Total</b></td><td></td><td><b>'
      + f2(out.matched_umol_g) + '</b></td><td><b>100%</b></td></tr>';
    el.smTable.innerHTML = h;

    if (out.spacing) {
      el.smSpacing.textContent = 'Surface pair statistics at θs = '
        + out.spacing.at_theta_s.toFixed(3) + ': modal spacing '
        + out.spacing.modal_gap_sites + ' bridging sites, P(spacing ≤ 2) '
        + '= ' + (out.spacing.P_gap_le_2 * 100).toFixed(1) + '%. '
        + (out.surface_reconstruction_regime
          ? 'The 4–5-site literature comparison applies near the 17–20% '
            + 'reconstruction-onset interval; this higher-coverage (1×1) '
            + 'snapshot is outside that validation range.'
          : 'Near the reconstruction onset, the literature reports modal '
            + 'spacing 4–5 with close pairs suppressed.');
    } else {
      el.smSpacing.textContent = '';
    }

    var notes = '';
    if (env && env.v && env.T_C === out.T_C) {
      if (env.v.tier === 'invalid') {
        notes += '<div class="note" style="border-left-color:var(--red)">'
          + '<b>Phase validity.</b> Under the equilibrium workspace&rsquo;s '
          + 'charge, thermodynamics predicts '
          + env.v.reduced_active.map(subHtml).join(' + ')
          + '. Rutile is absent at ' + env.T_C + ' °C. The single-phase rutile '
          + 'defect model does not hold at equilibrium there; a fixed '
          + 'measured inventory in rutile is metastable at best (Oxygen '
          + 'environment card).</div>';
      } else if (env.v.tier === 'boundary') {
        notes += '<div class="note"><b>Phase validity.</b> The equilibrium '
          + 'charge sits on the TiO<sub>2</sub> / '
          + env.v.reduced_active.map(subHtml).join(' + ')
          + ' boundary at ' + env.T_C + ' °C. Reduction begins at this '
          + 'condition (Oxygen environment card).</div>';
      }
    }
    out.warnings.forEach(function (w) {
      var calm = w.kind === 'surface_reconstruction' || w.kind === 'x_near';
      notes += '<div class="note' + (calm ? ' calm' : '') + '">'
        + (w.kind === 'x_max' ? '<b>Model range.</b> ' : '')
        + esc(w.text) + '</div>';
    });
    el.smNotes.innerHTML = notes;

    renderAccess(out, fresh, ctx, acc);
  }

  /* ------------------------------------------- CO2 accessibility */

  function renderAccess(out, fresh, ctx, acc) {
    var di = SM.dielectric(DATA, out.VO_total_umol_g);
    var k = '<div class="kpi"><div class="k">CO<sub>2</sub>-accessible '
      + 'V<sub>O</sub></div><div class="v">' + f2(acc.accessible_umol_g)
      + '</div><div class="u">μmol-O/g of ' + f2(out.VO_total_umol_g)
      + ' total</div></div>';
    k += '<div class="kpi"><div class="k">Recoverable fraction</div>'
      + '<div class="v">' + pctf(acc.f_recoverable) + '</div>'
      + '<div class="u">V<sub>O,acc</sub> / V<sub>O,tot</sub> at '
      + acc.params.T_reox_C + ' °C, ' + acc.params.exposure_s + ' s</div></div>';
    var dc = DATA.dielectric_correlation;
    if (di) {
      k += '<div class="kpi"><div class="k">ε″ (this work)</div>'
        + '<div class="v">'
        + (di.eps2 >= 10 ? di.eps2.toFixed(1) : di.eps2.toFixed(3))
        + '</div><div class="u">'
        + (di.tan_delta != null
           ? 'tan δ = ' + di.tan_delta.toFixed(3) + ' · ' : '')
        + (di.extrapolated
           ? 'outside the ' + dc.fit_range_umol_g[0] + '–'
             + dc.fit_range_umol_g[1] + ' fit range'
           : 'power-law fit, R² = ' + dc.fit_stats.r_squared.toFixed(2))
        + '</div></div>';
    } else {
      k += '<div class="kpi"><div class="k">ε″ (this work)</div>'
        + '<div class="v">N/A</div><div class="u">needs V<sub>O</sub> &gt; 0'
        + '</div></div>';
    }
    el.smAccKpis.innerHTML = k;

    el.smAccNote.textContent = 'Refill probabilities at '
      + acc.params.T_reox_C + ' °C, ' + acc.params.exposure_s + ' s, '
      + acc.params.p_CO2_atm + ' atm: surface ' + probf(acc.P.surface)
      + ', subsurface ' + probf(acc.P.subsurface) + ', bulk '
      + probf(acc.P.bulk) + '. Placeholder kinetics. The prefactors are '
      + 'effective values, not literature numbers.';

    el.smCoupledNote.innerHTML = '<div class="readout">Dielectric '
      + 'correlation (measured in this work): ε″ = 10<sup>'
      + dc.log10_intercept.toFixed(2) + '</sup> · V<sub>O</sub><sup>'
      + dc.exponent.toFixed(2) + '</sup>. Linear fit of log ε″ against '
      + 'log V<sub>O</sub>, N = ' + dc.fit_stats.n_points + ', Pearson r = '
      + dc.fit_stats.pearson_r.toFixed(3) + ', R² = '
      + dc.fit_stats.r_squared.toFixed(3) + ' (OriginLab, no weighting), '
      + 'fitted over ~' + dc.fit_range_umol_g[0] + '–'
      + dc.fit_range_umol_g[1] + ' μmol-O/g; below that range the value is '
      + 'an extrapolation. tan δ appears when ε′ is set in '
      + '<code>data/rutile_dft.json</code>.</div>';

    var vo = out.VO_total_umol_g;
    var xmax = 150;
    while (xmax < vo * 1.25) xmax += 50;
    var voList = [];
    for (var i = 0; i <= 60; i++) voList.push(xmax * i / 60);
    var rows = SM.sweep(DATA, ctx.tC, ctx.geom, voList,
                        { theta_s_fn: ctx.fn, eps: ctx.eps, kin: ctx.kin });
    drawSweep(rows, vo, xmax);

    renderParticle(out, acc, fresh);
  }

  function probf(v) {
    if (v >= 0.9995) return '1.000';
    if (v >= 0.001) return v.toFixed(3);
    if (v > 0) return v.toExponential(1);
    return '0';
  }

  /* ------------------------------------------------------- chart */

  var SVG = 'http://www.w3.org/2000/svg';

  function sv(tag, attrs, parent) {
    var n = document.createElementNS(SVG, tag);
    for (var a in attrs) n.setAttribute(a, attrs[a]);
    if (parent) parent.appendChild(n);
    return n;
  }

  function drawChart(out) {
    var svg = el.smChart;
    while (svg.firstChild) svg.removeChild(svg.firstChild);
    var W = 560, H = 300, L = 52, R = 12, T = 26, B = 48;
    var pw = W - L - R, ph = H - T - B;
    var fmax = 0;
    CLASSES.forEach(function (c) {
      if (out.fractions[c.key] > fmax) fmax = out.fractions[c.key];
    });
    // round tick steps: the axis reads 0/20/40/60/80%, never 17.5%
    var steps = [0.025, 0.05, 0.1, 0.2, 0.25];
    var step = steps[steps.length - 1];
    for (var si = 0; si < steps.length; si++) {
      if (4 * steps[si] >= Math.min(1, fmax * 1.1)) { step = steps[si]; break; }
    }
    var ymax = 4 * step;
    for (var i = 0; i <= 4; i++) {
      var fv = step * i;
      var y = T + ph - ph * fv / ymax;
      sv('line', { x1: L, x2: L + pw, y1: y, y2: y,
                   style: 'stroke:var(--line);stroke-width:1' }, svg);
      var lab = sv('text', { x: L - 8, y: y + 4, 'text-anchor': 'end',
        style: 'fill:var(--muted);font-size:15px' }, svg);
      lab.textContent = (Math.round(fv * 1000) / 10) + '%';
    }
    var bw = pw / 3;
    CLASSES.forEach(function (c, idx) {
      var f = out.fractions[c.key];
      var barW = Math.min(96, bw * 0.52);
      var x = L + bw * idx + (bw - barW) / 2;
      var hgt = ph * f / ymax;
      var y = T + ph - hgt;
      var r = Math.min(4, hgt);
      if (hgt > 0) {
        var d = 'M' + x + ',' + (T + ph)
          + ' L' + x + ',' + (y + r)
          + ' Q' + x + ',' + y + ' ' + (x + r) + ',' + y
          + ' L' + (x + barW - r) + ',' + y
          + ' Q' + (x + barW) + ',' + y + ' ' + (x + barW) + ',' + (y + r)
          + ' L' + (x + barW) + ',' + (T + ph) + ' Z';
        var p = sv('path', { d: d, style: 'fill:var(' + c.varname + ')' }, svg);
        var tt = document.createElementNS(SVG, 'title');
        tt.textContent = c.label + ': ' + f2(out.umol_g[c.key])
          + ' μmol-O/g (' + pctf(f) + ' of inventory), θ = '
          + thf(out.theta[c.key]);
        p.appendChild(tt);
      }
      var vt = sv('text', { x: x + barW / 2, y: y - 7, 'text-anchor': 'middle',
        style: 'fill:var(--ink);font-size:15px;font-weight:700' }, svg);
      vt.textContent = pctf(f);
      var nt = sv('text', { x: x + barW / 2, y: T + ph + 17,
        'text-anchor': 'middle',
        style: 'fill:var(--navy);font-size:14px;font-weight:700' }, svg);
      nt.textContent = c.label;
      var ut = sv('text', { x: x + barW / 2, y: T + ph + 33,
        'text-anchor': 'middle',
        style: 'fill:var(--muted);font-size:12.5px' }, svg);
      ut.textContent = f2(out.umol_g[c.key]) + ' μmol-O/g';
    });
    sv('line', { x1: L, x2: L + pw, y1: T + ph, y2: T + ph,
                 style: 'stroke:var(--muted);stroke-width:1' }, svg);
  }

  /* ----------------------------------------------------- sweep chart */

  var SWEEP = [
    { key: 'surface', label: 'Surface', varname: '--cSurf', dash: null },
    { key: 'subsurface', label: 'Subsurface', varname: '--cSub', dash: null },
    { key: 'bulk', label: 'Bulk', varname: '--cBulk', dash: null },
    { key: 'accessible', label: 'CO₂-accessible', varname: '--ink',
      dash: '7 5' }];

  var sweepState = null;

  function drawSweep(rows, currentVO, xmax) {
    var svg = el.smSweep;
    while (svg.firstChild) svg.removeChild(svg.firstChild);
    var W = 720, H = 380, L = 56, R = 116, T = 26, B = 52;
    var pw = W - L - R, ph = H - T - B;
    var vmax = 0;
    rows.forEach(function (r) {
      SWEEP.forEach(function (s) { if (r[s.key] > vmax) vmax = r[s.key]; });
    });
    var steps = [1, 2, 5, 10, 20, 25, 50];
    var step = steps[steps.length - 1];
    for (var si = 0; si < steps.length; si++) {
      if (4 * steps[si] >= vmax * 1.06) { step = steps[si]; break; }
    }
    var ymax = 4 * step;
    function X(v) { return L + pw * v / xmax; }
    function Y(v) { return T + ph - ph * v / ymax; }
    var i, yv, y;
    for (i = 0; i <= 4; i++) {
      yv = step * i;
      y = Y(yv);
      sv('line', { x1: L, x2: L + pw, y1: y, y2: y,
                   style: 'stroke:var(--line);stroke-width:1' }, svg);
      var lb = sv('text', { x: L - 8, y: y + 4, 'text-anchor': 'end',
        style: 'fill:var(--muted);font-size:15px' }, svg);
      lb.textContent = yv;
    }
    var yt = sv('text', { x: 6, y: 14,
      style: 'fill:var(--muted);font-size:15px' }, svg);
    yt.textContent = 'μmol-O/g';
    for (i = 0; i <= 5; i++) {
      var xv = xmax * i / 5;
      var xl = sv('text', { x: X(xv), y: T + ph + 18, 'text-anchor': 'middle',
        style: 'fill:var(--muted);font-size:15px' }, svg);
      xl.textContent = Math.round(xv);
    }
    var xt = sv('text', { x: L + pw / 2, y: T + ph + 38,
      'text-anchor': 'middle',
      style: 'fill:var(--muted);font-size:13.5px' }, svg);
    xt.textContent = 'total Vₒ (μmol-O/g)';

    SWEEP.forEach(function (s) {
      var d = '';
      rows.forEach(function (r, j) {
        d += (j ? 'L' : 'M') + X(r.VO_total_umol_g).toFixed(2) + ','
          + Y(r[s.key]).toFixed(2);
      });
      var p = sv('path', { d: d, fill: 'none',
        style: 'stroke:var(' + s.varname + ');stroke-width:2' }, svg);
      if (s.dash) p.setAttribute('stroke-dasharray', s.dash);
    });

    var last = rows[rows.length - 1];
    var ends = SWEEP.map(function (s) {
      return { label: s.label, y: Y(last[s.key]) };
    }).sort(function (a, b) { return a.y - b.y; });
    for (i = 1; i < ends.length; i++) {
      if (ends[i].y - ends[i - 1].y < 14) ends[i].y = ends[i - 1].y + 14;
    }
    ends.forEach(function (e) {
      var t2 = sv('text', { x: L + pw + 8, y: e.y + 4,
        style: 'fill:var(--ink);font-size:13.5px' }, svg);
      t2.textContent = e.label;
    });

    var gx = X(currentVO);
    sv('line', { x1: gx, x2: gx, y1: T, y2: T + ph,
                 'stroke-dasharray': '3 4',
                 style: 'stroke:var(--muted);stroke-width:1' }, svg);
    var tpos = currentVO / xmax * (rows.length - 1);
    var j0 = Math.min(rows.length - 2, Math.floor(tpos));
    var wgt = tpos - j0;
    SWEEP.forEach(function (s) {
      var val = rows[j0][s.key] * (1 - wgt) + rows[j0 + 1][s.key] * wgt;
      sv('circle', { cx: gx, cy: Y(val), r: 4,
        style: 'fill:var(' + s.varname + ');stroke:var(--panel);'
          + 'stroke-width:2' }, svg);
    });
    sv('line', { x1: L, x2: L + pw, y1: T + ph, y2: T + ph,
                 style: 'stroke:var(--muted);stroke-width:1' }, svg);

    var expPts = (DATA.experimental_recovery
                  && DATA.experimental_recovery.points) || [];
    expPts.forEach(function (q) {
      if (q.VO_total_umol_g <= xmax && q.recovered_umol_g <= ymax) {
        sv('circle', { cx: X(q.VO_total_umol_g), cy: Y(q.recovered_umol_g),
          r: 4.5, style: 'fill:none;stroke:var(--ink);stroke-width:2' }, svg);
      }
    });

    var ch = sv('line', { x1: L, x2: L, y1: T, y2: T + ph, opacity: 0,
                 style: 'stroke:var(--muted);stroke-width:1' }, svg);
    sweepState = { rows: rows, xmax: xmax, L: L, T: T, pw: pw, ph: ph,
                   W: W, ch: ch };
    if (!svg.dataset.wired) {
      svg.dataset.wired = '1';
      svg.addEventListener('mousemove', sweepMove);
      svg.addEventListener('mouseleave', sweepLeave);
    }

    var lg = '';
    SWEEP.forEach(function (s) {
      lg += '<span>' + (s.dash
        ? '<svg width="20" height="6" style="vertical-align:middle;'
          + 'margin-right:5px"><line x1="0" x2="20" y1="3" y2="3" '
          + 'stroke-dasharray="4 3" style="stroke:var(--ink);'
          + 'stroke-width:2"/></svg>'
        : '<span class="swatch" style="background:var(' + s.varname
          + ')"></span>') + s.label + '</span>';
    });
    if (expPts.length) {
      lg += '<span><svg width="12" height="12" style="vertical-align:middle;'
        + 'margin-right:5px"><circle cx="6" cy="6" r="4" fill="none" '
        + 'style="stroke:var(--ink);stroke-width:2"/></svg>'
        + 'measured CO₂ recovery</span>';
    }
    el.smSweepLegend.innerHTML = lg;
  }

  function sweepMove(ev) {
    var st = sweepState;
    if (!st) return;
    var svg = el.smSweep;
    var rect = svg.getBoundingClientRect();
    var sx = (ev.clientX - rect.left) * (st.W / rect.width);
    if (sx < st.L || sx > st.L + st.pw) { sweepLeave(); return; }
    var j = Math.round((sx - st.L) / st.pw * (st.rows.length - 1));
    j = Math.max(0, Math.min(st.rows.length - 1, j));
    var row = st.rows[j];
    var gx = st.L + st.pw * row.VO_total_umol_g / st.xmax;
    st.ch.setAttribute('x1', gx);
    st.ch.setAttribute('x2', gx);
    st.ch.setAttribute('opacity', 0.6);
    var tip = el.smTip;
    tip.style.display = 'block';
    tip.innerHTML = '<b>' + f2(row.VO_total_umol_g)
      + ' μmol-O/g total</b><br>surface ' + f2(row.surface)
      + ' · subsurface ' + f2(row.subsurface) + ' · bulk ' + f2(row.bulk)
      + '<br>accessible ' + f2(row.accessible) + ' ('
      + pctf(row.f_recoverable) + ')';
    var wrap = svg.parentNode;
    var wr = wrap.getBoundingClientRect();
    var px = (gx / st.W) * rect.width + (rect.left - wr.left);
    var flip = px > wr.width * 0.62;
    tip.style.left = (flip ? px - tip.offsetWidth - 14 : px + 12) + 'px';
    tip.style.top = (ev.clientY - wr.top - 14) + 'px';
  }

  function sweepLeave() {
    if (sweepState && sweepState.ch) sweepState.ch.setAttribute('opacity', 0);
    el.smTip.style.display = 'none';
  }

  /* ------------------------------------------ particle cross-section */

  var lastPv = null;

  function tx(x, y, str, style, anchor, parent) {
    var t = sv('text', { x: x, y: y,
      'text-anchor': anchor || 'start', style: style }, parent);
    t.textContent = str;
    return t;
  }

  function vacDot(x, y, cls, prob, accMode, parent) {
    var varn = cls === 'surface' ? '--cSurf'
      : cls === 'subsurface' ? '--cSub' : '--cBulk';
    var attrs = { cx: x, cy: y, r: 2.6, 'data-cls': cls, 'data-kind': 'vac' };
    if (!accMode) {
      attrs.style = 'fill:var(' + varn + ')';
    } else if (prob >= 0.05) {
      attrs.style = 'fill:var(' + varn + ');fill-opacity:'
        + Math.max(0.3, prob).toFixed(2);
    } else {
      attrs.style = 'fill:none;stroke:var(' + varn + ');stroke-width:1.3';
      attrs['data-locked'] = '1';
    }
    return sv('circle', attrs, parent);
  }

  function surfacePattern(nCols, thS, fresh) {
    if (fresh && iso) {
      var best = null;
      iso.forEach(function (q) {
        if (q.mu_eV == null || !q.surf_vac) return;
        var d = Math.abs(q.theta_s - thS);
        if (!best || d < best[0]) best = [d, q];
      });
      if (best) {
        var q = best[1];
        var occ = new Uint8Array(q.rows * q.row_sites);
        q.surf_vac.forEach(function (v) { occ[v] = 1; });
        var arr = [];
        for (var r = 0; r < q.rows && arr.length < nCols; r++) {
          for (var k = 0; k < q.row_sites && arr.length < nCols; k++) {
            arr.push(occ[r * q.row_sites + k]);
          }
        }
        return { arr: arr, sampled: true };
      }
    }
    var prng = SM.makeRng(4242);
    var arr2 = [];
    for (var j = 0; j < nCols; j++) arr2.push(prng() < thS ? 1 : 0);
    return { arr: arr2, sampled: false };
  }

  function renderParticle(out, acc, fresh) {
    lastPv = { out: out, acc: acc, fresh: fresh };
    var svg = el.smParticle;
    if (!svg) return;
    while (svg.firstChild) svg.removeChild(svg.firstChild);
    var accMode = el.smPvAcc.checked;
    var geom = out.geometry;
    var aA = DATA.lattice_constants_A.a;
    var cA = DATA.lattice_constants_A.c;
    var t1 = aA / Math.sqrt(2) / 10;                 // d110, nm
    var tShell = 2 * t1;
    var rho = SM.densityGCm3(DATA);
    var dUm = 6 / (rho * geom.area_m2_g * 1e4) * 1e4;
    var Rnm = dUm * 1000 / 2;
    var prng = SM.makeRng(777);

    /* ------------------------------------------- true-scale disc */
    var cx = 182, cy = 208, Rpx = 158;
    sv('circle', { cx: cx, cy: cy, r: Rpx,
      style: 'fill:var(--code);stroke:none' }, svg);
    var coreOpacity = Math.min(0.22, 0.035
      + 0.18 * Math.sqrt(Math.max(0, out.theta.bulk)));
    sv('circle', { cx: cx, cy: cy, r: Rpx - 3,
      style: 'fill:var(--cBulk);opacity:' + coreOpacity.toFixed(3)
        + ';stroke:none' }, svg);
    sv('circle', { cx: cx, cy: cy, r: Rpx,
      style: 'fill:none;stroke:var(--cSurf);stroke-width:1.3' }, svg);
    sv('circle', { cx: cx, cy: cy, r: Rpx - 1.6,
      style: 'fill:none;stroke:var(--cSub);stroke-width:1.1;opacity:.9' },
       svg);

    tx(cx, cy - 4, 'bulk occupation',
       'fill:var(--navy);font-size:13px;font-weight:bold', 'middle', svg);
    tx(cx, cy + 13, 'θb = ' + thf(out.theta.bulk),
       'fill:var(--muted);font-size:12px', 'middle', svg);

    /* zoom marker + leaders to the window frame */
    var x0 = 386, y0 = 24, fw = 331, fh = 250;
    var boxW = 16, boxH = 9;
    sv('rect', { x: cx - boxW / 2, y: cy - Rpx - 2, width: boxW,
      height: boxH, style: 'fill:none;stroke:var(--ink);stroke-width:1.1' },
       svg);
    ['M' + (cx + boxW / 2) + ',' + (cy - Rpx - 2) + ' L' + x0 + ',' + y0,
     'M' + (cx + boxW / 2) + ',' + (cy - Rpx + boxH - 2) + ' L' + x0 + ','
       + (y0 + fh)].forEach(function (d) {
      sv('path', { d: d, 'stroke-dasharray': '4 4',
        style: 'stroke:var(--muted);stroke-width:1;fill:none;opacity:.7' },
         svg);
    });

    var sb = Math.max(20, Math.round(Rnm / 4 / 50) * 50) || 50;
    var sbPx = sb * Rpx / Rnm;
    sv('line', { x1: cx - Rpx, x2: cx - Rpx + sbPx, y1: cy + Rpx + 22,
      y2: cy + Rpx + 22, style: 'stroke:var(--ink);stroke-width:2' }, svg);
    tx(cx - Rpx + sbPx / 2, cy + Rpx + 36, sb + ' nm',
       'fill:var(--muted);font-size:15px', 'middle', svg);
    tx(cx + 24, cy + Rpx + 36, 'd = ' + dUm.toFixed(2) + ' μm',
       'fill:var(--muted);font-size:15px', 'middle', svg);
    tx(cx, y0 + 6, 'shell 0.65 nm, not resolvable at this scale',
       'fill:var(--muted);font-size:12.5px', 'middle', svg);

    /* -------------------------------------------- magnified window */
    var nCols = 46;
    var sitePx = fw / nCols;
    var gasH = 24;
    var pxNm = 27.5;
    sv('rect', { x: x0, y: y0, width: fw, height: fh,
      style: 'fill:var(--panel);stroke:var(--line)' }, svg);
    sv('rect', { x: x0, y: y0, width: fw, height: gasH,
      style: 'fill:var(--code);opacity:.7' }, svg);
    tx(x0 + fw - 6, y0 + 16, 'gas', 'fill:var(--muted);font-size:12.5px',
       'end', svg);

    var surfY = y0 + gasH + 10;
    var rowY = function (k) { return surfY + k * t1 * pxNm; };
    var nRows = Math.floor((y0 + fh - 8 - surfY) / (t1 * pxNm)) + 1;
    for (var k = 0; k < nRows; k++) {
      sv('line', { x1: x0 + sitePx / 2 - 0.8, x2: x0 + fw - 2,
        y1: rowY(k), y2: rowY(k),
        'stroke-dasharray': '1.6 ' + (sitePx - 1.6).toFixed(2),
        style: 'stroke:var(--muted);stroke-width:1.6;opacity:.25' }, svg);
    }
    var siteX = function (j) { return x0 + (j + 0.5) * sitePx; };

    /* Surface: sampled configuration. Subsurface: expected count at seeded
       sites. Bulk: continuous class-average shading. */
    function pickCells(n, kk, rng2) {
      var idx = [];
      var t;
      for (t = 0; t < n; t++) idx.push(t);
      for (t = 0; t < kk && t < n; t++) {
        var swap = t + Math.floor(rng2() * (n - t));
        var tmp = idx[t];
        idx[t] = idx[swap];
        idx[swap] = tmp;
      }
      return idx.slice(0, Math.min(kk, n));
    }
    var pat = surfacePattern(nCols, out.theta.surface, fresh);
    for (var j = 0; j < nCols; j++) {
      if (pat.arr[j]) {
        vacDot(siteX(j), rowY(0), 'surface', acc.P.surface, accMode, svg);
      }
    }
    var nss = Math.round(nCols * out.theta.subsurface);
    pickCells(nCols, nss, prng).forEach(function (jj) {
      vacDot(siteX(jj), rowY(1), 'subsurface', acc.P.subsurface, accMode,
             svg);
    });
    var bulkTop = rowY(2) - 0.42 * t1 * pxNm;
    var bulkBottom = Math.min(y0 + fh - 8,
                              rowY(Math.max(2, nRows - 1)) + 0.42 * t1 * pxNm);
    sv('rect', { x: x0 + 2, y: bulkTop, width: fw - 4,
      height: Math.max(0, bulkBottom - bulkTop),
      style: 'fill:var(--cBulk);opacity:' + coreOpacity.toFixed(3)
        + ';stroke:none' }, svg);

    sv('line', { x1: x0 + 10, x2: x0 + 10, y1: y0 + fh - 10 - pxNm,
      y2: y0 + fh - 10, style: 'stroke:var(--ink);stroke-width:2' }, svg);
    tx(x0 + 16, y0 + fh - 12, '1 nm', 'fill:var(--muted);font-size:12.5px',
       'start', svg);

    var lx = x0 + fw + 10;
    var lab = function (y, cls, name, u, th, prob) {
      sv('rect', { x: lx, y: y - 8, width: 9, height: 9, rx: 2,
        style: 'fill:var(--c' + cls + ')' }, svg);
      tx(lx + 14, y, name, 'fill:var(--navy);font-size:13.5px;'
         + 'font-weight:bold', 'start', svg);
      tx(lx + 14, y + 13, f2(u) + ' μmol/g · θ ' + thf(th),
         'fill:var(--muted);font-size:12.5px', 'start', svg);
      tx(lx + 14, y + 25, 'P(CO₂) ' + probf(prob),
         'fill:var(--muted);font-size:12.5px', 'start', svg);
    };
    lab(rowY(0), 'Surf', 'Surface', out.umol_g.surface,
        out.theta.surface, acc.P.surface);
    lab(rowY(1) + 34, 'Sub', 'Subsurface', out.umol_g.subsurface,
        out.theta.subsurface, acc.P.subsurface);
    lab(rowY(1) + 108, 'Bulk', 'Bulk', out.umol_g.bulk,
        out.theta.bulk, acc.P.bulk);

    /* ------------------------------------------------ depth profile */
    var px0 = 386, py0 = 316, pw2 = 474, ph2 = 128;
    var papTop = py0 + 4, papBot = py0 + ph2 - 24;
    var linW = 0.5 * pw2, logW = 0.42 * pw2;
    var lgMax = Math.log10(Rnm);
    var XD = function (dep) {
      if (dep <= 1) return px0 + dep * linW;
      return px0 + linW + 14 + Math.log10(dep) / lgMax * logW;
    };
    var YT = function (th) {
      var lg = Math.log10(Math.max(th, 1e-4));
      return papBot - (lg + 4) / 4 * (papBot - papTop);
    };
    [[1, '100%'], [0.01, '1%'], [0.0001, '0.01%']].forEach(function (g) {
      var y = YT(g[0]);
      sv('line', { x1: px0, x2: px0 + pw2, y1: y, y2: y,
        style: 'stroke:var(--line);stroke-width:1' }, svg);
      tx(px0 - 6, y + 3.5, g[1], 'fill:var(--muted);font-size:14px', 'end',
         svg);
    });
    // axis break marks
    var bx = px0 + linW + 7;
    ['M' + (bx - 3) + ',' + (papBot + 5) + ' L' + (bx + 1) + ','
     + (papBot - 3), 'M' + (bx + 2) + ',' + (papBot + 5) + ' L' + (bx + 6)
     + ',' + (papBot - 3)].forEach(function (d) {
      sv('path', { d: d, style: 'stroke:var(--muted);stroke-width:1.2;'
        + 'fill:none' }, svg);
    });
    var segs = [[0, t1, out.theta.surface, '--cSurf'],
                [t1, tShell, out.theta.subsurface, '--cSub'],
                [tShell, Rnm, out.theta.bulk, '--cBulk']];
    segs.forEach(function (s, si) {
      var y = YT(s[2]);
      if (si === 2) {
        // the bulk segment spans the axis break: draw around the gap
        sv('line', { x1: XD(s[0]), x2: px0 + linW, y1: y, y2: y,
          style: 'stroke:var(' + s[3] + ');stroke-width:3' }, svg);
        sv('line', { x1: px0 + linW + 14, x2: XD(s[1]), y1: y, y2: y,
          style: 'stroke:var(' + s[3] + ');stroke-width:3' }, svg);
      } else {
        sv('line', { x1: XD(s[0]), x2: XD(s[1]), y1: y, y2: y,
          style: 'stroke:var(' + s[3] + ');stroke-width:3' }, svg);
      }
      if (si) {
        sv('line', { x1: XD(s[0]), x2: XD(s[0]), y1: YT(segs[si - 1][2]),
          y2: y, style: 'stroke:var(--muted);stroke-width:1;opacity:.6' },
           svg);
      }
    });
    sv('line', { x1: px0, x2: px0 + pw2, y1: papBot, y2: papBot,
      style: 'stroke:var(--muted);stroke-width:1' }, svg);
    [[0, '0'], [t1, ''], [tShell, '0.65'], [1, '1']].forEach(function (g) {
      tx(XD(g[0]), papBot + 13, g[1], 'fill:var(--muted);font-size:14px',
         'middle', svg);
    });
    [[10, '10'], [100, '100'], [Rnm, Math.round(Rnm) + '']].forEach(
      function (g) {
        if (g[0] <= Rnm) {
          tx(XD(g[0]), papBot + 13, g[1],
             'fill:var(--muted);font-size:14px', 'middle', svg);
        }
      });
    tx(px0 + pw2 / 2, papBot + 26, 'depth from surface (nm)',
       'fill:var(--muted);font-size:12.5px', 'middle', svg);
    var shellU = out.umol_g.surface + out.umol_g.subsurface;
    var shellPct = out.matched_umol_g
      ? shellU / out.matched_umol_g * 100 : 0;
    tx(XD(tShell) + 8, papTop + 10, shellPct.toFixed(1)
       + '% of the inventory within ' + tShell.toFixed(2) + ' nm',
       'fill:var(--ink);font-size:15px', 'start', svg);

    /* -------------------------------------------- legend + caption */
    var lg2 = '';
    [['Surface', '--cSurf'], ['Subsurface', '--cSub'],
     ['Bulk', '--cBulk']].forEach(function (s) {
      lg2 += '<span><span class="swatch" style="background:var(' + s[1]
        + ')"></span>' + s[0] + '</span>';
    });
    lg2 += '<span>' + (accMode
      ? 'fill strength = P(CO₂ refill); hollow = locked'
      : 'dot = one vacancy (window)') + '</span>';
    lg2 += '<span>disc core = bulk occupation</span>';
    el.smPvLegend.innerHTML = lg2;

    var shellTheta = (geom.N_s + geom.N_ss)
      ? shellU / (geom.N_s + geom.N_ss) : 0;
    el.smPvNote.textContent = 'Shell (' + tShell.toFixed(2) + ' nm): '
      + f2(shellU) + ' μmol-O/g = ' + shellPct.toFixed(1)
      + '% of the inventory at ' + (shellTheta * 100).toFixed(0)
      + '% local depletion; core (' + Math.round(Rnm) + ' nm): '
      + f2(out.umol_g.bulk) + ' μmol-O/g at θ = ' + thf(out.theta.bulk)
      + '. Surface row: ' + (pat.sampled
        ? 'sampled Monte Carlo configuration.'
        : 'ideal-preview draw. Run the model for the sampled arrangement.');
  }

  function exportParticleSvg() {
    var src = new XMLSerializer().serializeToString(el.smParticle);
    var cs = getComputedStyle(document.documentElement);
    ['--cSurf', '--cSub', '--cBulk', '--ink', '--muted', '--line',
     '--navy', '--panel', '--code', '--bg'].forEach(function (v) {
      var val = cs.getPropertyValue(v).trim();
      if (val) src = src.split('var(' + v + ')').join(val);
    });
    var bg = cs.getPropertyValue('--bg').trim() || '#ffffff';
    src = src.replace(/(<svg[^>]*>)/,
                      '$1<rect width="880" height="470" fill="' + bg + '"/>');
    var blob = new Blob([src], { type: 'image/svg+xml' });
    var a = document.createElement('a');
    a.href = URL.createObjectURL(blob);
    a.download = 'particle_cross_section.svg';
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
  }

  /* --------------------------------------------------------- wiring */

  function fillPreset(key) {
    var q = DATA.vacancy_energetics.presets[key];
    if (!q) return;
    el.smEpsSS.value = q.subsurface_eV;
    el.smEpsB.value = q.bulk_eV;
  }

  el.smPreset.addEventListener('change', function () {
    if (el.smPreset.value !== 'custom') fillPreset(el.smPreset.value);
    solve();
    scheduleRun();
  });
  [el.smEpsSS, el.smEpsB, el.smV1, el.smV2, el.smV3].forEach(function (n) {
    n.addEventListener('input', function () {
      var key = el.smPreset.value;
      var q = DATA.vacancy_energetics.presets[key];
      if (q && (num(el.smEpsSS, -1) !== q.subsurface_eV
                || num(el.smEpsB, -1) !== q.bulk_eV)) {
        el.smPreset.value = 'custom';
      }
      solve();
      scheduleRun();
    });
  });
  [el.smT, el.smQuality, el.smSeed].forEach(function (n) {
    n.addEventListener(n === el.smQuality ? 'change' : 'input', function () {
      solve();
      scheduleRun();
    });
  });
  el.smT.addEventListener('input', scheduleEnv);
  [el.smVO, el.smGeoVal, el.smEaS, el.smEaSS, el.smEaB, el.smAS, el.smASS,
   el.smAB, el.smKinT, el.smExp, el.smPCO2].forEach(function (n) {
    n.addEventListener('input', solve);
  });
  el.smPvAcc.addEventListener('change', function () {
    if (lastPv) renderParticle(lastPv.out, lastPv.acc, lastPv.fresh);
  });
  el.smPvDl.addEventListener('click', exportParticleSvg);
  el.smGeoMode.addEventListener('change', function () {
    if (el.smGeoMode.value === 'bet') {
      var g = SM.geometry(DATA, { d_um: num(el.smGeoVal, 0.9) });
      el.smGeoVal.value = g.area_m2_g.toFixed(2);
      el.smGeoUnit.textContent = 'm²/g (BET)';
    } else {
      el.smGeoVal.value = DATA.defaults.particle_diameter_um;
      el.smGeoUnit.textContent = 'μm (spherical particle)';
    }
    solve();
  });

  /* ---------------------------------------------------------- init */

  el.smT.value = DATA.defaults.T_C;
  el.smVO.value = DATA.defaults.VO_total_umol_g;
  el.smGeoVal.value = DATA.defaults.particle_diameter_um;
  el.smSeed.value = DATA.defaults.mc.seed;
  fillPreset(DATA.vacancy_energetics.default_preset);
  var ar = DATA.surface_ordering.pair_eV.along_row;
  el.smV1.value = ar[0];
  el.smV2.value = ar[1];
  el.smV3.value = ar[2];
  var kj = DATA.co2_kinetics;
  el.smEaS.value = kj.Ea_eV.surface;
  el.smEaSS.value = kj.Ea_eV.subsurface;
  el.smEaB.value = kj.Ea_eV.bulk;
  el.smAS.value = kj.A_eff_per_s_atm.surface;
  el.smASS.value = kj.A_eff_per_s_atm.subsurface;
  el.smAB.value = kj.A_eff_per_s_atm.bulk;
  el.smKinT.value = kj.defaults.T_reox_C;
  el.smExp.value = kj.defaults.exposure_s;
  el.smPCO2.value = kj.defaults.p_CO2_atm;
})();
