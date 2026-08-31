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
  ['smT', 'smVO', 'smGeoMode', 'smGeoVal', 'smGeoUnit', 'smSS', 'smDerived',
   'smPreset', 'smEpsSS', 'smEpsB',
   'smV1', 'smV2', 'smV3', 'smQuality',
   'smSeed', 'smRun', 'smStatus', 'smSummary', 'smKpis',
   'smTable', 'smSpacing', 'smNotes',
   'smEaS', 'smEaSS', 'smEaB', 'smAS', 'smASS', 'smAB', 'smKinT', 'smExp',
   'smPCO2', 'smAccKpis', 'smAccNote',
   'smCoupledNote', 'smEnv', 'smPvAcc', 'smPvNote',
   'cgEm', 'cgCells', 'cgTarget', 'cgFit',
   'cgFitOut', 'cgKpis', 'cgNote'].forEach(function (id) {
     el[id] = $(id);
   });

  /* ------------------------------------------------- published figures */

  /* Every quantitative plot on this page is a publication figure built by
     web/figures_defect.js through the shared kit. The page holds the SVG
     string it was given and hands the same string to the download, so the
     preview and the file cannot differ. */
  var FIG = window.DefectFigures;
  var KIT = window.FigKit;
  var figState = {};

  function mountFigure(hostId, name, svg) {
    var host = $(hostId);
    if (!host) return;
    figState[hostId] = { svg: svg, name: name };
    var box = host.querySelector('.figbox');
    box.innerHTML = svg;
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
    if (host) {
      host.querySelector('.figbox').textContent = why;
    }
  }

  /* ------------------------------------------------- result tabs */

  /* Three views of one solve. The panels are hidden rather than unbuilt,
     so switching tabs never re-runs anything and the figures keep the SVG
     they were mounted with. */
  var dfxTabs = document.querySelectorAll('#dfxTabs button');
  Array.prototype.forEach.call(dfxTabs, function (btn) {
    btn.addEventListener('click', function () {
      Array.prototype.forEach.call(dfxTabs, function (b) {
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
      Array.prototype.forEach.call(tabs, function (b) {
        b.classList.toggle('active', b === btn);
      });
      ['ws-thermo', 'ws-defect', 'ws-redox'].forEach(function (ws) {
        var host = $(ws);
        if (host) host.style.display = ws === btn.dataset.ws ? '' : 'none';
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
  function currentSS() {
    var v = Math.round(num(el.smSS, DATA.defaults.subsurface_layers));
    return Math.max(1, Math.min(4, v));
  }
  function currentGeom() {
    var n = currentSS();
    if (el.smGeoMode.value === 'bet') {
      return SM.geometry(DATA, { bet_m2_g: num(el.smGeoVal, 1.57),
                                 ss_layers: n });
    }
    return SM.geometry(DATA, {
      d_um: num(el.smGeoVal, DATA.defaults.particle_diameter_um),
      ss_layers: n });
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
    /* statmech.js borrows the free energy from particle.js, so the blob
       has to carry it first - the worker has no page globals. */
    var src = document.getElementById('px-engine').textContent + '\n'
      + document.getElementById('sm-engine').textContent + '\n'
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
      el.smEnv.className = 'dfx-status';
      el.smEnv.innerHTML = '<span class="more">The equilibrium workspace '
        + 'is unavailable on this page.</span>';
      return;
    }
    var desc = TB.describe();
    if (!desc) {
      env = null;
      el.smEnv.className = 'dfx-status';
      el.smEnv.innerHTML = '<span class="more">The equilibrium '
        + "workspace's gas charge sums to zero.</span>";
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

  /* One line, because the equilibrium verdict bounds this workspace
     rather than being part of it: the reader needs to know whether rutile
     is the phase to be in, not the whole charge. The detail is one click
     away in the workspace that owns it. */
  function renderEnv() {
    var box = el.smEnv;
    if (!env) return;
    box.className = 'dfx-status';
    var h;
    if (env.outOfRange) {
      h = env.T_C + ' \u00B0C is outside the equilibrium data range '
        + '(300\u20131650 \u00B0C); no phase verdict at this temperature.';
      box.className += ' warn';
    } else if (env.error) {
      h = 'Equilibrium solve failed: ' + esc(env.error);
      box.className += ' bad';
    } else {
      var v = env.v;
      if (v.tier === 'stable') {
        h = '<b>Rutile stable</b> \u00B7 margin to ' + subHtml(v.nearest_phase)
          + ': +' + v.margin_per_mol_O_kJ.toFixed(1) + ' kJ\u00A0mol\u207B\u00B9-O';
      } else if (v.tier === 'boundary') {
        h = '<b>On the reduction boundary</b> \u00B7 TiO<sub>2</sub> coexists '
          + 'with ' + v.reduced_active.map(subHtml).join(' + ');
        box.className += ' warn';
      } else {
        h = '<b>Rutile absent at equilibrium</b> \u00B7 thermodynamics '
          + 'predicts ' + v.reduced_active.map(subHtml).join(' + ')
          + '; the single-phase model does not apply here';
        box.className += ' bad';
      }
      h += ' <span class="more">\u03BC<sub>O</sub> = ' + env.mu_O.toFixed(2)
        + ' kJ/mol-O at ' + env.T_C + ' \u00B0C \u00B7 ' + env.desc + '</span>';
    }
    h += '<span class="more"><button type="button" id="smEnvOpen">Edit '
      + 'charge</button>';
    if (env.thermoT != null && Math.abs(env.thermoT - env.T_C) > 1e-9) {
      h += ' \u00B7 <button type="button" id="smEnvMatch">Use equilibrium T ('
        + env.thermoT + ' \u00B0C)</button>';
    }
    h += '</span>';
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
    /* decoupled by design: the layer partition is the analytic
       phase-aware solve on the formation energies; the sampled isotherm
       feeds only the ordering exhibit, at the coverage of the surviving
       (1x1) patches once coexistence is reached */
    var fresh = !!iso && isoKey === isoSignature();
    var out = SM.distribute(DATA, tC, vo, geom, { eps: eps });
    out.geometry = geom;
    out.spacing = fresh ? SM.spacingSummary(iso,
      Math.min(out.theta.surface, DATA.surface_phases.theta_transition))
      : null;
    render(out, fresh, { tC: tC, eps: eps, geom: geom,
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

  var REGIME_LABEL = {
    dilute: '(1×1) point-defect regime, below the reconstruction onset.',
    surface_coexistence: 'Surface in (1×1)/(1×2) two-phase coexistence '
      + '(lever rule at pinned μ<sub>V</sub>).',
    reconstructed: 'Surface fully reconstructed to the (1×2) added-row '
      + 'Ti₂O₃ phase; deeper classes on their isotherms.',
    cs_coexistence: 'μ<sub>V</sub> pinned at the rutile/CS coexistence '
      + 'estimate; the excess precipitates as extended defects.'
  };

  function render(out, fresh, ctx) {
    var geom = out.geometry;
    var acc = SM.accessibility(DATA, out.umol_g, ctx.kin,
                               out.VO_total_umol_g);
    var b = out.phase_boundaries;
    var slabSum = geom.ss_layers;
    el.smDerived.textContent = 'x in TiO2-x = ' + out.x_TiO2mx.toFixed(4)
      + ' · Ti³⁺/Ti = ' + (out.Ti3_frac * 100).toFixed(2)
      + '% · area ' + geom.area_m2_g.toFixed(2)
      + ' m²/g · shell = bridging plane + ' + slabSum
      + ' × 0.325 nm oxygen layer' + (slabSum > 1 ? 's' : '')
      + ' (N_s ' + geom.N_s.toFixed(2) + ' + N_ss '
      + geom.N_ss.toFixed(2) + ' μmol-O/g) · phase ladder at ' + out.T_C
      + ' °C: (1×2) onset '
      + b.VO_onset_umol_g.toFixed(2) + ' → reconstruction complete '
      + b.VO_recon_complete_umol_g.toFixed(2) + ' → CS ceiling ~'
      + b.VO_cs_umol_g.toFixed(0) + ' μmol-O/g (est.)';

    var bad = out.warnings.some(function (w) {
      return w.kind === 'unresolved';
    });
    var parts = CLASSES.map(function (c) {
      return c.label.toLowerCase() + ' ' + f2(out.umol_g[c.key]) + ' ('
        + pctf(out.fractions[c.key]) + ')';
    });
    if (out.extended_defects_umol_g > 0) {
      parts.push('extended defects ' + f2(out.extended_defects_umol_g)
        + ' (' + pctf(out.fractions.extended) + ')');
    }
    el.smSummary.className = 'verdict ' + (bad ? 'red' : 'none');
    var sp0 = out.surface_phase;
    var phaseTxt = sp0.phase === '1x1' ? '(1×1) vacancy lattice gas'
      : sp0.phase === '1x2' ? '(1×2) added-row Ti₂O₃, θ_eff = 0.5'
      : '(1×1)+(1×2) coexistence, '
        + (sp0.phi_reconstructed * 100).toFixed(0) + '% of area reconstructed';
    el.smSummary.innerHTML = '<b>' + f2(out.VO_total_umol_g)
      + ' μmol-O/g at ' + out.T_C
      + ' °C (conditional equilibrium):</b> ' + parts.join(', ')
      + ' μmol-O/g. ' + REGIME_LABEL[out.regime]
      + '<span class="more">Surface phase: ' + phaseTxt
      + ' · μ<sub>V</sub> = ' + out.mu_V_eV.toFixed(3) + ' eV'
      + (out.extended_defects_umol_g > 0
         ? ' · the excess is CS planes by the lever rule, an estimate' : '')
      + '</span>';

    /* Four, and no more: the three site classes and the potential that
       makes them one solve. The surface phase, the extended-defect count
       and the recoverable fraction all have a home of their own - the
       summary line, the table, and the accessibility tab - and repeating
       them here only makes the row long enough to be scrolled past. */
    var k = '';
    CLASSES.forEach(function (c) {
      k += '<div class="kpi"><div class="k"><span class="swatch" style="'
        + 'background:var(' + c.varname + ')"></span>' + c.label
        + ' V<sub>O</sub></div><div class="v">' + f2(out.umol_g[c.key])
        + '</div><div class="u">μmol-O/g · '
        + pctf(out.fractions[c.key]) + ' of inventory</div></div>';
    });
    k += '<div class="kpi"><div class="k">μ<sub>V</sub> (vacancy)</div>'
      + '<div class="v">' + out.mu_V_eV.toFixed(3) + '</div>'
      + '<div class="u">eV'
      + (out.regime === 'surface_coexistence'
         || out.regime === 'cs_coexistence'
         ? ', pinned at coexistence' : ', common to all classes')
      + '</div></div>';
    el.smKpis.innerHTML = k;

    try {
      mountFigure('figDistribution', 'vacancy-distribution',
        FIG.distribution({ out: out }));
    } catch (e) { figFail('figDistribution', e.message); }

    var h = '<tr><th>Site class</th><th>θ occupancy</th>'
      + '<th>μmol-O/g</th><th>share of inventory</th></tr>';
    CLASSES.forEach(function (c) {
      var th = thf(out.theta[c.key]);
      if (c.key === 'surface' && sp0.phase !== '1x1') {
        th += sp0.phase === '1x2' ? ' · (1×2)' : ' · mixed';
      }
      h += '<tr><td><span class="swatch" style="background:var(' + c.varname
        + ')"></span>' + c.label + '</td><td>' + th
        + '</td><td>' + f2(out.umol_g[c.key]) + '</td><td>'
        + pctf(out.fractions[c.key]) + '</td></tr>';
    });
    if (out.extended_defects_umol_g > 0) {
      parts.push('extended defects ' + f2(out.extended_defects_umol_g)
        + ' (' + pctf(out.fractions.extended) + ')');
    }
    el.smSummary.className = 'verdict ' + (bad ? 'red' : 'none');
    var sp0 = out.surface_phase;
    var phaseTxt = sp0.phase === '1x1' ? '(1×1) vacancy lattice gas'
      : sp0.phase === '1x2' ? '(1×2) added-row Ti₂O₃, θ_eff = 0.5'
      : '(1×1)+(1×2) coexistence, '
        + (sp0.phi_reconstructed * 100).toFixed(0) + '% of area reconstructed';
    el.smSummary.innerHTML = '<b>' + f2(out.VO_total_umol_g)
      + ' μmol-O/g at ' + out.T_C
      + ' °C (conditional equilibrium):</b> ' + parts.join(', ')
      + ' μmol-O/g. ' + REGIME_LABEL[out.regime]
      + '<span class="more">Surface phase: ' + phaseTxt
      + ' · μ<sub>V</sub> = ' + out.mu_V_eV.toFixed(3) + ' eV'
      + (out.extended_defects_umol_g > 0
         ? ' · the excess is CS planes by the lever rule, an estimate' : '')
      + '</span>';

    /* Four, and no more: the three site classes and the potential that
       makes them one solve. The surface phase, the extended-defect count
       and the recoverable fraction all have a home of their own - the
       summary line, the table, and the accessibility tab - and repeating
       them here only makes the row long enough to be scrolled past. */
    var k = '';
    CLASSES.forEach(function (c) {
      k += '<div class="kpi"><div class="k"><span class="swatch" style="'
        + 'background:var(' + c.varname + ')"></span>' + c.label
        + ' V<sub>O</sub></div><div class="v">' + f2(out.umol_g[c.key])
        + '</div><div class="u">μmol-O/g · '
        + pctf(out.fractions[c.key]) + ' of inventory</div></div>';
    });
    k += '<div class="kpi"><div class="k">μ<sub>V</sub> (vacancy)</div>'
      + '<div class="v">' + out.mu_V_eV.toFixed(3) + '</div>'
      + '<div class="u">eV'
      + (out.regime === 'surface_coexistence'
         || out.regime === 'cs_coexistence'
         ? ', pinned at coexistence' : ', common to all classes')
      + '</div></div>';
    el.smKpis.innerHTML = k;

    try {
      mountFigure('figDistribution', 'vacancy-distribution',
        FIG.distribution({ out: out }));
    } catch (e) { figFail('figDistribution', e.message); }

    var h = '<tr><th>Site class</th><th>θ occupancy</th>'
      + '<th>μmol-O/g</th><th>share of inventory</th></tr>';
    CLASSES.forEach(function (c) {
      var th = thf(out.theta[c.key]);
      if (c.key === 'surface' && sp0.phase !== '1x1') {
        th += sp0.phase === '1x2' ? ' · (1×2)' : ' · mixed';
      }
      h += '<tr><td><span class="swatch" style="background:var(' + c.varname
        + ')"></span>' + c.label + '</td><td>' + th
        + '</td><td>' + f2(out.umol_g[c.key]) + '</td><td>'
        + pctf(out.fractions[c.key]) + '</td></tr>';
    });
    if (out.extended_defects_umol_g > 0) {
      h += '<tr><td>Extended defects (CS)</td><td>—</td><td>'
        + f2(out.extended_defects_umol_g) + '</td><td>'
        + pctf(out.fractions.extended) + '</td></tr>';
    }
    h += '<tr><td><b>Total</b></td><td></td><td><b>'
      + f2(out.matched_umol_g) + '</b></td><td><b>100%</b></td></tr>';
    el.smTable.innerHTML = h;

    if (out.spacing) {
      el.smSpacing.textContent = 'Surface ordering (separate calculation; '
        + 'effective pair interactions constrained to the Birschitzky '
        + 'et al. statistics): at θs = '
        + out.spacing.at_theta_s.toFixed(3) + ', modal spacing '
        + out.spacing.modal_gap_sites + ' bridging sites, P(spacing ≤ 2) '
        + '= ' + (out.spacing.P_gap_le_2 * 100).toFixed(1) + '%. '
        + (out.surface_reconstruction_regime
          ? 'Sampled at the coverage of the surviving (1×1) patches '
            + '(θ_t) — in coexistence the (1×1) domains sit there.'
          : 'Literature window: modal 4–5, close pairs suppressed.');
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
      var calm = w.kind === 'surface_phase' || w.kind === 'x_near';
      notes += '<div class="note' + (calm ? ' calm' : '') + '">'
        + (w.kind === 'cs_precipitation' ? '<b>Phase change.</b> '
           : w.kind === 'subsurface_dense' ? '<b>Reduced shell.</b> ' : '')
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
      + ' total · illustrative</div></div>';
    k += '<div class="kpi"><div class="k">Recoverable fraction</div>'
      + '<div class="v">' + pctf(acc.f_recoverable) + '</div>'
      + '<div class="u">illustrative until fitted · '
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
    var b = out.phase_boundaries;
    var xmax = Math.max(220, b.VO_cs_umol_g * 1.3, vo * 1.25);
    var voList = [];
    for (var i = 0; i <= 72; i++) voList.push(xmax * i / 72);
    var rows = SM.sweep(DATA, ctx.tC, ctx.geom, voList,
                        { eps: ctx.eps, kin: ctx.kin });
    try {
      mountFigure('figAccessible', 'accessible-inventory', FIG.accessible({
        vo: vo, bounds: b, acc: acc,
        sweep: {
          xmax: xmax,
          rows: rows.map(function (r) {
            return { VO_total_umol_g: r.VO_total_umol_g,
                     surface: r.surface, subsurface: r.subsurface,
                     bulk: r.bulk, extended: r.extended,
                     accessible: r.accessible };
          })
        }
      }));
    } catch (e) { figFail('figAccessible', e.message); }

    renderParticle(out, acc, fresh);
    renderCg(out, acc, ctx);
  }

  function probf(v) {
    if (v >= 0.9995) return '1.000';
    if (v >= 0.001) return v.toFixed(3);
    if (v > 0) return v.toExponential(1);
    return '0';
  }

  /* ------------------------------------------------------- chart */

  var SVG = 'http://www.w3.org/2000/svg';

  /* ------------------------------------------- particle cross-section */

  /* The drawing lives in web/figures_defect.js with the other four
     figures; this only assembles what it needs. The surface arrangement
     is the sampled Monte Carlo configuration when one exists, so the
     bridging row shows a real snapshot rather than a coin toss. */

  var lastPv = null;
  var cgFitVal = null;

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
    var geom = out.geometry;
    var aA = DATA.lattice_constants_A.a;
    var cA = DATA.lattice_constants_A.c;
    var t1 = aA / Math.sqrt(2) / 10;                 // d110 slab, nm
    var nLay = geom.ss_layers;                       // slabs in the shell
    var tShell = geom.shell_nm;
    var rho = SM.densityGCm3(DATA);
    var dUm = 6 / (rho * geom.area_m2_g * 1e4) * 1e4;
    var pat = surfacePattern(32, out.theta.surface, fresh);
    var shellU = out.umol_g.surface + out.umol_g.subsurface;
    var shellPct = out.matched_umol_g
      ? shellU / out.matched_umol_g * 100 : 0;

    var pv = {
      theta: out.theta,
      theta_bulk_txt: thf(out.theta.bulk),
      P: acc.P,
      accMode: el.smPvAcc.checked,
      n_layers: nLay,
      layer_nm: t1,
      cell_nm: cA / 10,
      shell_nm: tShell,
      d_um: dUm,
      R_nm: dUm * 1000 / 2,
      surf_occ: pat.arr,
      shell_pct: shellPct
    };
    try {
      mountFigure('figParticle', 'particle-cross-section',
        FIG.particle({ pv: pv }));
    } catch (e) { figFail('figParticle', e.message); }

    el.smPvNote.textContent = 'Shell (bridging plane + ' + nLay
      + ' × ' + t1.toFixed(3) + ' nm oxygen layer'
      + (nLay > 1 ? 's' : '') + ' = 1 + ' + (4 * nLay)
      + ' O per (1×1) cell, ' + tShell.toFixed(2) + ' nm): '
      + f2(shellU) + ' μmol-O/g = ' + shellPct.toFixed(1)
      + '% of the inventory. Surface arrangement '
      + (pat.sampled ? 'from the sampled Monte Carlo configuration'
         : 'is an uncorrelated preview until a run finishes')
      + '. The bulk class is drawn as a tint, not as placed vacancies: '
      + 'tint = bulk occupation, because the canonical partition carries '
      + 'no depth structure inside that class.';
  }

  function tfmt(t) {
    if (t == null) return '—';
    if (t < 1e-3) return (t * 1e6).toPrecision(2) + ' μs';
    if (t < 1.0) return (t * 1e3).toPrecision(2) + ' ms';
    if (t < 600.0) return t.toPrecision(3) + ' s';
    return Math.round(t) + ' s';
  }

  function effDiameterUm(geom) {
    var rho = SM.densityGCm3(DATA);
    return 6 / (rho * geom.area_m2_g * 1e4) * 1e4;
  }

  function renderCg(out, acc, ctx) {
    if (!window.CGMC || !el.cgKpis) return;
    var tEnd = Math.max(num(el.smExp, 300), 1e-3);
    var sys = CGMC.build(DATA, {
      vo_total: out.VO_total_umol_g, d_um: effDiameterUm(out.geometry),
      kin: ctx.kin, eps: ctx.eps,
      cg: { E_m_bulk_eV: num(el.cgEm, 0.65),
            n_cells: Math.max(6, Math.round(num(el.cgCells, 28))) },
      dist_fractions: out.fractions });
    var r = CGMC.run(sys, tEnd, {});
    var t50 = null;
    for (var i = 0; i < r.t_s.length; i++) {
      if (r.recovery[i] >= 0.5) { t50 = r.t_s[i]; break; }
    }
    var k = '<div class="kpi"><div class="k">CGMC recovery</div><div '
      + 'class="v">' + pctf(r.recovered_fraction) + '</div><div class="u">'
      + 'at ' + tfmt(tEnd) + ', E<sub>m</sub> = '
      + num(el.cgEm, 0.65).toFixed(2) + ' eV</div></div>';
    k += '<div class="kpi"><div class="k">t(50% recovered)</div><div '
      + 'class="v">' + tfmt(t50) + '</div><div class="u">'
      + (t50 == null ? 'not reached in this exposure'
         : 'transport + surface exchange') + '</div></div>';
    k += '<div class="kpi"><div class="k">First-order layer</div><div '
      + 'class="v">' + pctf(acc.f_recoverable) + '</div><div class="u">'
      + 'per-class refill, no transport</div></div>';
    k += '<div class="kpi"><div class="k">Fitted E<sub>m,eff</sub></div>'
      + '<div class="v">' + (cgFitVal == null ? '—'
        : cgFitVal.toFixed(2)) + '</div><div class="u">'
      + (cgFitVal == null ? 'enter a measured recovery and fit'
        : 'eV, from the measured recovery') + '</div></div>';
    el.cgKpis.innerHTML = k;
    var cg = { r: r, sys: sys, tEnd: tEnd,
               E_m_eV: num(el.cgEm, 0.65),
               target_s: tEnd, target_pct: r.recovered_fraction * 100 };
    try {
      mountFigure('figRecovery', 'cgmc-recovery', FIG.recovery({ cg: cg }));
    } catch (e) { figFail('figRecovery', e.message); }
    try {
      mountFigure('figRadial', 'radial-profile', FIG.radial({ cg: cg }));
    } catch (e) { figFail('figRadial', e.message); }
    el.cgNote.textContent = 'Closed-system check: with the gas channel '
      + 'shut, this transport network relaxes to the exact product-'
      + 'binomial equilibrium of its own (1×1) site-class Hamiltonian '
      + '(enforced by tests/test_cgmc.py). Reconstruction kinetics are '
      + 'not modelled; the initial state is the phase-aware partition. '
      + 'Deterministic mean of the CGMC process; fluctuations scale as '
      + '1/√N ≈ 0.01% at the ~10⁸ vacancies of one particle. '
      + 'E_m = 0.65 eV is the neutral-vacancy bulk barrier (Iddir et al.); '
      + 'an effective fitted barrier above it points at charged migration, '
      + 'trapping or defect aggregation rather than neutral hops.';
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
  [el.smVO, el.smGeoVal, el.smSS, el.smEaS, el.smEaSS, el.smEaB, el.smAS,
   el.smASS, el.smAB, el.smKinT, el.smExp, el.smPCO2].forEach(function (n) {
    n.addEventListener('input', solve);
  });
  el.smPvAcc.addEventListener('change', function () {
    if (lastPv) renderParticle(lastPv.out, lastPv.acc, lastPv.fresh);
  });
  [el.cgEm, el.cgCells, el.cgTarget].forEach(function (n) {
    if (n) n.addEventListener('input', solve);
  });
  if (el.cgFit) {
    el.cgFit.addEventListener('click', function () {
      var tgt = num(el.cgTarget, NaN);
      if (!(tgt > 0 && tgt < 100)) {
        el.cgFitOut.textContent = 'enter a measured % first';
        return;
      }
      el.cgFitOut.textContent = 'fitting…';
      setTimeout(function () {
        var geom = currentGeom();
        var em = CGMC.effectiveBarrier(DATA, tgt / 100,
          Math.max(num(el.smExp, 300), 1e-3),
          { vo_total: num(el.smVO, DATA.defaults.VO_total_umol_g),
            d_um: effDiameterUm(geom), kin: currentKin(),
            eps: currentEps(),
            cg: { n_cells: Math.max(6, Math.round(num(el.cgCells, 28))) } });
        cgFitVal = em;
        el.cgEm.value = em.toFixed(3);
        el.cgFitOut.textContent = 'effective Eₘ = ' + em.toFixed(2)
          + ' eV';
        solve();
      }, 30);
    });
  }
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
  if (el.cgEm) {
    el.cgEm.value = DATA.kinetics_cgmc.E_m_bulk_eV;
    el.cgCells.value = DATA.kinetics_cgmc.n_cells;
  }
})();
