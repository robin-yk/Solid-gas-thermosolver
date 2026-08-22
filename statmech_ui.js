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
   'smTable', 'smSpacing', 'smNotes'].forEach(function (id) { el[id] = $(id); });

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
      if (btn.dataset.ws === 'ws-defect' && !started) {
        started = true;
        solve();
        launch();
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

  /* ------------------------------------------------ solve + render */

  function solve() {
    var geom = currentGeom();
    var tC = currentT();
    var vo = num(el.smVO, DATA.defaults.VO_total_umol_g);
    var eps = currentEps();
    var kt = SM.KB_EV * (tC + 273.15);
    var fresh = !!iso && isoKey === isoSignature();
    var opts = { eps: eps };
    if (fresh) opts.theta_s_fn = SM.thetaSInterp(iso, eps[0], kt);
    var out = SM.distribute(DATA, tC, vo, geom, opts);
    out.geometry = geom;
    out.spacing = fresh ? SM.spacingSummary(iso, out.theta.surface) : null;
    render(out, fresh);
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

  function render(out, fresh) {
    var geom = out.geometry;
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
         + 'preview — the interacting isotherm is being sampled.</span>');

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
        + '= ' + (out.spacing.P_gap_le_2 * 100).toFixed(1)
        + '%. Literature window: modal 4–5, close pairs suppressed.';
    } else {
      el.smSpacing.textContent = '';
    }

    var notes = '';
    out.warnings.forEach(function (w) {
      var calm = w.kind === 'surface_cap' || w.kind === 'x_near';
      notes += '<div class="note' + (calm ? ' calm' : '') + '">'
        + (w.kind === 'x_max' ? '<b>Model range.</b> ' : '')
        + esc(w.text) + '</div>';
    });
    el.smNotes.innerHTML = notes;
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
        style: 'fill:var(--muted);font-size:11px' }, svg);
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
        style: 'fill:var(--ink);font-size:13px;font-weight:700' }, svg);
      vt.textContent = pctf(f);
      var nt = sv('text', { x: x + barW / 2, y: T + ph + 17,
        'text-anchor': 'middle',
        style: 'fill:var(--navy);font-size:12px;font-weight:700' }, svg);
      nt.textContent = c.label;
      var ut = sv('text', { x: x + barW / 2, y: T + ph + 33,
        'text-anchor': 'middle',
        style: 'fill:var(--muted);font-size:10.5px' }, svg);
      ut.textContent = f2(out.umol_g[c.key]) + ' μmol-O/g';
    });
    sv('line', { x1: L, x2: L + pw, y1: T + ph, y2: T + ph,
                 style: 'stroke:var(--muted);stroke-width:1' }, svg);
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
  [el.smVO, el.smGeoVal].forEach(function (n) {
    n.addEventListener('input', solve);
  });
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
})();
