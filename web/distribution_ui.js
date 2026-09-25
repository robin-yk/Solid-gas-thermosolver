/* The vacancy-distribution workspace. Every number shown comes from
   paper_outputs/tof_range.json, which the one-model run wrote; the page
   only picks a parameter point and divides the measured rate by the
   site count stored for it. The two 3D figures draw those numbers. */
(function () {
  'use strict';
  var $ = function (id) { return document.getElementById(id); };
  var el = $('tof-data');
  if (!el || !$('ws-distribution')) return;
  var D = JSON.parse(el.textContent);
  var K = window.FigKit, F = window.DistributionFigures, V = window.Vacancy3D;
  var AX = D.axes, SAMPLES = D.samples;
  var NPTS = AX.energy_map.length * AX.cutoff_nm.length * AX.dG_eV.length * AX.f110.length * AX.eps.length;
  var START = { sample: 'R600', map: 'PAB', cutoff: 0.34, dG: 0, f110: 0.75, eps: 'a_axis', reactive: 'BRI' };
  var REACT = {
    BRI: 'All BRI vacancies',
    ISO_z2: 'BRI vacancies, z = 2 intact neighbours',
    ISO_z4: 'BRI vacancies, z = 4 intact neighbours',
    ISO_z8: 'BRI vacancies, z = 8 intact neighbours',
    'BRI+BASAL': 'BRI + layer-1 IPL vacancies'
  };
  var MAPS = { PAB: 'PAB', HAM: 'HAM', LI_SBR1: 'LI_SBR1', LI_L2: 'LI_L2' };

  function sig(v, d) {
    if (v == null || !isFinite(v)) return '—';
    if (v === 0) return '0';
    var a = Math.abs(v);
    if (a >= 1e4 || a < 1e-3) {
      var p = v.toExponential((d || 3) - 1).split('e');
      return p[0] + '×10<sup>' + Number(p[1]) + '</sup>';
    }
    return Number(v.toPrecision(d || 3)).toString();
  }
  function pct(v) { return (100 * v).toFixed(v > 0.999 ? 2 : 1) + '%'; }
  function kv(rows) {
    return '<tbody>' + rows.map(function (r) {
      return '<tr><th>' + r[0] + '</th><td>' + r[1] + '</td></tr>';
    }).join('') + '</tbody>';
  }
  function opt(sel, items, label) {
    sel.innerHTML = items.map(function (v) {
      return '<option value="' + v + '">' + (label ? label(v) : v) + '</option>';
    }).join('');
  }

  opt($('vdSample'), SAMPLES.map(function (s) { return s.sample; }));
  opt($('vdMap'), AX.energy_map, function (v) { return MAPS[v] || v; });
  opt($('vdCutoff'), AX.cutoff_nm, function (v) { return v + ' nm'; });
  opt($('vdDG'), AX.dG_eV, function (v) { return (v > 0 ? '+' : v < 0 ? '−' : '') + Math.abs(v).toFixed(1) + ' eV'; });
  opt($('vdF110'), AX.f110, function (v) { return Math.round(100 * v) + '%'; });
  opt($('vdEps'), AX.eps.map(function (e) { return e.key; }), function (k) {
    var e = AX.eps.filter(function (x) { return x.key === k; })[0];
    return e.value + ' (' + k.replace('_axis', ' axis') + ')';
  });
  opt($('vdReactive'), D.reactive, function (v) { return REACT[v] || v; });

  function state() {
    return { sample: $('vdSample').value, map: $('vdMap').value, cutoff: +$('vdCutoff').value,
             dG: +$('vdDG').value, f110: +$('vdF110').value, eps: $('vdEps').value,
             reactive: $('vdReactive').value };
  }
  function setState(s) {
    $('vdSample').value = s.sample; $('vdMap').value = s.map; $('vdCutoff').value = String(s.cutoff);
    $('vdDG').value = String(s.dG); $('vdF110').value = String(s.f110); $('vdEps').value = s.eps;
    $('vdReactive').value = s.reactive;
  }
  /* the run's grid order: map, cutoff, dG, f110, eps */
  function point(s) {
    var i = AX.energy_map.indexOf(s.map);
    i = i * AX.cutoff_nm.length + AX.cutoff_nm.indexOf(s.cutoff);
    i = i * AX.dG_eV.length + AX.dG_eV.indexOf(s.dG);
    i = i * AX.f110.length + AX.f110.indexOf(s.f110);
    return i * AX.eps.length + AX.eps.map(function (e) { return e.key; }).indexOf(s.eps);
  }
  function sampleOf(name) { return SAMPLES.filter(function (x) { return x.sample === name; })[0]; }
  function caseOf(smp, p, r) {
    var c = smp.cases, j = D.reactive.indexOf(r), pools = {};
    D.pools.forEach(function (q, k) { pools[q] = c.pools[p][k]; });
    var n = c.sites[p][j];
    return { theta: c.theta[p], f_rec: c.f_rec[p], pools: pools, sites: n,
             below: !!c.below[p][j], tof: n > 0 ? smp.rate_co_umol_g_s / n : Infinity };
  }

  /* ------------------------------------------------------------ 3D */
  var P3 = null, S3 = null, tried = false;
  function mount3d() {
    if (tried) return;
    var hp = $('figParticle').querySelector('.figbox'), hs = $('figSlab').querySelector('.figbox');
    if (!hp.clientWidth) return;
    tried = true;
    if (!V || !V.available()) {
      hp.textContent = hs.textContent = 'This figure needs WebGL, which this browser does not provide.';
      return;
    }
    P3 = V.particle(hp);
    S3 = V.slab(hs);
    draw(false);
  }
  (function wait() { if (!tried) { mount3d(); requestAnimationFrame(wait); } })();

  /* ------------------------------------------------------- figures */
  var figState = {};
  function wireDownload(hostId) {
    var dl = $(hostId).querySelector('.figdl');
    if (!dl || dl.dataset.wired) return;
    dl.dataset.wired = '1';
    dl.addEventListener('click', function (ev) {
      var fmt = ev.target && ev.target.getAttribute('data-fmt');
      var st = figState[hostId];
      if (!fmt || !st) return;
      if (st.canvas) {
        st.canvas().toBlob(function (b) { K.saveBlob(b, st.name + '.png'); }, 'image/png');
      } else if (fmt === 'svg') {
        K.downloadSVG(st.svg, st.name);
      } else if (fmt === 'png') {
        K.downloadPNG(st.svg, st.name);
      } else if (fmt === 'csv') {
        K.saveBlob(new Blob([st.csv], { type: 'text/csv;charset=utf-8' }), st.name + '.csv');
      }
      dl.open = false;
    });
  }
  ['figParticle', 'figSlab', 'figTofRange'].forEach(wireDownload);
  figState.figParticle = { name: 'vacancy-particle', canvas: function () { return P3.stage.renderer.domElement; } };
  figState.figSlab = { name: 'vacancy-surface', canvas: function () { return S3.stage.renderer.domElement; } };

  function num(x) { return +x; }
  function draw(animate) {
    var s = state(), p = point(s), smp = sampleOf(s.sample), c = caseOf(smp, p, s.reactive);
    var ref = caseOf(sampleOf('R600'), p, s.reactive);
    var cap = D.capacity_umol_g[String(s.f110)];
    var pools = c.pools, inv = smp.inventory_umol_g;
    var surf = pools.bridging + pools.reconstructed_row + pools.basal_L1;
    var sub = pools.L1_subbridging + pools.subsurface_L2_4;
    var sm = smp.summary;
    var strong = smp.treatment_T_C >= D.strong_reduction_C;

    $('vdSampleNote').innerHTML = 'Measured: inventory ' + sig(inv, 4) + ' µmol O g⁻¹, r<sub>CO</sub> '
      + sig(smp.rate_co_umol_g_s) + ' µmol g⁻¹ s⁻¹. Treatment ' + smp.treatment_T_C + ' °C.';
    $('vdReactiveNote').innerHTML = s.reactive.indexOf('ISO') === 0
      ? 'N<sub>react</sub> = cθ(1−θ)<sup>' + s.reactive.slice(5) + '</sup>'
      : s.reactive === 'BRI' ? 'N<sub>react</sub> = cθ' : 'N<sub>react</sub> = cθ + N<sub>IPL,1</sub>';

    $('vdBasis').innerHTML = kv([
      ['Inventory, measured', sig(inv, 4) + ' µmol O g⁻¹'],
      ['Bulk, calculated', sig(pools.bulk) + ' µmol O g⁻¹ (' + pct(pools.bulk / inv) + ' of inventory)'],
      ['Subsurface (SBR layer 1, layers 2–4)', sig(sub) + ' µmol O g⁻¹'],
      ['Surface (BRI, reconstructed, IPL layer 1)', sig(surf) + ' µmol O g⁻¹'],
      ['θ, vacant fraction of BRI sites', sig(c.theta)],
      ['Reconstructed fraction of (1×2) cells', pct(c.f_rec)]
    ]);
    $('vdVerdict').innerHTML = c.below
      ? '<div class="verdict red">N<sub>react</sub> below threshold (&lt; 0.01 µmol g⁻¹ or &lt; 1% of BRI sites); '
        + 'case excluded from the TOF range.</div>' : '';
    $('vdKpis').innerHTML = kv([
      ['N<sub>react</sub>, calculated', sig(c.sites) + ' µmol g⁻¹'],
      ['TOF = r<sub>CO</sub>/N<sub>react</sub>', sig(c.tof) + ' s⁻¹'],
      ['TOF/TOF(R600), same parameters', ref.below ? '— (R600 below threshold)' : sig(c.tof / ref.tof)],
      ['TOF range, ' + sm.n_ok + ' of ' + sm.n_cases + ' cases', sig(num(sm.TOF_min_s_1)) + ' to ' + sig(num(sm.TOF_max_s_1)) + ' s⁻¹'],
      ['TOF, N<sub>react</sub> = 2.31 µmol g⁻¹', sig(num(sm.fixed_TOF_SI2a_s_1)) + ' s⁻¹']
    ].concat(strong ? [['Treatment ≥ ' + D.strong_reduction_C + ' °C', 'Ti₂O₃-(1×2) regime (Yuan 2024)']] : []));
    $('vdStatus').textContent = 'Parameter point ' + (p + 1) + ' of ' + NPTS;

    /* figure 3 */
    var pts = SAMPLES.map(function (q) {
      var cc = caseOf(q, p, s.reactive), m = q.summary;
      return { sample: q.sample, min: num(m.TOF_min_s_1), median: num(m.TOF_median_s_1),
               max: num(m.TOF_max_s_1), fixed: num(m.fixed_TOF_SI2a_s_1),
               nBelow: num(m.n_below_threshold), cur: cc.tof, curBelow: cc.below };
    });
    var svg = F.tofRange({ points: pts, nPoints: NPTS, fixed: D.fixed_sites_umol_g });
    $('figTofRange').querySelector('.figbox').innerHTML = svg;
    figState.figTofRange = { svg: svg, name: 'apparent-tof-range',
      csv: 'sample,TOF_min_s_1,TOF_median_s_1,TOF_max_s_1,TOF_fixed_s_1,below_threshold,TOF_this_scenario_s_1\n'
        + pts.map(function (q) { return [q.sample, q.min, q.median, q.max, q.fixed, q.nBelow, q.cur].join(','); }).join('\n') + '\n' };

    /* samples table */
    $('vdTable').innerHTML = SAMPLES.map(function (q, i) {
      var m = q.summary, x = pts[i];
      return '<tr' + (q.sample === s.sample ? ' class="hl"' : '') + '><td>' + q.sample + '</td><td>'
        + q.treatment_T_C + ' °C</td><td class="n">' + sig(q.inventory_umol_g, 4) + '</td><td class="n">'
        + sig(q.rate_co_umol_g_s) + '</td><td class="n">' + sig(x.cur) + (x.curBelow ? '*' : '')
        + '</td><td class="n">' + sig(x.min) + '</td><td class="n">' + sig(x.median) + '</td><td class="n">'
        + sig(x.max) + '</td><td class="n">' + sig(x.fixed) + '</td><td class="n">' + m.n_below_threshold
        + ' of ' + m.n_cases + '</td></tr>';
    }).join('');

    /* the 3D figures */
    if (P3) {
      P3.set({ f110: s.f110, pools: pools, animate: animate });
      var shown = P3.want;
      var r = D.explicit_depth_nm, d = D.diameters_nm;
      $('vdParticleCaption').innerHTML = '1 bead = ' + P3.dotUmol + ' µmol O g⁻¹ (calculated). Purple: bulk, '
        + Math.round(shown.bulk) + ' beads; orange: subsurface, ' + Math.round(shown.subsurface)
        + '; blue: surface, ' + Math.round(shown.surface) + '; beads in the removed octant not drawn. '
        + 'Blue tiles: (110), ' + Math.round(100 * s.f110) + '% of the area. Surface band ('
        + r + ' nm) drawn ' + Math.round(0.05 / (2 * r / d[0])) + '–'
        + Math.round(0.05 / (2 * r / d[d.length - 1])) + '× thicker than scale for '
        + d[0] + '–' + d[d.length - 1] + ' nm particles.';
    }
    if (S3) {
      S3.set({ theta: c.theta, f_rec: c.f_rec, reactive: s.reactive,
               x_basal: pools.basal_L1 / cap.basal_L1, x_sbr: pools.L1_subbridging / cap.L1_subbridging,
               x_l24: pools.subsurface_L2_4 / cap.subsurface_L2_4 });
      var n = S3.counts, P = V.PATCH;
      $('vdSlabCaption').innerHTML = P.NX + ' × ' + P.NY + ' (1×1) cells, ' + (P.NX * 0.6497).toFixed(1)
        + ' × ' + (P.NY * 0.2959).toFixed(1) + ' nm; unrelaxed bulk positions (a = 0.4594 nm, c = 0.2959 nm, '
        + 'u = 0.305). Vacant sites at the calculated site fractions: BRI ' + n.BRI + ' (θ = ' + sig(c.theta)
        + '), IPL layer 1 ' + n.IPL + ', SBR layer 1 ' + n.SBR + ', layers 2–4 ' + n.L24
        + '. Blue: surface vacancy; orange: subsurface vacancy; pulsing: counted in N<sub>react</sub> ('
        + (REACT[s.reactive] || s.reactive) + '), ' + n.reactive + '. Light blue: (1×2) added rows, '
        + n.cells + ' of ' + (P.NX / 2 * P.NY) + ' cells, schematic positions.';
    }
  }

  ['vdSample', 'vdMap', 'vdCutoff', 'vdDG', 'vdF110', 'vdEps', 'vdReactive'].forEach(function (id) {
    $(id).addEventListener('change', function () { stopPlay(); draw(true); });
  });
  $('vdReset').addEventListener('click', function () { stopPlay(); setState(START); draw(true); });

  var timer = null;
  function stopPlay() {
    if (timer) { clearInterval(timer); timer = null; $('vdPlay').textContent = 'Play the series'; }
  }
  $('vdPlay').addEventListener('click', function () {
    if (timer) { stopPlay(); return; }
    var names = SAMPLES.map(function (q) { return q.sample; }), k = 0;
    $('vdPlay').textContent = 'Stop';
    var step = function () {
      $('vdSample').value = names[k];
      draw(true);
      k += 1;
      if (k >= names.length) stopPlay();
    };
    step();
    timer = setInterval(step, 2200);
  });

  var grid = [
    ['Vacancy formation-energy set', AX.energy_map.join(', ')],
    ['Aggregate cutoff distance', AX.cutoff_nm.join(', ') + ' nm'],
    ['(1×2) reconstruction ΔG', AX.dG_eV.map(function (v) { return (v > 0 ? '+' : '') + v; }).join(', ') + ' eV per cell'],
    ['(110) area fraction', AX.f110.join(', ')],
    ['Static dielectric constant ε', AX.eps.map(function (e) { return e.value; }).join(', ')],
    ['Parameter points', String(NPTS)],
    ['Reactive-site definitions', String(D.reactive.length)],
    ['Particle diameters', D.diameters_nm.join(', ') + ' nm, equal mass']
  ];
  $('vdGrid').innerHTML = kv(grid);
  $('vdChecks').innerHTML = kv([
    ['Cases per sample', String(NPTS * D.reactive.length)],
    ['N<sub>react</sub> threshold', D.site_threshold.min_sites_umol_g + ' µmol g⁻¹ and '
      + Math.round(100 * D.site_threshold.min_fraction_of_bridging) + '% of BRI sites']
  ]);

  setState(START);
  draw(false);
  window.VacancyDistribution = { draw: draw, state: state, setState: setState, caseOf: caseOf, point: point,
    three: function () { return { particle: P3, slab: S3 }; } };
})();
