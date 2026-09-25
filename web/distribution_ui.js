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
    BRI: 'Every bridging vacancy',
    ISO_z2: 'Isolated bridging, 2 intact neighbours',
    ISO_z4: 'Isolated bridging, 4 intact neighbours',
    ISO_z8: 'Isolated bridging, 8 intact neighbours',
    'BRI+BASAL': 'Bridging plus layer-1 in-plane'
  };
  var MAPS = { PAB: 'PAB', HAM: 'HAM', LI_SBR1: 'LI, sub-bridging layer 1', LI_L2: 'LI, layer 2' };

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

    $('vdSampleNote').innerHTML = 'Inventory ' + sig(inv) + ' µmol-O g⁻¹ · r<sub>CO</sub> '
      + sig(smp.rate_co_umol_g_s) + ' µmol g⁻¹ s⁻¹ · treated at ' + smp.treatment_T_C + ' °C';
    $('vdReactiveNote').textContent = s.reactive.indexOf('ISO') === 0
      ? 'N = cθ(1−θ)^z over the bridging sites, z = ' + s.reactive.slice(5) + '.'
      : s.reactive === 'BRI' ? 'N = cθ over the bridging sites.' : 'Bridging plus first-layer in-plane vacancies.';

    $('vdBasis').innerHTML = kv([
      ['Measured inventory', sig(inv, 4) + ' µmol-O g⁻¹'],
      ['Bulk', sig(pools.bulk) + ' µmol g⁻¹ (' + pct(pools.bulk / inv) + ')'],
      ['Subsurface', sig(sub) + ' µmol g⁻¹'],
      ['Surface', sig(surf) + ' µmol g⁻¹'],
      ['Bridging coverage θ', sig(c.theta)],
      ['Reconstructed cells', pct(c.f_rec)]
    ]);
    var ratio = c.tof / num(sm.fixed_TOF_SI2a_s_1);
    var verdict;
    if (c.below) {
      verdict = '<div class="verdict red">At this scenario ' + s.sample + ' has <b>'
        + sig(c.sites) + ' µmol g⁻¹</b> reactive sites, below the counting threshold. '
        + 'Its TOF is kept in the results but left out of the range.</div>';
    } else {
      verdict = '<div class="verdict">At this scenario ' + s.sample + ' turns over <b>'
        + sig(c.tof) + ' s⁻¹</b> on ' + sig(c.sites) + ' µmol g⁻¹ of reactive sites, '
        + sig(ratio, 2) + '× the fixed-denominator value.</div>';
    }
    $('vdVerdict').innerHTML = verdict;
    $('vdKpis').innerHTML = kv([
      ['Reactive sites', sig(c.sites) + ' µmol g⁻¹'],
      ['TOF, this scenario', sig(c.tof) + ' s⁻¹'],
      ['TOF / TOF(R600)', ref.below ? '— (R600 below threshold)' : sig(c.tof / ref.tof)],
      ['Range, ' + sm.n_ok + ' of ' + sm.n_cases + ' cases', sig(num(sm.TOF_min_s_1)) + ' to ' + sig(num(sm.TOF_max_s_1)) + ' s⁻¹'],
      ['Fixed 2.31 µmol g⁻¹', sig(num(sm.fixed_TOF_SI2a_s_1)) + ' s⁻¹']
    ].concat(strong ? [['Flag', 'treated above ' + D.strong_reduction_C + ' °C (Ti₂O₃-(1×2) expected)']] : []));
    $('vdStatus').textContent = 'Parameter point ' + (p + 1) + ' of ' + NPTS + '.';

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
      $('vdParticleCaption').innerHTML = 'One bead is ' + P3.dotUmol + ' µmol g⁻¹ of vacancies: '
        + Math.round(shown.bulk) + ' bulk (purple), ' + Math.round(shown.subsurface) + ' subsurface (orange), '
        + Math.round(shown.surface) + ' surface (blue) beads, less the cut octant. Blue tiles are the '
        + Math.round(100 * s.f110) + '% (110) share. The surface band holds the ' + r
        + ' nm of explicit trilayers and is drawn about ' + Math.round(0.05 / (2 * r / d[0]))
        + ' to ' + Math.round(0.05 / (2 * r / d[d.length - 1])) + ' times thicker than scale. Drag to turn.';
    }
    if (S3) {
      S3.set({ theta: c.theta, f_rec: c.f_rec, reactive: s.reactive,
               x_basal: pools.basal_L1 / cap.basal_L1, x_sbr: pools.L1_subbridging / cap.L1_subbridging,
               x_l24: pools.subsurface_L2_4 / cap.subsurface_L2_4 });
      var n = S3.counts, P = V.PATCH;
      $('vdSlabCaption').innerHTML = P.NX + ' × ' + P.NY + ' surface cells (' + (P.NX * 0.6497).toFixed(1)
        + ' × ' + (P.NY * 0.2959).toFixed(1) + ' nm), four trilayers. Each site is missing when its fixed random '
        + 'rank is below the computed fraction: here ' + n.BRI + ' bridging (θ = ' + sig(c.theta) + '), '
        + n.IPL + ' in-plane, ' + n.SBR + ' sub-bridging and ' + n.L24 + ' deeper; '
        + n.reactive + ' pulse as reactive under "' + (REACT[s.reactive] || s.reactive).toLowerCase() + '". '
        + n.cells + ' of ' + (P.NX / 2 * P.NY) + ' (1×2) cells carry an added row (light blue, schematic). '
        + 'Unrelaxed bulk positions.';
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
    ['Vacancy energy map', AX.energy_map.join(', ')],
    ['Aggregate cutoff', AX.cutoff_nm.join(', ') + ' nm'],
    ['Reconstruction ΔG', AX.dG_eV.map(function (v) { return (v > 0 ? '+' : '') + v; }).join(', ') + ' eV'],
    ['(110) share', AX.f110.join(', ')],
    ['Dielectric constant', AX.eps.map(function (e) { return e.value; }).join(', ')],
    ['Parameter points', String(NPTS)],
    ['Reactive-site definitions', String(D.reactive.length)],
    ['Diameters averaged', D.diameters_nm.join(', ') + ' nm']
  ];
  $('vdGrid').innerHTML = kv(grid);
  $('vdChecks').innerHTML = kv([
    ['Cases per sample', String(NPTS * D.reactive.length)],
    ['Site threshold', D.site_threshold.min_sites_umol_g + ' µmol g⁻¹ and '
      + Math.round(100 * D.site_threshold.min_fraction_of_bridging) + '% of bridging sites'],
    ['Model tests', 'pilot/titania-super-multiscale/tests']
  ]);

  setState(START);
  draw(false);
  window.VacancyDistribution = { draw: draw, state: state, setState: setState, caseOf: caseOf, point: point,
    three: function () { return { particle: P3, slab: S3 }; } };
})();
