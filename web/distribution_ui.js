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
  var K = window.FigKit, F = window.DistributionFigures;
  var AX = D.axes, SAMPLES = D.samples;
  var NPTS = AX.energy_map.length * AX.cutoff_nm.length * AX.dG_eV.length * AX.f110.length * AX.eps.length;
  var START = { sample: 'R600', map: 'LI_SBR1', cutoff: 0.28, dG: 0, f110: 0.75, eps: 'a_axis', reactive: 'BRI' };
  var REACT = {
    BRI: 'All bridging-oxygen vacancies',
    ISO_z2: 'Bridging vacancies with 2 intact neighbours',
    ISO_z4: 'Bridging vacancies with 4 intact neighbours',
    ISO_z8: 'Bridging vacancies with 8 intact neighbours',
    'BRI+BASAL': 'Bridging + layer-1 in-plane vacancies'
  };
  var MAPS = { PAB: 'Pabisiak 2007 · GGA', HAM: 'Hameeuw 2006 · LDA',
    LI_SBR1: 'Li 2015 · sX · subsurface → SBR layer 1',
    LI_L2: 'Li 2015 · sX · subsurface → layer 2' };
  var SCENARIOS = [
    {id:'S1', map:'LI_SBR1', cutoff:0.28, dG:0, label:'Highest surface enrichment'},
    {id:'S2', map:'LI_SBR1', cutoff:0.28, dG:0.4, label:'Higher reconstruction cost'},
    {id:'S3', map:'HAM', cutoff:0.28, dG:0, label:'Intermediate surface enrichment'},
    {id:'S4', map:'HAM', cutoff:0.34, dG:0, label:'Lower surface enrichment'},
    {id:'S5', map:'PAB', cutoff:0.34, dG:0, label:'Lowest surface enrichment'}
  ];

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
  opt($('vdScenario'), SCENARIOS.map(function(q) { return q.id; }).concat(['custom']), function(id) {
    return id === 'custom' ? 'Custom parameters' : id + ': ' + SCENARIOS.filter(function(q) { return q.id === id; })[0].label;
  });
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

  /* ----------------------------------------------- rendered panels */
  /* Blender renders: the grains (four hover states and a mask whose red,
     green and blue pixels are the surface, subsurface and bulk sections)
     and one (110) surface per sample at the starting parameters. */
  var RD = JSON.parse($('render-data').textContent);
  var LAYER = { 1: 'surface', 2: 'subsurface', 3: 'bulk' };
  var mask = null;
  (function () {
    var im = new Image();
    im.onload = function () {
      var c = document.createElement('canvas');
      c.width = im.width; c.height = im.height;
      var g = c.getContext('2d');
      g.drawImage(im, 0, 0);
      mask = { w: im.width, h: im.height, d: g.getImageData(0, 0, im.width, im.height).data };
    };
    im.src = RD.mask;
  })();
  var defaultParticle = $('vdParticleImg').src;
  function layerAt(ev) {
    if (!mask) return null;
    var r = $('vdParticleImg').getBoundingClientRect();
    var x = Math.floor((ev.clientX - r.left) / r.width * mask.w);
    var y = Math.floor((ev.clientY - r.top) / r.height * mask.h);
    if (x < 0 || y < 0 || x >= mask.w || y >= mask.h) return null;
    var i = 4 * (y * mask.w + x), d = mask.d;
    if (d[i + 3] < 128) return null;
    var k = d[i] > 128 ? 1 : d[i + 1] > 128 ? 2 : d[i + 2] > 128 ? 3 : 0;
    return LAYER[k] || null;
  }
  var shownLayer = null, current = null;
  function showLayer(name, ev) {
    var tip = $('vdTip');
    if (name !== shownLayer) {
      $('vdParticleImg').src = name ? RD.particle[name] : defaultParticle;
      shownLayer = name;
      slab.highlight(name);
    }
    if (!name || !current) { tip.hidden = true; return; }
    var p = current.pools, inv = current.inv, v, what;
    if (name === 'surface') { v = p.bridging + p.reconstructed_row + p.basal_L1; what = 'BRI, reconstructed and layer-1 IPL sites'; }
    else if (name === 'subsurface') { v = p.L1_subbridging + p.subsurface_L2_4; what = 'layer-1 SBR and layers 2–4'; }
    else { v = p.bulk; what = 'below 1.30 nm'; }
    tip.innerHTML = '<b>' + name.charAt(0).toUpperCase() + name.slice(1) + ', calculated</b>'
      + sig(v) + ' µmol O g⁻¹ (' + pct(v / inv) + ' of ' + current.sample + ' inventory)<br>' + what;
    var r = $('vdParticleStage').getBoundingClientRect();
    tip.hidden = false;
    var x = ev.clientX - r.left + 14, y = ev.clientY - r.top + 14;
    tip.style.left = Math.max(4, Math.min(x, r.width - tip.offsetWidth - 4)) + 'px';
    tip.style.top = Math.max(4, Math.min(y, r.height - tip.offsetHeight - 4)) + 'px';
  }
  $('vdParticleStage').addEventListener('mousemove', function (ev) { showLayer(layerAt(ev), ev); });
  $('vdParticleStage').addEventListener('mouseleave', function () { showLayer(null); });
  $('vdParticleStage').addEventListener('click', function (ev) { showLayer(layerAt(ev), ev); });

  /* The (110) surface, drawn live from the CIF atom list at the site
     fractions of the case shown. */
  var SC = window.SlabCanvas, slab = SC.mount($('vdSlabCanvas'), RD.slab), slabShown = null;
  $('vdSlabPause').hidden = slab.reduced;
  $('vdSlabPause').addEventListener('click', function () {
    var p = slab.pause(!slab.paused());
    this.textContent = p ? 'Play' : 'Pause';
    this.setAttribute('aria-pressed', String(p));
  });

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
      if (fmt === 'svg') {
        K.downloadSVG(st.svg, st.name);
      } else if (fmt === 'png') {
        K.downloadPNG(st.svg, st.name);
      } else if (fmt === 'csv') {
        K.saveBlob(new Blob([st.csv], { type: 'text/csv;charset=utf-8' }), st.name + '.csv');
      }
      dl.open = false;
    });
  }
  wireDownload('figTofRange');

  function num(x) { return +x; }
  function draw(animate) {
    var s = state(), p = point(s), smp = sampleOf(s.sample), c = caseOf(smp, p, s.reactive);
    var preset = SCENARIOS.filter(function(q) {
      return q.map === s.map && q.cutoff === s.cutoff && q.dG === s.dG && s.f110 === 0.75 && s.eps === 'a_axis' && s.reactive === 'BRI';
    })[0];
    $('vdScenario').value = preset ? preset.id : 'custom';
    var ref = caseOf(sampleOf('R600'), p, s.reactive);
    var cap = D.capacity_umol_g[String(s.f110)];
    var pools = c.pools, inv = smp.inventory_umol_g;
    var surf = pools.bridging + pools.reconstructed_row + pools.basal_L1;
    var sub = pools.L1_subbridging + pools.subsurface_L2_4;
    var strong = smp.treatment_T_C >= D.strong_reduction_C;

    $('vdSampleNote').innerHTML = 'Measured oxygen removal: ' + sig(inv, 4) + ' µmol O g⁻¹; CO formation rate: '
      + sig(smp.rate_co_umol_g_s) + ' µmol g⁻¹ s⁻¹. Treatment ' + smp.treatment_T_C + ' °C.';
    $('vdReactiveNote').innerHTML = s.reactive.indexOf('ISO') === 0
      ? 'N<sub>react</sub> = cθ(1−θ)<sup>' + s.reactive.slice(5) + '</sup>. Isolated-vacancy estimate with independent site occupancy.'
      : s.reactive === 'BRI' ? 'N<sub>react</sub> = cθ. All unreconstructed bridging-oxygen vacancies are counted as reactive sites.' : 'N<sub>react</sub> = cθ + N<sub>IPL,1</sub>';

    $('vdBasis').innerHTML = kv([
      ['Total oxygen removal, measured', sig(inv, 4) + ' µmol O g⁻¹'],
      ['Bulk, calculated', sig(pools.bulk) + ' µmol O g⁻¹ (' + pct(pools.bulk / inv) + ' of inventory)'],
      ['Subsurface (SBR layer 1, layers 2–4)', sig(sub) + ' µmol O g⁻¹'],
      ['Surface (BRI, reconstructed, IPL layer 1)', sig(surf) + ' µmol O g⁻¹'],
      ['Bridging-oxygen vacancy fraction θ (unreconstructed sites)', sig(c.theta)],
      ['Reconstructed fraction of (1×2) cells', pct(c.f_rec)]
    ]);
    $('vdVerdict').innerHTML = c.below
      ? '<div class="verdict red">N<sub>react</sub> below threshold (&lt; 0.01 µmol g⁻¹ or &lt; 1% of BRI sites); '
        + 'TOF uses a small calculated site concentration.</div>' : '';
    $('vdKpis').innerHTML = kv([
      ['Site concentration under selected definition', sig(c.sites) + ' µmol g⁻¹'],
      ['Apparent CO TOF = r<sub>CO</sub>/N<sub>react</sub>', sig(c.tof) + ' s⁻¹'],
      ['TOF/TOF(R600), same parameters', ref.below ? '— (R600 below threshold)' : sig(c.tof / ref.tof)]
    ].concat(strong ? [['Treatment ≥ ' + D.strong_reduction_C + ' °C', 'Ti₂O₃-(1×2) regime (Yuan et al. 2024)']] : []));
    $('vdStatus').textContent = 'Parameter point ' + (p + 1) + ' of ' + NPTS;

    /* figure 3 */
    var pts = SAMPLES.map(function (q) {
      return { sample: q.sample, scenarios: SCENARIOS.map(function (sc) {
        var ps = Object.assign({}, sc, {f110:0.75, eps:'a_axis', reactive:'BRI'});
        var cc = caseOf(q, point(ps), 'BRI');
        return {id:sc.id, tof:cc.tof, sites:cc.sites, below:cc.below};
      }) };
    });
    var svg = F.tofRange({ points: pts });
    $('figTofRange').querySelector('.figbox').innerHTML = svg;
    figState.figTofRange = { svg: svg, name: 'apparent-tof-five-scenarios',
      csv: 'sample,scenario,energy_set,cutoff_nm,dG_eV,f110,epsilon,reactive_definition,Nreact_umol_g,TOF_s_1,below_threshold\n'
        + pts.map(function (q) { return q.scenarios.map(function (cc, j) {
          var sc = SCENARIOS[j];
          return [q.sample,sc.id,sc.map,sc.cutoff,sc.dG,0.75,64,'BRI',cc.sites,cc.tof,cc.below].join(',');
        }).join('\n'); }).join('\n') + '\n' };

    /* samples table */
    $('vdTable').innerHTML = SAMPLES.map(function (q, i) {
      var x = pts[i];
      return '<tr' + (q.sample === s.sample ? ' class="hl"' : '') + '><td>' + q.sample + '</td><td>'
        + q.treatment_T_C + ' °C</td><td class="n">' + sig(q.inventory_umol_g, 4) + '</td><td class="n">'
        + sig(q.rate_co_umol_g_s) + '</td>' + x.scenarios.map(function (cc) {
          return '<td class="n">' + sig(cc.tof) + (cc.below ? '*' : '') + '</td>';
        }).join('') + '</tr>';
    }).join('');

    /* the rendered panels */
    current = { pools: pools, inv: inv, sample: s.sample };
    var nv = slab.set(SC.fractions(D, smp, p, s.f110)), cells = RD.slab.cells;
    slabShown = s.sample;
    $('vdSlabCaption').innerHTML = '<b>Rutile (110), ' + s.sample + ', selected parameters.</b> ' + cells[0]
      + ' × ' + cells[1] + ' cells, four trilayers; atoms from the rutile CIF (P4₂/mnm, a = 0.4594 nm, '
      + 'c = 0.2959 nm, x(O) = 0.3048), unrelaxed. Vacancies shown at the calculated site fractions: '
      + (nv.BRI + nv.IPL) + ' in the top layer (blue rings), ' + (nv.SBR + nv.L24) + ' below it (orange). '
      + 'Pink: Ti³⁺, placed on the two Ti nearest each vacancy for illustration. The calculation gives Ti³⁺ populations by layer. '
      + 'Rotation and atomic motion are visual effects.';
  }

  ['vdSample', 'vdMap', 'vdCutoff', 'vdDG', 'vdF110', 'vdEps', 'vdReactive'].forEach(function (id) {
    $(id).addEventListener('change', function () { stopPlay(); draw(true); });
  });
  $('vdScenario').addEventListener('change', function() {
    var q = SCENARIOS.filter(function(q) { return q.id === $('vdScenario').value; })[0];
    if (!q) return;
    stopPlay();
    setState({sample:state().sample, map:q.map, cutoff:q.cutoff, dG:q.dG, f110:0.75, eps:'a_axis', reactive:'BRI'});
    draw(true);
  });
  $('vdReset').addEventListener('click', function () { stopPlay(); setState(START); draw(true); });

  var timer = null;
  function stopPlay() {
    if (timer) { clearInterval(timer); timer = null; $('vdPlay').textContent = 'Step through samples'; }
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
    ['Vacancy segregation energies', AX.energy_map.map(function(k) { return MAPS[k]; }).join('; ')],
    ['Vacancy–vacancy interaction cutoff', AX.cutoff_nm.join(', ') + ' nm'],
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
    scenarios: SCENARIOS, renders: RD, layerAt: layerAt, slab: slab, slabShown: function () { return slabShown; } };
})();
