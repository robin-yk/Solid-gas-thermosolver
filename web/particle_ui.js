/* The depth-resolved workspace card.

   Runs web/particle.js on whatever the defect panel is currently set to
   - inventory, temperature, particle size, energy preset - and draws the
   two figures plus the three readings the class model could not make:
   where the oxygen is, whether the equilibrium is reachable inside the
   exposure, and whether local electroneutrality is allowed.

   It listens to the defect panel rather than owning inputs of its own,
   so the two engines are always answering the same question and the
   difference between them is visible instead of arranged. */

(function () {
  'use strict';

  var PA = window.Particle, FIG = window.ParticleFigures, KIT = window.FigKit;
  if (!PA || !FIG || !KIT) return;

  var $ = function (id) { return document.getElementById(id); };
  if (!$('figDepth') || !$('figWhere')) return;

  var DATA = null;
  try {
    DATA = JSON.parse($('statmech-data').textContent);
  } catch (e) { return; }

  var figState = {};
  function mount(hostId, name, svg) {
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

  function num(id, fallback) {
    var el = $(id);
    var v = el ? parseFloat(el.value) : NaN;
    return isFinite(v) ? v : fallback;
  }

  function fmt(v, n) {
    if (v == null || !isFinite(v)) return '—';
    if (v !== 0 && (Math.abs(v) < 1e-3 || Math.abs(v) >= 1e5)) {
      return v.toExponential(n == null ? 2 : n);
    }
    return v.toFixed(n == null ? 3 : n);
  }

  function rows(pairs) {
    return '<table class="kv">' + pairs.map(function (q) {
      return '<tr><th>' + q[0] + '</th><td>' + q[1] + '</td></tr>';
    }).join('') + '</table>';
  }

  function draw() {
    var d = DATA.defaults;
    var vo = num('smVO', d.VO_total_umol_g);
    var tC = num('smT', d.T_C);
    var modeEl = $('smGeoMode');
    var bet = modeEl && modeEl.value === 'bet';
    var geo = num('smGeoVal', d.particle_diameter_um);
    var presetEl = $('smPreset');
    var preset = presetEl ? presetEl.value : d.energetics_preset;
    /* a custom energy triple belongs to the class model's own controls;
       the depth profile needs a triple that decays, so it stays on a
       named preset and says so rather than inventing a decay length */
    var custom = preset === 'custom';
    if (custom) preset = d.energetics_preset;
    var exposure = num('smExp',
                       (DATA.co2_kinetics.defaults || {}).exposure_s || 300);

    var part, res;
    try {
      part = PA.particle(DATA, bet ? { bet_m2_g: geo, preset: preset }
                                   : { d_um: geo, preset: preset });
      res = PA.equilibrate(vo, tC, part);
    } catch (err) {
      $('pxSummary').innerHTML = '<p class="warn">' + err.message + '</p>';
      return;
    }

    var seg = part.segregation, kt = res.kt_eV;
    var layers = [];
    for (var k = 1; k <= 60; k++) {
      var z = k * seg.d110_nm;
      layers.push({ k: k, z_nm: z,
                    theta: PA.thetaOfMu(res.mu_eV,
                                        PA.formationEnergy(z, seg), kt,
                                        res.compensated, res.polaron_ratio) });
    }
    var depths = [], n = 44;
    var hi = Math.log(part.radius_nm * 0.999) / Math.LN10;
    for (var i = 0; i <= n; i++) {
      depths.push(Math.pow(10, -1 + i * (hi + 1) / n));
    }
    /* the figure quotes the one-nanometre shell, so put it on the grid
       exactly rather than reporting whichever sample landed near it */
    depths.push(1.0);
    depths.sort(function (a, b) { return a - b; });
    var shells = PA.shellSplit(res, part, depths);

    mount('figDepth', 'particle-depth',
          FIG.depth(FIG.depthData(res, part, layers)));
    mount('figWhere', 'particle-where',
          FIG.where(FIG.whereData(shells, res)));

    var ceiling = res.point_defect_ceiling_umol_g;
    var pin = res.limited_by
      ? 'pinned by ' + res.limited_by
      : 'not pinned: both competing phases are still off';
    $('pxSummary').innerHTML = '<p class="sub">The same inventory, resolved '
      + 'by depth instead of by site class. E(z) = ' + fmt(seg.e_bulk_eV, 2)
      + ' &minus; ' + fmt(seg.dE_seg_eV, 2) + ' exp(&minus;z/'
      + fmt(seg.xi_nm, 3) + ' nm), with the decay length derived from the '
      + seg.label + ' triple rather than assumed. ' + pin + '.</p>';

    $('pxState').innerHTML = rows([
      ['&mu;<sub>V</sub>', fmt(res.mu_eV, 6) + ' eV'],
      ['bridging row', fmt(res.surface.theta, 4) + ' &nbsp;<span class="dim">'
        + res.surface.phase + '</span>'],
      ['first subsurface layer', fmt(res.theta_first_layer, 5)],
      ['bulk plateau', fmt(res.theta_bulk_plateau, 5)],
      ['enrichment, layer 1 / bulk',
        fmt(res.theta_first_layer / res.theta_bulk_plateau, 1) + '&times;'],
      ['on the row', fmt(res.umol_surface, 2) + ' &micro;mol-O/g'],
      ['below it', fmt(res.umol_interior, 2) + ' &micro;mol-O/g'],
      ['extended defects', fmt(res.umol_extended, 2) + ' &micro;mol-O/g'],
      ['point-defect ceiling', ceiling == null ? '&mdash;'
        : fmt(ceiling, 1) + ' &micro;mol-O/g &nbsp;<span class="dim">this '
          + 'loading is ' + (100 * vo / ceiling).toFixed(0) + '% of it</span>']
    ]);

    var v = PA.verdict(tC, part, exposure, DATA);
    $('pxReach').innerHTML = rows([
      ['D', fmt(v.D_cm2_per_s, 3) + ' cm<sup>2</sup> s<sup>&minus;1</sup>'],
      ['99% homogenised in', fmt(v.t_homogenise_s, 3) + ' s'],
      ['exposure', fmt(exposure, 0) + ' s'],
      ['reachable', v.equilibrium_reachable ? 'yes, with '
        + (1 / v.ratio).toExponential(1) + ' to spare' : 'NO'],
      ['barrier that would break it',
        fmt(v.E_m_that_would_break_it_eV, 3) + ' eV &nbsp;<span class="dim">'
        + 'literature ' + fmt(v.E_m_eV, 2) + ', margin '
        + fmt(v.barrier_margin_eV, 3) + '</span>']
    ]);

    var a = PA.assess(res, part, DATA);
    $('pxCharge').innerHTML = rows([
      ['Debye length at the plateau', fmt(a.debye_nm_at_bulk, 3) + ' nm'],
      ['segregation depth &xi;', fmt(a.xi_nm, 3) + ' nm'],
      ['ratio', fmt(a.ratio_debye_over_xi, 2)],
      ['&epsilon;<sub>r</sub> assumed', fmt(a.eps_r, 0) + ' &nbsp;<span '
        + 'class="dim">50 &rarr; 170 moves the ratio to '
        + fmt(PA.assess(res, part, DATA, 50).ratio_debye_over_xi, 2) + ' &rarr; '
        + fmt(PA.assess(res, part, DATA, 170).ratio_debye_over_xi, 2)
        + '</span>']
    ]) + '<p class="readout" style="margin-top:8px">' + a.verdict + '</p>';
  }

  ['smVO', 'smT', 'smGeoVal', 'smGeoMode', 'smPreset', 'smExp']
    .forEach(function (id) {
      var el = $(id);
      if (el) {
        el.addEventListener('change', draw);
        el.addEventListener('input', draw);
      }
    });
  var run = $('smRun');
  if (run) run.addEventListener('click', draw);
  draw();
})();
