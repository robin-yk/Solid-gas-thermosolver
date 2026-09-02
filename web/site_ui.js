/* The surface workspace, and the switch between the two.

   The equilibrium half is driven by ti_solver_page.js and thermo_ui.js,
   which solve in the browser. This half needs no solver: the capacity is
   lattice geometry and the bound is one division, both mirrored in
   vacancy.js and gated against the Python module. */

(function () {
  'use strict';
  var D = window.SITE_DATA;
  var SF = window.SurfaceFigures, VI = window.Vacancy;
  var PF = window.PopulationFigures, PM = window.Population;
  function $(id) { return document.getElementById(id); }
  function n(v, d) {
    return (v == null || !isFinite(v)) ? '—' : Number(v).toFixed(d == null ? 2 : d);
  }
  /* Inventories run from 8e-4 to 766 umol-O/g. Two fixed decimals print
     the smallest of them as 0.00, which reads as "there are none" and
     contradicts the figure directly above the table. */
  function amt(v) {
    if (v == null || !isFinite(v)) return '—';
    if (v === 0) return '0';
    return Math.abs(v) >= 0.01 ? v.toFixed(2) : v.toExponential(1);
  }
  function kpi(label, value) {
    return '<div class="kpi"><div class="sub">' + label + '</div><b>'
      + value + '</b></div>';
  }

  function showWorkspace(id) {
    ['ws-thermo', 'ws-surface', 'ws-population'].forEach(function (w) {
      var el = $(w);
      if (el) el.style.display = (w === id) ? '' : 'none';
    });
    Array.prototype.forEach.call(document.querySelectorAll('.wstab'), function (b) {
      b.classList.toggle('active', b.getAttribute('data-ws') === id);
    });
    Array.prototype.forEach.call(document.querySelectorAll('.workspace-head'), function (h) {
      h.classList.toggle('active', h.id === 'head-' + id.replace('ws-', ''));
    });
    if (location.hash.slice(1) !== id) history.replaceState(null, '', '#' + id);
  }

  /* --------------------------------------------------------- surface */

  function inputs() {
    var d = parseFloat($('sfD').value);
    var bet = parseFloat($('sfBET').value);
    var inv = parseFloat($('sfV').value);
    return {
      d_um: isFinite(d) && d > 0 ? d : null,
      bet: isFinite(bet) && bet > 0 ? bet : null,
      inv: isFinite(inv) && inv > 0 ? inv : null
    };
  }

  function series(g, o) {
    /* the bar figure wants a ladder of inventories; keep the one the user
       typed in it so the marked bar is always the one being read out */
    var base = [30, 40, 60, 95, 130, 150, 190, 220];
    if (o.inv && base.indexOf(o.inv) < 0) base.push(o.inv);
    base.sort(function (a, b) { return a - b; });
    return base.map(function (v) {
      var b = VI.bound(g, v, o.d_um, o.bet);
      return {
        diameter_um: o.d_um || 0, area_m2_g: b.area_m2_g,
        measured_umol_O_g: v,
        bridging_capacity_umol_O_g: b.bridging_capacity_umol_o_g,
        reconstruction_deficiency_umol_O_g: b.reconstruction_deficiency_umol_o_g,
        inventory_over_capacity: b.inventory_over_capacity,
        surface_fraction_max_pct: b.surface_fraction_max * 100,
        below_surface_fraction_min_pct: b.below_surface_fraction_min * 100
      };
    });
  }

  function drawSurface() {
    var g = D.surface.geometry_card, o = inputs();
    var st = $('sfStatus');
    if (!o.inv || (!o.d_um && !o.bet)) {
      st.textContent = 'give a positive oxygen removal and either a diameter '
        + 'or a BET area';
      return;
    }
    st.textContent = o.bet
      ? 'area from the measured BET value; the diameter is ignored'
      : 'area from the smooth-sphere estimate 6/(ρd), a lower bound on a rough particle';

    var rows = series(g, o);
    /* every row shares one area here, so the figures key off it directly */
    rows.forEach(function (r) { r.diameter_um = rows[0].diameter_um; });
    var d0 = rows[0].diameter_um;
    $('figCapacity').innerHTML = SF.capacity(SF.capacityData(rows, d0, o.inv));
    $('figBound').innerHTML = SF.bound(SF.boundData(rows, d0, o.inv));

    var b = VI.bound(g, o.inv, o.d_um, o.bet);
    $('sfKpis').innerHTML = [
      kpi('Surface area', n(b.area_m2_g, 3) + ' m² g⁻¹'),
      kpi('Bridging-O capacity', n(b.bridging_capacity_umol_o_g, 2) + ' µmol-O g⁻¹'),
      kpi('(1×2) deficiency', n(b.reconstruction_deficiency_umol_o_g, 2) + ' µmol-O g⁻¹'),
      kpi('Measured / capacity', n(b.inventory_over_capacity, 1) + '×'),
      kpi('Surface share, at most', n(b.surface_fraction_max * 100, 1) + '%'),
      kpi('Below surface, at least', n(b.below_surface_fraction_min * 100, 1) + '%')
    ].join('');
    $('sfDerived').innerHTML =
      '<tr><td>x in TiO<sub>2−x</sub></td><td class="k2">'
      + n(b.x_in_TiO2mx, 5) + '</td></tr>'
      + '<tr><td>Ti reduced to Ti³⁺</td><td class="k2">'
      + n(b.Ti3_fraction * 100, 2) + '%</td></tr>'
      + '<tr><td>Whole bridging row stripped, at most</td><td class="k2">'
      + n(b.surface_fraction_max_full_row * 100, 1) + '%</td></tr>';
  }

  /* ------------------------------------------------------ population */

  var REPORTED = { ns: PM.REPORTED.ns, E: PM.REPORTED.E_loss,
                   lognu: PM.REPORTED.log10_nu };

  function popKind() {
    var el = $('pnKind');
    return el ? el.value : 'log';
  }

  function popParams() {
    return { ns: parseFloat($('pnNs').value),
             E: parseFloat($('pnE').value),
             lognu: parseFloat($('pnNu').value) };
  }

  function drawPopulation() {
    var p = popParams(), st = $('pnStatus');
    var series = D.population.series;
    if (!(p.ns > 0) || !(p.E >= 0) || !isFinite(p.lognu)) {
      st.textContent = 'give a positive capacity, a non-negative barrier and a finite log prefactor';
      return;
    }
    var nu = Math.pow(10, p.lognu);
    var rows = PM.partition(series, p.ns, p.E, nu);
    $('figPopulation').innerHTML = PF.population(PF.partitionData(rows));
    $('figRates').innerHTML = PF.rates(PF.ratesData(rows));

    var ssq = PM.objective(series, p.ns, p.E, nu, popKind());
    var ratio = PM.rateRatio(series, p.ns, p.E, nu);
    var meas = PM.measuredRatio(series);
    var isoMax = Math.max.apply(null, rows.map(function (q) { return q.n_iso_umol_g; }));
    var isoMin = Math.min.apply(null, rows.map(function (q) { return q.n_iso_umol_g; }));
    var closure = Math.max.apply(null, rows.map(function (q) {
      return Math.abs(q.closure_umol_g); }));

    var lo = D.population.ns_bounds[0], hi = D.population.ns_bounds[1];
    st.textContent = (p.ns < lo || p.ns > hi)
      ? 'n_s is outside the geometric bounds ' + lo + '–' + hi + ' µmol-O g⁻¹'
      : 'n_s sits between the bridging row and the bridging plus first subsurface layer';

    $('pnKpis').innerHTML = [
      kpi('Objective (' + popKind() + ')', ssq.toExponential(3)),
      kpi('R800 / R600, fitted', n(ratio, 3)),
      kpi('R800 / R600, measured', n(meas, 3)),
      kpi('Isolated, highest', n(isoMax, 2) + ' µmol-O g⁻¹'),
      kpi('Isolated, lowest', n(isoMin, 3) + ' µmol-O g⁻¹'),
      kpi('Mass balance', closure === 0 ? 'exact by construction'
          : n(closure, 9) + ' µmol-O g⁻¹')
    ].join('');

    $('pnTable').innerHTML = rows.map(function (q) {
      return '<tr><td>' + q.sample + '</td><td>' + q.treatment + '</td>'
        + '<td class="k2">' + n(q.nV_total_umol_g, 2) + '</td>'
        + '<td class="k2">' + n(q.r_CO_measured, 4) + '</td>'
        + '<td class="k2">' + n(q.r_CO_fitted, 4) + '</td>'
        + '<td class="k2">' + amt(q.n_iso_umol_g) + '</td>'
        + '<td class="k2">' + amt(q.n_assoc_umol_g) + '</td>'
        + '<td class="k2">' + amt(q.n_below_umol_g) + '</td></tr>';
    }).join('');

    var below = rows.every(function (q) { return q.below_fraction > 0.5; });
    $('pnIdent').innerHTML =
      '<tr><td>ν<sub>eff</sub></td><td class="k2">'
      + nu.toExponential(2) + ' s⁻¹</td></tr>'
      + '<tr><td>Most of the inventory below the active region</td>'
      + '<td class="k2">' + (below ? 'yes, every sample' : 'not every sample') + '</td></tr>'
      + '<tr><td>Isolated population, R1000 over its maximum</td><td class="k2">'
      + n(isoMin / isoMax * 100, 2) + '%</td></tr>'
      + '<tr><td>Measured inventory, R1000 over A600</td><td class="k2">'
      + n(series[series.length - 1].nV_total_umol_g / series[0].nV_total_umol_g, 1)
      + '×</td></tr>';
  }

  function fitPopulation() {
    var st = $('pnStatus');
    st.textContent = 'fitting…';
    setTimeout(function () {
      var f = PM.fit(D.population.series, popKind(),
                     D.population.ns_bounds);
      $('pnNs').value = f.ns.toFixed(3);
      $('pnE').value = f.E_loss.toFixed(4);
      $('pnNu').value = f.log10_nu.toFixed(4);
      drawPopulation();
    }, 10);
  }

  function init() {
    ['pnNs', 'pnE', 'pnNu', 'pnKind'].forEach(function (id) {
      var el = $(id);
      if (el) el.addEventListener('input', drawPopulation);
      if (el && el.tagName === 'SELECT') el.addEventListener('change', drawPopulation);
    });
    if ($('pnFit')) $('pnFit').addEventListener('click', fitPopulation);
    if ($('pnReset')) $('pnReset').addEventListener('click', function () {
      $('pnNs').value = REPORTED.ns.toFixed(3);
      $('pnE').value = REPORTED.E.toFixed(4);
      $('pnNu').value = REPORTED.lognu.toFixed(4);
      drawPopulation();
    });
    ['sfD', 'sfBET', 'sfV'].forEach(function (id) {
      var el = $(id);
      if (el) el.addEventListener('input', drawSurface);
    });
    Array.prototype.forEach.call(document.querySelectorAll('.wstab'), function (b) {
      b.addEventListener('click', function () {
        showWorkspace(b.getAttribute('data-ws'));
      });
    });
    drawSurface();
    drawPopulation();
    var want = location.hash.slice(1);
    showWorkspace(['ws-surface', 'ws-population'].indexOf(want) >= 0
                  ? want : 'ws-thermo');
  }

  if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', init);
  } else {
    init();
  }
})();
