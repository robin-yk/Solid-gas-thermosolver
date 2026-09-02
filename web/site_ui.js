/* The surface workspace, and the switch between the two.

   The equilibrium half is driven by ti_solver_page.js and thermo_ui.js,
   which solve in the browser. This half needs no solver: the capacity is
   lattice geometry and the bound is one division, both mirrored in
   vacancy.js and gated against the Python module. */

(function () {
  'use strict';
  var D = window.SITE_DATA;
  var SF = window.SurfaceFigures, VI = window.Vacancy;
  function $(id) { return document.getElementById(id); }
  function n(v, d) {
    return (v == null || !isFinite(v)) ? '—' : Number(v).toFixed(d == null ? 2 : d);
  }
  function kpi(label, value) {
    return '<div class="kpi"><div class="sub">' + label + '</div><b>'
      + value + '</b></div>';
  }

  function showWorkspace(id) {
    ['ws-thermo', 'ws-surface'].forEach(function (w) {
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

  function init() {
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
    showWorkspace(location.hash.slice(1) === 'ws-surface' ? 'ws-surface' : 'ws-thermo');
  }

  if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', init);
  } else {
    init();
  }
})();
