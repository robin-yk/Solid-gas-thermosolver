/* The page. It displays; it does not calculate.

   Every number here was produced by scripts/reproduce_paper.py and is
   inlined as SITE_DATA at build time. There is no solver in the browser
   and no second implementation of any physics in this repository: if a
   value is on this page, it came out of the Python package, and running
   the script again is what changes it. The only arithmetic below is
   formatting and picking which row to show. */

(function () {
  'use strict';
  var D = window.SITE_DATA;
  var TF = window.ThermoFigures, SF = window.SurfaceFigures;
  function $(id) { return document.getElementById(id); }
  function chem(s) { return window.FigKit.chem(String(s)); }
  function n(v, d) {
    return (v == null || !isFinite(v)) ? '—' : Number(v).toFixed(d == null ? 2 : d);
  }

  /* ------------------------------------------------------- workspaces */

  function showWorkspace(id) {
    ['ws-equilibrium', 'ws-surface'].forEach(function (w) {
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

  /* ----------------------------------------------------- equilibrium */

  function equilibriumTemperature() {
    var sel = $('eqT');
    return sel ? Number(sel.value) : 600;
  }

  function rowAt(T_C) {
    var rows = D.equilibrium.rows;
    return rows.reduce(function (a, b) {
      return Math.abs(b.T_C - T_C) < Math.abs(a.T_C - T_C) ? b : a;
    });
  }

  function drawEquilibrium() {
    var T = equilibriumTemperature();
    var r = rowAt(T);
    var rows = D.equilibrium.rows;

    /* the drawers take the whole payload object, keyed by figure name */
    $('figConversion').innerHTML = TF.conversion(
      { conversion: TF.conversionData(rows, T, 'TiO2', r.conversion_CO2_pct) });
    $('figMargin').innerHTML = TF.margin({ margin: TF.marginData(rows, T) });
    var curve = D.equilibrium.boundary;
    var at = curve.reduce(function (a, b) {
      return Math.abs(b.T_C - T) < Math.abs(a.T_C - T) ? b : a;
    });
    $('figBoundary').innerHTML = TF.boundary(
      { boundary: TF.boundaryData(curve, null, T, at) });

    var costs = r.inactive_phase_reduced_costs || {};
    var near = Object.keys(costs).reduce(function (a, b) {
      return (a === null || costs[b].per_formula_kJ < costs[a].per_formula_kJ) ? b : a;
    }, null);

    $('eqKpis').innerHTML = [
      kpi('CO₂ conversion', n(r.conversion_CO2_pct, 1) + '%'),
      kpi('Stable solid', chem(r.active_condensed_phases.join(' + '))),
      kpi('Nearest reduced phase', chem(near || '—')),
      kpi('Margin to it', near ? n(costs[near].per_formula_kJ, 1) + ' kJ' : '—'),
      kpi('CO₂ needed to hold rutile', n(at.y_CO2 * 100, 4) + ' mol%')
    ].join('');

    var g = r.gas_fractions;
    $('eqGas').innerHTML = ['CO2', 'H2', 'CO', 'H2O', 'N2'].map(function (s) {
      return '<tr><td>' + chem(s) + '</td><td class="k2">'
        + n(g[s] * 100, 3) + '</td></tr>';
    }).join('');
  }

  function kpi(label, value) {
    return '<div class="kpi"><div class="sub">' + label + '</div><b>'
      + value + '</b></div>';
  }

  /* --------------------------------------------------------- surface */

  function surfaceDiameter() {
    var sel = $('sfD');
    return sel ? Number(sel.value) : 0.9;
  }
  function surfaceInventory() {
    var sel = $('sfV');
    return sel ? Number(sel.value) : 95;
  }

  function drawSurface() {
    var d = surfaceDiameter(), inv = surfaceInventory();
    var rows = D.surface.rows;
    $('figCapacity').innerHTML = SF.capacity(SF.capacityData(rows, d, inv));
    $('figBound').innerHTML = SF.bound(SF.boundData(rows, d, inv));

    var g = D.surface.by_diameter[String(d)];
    var row = rows.filter(function (r) {
      return Math.abs(r.diameter_um - d) < 1e-9
        && Math.abs(r.measured_umol_O_g - inv) < 1e-9;
    })[0];
    if (!row) return;

    $('sfKpis').innerHTML = [
      kpi('Surface area', n(g.area_m2_g, 3) + ' m² g⁻¹'),
      kpi('Bridging-O capacity', n(g.bridging_capacity_umol_O_g, 2) + ' µmol-O g⁻¹'),
      kpi('(1×2) deficiency', n(g.reconstruction_deficiency_umol_O_g, 2) + ' µmol-O g⁻¹'),
      kpi('Measured / capacity', n(row.inventory_over_capacity, 1) + '×'),
      kpi('Surface share, at most', n(row.surface_fraction_max_pct, 1) + '%'),
      kpi('Below surface, at least', n(row.below_surface_fraction_min_pct, 1) + '%')
    ].join('');

    $('sfDerived').innerHTML =
      '<tr><td>x in TiO<sub>2−x</sub></td><td class="k2">'
      + n(row.x_in_TiO2mx, 5) + '</td></tr>'
      + '<tr><td>Ti reduced to Ti³⁺</td><td class="k2">'
      + n(row.Ti3_fraction_pct, 2) + '%</td></tr>';
  }

  /* ------------------------------------------------------------ wire */

  function init() {
    if (!D) return;
    var eqT = $('eqT');
    D.equilibrium.rows.forEach(function (r) {
      if (r.T_C % 100 !== 0) return;
      var o = document.createElement('option');
      o.value = r.T_C; o.textContent = r.T_C + ' °C';
      if (r.T_C === 600) o.selected = true;
      eqT.appendChild(o);
    });
    var sfV = $('sfV');
    D.surface.inventories.forEach(function (v) {
      var o = document.createElement('option');
      o.value = v; o.textContent = v + ' µmol-O g⁻¹';
      if (v === 95) o.selected = true;
      sfV.appendChild(o);
    });

    eqT.addEventListener('change', drawEquilibrium);
    $('sfD').addEventListener('change', drawSurface);
    sfV.addEventListener('change', drawSurface);
    Array.prototype.forEach.call(document.querySelectorAll('.wstab'), function (b) {
      b.addEventListener('click', function () {
        showWorkspace(b.getAttribute('data-ws'));
      });
    });

    drawEquilibrium();
    drawSurface();
    var want = location.hash.slice(1);
    showWorkspace(want === 'ws-surface' ? 'ws-surface' : 'ws-equilibrium');
  }

  if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', init);
  } else {
    init();
  }
})();
