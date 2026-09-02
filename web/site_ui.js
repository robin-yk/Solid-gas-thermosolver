/* The home screen, the switch between the two workspaces, and the
   population workspace, which integrates one rate equation mirrored in
   population.js and gated against the Python module. The equilibrium
   workspace is driven by ti_solver_page.js and thermo_ui.js. */

(function () {
  'use strict';
  var D = window.SITE_DATA;
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
  function kv(rows) {
    return '<tbody>' + rows.map(function (r) {
      return '<tr><th>' + r[0] + '</th><td>' + r[1] + '</td></tr>';
    }).join('') + '</tbody>';
  }
  function row(label, value) { return [label, value]; }

  /* Three workspaces behind three hashes; no hash is the home screen. */
  var HASH = { equilibrium: 'ws-thermo', population: 'ws-population' };
  var WS = ['ws-thermo', 'ws-population'];

  function hashOf(id) {
    for (var k in HASH) if (HASH[k] === id) return k;
    return 'home';
  }

  function showWorkspace(id) {
    var home = WS.indexOf(id) < 0;
    WS.forEach(function (w) {
      var el = $(w);
      if (el) el.style.display = (w === id) ? '' : 'none';
    });
    $('home').style.display = home ? '' : 'none';
    $('topbar').style.display = home ? 'none' : '';
    Array.prototype.forEach.call(document.querySelectorAll('.wstab'), function (b) {
      b.classList.toggle('active', b.getAttribute('data-ws') === id);
    });
    var want = home ? '' : '#' + hashOf(id);
    if (location.hash !== want) history.replaceState(null, '', location.pathname + want);
    window.scrollTo(0, 0);
  }

  function route() {
    var h = location.hash.replace(/^#/, '');
    showWorkspace(HASH[h] || 'home');
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

    var hiRow = rows.filter(function (q) { return q.n_iso_umol_g === isoMax; })[0];
    var loRow = rows.filter(function (q) { return q.n_iso_umol_g === isoMin; })[0];
    var rise = series[series.length - 1].nV_total_umol_g / series[0].nV_total_umol_g;
    $('pnBasis').innerHTML = kv([
      row('Active-region capacity n<sub>s</sub>', n(p.ns, 3) + ' µmol-O g⁻¹'),
      row('E<sub>loss</sub>', n(p.E, 4) + ' eV'),
      row('ν<sub>eff</sub>', nu.toExponential(2) + ' s⁻¹'),
      row('Residual', popKind() + ' ratio'),
      row('Samples', series.length + ' (' + series[0].sample + ' to '
          + series[series.length - 1].sample + ')')
    ]);
    $('pnVerdict').innerHTML = '<div class="verdict">The isolated population peaks at <b>'
      + hiRow.sample + '</b> and falls to ' + n(isoMin / isoMax * 100, 1) + '% of that at '
      + loRow.sample + ', while the total inventory rises ' + n(rise, 0) + '×.</div>';
    $('pnKpis').innerHTML = kv([
      row('Isolated, highest', n(isoMax, 2) + ' µmol-O g⁻¹ (' + hiRow.sample + ')'),
      row('Isolated, lowest', n(isoMin, 3) + ' µmol-O g⁻¹ (' + loRow.sample + ')'),
      row('R800 / R600, fitted', n(ratio, 3)),
      row('R800 / R600, measured', n(meas, 3))
    ]);

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
      '<tr><td>Objective (' + popKind() + ' residual)</td><td class="k2">'
      + ssq.toExponential(3) + '</td></tr>'
      + '<tr><td>Mass balance</td><td class="k2">'
      + (closure === 0 ? 'closes exactly, by construction'
         : 'closes to ' + n(closure, 9) + ' µmol-O g⁻¹') + '</td></tr>'
      + '<tr><td>Parameters fitted</td><td class="k2">3 (n<sub>s</sub>, E<sub>loss</sub>, ν<sub>eff</sub>)</td></tr>'
      + '<tr><td>Rates fitted</td><td class="k2">' + series.length + '</td></tr>'
      + '<tr><td>Inventory mostly below the active region</td>'
      + '<td class="k2">' + (below ? 'every sample' : 'not every sample') + '</td></tr>';
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
    Array.prototype.forEach.call(document.querySelectorAll('.wstab, .card'), function (b) {
      b.addEventListener('click', function () {
        showWorkspace(b.getAttribute('data-ws'));
      });
    });
    window.addEventListener('hashchange', route);
    drawPopulation();
    route();
  }

  if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', init);
  } else {
    init();
  }
})();
