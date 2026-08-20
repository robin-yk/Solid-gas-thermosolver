/* UI for the Ti-O active-set solver. Every number comes from
   ActiveSet.Solver.solve; the stability verdict is read off the KKT reduced
   costs, never off the size of a mole number. */
(function () {
  'use strict';
  var $ = function (id) { return document.getElementById(id); };
  var D = JSON.parse($('activeset-data').textContent);
  var REF = JSON.parse($('reference-data').textContent);
  var solver = new ActiveSet.Solver(D);
  var LAST = null;

  var sub = function (s) {
    return String(s).replace(/([A-Za-z])(\d+)/g, '$1<sub>$2</sub>');
  };
  function sci(v, digits) {
    if (v === 0) return '0';
    if (v === null || v === undefined || !isFinite(v)) return String(v);
    var a = Math.abs(v);
    if (a >= 0.001 && a < 100000) {
      return Number(v.toPrecision(digits || 4)).toString();
    }
    var p = v.toExponential((digits || 4) - 1).split('e');
    return p[0] + '×10' + sup(+p[1]);
  }
  function sup(n) {
    return String(n).replace(/-/g, '⁻').replace(/[0-9]/g, function (d) {
      return '⁰¹²³⁴⁵⁶⁷⁸⁹'[+d];
    });
  }

  /* ------------------------------------------------------------- inputs */

  /* A preset always sets the gas composition; the condition presets also set
     temperature, mass, volume, and the initial solid. */
  var PRESETS = {
    rwgs: { gas: { CO2: 1, H2: 1, CO: 0, H2O: 0, N2: 0, O2: 0 } },
    rwgsdil: { gas: { CO2: 10, H2: 10, CO: 0, H2O: 0, N2: 80, O2: 0 } },
    h2rich: { gas: { CO2: 1, H2: 3, CO: 0, H2O: 0, N2: 0, O2: 0 } },
    product: { gas: { CO2: 3, H2: 3, CO: 2, H2O: 2, N2: 0, O2: 0 } },
    h2: { gas: { CO2: 0, H2: 1, CO: 0, H2O: 0, N2: 0, O2: 0 } },
    h2co2: { gas: { CO2: 1, H2: 99, CO: 0, H2O: 0, N2: 0, O2: 0 } },
    n2: { gas: { CO2: 0, H2: 0, CO: 0, H2O: 0, N2: 1, O2: 0 } },
    n2o2: { gas: { CO2: 0, H2: 0, CO: 0, H2O: 0, N2: 999990, O2: 10 } },
    co2: { gas: { CO2: 1, H2: 0, CO: 0, H2O: 0, N2: 0, O2: 0 } },
    magneli: { gas: { CO2: 0, H2: 1, CO: 0, H2O: 0, N2: 0, O2: 0 },
               T: 1200, mass: 2, V: 25, solid: 'TiO2' },
    reox: { gas: { CO2: 1, H2: 0, CO: 0, H2O: 0, N2: 0, O2: 0 },
            T: 600, mass: 100, V: 25, solid: 'Ti4O7' },
  };

  D.gases.forEach(function (g) {
    var div = document.createElement('div');
    div.className = 'g';
    div.innerHTML = '<span>' + sub(g) + '</span>'
      + '<input type="number" step="any" min="0" id="gas_' + g + '" value="0">';
    $('gasInputs').appendChild(div);
  });
  D.solids.forEach(function (s) {
    var o = document.createElement('option');
    o.value = s;
    o.innerHTML = sub(s) + (s === 'TiO2' ? ' (rutile)' : '');
    $('solid0').appendChild(o);
  });

  function feedFromInputs() {
    var feed = {}, tot = 0;
    D.gases.forEach(function (g) {
      var v = parseFloat($('gas_' + g).value);
      if (!(v >= 0)) v = 0;
      feed[g] = v;
      tot += v;
    });
    return { feed: feed, total: tot };
  }
  function setFeed(feed) {
    D.gases.forEach(function (g) {
      $('gas_' + g).value = feed[g] !== undefined ? feed[g] : 0;
    });
    showNormalized();
  }
  function showNormalized() {
    var f = feedFromInputs();
    if (!(f.total > 0)) {
      $('normRead').innerHTML = '<span class="bad">gas feed sums to zero - '
        + 'enter at least one species</span>';
      return;
    }
    var parts = [];
    D.gases.forEach(function (g) {
      var y = f.feed[g] / f.total;
      if (y > 0) {
        parts.push(sub(g) + ' ' + (y >= 1e-4
          ? (100 * y).toFixed(3) + '%' : sci(y, 3)));
      }
    });
    $('normRead').innerHTML = 'used by the solver (normalised): ' + parts.join(' &middot; ');
  }
  D.gases.forEach(function (g) {
    $('gas_' + g).addEventListener('input', function () {
      $('preset').value = '';
      showNormalized();
    });
  });
  $('normBtn').addEventListener('click', function () {
    var f = feedFromInputs();
    if (!(f.total > 0)) return;
    D.gases.forEach(function (g) {
      var y = f.feed[g] / f.total;
      $('gas_' + g).value = y === 0 ? 0 : Number((100 * y).toPrecision(10));
    });
    showNormalized();
  });
  $('o2ppmSet').addEventListener('click', function () {
    var ppm = parseFloat($('o2ppm').value);
    if (!(ppm >= 0)) return;
    var f = feedFromInputs();
    var others = f.total - f.feed.O2;
    if (!(others > 0)) return;
    var frac = ppm * 1e-6;
    $('gas_O2').value = frac >= 1 ? 0 : others * frac / (1 - frac);
    $('preset').value = '';
    showNormalized();
  });
  $('preset').addEventListener('change', function () {
    var p = PRESETS[$('preset').value];
    if (!p) return;
    setFeed(p.gas);
    if (p.T !== undefined) $('T').value = p.T;
    if (p.mass !== undefined) $('mass').value = p.mass;
    if (p.V !== undefined) $('V').value = p.V;
    if (p.solid !== undefined) $('solid0').value = p.solid;
  });
  $('chargeAtT').addEventListener('change', function () {
    $('chargeTWrap').style.display = $('chargeAtT').checked ? 'none' : 'grid';
  });

  /* --------------------------------------------------------------- run */

  function currentOpts() {
    var f = feedFromInputs();
    if (!(f.total > 0)) return null;
    var opts = {
      feed: f.feed,
      T_K: (parseFloat($('T').value) || 900) + 273.15,
      P_atm: parseFloat($('P').value) || 1,
      V_cm3: parseFloat($('V').value) || 25,
      mass_g: (parseFloat($('mass').value) || 100) / 1000,
      initial_solid: $('solid0').value,
    };
    if (!$('chargeAtT').checked) {
      opts.T_charge_K = (parseFloat($('Tcharge').value) || 25) + 273.15;
    }
    return opts;
  }

  function run() {
    var opts = currentOpts();
    if (!opts) {
      $('status').innerHTML = '<span class="bad">gas feed sums to zero</span>';
      return;
    }
    var t0 = performance.now();
    var R;
    try {
      R = solver.solve(opts);
    } catch (e) {
      $('status').innerHTML = '<span class="bad">' + e.message + '</span>';
      return;
    }
    LAST = R;
    $('status').innerHTML = 'solved in ' + (performance.now() - t0).toFixed(1)
      + ' ms &middot; ' + R.iterations + ' Newton steps &middot; residual '
      + sci(R.newton_residual, 2) + ' &middot; method <code>gibbs_min</code>';
    render(R);
  }

  /* ------------------------------------------------------------ render */

  function nearestBoundary(R) {
    var best = null, bestV = Infinity, s;
    for (s in R.inactive_phase_reduced_costs) {
      var v = R.inactive_phase_reduced_costs[s].per_mol_O_kJ;
      if (v === null || !isFinite(v) || v <= 0) continue;
      if (v < bestV) { bestV = v; best = s; }
    }
    return best ? { phase: best, r: bestV } : null;
  }

  function render(R) {
    var reduced = R.active_condensed_phases.filter(function (s) {
      return D.ti3[s] > 0;
    });
    var v;
    if (reduced.length) {
      v = '<div class="verdict red"><b>Reduction is thermodynamically favoured.</b> '
        + 'Stable assemblage: <b>' + R.active_condensed_phases.map(sub).join(' + ')
        + '</b>; Ti³⁺ = ' + sci(R.reduced_pct, 4) + '% of the titanium'
        + (R.reduced_pct < 1e-6
           ? ' (trace amount: log₁₀ n = '
             + (R.log10_trace_amounts[reduced[0]] || 0).toFixed(1)
             + '; two-phase buffer under an oxygen-free gas)' : '')
        + '.</div>';
    } else {
      var nb = nearestBoundary(R);
      v = '<div class="verdict none"><b>No reduced phase is stable at this '
        + 'condition.</b> All reduced phases are held at zero moles with '
        + 'positive KKT reduced costs' + (nb
          ? '; the nearest boundary is <b>' + sub(nb.phase) + '</b> at r = +'
            + sci(nb.r, 4) + ' kJ per mol O' : '')
        + '.</div>';
    }
    if (R.boundary_or_degenerate) {
      v += '<div class="note">Degenerate boundary: '
        + (R.degenerate_phases.map(sub).join(', ') || 'triple-point probe')
        + '. A reduced cost is zero within '
        + sci(R.solver_tolerance.degeneracy_kJ_per_mol_O, 1)
        + ' kJ/mol O; two assemblages coexist.</div>';
    }
    $('verdict').innerHTML = v;

    var kpis = [
      ['Assemblage', R.active_condensed_phases.map(sub).join(' + '), 'active condensed phases', reduced.length ? 'warn' : 'good'],
      ['Ti³⁺', sci(R.reduced_pct, 4) + '%', 'of total titanium', ''],
      ['x in TiO₂₋ₓ', sci(R.x_in_TiO2x, 4), 'solid O deficit', ''],
      ['O moved to gas', sci(R.oxygen_to_gas_mol_O, 4), 'mol O (negative = into solid)', ''],
    ];
    if (R.conversion_CO2_pct !== null) {
      kpis.push(['CO₂ conversion', R.conversion_CO2_pct.toFixed(2) + '%', 'net, vs feed', '']);
    }
    kpis.push(['Gas charge', sci(R.n_gas_charge_mol, 4), 'mol at Tₑ = '
               + (R.T_charge_K - 273.15).toFixed(0) + ' °C', '']);
    $('kpis').innerHTML = kpis.map(function (k) {
      return '<div class="kpi ' + k[3] + '"><div class="k">' + k[0]
        + '</div><div class="v">' + k[1] + '</div><div class="u">' + k[2]
        + '</div></div>';
    }).join('');

    /* phase table */
    var h = '<thead><tr><th>phase</th><th>mol</th><th>% of Ti</th><th>status</th></tr></thead><tbody>';
    D.solids.forEach(function (s) {
      var n = R.solid_moles[s];
      var status;
      if (R.active_condensed_phases.indexOf(s) >= 0) {
        status = n < 1e-12
          ? 'active, trace: log₁₀ n = ' + R.log10_trace_amounts[s].toFixed(2)
          : 'active';
      } else {
        var rc = R.inactive_phase_reduced_costs[s];
        status = 'excluded, exactly 0; r = '
          + (rc.per_mol_O_kJ === null ? '—'
             : !isFinite(rc.per_mol_O_kJ) ? '−∞'
             : '+' + sci(rc.per_mol_O_kJ, 4)) + ' kJ/mol O';
      }
      h += '<tr><td>' + sub(s) + '</td><td>' + (n === 0 ? '0 (exact)' : sci(n, 5))
        + '</td><td>' + (R.phase_split_pct[s] === 0 ? '0'
                         : sci(R.phase_split_pct[s], 5)) + '</td><td>'
        + status + '</td></tr>';
    });
    $('tPhase').innerHTML = h + '</tbody>';

    /* gas table */
    h = '<thead><tr><th>species</th><th>y feed</th><th>n feed (mol)</th>'
      + '<th>y equilibrium</th><th>n equilibrium (mol)</th><th>Δn (mol)</th></tr></thead><tbody>';
    D.gases.forEach(function (g) {
      var y0 = R.feed_mole_fractions[g] || 0;
      var n0 = y0 * R.n_gas_charge_mol;
      h += '<tr><td>' + sub(g) + '</td><td>' + (y0 ? sci(y0, 5) : '0')
        + '</td><td>' + (n0 ? sci(n0, 5) : '0')
        + '</td><td>' + sci(R.gas_fractions[g], 5)
        + '</td><td>' + sci(R.species_moles[g], 5)
        + '</td><td>' + sci(R.delta_n_gas[g], 4) + '</td></tr>';
    });
    $('tGas').innerHTML = h + '</tbody>';

    $('resNotes').innerHTML =
      '<div class="note calm">G<sub>total</sub> = ' + sci(R.G_total_kJ, 8)
      + ' kJ; common reference (unreacted charge) = ' + sci(R.common_reference_G_kJ, 8)
      + ' kJ; G<sub>rel</sub> = '
      + (R.G_rel_below_double_resolution
         ? 'below double-precision resolution (|G_rel| &lt; ' + sci(R.reporting_detection_limit.G_rel_kJ_floor, 2)
           + ' kJ); the exact value is in the reference dataset'
         : sci(R.G_rel_kJ, 6) + ' kJ')
      + '. Mode: ' + R.mode + '.</div>';

    /* KKT tab */
    h = '<thead><tr><th>inactive phase</th><th>r per formula (kJ)</th>'
      + '<th>r per mol O (kJ)</th><th>O removed per extent</th><th>reading</th></tr></thead><tbody>';
    D.solids.forEach(function (s) {
      var rc = R.inactive_phase_reduced_costs[s];
      if (!rc) return;
      var pf = rc.per_formula_kJ, po = rc.per_mol_O_kJ;
      var reading = pf === null ? '—'
        : !isFinite(pf) ? 'oxygen-free gas; formation direction unbounded'
        : Math.abs(po) < R.solver_tolerance.degeneracy_kJ_per_mol_O
          ? 'on the boundary (degenerate)'
        : pf > 0 ? 'absent (exactly zero)' : 'negative reduced cost (not expected)';
      h += '<tr><td>' + sub(s) + '</td><td>'
        + (pf === null ? '—' : !isFinite(pf) ? '−∞' : sci(pf, 5))
        + '</td><td>'
        + (po === null ? '—' : !isFinite(po) ? '−∞' : sci(po, 5))
        + '</td><td>' + rc.oxygen_removed_per_extent + '</td><td>' + reading
        + '</td></tr>';
    });
    $('tKKT').innerHTML = h + '</tbody>';

    h = '<thead><tr><th>element</th><th>λ (kJ/mol)</th></tr></thead><tbody>';
    Object.keys(R.lambda_kJ_per_mol).sort().forEach(function (e) {
      h += '<tr><td>' + e + '</td><td>' + sci(R.lambda_kJ_per_mol[e], 6) + '</td></tr>';
    });
    $('tLam').innerHTML = h + '</tbody>';

    h = '<thead><tr><th>candidate</th><th>feasible</th><th>KKT</th>'
      + '<th>worst inactive r (kJ)</th><th>note</th></tr></thead><tbody>';
    R.candidates.forEach(function (c) {
      var isWin = c.active.join() === R.active_condensed_phases.join();
      h += '<tr' + (isWin ? ' style="font-weight:700"' : '') + '><td>'
        + c.active.map(sub).join(' + ') + (isWin ? ' ← selected' : '')
        + '</td><td>' + (c.feasible ? '<span class="ok">yes</span>' : 'no')
        + '</td><td>' + (c.kkt_ok ? '<span class="ok">clean</span>' : 'fails')
        + '</td><td>' + (c.worst_inactive_r_kJ === null ? '−∞ / —'
                         : sci(c.worst_inactive_r_kJ, 4))
        + '</td><td>' + c.reason + '</td></tr>';
    });
    $('tCand').innerHTML = h + '</tbody>';

    /* balance tab */
    h = '<thead><tr><th>element</th><th>charge (mol)</th><th>solution (mol)</th>'
      + '<th>residual (mol)</th></tr></thead><tbody>';
    ['C', 'H', 'N', 'O', 'Ti'].forEach(function (e) {
      var target = 0;
      D.gases.forEach(function (g) {
        target += (D.gas_formula[g][e] || 0) * (R.feed_mole_fractions[g] || 0)
                  * R.n_gas_charge_mol;
      });
      var st = D.stoich[R.initial_solid];
      var nS0 = R.mass_g / (st[0] * D.atomic_weight.Ti + st[1] * D.atomic_weight.O);
      if (e === 'Ti') target += st[0] * nS0;
      if (e === 'O') target += st[1] * nS0;
      h += '<tr><td>' + e + '</td><td>' + sci(target, 6) + '</td><td>'
        + sci(target + R.balance_by_element[e], 6) + '</td><td>'
        + sci(R.balance_by_element[e], 3) + '</td></tr>';
    });
    $('tBal').innerHTML = h + '</tbody>';
    $('balNote').innerHTML = '<div class="note calm">Solver tolerance (scaled): '
      + sci(R.solver_tolerance.newton_scaled_residual, 2)
      + '; KKT tolerance: ' + sci(R.solver_tolerance.kkt_kJ, 2) + ' kJ.</div>';
  }

  /* -------------------------------------------------------------- CSV */

  $('csv').addEventListener('click', function () {
    if (!LAST) return;
    var R = LAST, lines = [], k;
    lines.push('field,value');
    ['system', 'method', 'solver', 'mode', 'T_C', 'P_atm', 'V_cm3',
     'T_charge_K', 'mass_g', 'initial_solid', 'n_gas_charge_mol', 'n_Ti_mol',
     'reduced_pct', 'x_in_TiO2x', 'oxygen_to_gas_mol_O', 'conversion_CO2_pct',
     'G_total_kJ', 'common_reference_G_kJ', 'G_rel_kJ',
     'G_rel_below_double_resolution', 'boundary_or_degenerate', 'converged']
      .forEach(function (f) {
        lines.push(f + ',"' + String(R[f]).replace(/"/g, '""') + '"');
      });
    lines.push('active_condensed_phases,"' + R.active_condensed_phases.join('+') + '"');
    for (k in R.solid_moles) lines.push('solid_mol_' + k + ',' + R.solid_moles[k]);
    for (k in R.species_moles) {
      lines.push('gas_mol_' + k + ',' + R.species_moles[k]);
      lines.push('gas_frac_' + k + ',' + R.gas_fractions[k]);
      lines.push('delta_n_' + k + ',' + R.delta_n_gas[k]);
    }
    for (k in R.inactive_phase_reduced_costs) {
      lines.push('r_per_mol_O_' + k + ','
        + R.inactive_phase_reduced_costs[k].per_mol_O_kJ);
    }
    for (k in R.balance_by_element) {
      lines.push('balance_' + k + ',' + R.balance_by_element[k]);
    }
    var blob = new Blob([lines.join('\n')], { type: 'text/csv' });
    var a = document.createElement('a');
    a.href = URL.createObjectURL(blob);
    a.download = 'ti_o_equilibrium.csv';
    a.click();
    URL.revokeObjectURL(a.href);
  });

  /* -------------------------------------------------------- validation */

  $('valRun').addEventListener('click', function () {
    var cases = REF.cases_feeds;
    var out = '<table><thead><tr><th>case</th><th>T (°C)</th>'
      + '<th>assemblage</th><th>oracle agrees</th><th>worst |Δy/y|</th>'
      + '<th>worst |Δr| kJ</th><th>trace log₁₀ Δ</th></tr></thead><tbody>';
    var allOk = true;
    REF.rows.forEach(function (row) {
      var feed = cases[row.case];
      if (!feed) return;
      var R = solver.solve({ feed: feed, T_K: row.T_C + 273.15 });
      var okSet = R.active_condensed_phases.join() === row.active.join();
      var wy = 0, wr = 0, wt = 0, g, s;
      for (g in row.gas_fractions) {
        var ov = Number(row.gas_fractions[g]);
        if (ov > 1e-25) {
          wy = Math.max(wy, Math.abs(R.gas_fractions[g] - ov) / ov);
        }
      }
      for (s in row.r_per_mol_O) {
        var rv = row.r_per_mol_O[s];
        var mine = R.inactive_phase_reduced_costs[s];
        if (rv !== null && mine && mine.per_mol_O_kJ !== null
            && isFinite(mine.per_mol_O_kJ)) {
          wr = Math.max(wr, Math.abs(mine.per_mol_O_kJ - Number(rv)));
        }
      }
      for (s in row.log10_traces) {
        if (R.log10_trace_amounts[s] !== undefined) {
          wt = Math.max(wt, Math.abs(R.log10_trace_amounts[s] - row.log10_traces[s]));
        }
      }
      var ok = okSet && wy < 1e-9 && wr < 1e-8 && wt < 1e-6;
      allOk = allOk && ok;
      out += '<tr><td>' + row.case + '</td><td>' + row.T_C + '</td><td>'
        + R.active_condensed_phases.map(sub).join('+') + '</td><td>'
        + (ok ? '<span class="ok">yes</span>' : '<span class="bad">NO</span>')
        + '</td><td>' + sci(wy, 2) + '</td><td>' + sci(wr, 2) + '</td><td>'
        + sci(wt, 2) + '</td></tr>';
    });
    out += '</tbody></table>';
    out = '<div class="note ' + (allOk ? 'calm' : '') + '"><b>'
      + (allOk ? 'Every reference point reproduced within tolerance.'
               : 'Disagreement with the reference; this build is not validated.')
      + '</b> Reference: ' + REF.precision + '.</div>' + out;
    $('valBody').innerHTML = out;
  });

  /* ------------------------------------------------------ coefficients */

  function renderCoefs() {
    var h = '<h3>Gases - NIST-JANAF Shomate</h3>'
      + '<div class="tablewrap"><table><thead><tr><th>species</th><th>set</th>'
      + '<th>A</th><th>B</th><th>C</th><th>D</th><th>E</th><th>F</th><th>G</th>'
      + '<th>H</th><th>ΔH°f (kJ/mol)</th><th>break (K)</th></tr></thead><tbody>';
    D.gases.forEach(function (g) {
      var e = D.shomate[g];
      ['lo', 'hi'].forEach(function (band) {
        h += '<tr><td>' + (band === 'lo' ? sub(g) : '') + '</td><td>' + band + '</td>'
          + e[band].map(function (c) { return '<td>' + c + '</td>'; }).join('')
          + '<td>' + e.dHf + '</td><td>' + (band === 'lo' ? e.brk : '') + '</td></tr>';
      });
    });
    h += '</tbody></table></div>'
      + '<h3>Condensed phases - Waldner &amp; Eriksson CALPHAD</h3>'
      + '<div class="tablewrap"><table><thead><tr><th>phase</th>'
      + '<th>ΔH°f,298 (J/mol)</th><th>S°298 (J/mol K)</th>'
      + '<th>range (K)</th><th>form</th><th>Cp coefficients</th>'
      + '<th>ΔH trans (J)</th></tr></thead><tbody>';
    D.solids.forEach(function (s) {
      var e = D.waldner[s];
      e.ranges.forEach(function (r, i) {
        h += '<tr><td>' + (i ? '' : sub(s)) + '</td><td>' + (i ? '' : e.dHf298)
          + '</td><td>' + (i ? '' : e.S298) + '</td><td>' + r.T0 + '–' + r.T1
          + '</td><td>' + r.form + '</td><td>' + r.c.join(', ') + '</td><td>'
          + r.dHt + '</td></tr>';
      });
    });
    h += '</tbody></table></div>'
      + '<div class="note calm">' + D.provenance.gases + ' &middot; '
      + D.provenance.solids + '</div>';
    $('coefBody').innerHTML = h;
  }

  /* ------------------------------------------------------------- shell */

  document.querySelectorAll('.tab').forEach(function (b) {
    b.addEventListener('click', function () {
      document.querySelectorAll('.tab').forEach(function (x) {
        x.classList.toggle('active', x === b);
      });
      document.querySelectorAll('.page').forEach(function (p) {
        p.classList.toggle('active', p.id === b.dataset.page);
      });
    });
  });

  $('themeBtn').addEventListener('click', function () {
    var r = document.documentElement;
    var dark = r.getAttribute('data-theme') === 'dark';
    if (dark) r.removeAttribute('data-theme');
    else r.setAttribute('data-theme', 'dark');
    try { localStorage.setItem('tio_theme', dark ? 'light' : 'dark'); } catch (e) {}
  });
  try {
    if (localStorage.getItem('tio_theme') === 'dark') {
      document.documentElement.setAttribute('data-theme', 'dark');
    }
  } catch (e) {}

  $('rawSrc').textContent = [
    'mu0Gas', 'mu0Solid', 'buildCharge', 'solveCandidate',
    'enumerateCandidates', 'tripleProbes', 'solve',
  ].map(function (f) {
    return '// Solver.prototype.' + f + '\n'
      + ActiveSet.Solver.prototype[f].toString();
  }).join('\n\n');

  $('run').addEventListener('click', run);
  renderCoefs();
  setFeed(PRESETS.rwgs.gas);
  $('preset').value = 'rwgs';
  run();
})();
