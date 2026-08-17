/* UI for the gas-solid Gibbs minimisation. All numbers come from Oxide.minimise;
   the verification tab is the only place the second route is allowed to appear,
   and it is labelled there. */
(function () {
  'use strict';
  const $ = id => document.getElementById(id);
  const D = Oxide.D;
  const sub = s => String(s).replace(/([A-Za-z])(\d+)/g, '$1<sub>$2</sub>');
  const sci = v => { const [m, e] = v.toExponential(2).split('e'); return m + '×10' + sup(+e); };
  const sup = n => String(n).replace(/-/g, '\u207b').replace(/[0-9]/g,
    d => '⁰¹²³⁴⁵⁶⁷⁸⁹'[+d]);

  /* Below this share of the metal, a phase is the optimiser's floor rather than
     a phase. The present build cannot represent an absent phase as exactly zero
     - see the defect note on the How tab - so the limit is stated and shown
     rather than quietly rounded away. */
  const DETECT_PCT = 1e-6;
  const LIMIT = '&lt; 10\u207b\u2076 %';
  const GAS = D.gases;
  const DEPTHS = {
    full: ['TiO2', 'Ti10O19', 'Ti8O15', 'Ti6O11', 'Ti5O9', 'Ti4O7', 'Ti3O5', 'Ti2O3'],
    ti3o5: ['TiO2', 'Ti4O7', 'Ti3O5'],
    ti4o7: ['TiO2', 'Ti4O7'],
  };
  let LAST = null, SWEEP = null;

  function cfg() {
    const name = $('sys').value, sy = D.systems[name];
    const solids = name === 'Ti-O' ? DEPTHS[$('depth').value] : sy.solids;
    return {
      name, host: sy.host, metal: sy.metal, mw: sy.mw, source: sy.source, solids,
      P: +$('P').value || 1, V: +$('V').value || 25,
      mass: (+$('mass').value || 100) / 1000, T_C: +$('T').value,
    };
  }

  function solveAt(c, T_C) {
    const T = T_C + 273.15, ng = c.P * c.V / (D.R_ATM * T), n0 = ng / 2;
    const nM = c.mass / c.mw;
    const S = new Oxide.GibbsSystem(GAS, c.solids, ['CO2', 'H2', 'H2O', c.host]);
    const b = S.bFrom({ CO2: n0, H2: n0, [c.host]: nM });
    const t0 = performance.now();
    const r = Oxide.minimise(S, b, T, Oxide.buildSeeds(S, ng, nM), c.P);
    if (!r) return null;
    const ms = performance.now() - t0;
    const y = S.gasFractions(r.n);
    const split = S.phaseSplit(r.n, nM);
    const present = Object.keys(split).filter(p => split[p] > DETECT_PCT);
    let O = 0, M = 0;
    Object.keys(split).forEach(p => {
      const moles = r.n[S.species.indexOf(p)];
      O += moles * (D.formulas[p].O || 0); M += moles * (D.formulas[p][c.metal] || 0);
    });
    return {
      T_C, T, S, b, n: r.n, y, split, present, ms, solver: r, nM, n0, ng,
      reduced: S.reducedPercent(r.n, nM),
      x: M > 0 ? 2 - O / M : 0,
      conv: (n0 - r.n[S.species.indexOf('CO2')]) / n0 * 100,
      Kp: Math.exp(-(Oxide.mu0('CO', T) + Oxide.mu0('H2O', T)
                     - Oxide.mu0('CO2', T) - Oxide.mu0('H2', T)) / (Oxide.R * T)),
      Q: (y.CO * y.H2O) / (y.CO2 * y.H2),
      muO: Oxide.Verify.muOFromGas(y, T),
    };
  }

  function table(id, cols, rows) {
    const t = $(id); t.textContent = '';
    const th = document.createElement('thead'), tr = document.createElement('tr');
    cols.forEach(c => { const e = document.createElement('th'); e.innerHTML = c; tr.appendChild(e); });
    th.appendChild(tr); t.appendChild(th);
    const tb = document.createElement('tbody');
    rows.forEach(r => {
      const e = document.createElement('tr');
      r.forEach(c => { const d = document.createElement('td'); d.innerHTML = c; e.appendChild(d); });
      tb.appendChild(e);
    });
    t.appendChild(tb);
  }

  function render() {
    const c = cfg(), R = solveAt(c, c.T_C);
    if (!R) { $('status').innerHTML = '<span class="bad">no feasible solution</span>'; return; }
    LAST = { c, R };

    const reduced = R.reduced > DETECT_PCT;
    $('kpis').innerHTML = [
      ['Solid', R.present.map(sub).join(' + ') || sub(c.host), 'phases present', ''],
      [c.metal + '(3+)', reduced ? R.reduced.toFixed(R.reduced < 1 ? 4 : 2) + '%' : LIMIT,
       reduced ? 'of total ' + c.metal : 'below detection limit', reduced ? 'warn' : 'good'],
      ['x in ' + c.metal + 'O<sub>2−x</sub>', reduced ? R.x.toExponential(2) : LIMIT,
       'oxygen deficiency', reduced ? 'warn' : 'good'],
      ['μ<sub>O</sub>', R.muO.muO == null ? '—' : R.muO.muO.toFixed(1), 'kJ per mol O', ''],
      ['pO<sub>2</sub>', R.muO.muO == null ? '—' : sci(Oxide.Verify.muOToPO2(R.muO.muO, R.T)), 'atm', ''],
      ['CO<sub>2</sub> conversion', R.conv.toFixed(2) + '%', 'of the feed', ''],
    ].map(k => '<div class="kpi ' + k[3] + '"><div class="k">' + k[0] + '</div><div class="v">'
      + k[1] + '</div><div class="u">' + k[2] + '</div></div>').join('');

    table('tPhase', ['Phase', 'O / ' + c.metal, 'share of ' + c.metal + ' (%)', 'moles'],
      c.solids.map(p => {
        const f = D.formulas[p], moles = R.n[R.S.species.indexOf(p)];
        const on = R.split[p] > DETECT_PCT;
        return [sub(p), (f.O / f[c.metal]).toFixed(4),
          on ? R.split[p].toFixed(6) : '<span style="color:var(--muted)">' + LIMIT + '</span>',
          on ? moles.toExponential(3) : '0'];
      }));

    table('tGas', ['Species', 'mol %', 'moles', 'p<sub>i</sub> (atm)'],
      GAS.map(s => [sub(s), R.y[s] > 1e-12 ? (R.y[s] * 100).toFixed(6) : sci(R.y[s] * 100),
        R.n[R.S.species.indexOf(s)].toExponential(3), sci(R.y[s] * c.P)]));

    $('status').innerHTML = 'converged &middot; ' + R.solver.seeds + ' starts, '
      + R.solver.iterations + ' iterations, ' + R.ms.toFixed(0) + ' ms<br>'
      + 'balance residual ' + R.solver.residual.toExponential(1) + ' mol &middot; '
      + 'method <code>' + R.solver.method + '</code>';

    renderNotes(c, R); renderVerify(c, R); renderBalance(c, R); renderCoef(c);
    $('prov').textContent = 'Solved in this page: ' + R.solver.seeds + ' starts, '
      + R.solver.iterations + ' iterations, ' + R.ms.toFixed(0) + ' ms. ';
  }

  /* What this run cannot tell you, next to what it does. */
  function renderNotes(c, R) {
    const n = [];
    if (c.name === 'Ti-O') {
      n.push('<div class="note"><strong>The series is cut at ' + sub(c.solids[c.solids.length - 1])
        + '.</strong> Where it is cut changes the answer, because each shallower Magnéli member '
        + 'forms more easily than the last. Use the <em>Series cut</em> control and watch μ<sub>O</sub> '
        + 'against the boundary on the Verification tab: cutting at Ti<sub>4</sub>O<sub>7</sub> alone '
        + 'measures against too hard a step and reads about 14 kJ/mol-O too safe at 1500 °C.</div>');
    }
    n.push('<div class="note stop"><strong>Not modelled.</strong> No graphite, so carbon cannot '
      + 'deposit and must leave as CO or CH<sub>4</sub>. No oxygen-deficient rutile '
      + '(' + sub(c.host) + '<sub>2−x</sub> as a solution phase) &mdash; every condensed phase here is a '
      + 'line compound, present or absent at unit activity, which is exactly what the n → ∞ '
      + 'limit of the series is not. Equilibrium only: nothing here says how fast any of it happens.</div>');
    n.push('<div class="note stop"><strong>No experimental validation.</strong> Every check in this '
      + 'page is internal consistency — two independent routes agreeing, coexisting phases landing '
      + 'on their computed boundary, Q matching K<sub>p</sub>, atoms balancing. Nothing has been '
      + 'compared against a measured phase boundary.</div>');
    $('resNotes').innerHTML = n.join('');
  }

  function renderVerify(c, R) {
    const T = R.T, V = Oxide.Verify;
    let h = '<h3>1. Reaction quotient against the equilibrium constant</h3>'
      + '<p class="sub">The minimiser is not told about the RWGS reaction. If it found the '
      + 'equilibrium, Q must come out equal to K<sub>p</sub> anyway.</p>'
      + '<div class="tablewrap"><table><thead><tr><th>quantity</th><th>value</th></tr></thead><tbody>'
      + '<tr><td>K<sub>p</sub> from standard potentials</td><td>' + R.Kp.toFixed(6) + '</td></tr>'
      + '<tr><td>Q from the minimised composition</td><td>' + R.Q.toFixed(6) + '</td></tr>'
      + '<tr><td>relative difference</td><td class="' + (Math.abs(R.Q / R.Kp - 1) < 1e-4 ? 'ok' : 'bad')
      + '">' + Math.abs(R.Q / R.Kp - 1).toExponential(2) + '</td></tr></tbody></table></div>';

    h += '<h3>2. Oxygen potential against the phase boundaries</h3>'
      + '<p class="sub">This is the second route: the metal does not enter the gas, so each phase '
      + 'boundary is a plain function of temperature computed from the phase data alone — with no '
      + 'optimiser in it, and using nothing the minimiser produced. Where the minimisation leaves two '
      + 'phases together, coexistence pins the gas exactly onto that boundary.</p>';
    const lad = V.ladder(c.solids, T, c.metal);
    const rows = lad.map(l => {
      const d = R.muO.muO == null ? null : R.muO.muO - l.crit;
      const both = R.split[l.higher] > DETECT_PCT && R.split[l.lower] > DETECT_PCT;
      return '<tr><td>' + sub(l.higher) + ' → ' + sub(l.lower) + '</td><td>'
        + l.crit.toFixed(3) + '</td><td>' + (d == null ? '—'
          : (d > 0 ? '+' : '') + d.toFixed(3)) + '</td><td>'
        + (both ? '<span class="ok">coexisting — gas must sit here</span>'
          : (d > 0 ? 'gas above: this step does not happen' : 'gas below')) + '</td></tr>';
    }).join('');
    h += '<div class="tablewrap"><table><thead><tr><th>step</th>'
      + '<th>boundary μ<sub>O</sub> (kJ/mol-O)</th><th>gas − boundary</th><th></th></tr></thead>'
      + '<tbody>' + rows + '</tbody></table></div>';

    const coex = lad.find(l => R.split[l.higher] > DETECT_PCT && R.split[l.lower] > DETECT_PCT);
    if (coex) {
      const d = R.muO.muO - coex.crit;
      h += '<div class="note"><strong>Two-phase check: '
        + (Math.abs(d) < 0.01 ? '<span class="ok">passed</span>' : '<span class="bad">failed</span>')
        + '.</strong> ' + sub(coex.higher) + ' and ' + sub(coex.lower) + ' coexist, so the gas must sit '
        + 'on their boundary. It is off by ' + Math.abs(d).toExponential(2) + ' kJ per mol O.</div>';
    } else {
      h += '<div class="note">One phase holds the metal at this condition, so the gas is not pinned '
        + 'to any boundary and there is nothing to check here. The column above says how far it is '
        + 'from each step instead.</div>';
    }

    h += '<h3>3. Couple agreement</h3><p class="sub">H<sub>2</sub>/H<sub>2</sub>O and '
      + 'CO/CO<sub>2</sub> are two independent ways to read the same gas. At equilibrium they must '
      + 'give the same oxygen potential.</p><div class="eq">'
      + 'μ<sub>O</sub> = ' + (R.muO.muO == null ? '—' : R.muO.muO.toFixed(4)) + ' kJ/mol-O'
      + '<br/>disagreement between the two couples = '
      + (R.muO.spread == null ? '—' : R.muO.spread.toExponential(2)) + ' kJ/mol-O</div>';

    h += '<div class="note stop"><strong>What this tab is and is not.</strong> The oxygen-potential '
      + 'route above is a check, never the reported answer. It is legitimate and it is far faster, '
      + 'but it assumes the series is monotonic in n — it walks the ladder rung by rung — and in '
      + 'the Waldner assessment that is not true everywhere. Where the two disagree, the reported '
      + 'value is the minimisation.</div>';
    $('verBody').innerHTML = h;
  }

  function renderBalance(c, R) {
    const S = R.S;
    table('tBal', ['element', 'in the charge (mol)', 'at equilibrium (mol)', 'residual (mol)'],
      S.elements.map((e, j) => {
        let v = 0;
        for (let i = 0; i < S.species.length; i++) v += S.E[i][j] * R.n[i];
        return [e, R.b[j].toExponential(8), v.toExponential(8),
          '<span class="' + (Math.abs(v - R.b[j]) < 1e-12 ? 'ok' : 'bad') + '">'
          + (v - R.b[j]).toExponential(2) + '</span>'];
      }));
    $('balNote').innerHTML = '<div class="note">Elements are the constraint, not a result: the '
      + 'minimisation is carried out subject to <code>A·n = A·n<sub>feed</sub></code> exactly. '
      + 'A residual here would mean the optimiser left the feasible set, which is a different failure '
      + 'from finding the wrong minimum inside it. Both are possible; this catches one of them.</div>';
  }

  function renderCoef(c) {
    let h = '<h3>Gases &mdash; NIST-JANAF Shomate</h3><div class="tablewrap"><table>'
      + '<thead><tr><th>species</th><th>A</th><th>B</th><th>C</th><th>D</th><th>E</th>'
      + '<th>F′</th><th>G</th><th>range split (K)</th></tr></thead><tbody>';
    GAS.forEach(s => {
      const e = D.shomate[s];
      h += '<tr><td>' + sub(s) + '</td>' + e.lo.map(v => '<td>' + v.toPrecision(7) + '</td>').join('')
        + '<td>' + e.brk + '</td></tr>';
    });
    h += '</tbody></table></div><div class="note">Low-temperature set shown; the high set takes over '
      + 'above the split. F′ = F − H + Δ<sub>f</sub>H°<sub>298</sub>, folded at export so the '
      + 'seven-coefficient evaluator returns the absolute potential.</div>';

    if (c.name === 'Ti-O') {
      h += '<h3>Ti&ndash;O condensed phases &mdash; CALPHAD (Waldner &amp; Eriksson)</h3>'
        + '<p class="sub">A different C<sub>p</sub> functional form from the gases, evaluated by a '
        + 'separate function. Magnéli form: C<sub>p</sub> = k₀ + k₁T<sup>−0.5</sup> '
        + '+ k₂T<sup>−2</sup> + k₃T<sup>−3</sup>. Poly form: '
        + 'C<sub>p</sub> = a + bT + cT² + d T<sup>−2</sup>.</p>'
        + '<div class="tablewrap"><table><thead><tr><th>phase</th>'
        + '<th>Δ<sub>f</sub>H°<sub>298</sub> (J/mol)</th><th>S°<sub>298</sub></th>'
        + '<th>range (K)</th><th>form</th><th>coefficients</th>'
        + '<th>ΔH<sub>trans</sub> (J/mol)</th></tr></thead><tbody>';
      c.solids.forEach(p => {
        const e = D.calphad[p];
        if (!e) return;
        e.ranges.forEach((r, i) => {
          h += '<tr><td>' + (i === 0 ? sub(p) : '') + '</td><td>'
            + (i === 0 ? e.dHf298.toFixed(0) : '') + '</td><td>'
            + (i === 0 ? e.S298.toFixed(3) : '') + '</td><td>' + r.T0 + '–' + r.T1
            + '</td><td>' + r.form + '</td><td>' + r.c.map(v => v.toPrecision(6)).join(', ')
            + '</td><td>' + (r.dHt || 0) + '</td></tr>';
        });
      });
      h += '</tbody></table></div>';
    } else {
      h += '<h3>Ce&ndash;O condensed phases &mdash; Shomate</h3><div class="tablewrap"><table>'
        + '<thead><tr><th>species</th><th>A</th><th>B</th><th>C</th><th>D</th><th>E</th>'
        + '<th>F′</th><th>G</th></tr></thead><tbody>';
      c.solids.forEach(s => {
        const e = D.shomate[s];
        h += '<tr><td>' + sub(s) + '</td>' + e.lo.map(v => '<td>' + v.toPrecision(7) + '</td>').join('') + '</tr>';
      });
      h += '</tbody></table></div><div class="note stop"><strong>Provenance is thin here.</strong> '
        + 'Both ceria entries carry one coefficient set standing in for the whole range rather than '
        + 'the usual low/high pair, and no fit range is recorded with them. They arrive that way from '
        + 'the notebook this work started in. Treat Ce&ndash;O numbers as provisional until that is '
        + 'traced back to a source.</div>';
    }
    $('coefBody').innerHTML = h;
  }

  /* ---- sweep and CSV ---- */
  function sweep() {
    const c = cfg(), rows = [];
    for (let T_C = 500; T_C <= 1500.001; T_C += 100) {
      const R = solveAt(c, T_C);
      if (R) rows.push(R);
    }
    SWEEP = { c, rows };
    table('tSweep', ['T (&deg;C)', 'solid', c.metal + '(3+) %', 'x', 'CO<sub>2</sub> conv %',
      'CO %', 'H<sub>2</sub>O %', 'CH<sub>4</sub> %', 'μ<sub>O</sub>', 'Q/K<sub>p</sub>', 'residual'],
      rows.map(R => [R.T_C, R.present.map(sub).join(' + '),
        R.reduced > DETECT_PCT ? R.reduced.toFixed(4) : '&lt;10\u207b\u2076',
        R.reduced > DETECT_PCT ? R.x.toExponential(2) : '\u2014',
        R.conv.toFixed(2), (R.y.CO * 100).toFixed(3), (R.y.H2O * 100).toFixed(3),
        (R.y.CH4 * 100).toPrecision(4),
        R.muO.muO == null ? '—' : R.muO.muO.toFixed(2),
        (R.Q / R.Kp).toFixed(6), R.solver.residual.toExponential(1)]));
  }

  function csv() {
    if (!SWEEP) sweep();
    const head = ['system', 'T_C', 'method', 'phases', 'reduced_pct', 'x', 'conv_CO2_pct']
      .concat(GAS.map(g => 'y_' + g)).concat(['mu_O_kJ', 'Kp', 'Q', 'residual_mol']);
    const lines = [head.join(',')].concat(SWEEP.rows.map(R => [
      SWEEP.c.name, R.T_C, R.solver.method, '"' + R.present.join(' + ') + '"',
      R.reduced, R.x, R.conv].concat(GAS.map(g => R.y[g]))
      .concat([R.muO.muO, R.Kp, R.Q, R.solver.residual]).join(',')));
    const blob = new Blob([lines.join('\n')], { type: 'text/csv' });
    const a = document.createElement('a');
    a.href = URL.createObjectURL(blob);
    a.download = 'oxide_' + SWEEP.c.name + '_gibbs_min.csv';
    document.body.appendChild(a); a.click(); a.remove();
  }

  /* ---- saved conditions ---- */
  const KEY = 'oxideSolverConditions';
  const getSaved = () => { try { return JSON.parse(localStorage.getItem(KEY) || '{}'); } catch (e) { return {}; } };
  function setSaved(o) {
    try { localStorage.setItem(KEY, JSON.stringify(o)); return true; }
    catch (e) { $('favNote').innerHTML = '<span class="bad">Browser storage unavailable; not saved.</span>'; return false; }
  }
  const FIELDS = ['sys', 'depth', 'P', 'V', 'mass', 'T'];
  function refreshFav(sel) {
    const o = getSaved();
    $('favSel').innerHTML = Object.keys(o).sort().map(n =>
      '<option' + (n === sel ? ' selected' : '') + '>' + n + '</option>').join('')
      || '<option value="">(none saved)</option>';
  }

  /* ---- wiring ---- */
  $('tabs').addEventListener('click', e => {
    const b = e.target.closest('.tab'); if (!b) return;
    document.querySelectorAll('.tab').forEach(x => x.classList.toggle('active', x === b));
    document.querySelectorAll('.page').forEach(x => x.classList.toggle('active', x.id === b.dataset.page));
  });
  $('T').addEventListener('input', () => { $('Tval').innerHTML = $('T').value + ' &deg;C'; });
  $('sys').addEventListener('change', () => {
    $('depthWrap').hidden = $('sys').value !== 'Ti-O';
    $('sysNote').textContent = D.systems[$('sys').value].source;
    render();
  });
  ['depth', 'P', 'V', 'mass'].forEach(id => $(id).addEventListener('change', render));
  $('T').addEventListener('change', render);
  $('run').addEventListener('click', render);
  $('sweep').addEventListener('click', sweep);
  $('csv').addEventListener('click', csv);
  $('favSave').addEventListener('click', () => {
    const n = ($('favName').value || '').trim();
    if (!n) { $('favNote').textContent = 'Give it a name first.'; return; }
    const o = getSaved(); o[n] = {}; FIELDS.forEach(f => { o[n][f] = $(f).value; });
    if (setSaved(o)) { refreshFav(n); $('favNote').textContent = 'Saved "' + n + '".'; $('favName').value = ''; }
  });
  $('favLoad').addEventListener('click', () => {
    const o = getSaved(), n = $('favSel').value;
    if (!o[n]) return;
    FIELDS.forEach(f => { if (o[n][f] !== undefined) $(f).value = o[n][f]; });
    $('Tval').innerHTML = $('T').value + ' &deg;C';
    $('depthWrap').hidden = $('sys').value !== 'Ti-O';
    $('favNote').textContent = 'Loaded "' + n + '".';
    render();
  });
  $('favDel').addEventListener('click', () => {
    const o = getSaved(), n = $('favSel').value;
    if (!o[n]) return;
    delete o[n]; if (setSaved(o)) { refreshFav(); $('favNote').textContent = 'Deleted "' + n + '".'; }
  });
  $('themeBtn').addEventListener('click', () => {
    const cur = document.documentElement.getAttribute('data-theme');
    const next = cur === 'dark' ? 'light' : cur === 'light' ? null : 'dark';
    if (next) document.documentElement.setAttribute('data-theme', next);
    else document.documentElement.removeAttribute('data-theme');
    $('themeBtn').textContent = 'Theme: ' + (next || 'system');
  });

  $('rawSrc').textContent = [
    Oxide.shomate, Oxide.calphad, Oxide.mu0,
    Oxide.GibbsSystem.prototype.dGForm, Oxide.GibbsSystem.prototype.dmix,
    Oxide.GibbsSystem.prototype.GRel, Oxide.GibbsSystem.prototype.gradGRel,
    Oxide.minimise, Oxide.buildSeeds, Oxide.Verify.muOCritical,
  ].map(f => f.toString()).join('\n\n');

  $('depthWrap').hidden = false;
  $('sysNote').textContent = D.systems['Ti-O'].source;
  $('themeBtn').textContent = 'Theme: system';
  refreshFav();
  render();
})();
