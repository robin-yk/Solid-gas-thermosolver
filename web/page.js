/* Page behaviour. The engine in solver.js does the thermodynamics; this file
   only decides what to draw and how. */
'use strict';
const D = JSON.parse(document.getElementById('thermo').textContent);
const E = new Thermo.Engine(D);
const $ = id => document.getElementById(id);
const cv = n => getComputedStyle(document.documentElement).getPropertyValue(n).trim();
const NS = 'http://www.w3.org/2000/svg';
const el = (t, a) => { const e = document.createElementNS(NS, t); for (const k in a) e.setAttribute(k, a[k]); return e; };
/* Ti10O19 -> Ti<sub>10</sub>O<sub>19</sub>; digits after a letter only, so
   a bare number is never mangled. */
const sub = s => String(s).replace(/([A-Za-z])(\d+)/g, '$1<sub>$2</sub>');
const sci = v => { const [m, e] = v.toExponential(2).split('e');
  return m + ' \u00d7 10<sup>' + (+e) + '</sup>'; };
let SRC = 'waldner', LAST = null;

function phaseList() {
  const m = $('deep').value;
  if (m === 'ti4o7') return ['TiO2', 'Ti4O7'];
  if (SRC === 'nist') return D.active.filter(p => D.shomate[p]);
  return m === 'full' ? D.waldner_order : D.active;
}

/* ---------------------------------------------------------------- the chart

   One function for all three plots. Gridlines are hairline and solid, the fill
   under a curve is faint enough to read as weight rather than as a second
   series, the final point carries a dot with a surface-coloured ring, and a
   crosshair reads every series at one temperature. */
function chart(svgId, legId, tipId, series, opt) {
  const svg = $(svgId), tip = $(tipId);
  svg.textContent = '';
  const W = 900, H = svg.viewBox.baseVal.height;
  const M = { l: 64, r: opt.labelRoom || 74, t: 14, b: 40 };
  const iw = W - M.l - M.r, ih = H - M.t - M.b, log = !!opt.log;
  const vals = series.flatMap(s => s.y).filter(v => v != null && isFinite(v) && (!log || v > 0));
  if (!vals.length) return;
  const x0 = opt.x0, x1 = opt.x1;
  let y0 = opt.y0 != null ? opt.y0 : Math.min(...vals);
  let y1 = opt.y1 != null ? opt.y1 : Math.max(...vals);
  if (y0 === y1) { y0 -= 1; y1 += 1; }
  if (opt.y0 == null && opt.y1 == null) { const pad = (y1 - y0) * 0.08; y0 -= pad; y1 += pad; }
  const X = v => M.l + (v - x0) / (x1 - x0) * iw;
  const Y = log
    ? v => M.t + ih - (Math.log10(Math.max(v, y0)) - Math.log10(y0)) / (Math.log10(y1) - Math.log10(y0)) * ih
    : v => M.t + ih - (v - y0) / (y1 - y0) * ih;

  const rule = cv('--rule'), muted = cv('--muted'), surface = cv('--surface');
  const ticks = log
    ? (() => { const a = []; for (let k = Math.ceil(Math.log10(y0)); k <= Math.floor(Math.log10(y1)); k++) a.push(10 ** k); return a; })()
    : Array.from({ length: 6 }, (_, i) => y0 + (y1 - y0) * i / 5);
  ticks.forEach(t => {
    svg.appendChild(el('line', { x1: M.l, x2: M.l + iw, y1: Y(t), y2: Y(t), stroke: rule, 'stroke-width': 1 }));
    const lab = log ? '10' : (+t.toFixed(opt.yfix || 0)).toLocaleString();
    const tx = el('text', { x: M.l - 9, y: Y(t) + 4, 'text-anchor': 'end', fill: muted, 'font-size': 11,
                            'font-family': 'ui-monospace,monospace' });
    tx.textContent = lab;
    if (log) { const e = el('tspan', { 'font-size': 8, dy: -4 }); e.textContent = Math.round(Math.log10(t)); tx.appendChild(e); }
    svg.appendChild(tx);
  });
  for (let t = x0; t <= x1 + 1; t += opt.xstep || 200) {
    svg.appendChild(el('line', { x1: X(t), x2: X(t), y1: M.t, y2: M.t + ih, stroke: rule, 'stroke-width': 1 }));
    const tx = el('text', { x: X(t), y: H - 16, 'text-anchor': 'middle', fill: muted, 'font-size': 11,
                            'font-family': 'ui-monospace,monospace' });
    tx.textContent = t; svg.appendChild(tx);
  }
  const ax = el('text', { x: M.l + iw / 2, y: H - 2, 'text-anchor': 'middle', fill: muted, 'font-size': 11 });
  ax.textContent = 'Temperature (°C)'; svg.appendChild(ax);
  const ay = el('text', { x: 13, y: M.t + ih / 2, fill: muted, 'font-size': 11, 'text-anchor': 'middle',
                          transform: 'rotate(-90 13 ' + (M.t + ih / 2) + ')' });
  ay.textContent = opt.ylabel; svg.appendChild(ay);

  const clean = series.map(s => ({ s, pts: s.x.map((v, i) => [v, s.y[i]])
    .filter(p => p[1] != null && isFinite(p[1]) && (!log || p[1] > 0)) })).filter(o => o.pts.length);

  /* fill first, so every stroke sits above every fill */
  clean.forEach(({ s, pts }) => {
    if (!s.fill) return;
    const d = pts.map((p, i) => (i ? 'L' : 'M') + X(p[0]).toFixed(1) + ' ' + Y(p[1]).toFixed(1)).join(' ')
      + ' L' + X(pts[pts.length - 1][0]).toFixed(1) + ' ' + (M.t + ih) + ' L' + X(pts[0][0]).toFixed(1) + ' ' + (M.t + ih) + ' Z';
    svg.appendChild(el('path', { d, fill: cv(s.color), opacity: s.fill, stroke: 'none' }));
  });
  clean.forEach(({ s, pts }) => {
    const d = pts.map((p, i) => (i ? 'L' : 'M') + X(p[0]).toFixed(1) + ' ' + Y(p[1]).toFixed(1)).join(' ');
    svg.appendChild(el('path', { d, fill: 'none', stroke: cv(s.color), 'stroke-width': s.dash ? 1.8 : 2,
                                 'stroke-dasharray': s.dash || '', 'stroke-linejoin': 'round', 'stroke-linecap': 'round' }));
  });

  /* end dot and a direct label, but only where labels will not collide */
  const ends = clean.map(({ s, pts }) => ({ s, p: pts[pts.length - 1] }))
                    .sort((a, b) => Y(a.p[1]) - Y(b.p[1]));
  const placed = [];
  ends.forEach(({ s, p }) => {
    svg.appendChild(el('circle', { cx: X(p[0]), cy: Y(p[1]), r: 4.5, fill: cv(s.color),
                                   stroke: surface, 'stroke-width': 2 }));
    if (!opt.endLabels) return;
    const y = Y(p[1]);
    if (placed.some(q => Math.abs(q - y) < 13)) return;
    placed.push(y);
    const tx = el('text', { x: X(p[0]) + 9, y: y + 4, fill: muted, 'font-size': 11 });
    tx.textContent = s.short || s.label; svg.appendChild(tx);
  });

  /* crosshair layer */
  const hair = el('line', { x1: 0, x2: 0, y1: M.t, y2: M.t + ih, stroke: muted, 'stroke-width': 1, opacity: 0 });
  svg.appendChild(hair);
  const hit = el('rect', { x: M.l, y: M.t, width: iw, height: ih, fill: 'transparent' });
  svg.appendChild(hit);
  function move(ev) {
    const r = svg.getBoundingClientRect();
    const px = (ev.clientX - r.left) / r.width * W;
    if (px < M.l || px > M.l + iw) return hide();
    const t = x0 + (px - M.l) / iw * (x1 - x0);
    const idx = series[0].x.reduce((a, v, i) => Math.abs(v - t) < Math.abs(series[0].x[a] - t) ? i : a, 0);
    const tx = series[0].x[idx];
    hair.setAttribute('x1', X(tx)); hair.setAttribute('x2', X(tx)); hair.setAttribute('opacity', .45);
    const head = document.createElement('div'); head.className = 'th';
    head.textContent = tx.toFixed(0) + ' °C';
    tip.textContent = ''; tip.appendChild(head);
    series.forEach(s => {
      const v = s.y[idx];
      if (v == null || !isFinite(v)) return;
      const row = document.createElement('div'); row.className = 'r';
      const key = document.createElement('i'); key.style.background = cv(s.color);
      const val = document.createElement('b'); val.textContent = opt.fmt ? opt.fmt(v) : v.toFixed(1);
      const nm = document.createElement('span'); nm.textContent = s.label;
      row.append(key, val, nm); tip.appendChild(row);
    });
    tip.style.opacity = 1;
    const wrapW = svg.parentElement.clientWidth, tw = tip.offsetWidth;
    const left = X(tx) / W * wrapW + 14;
    tip.style.left = Math.min(left, wrapW - tw - 8) + 'px';
    tip.style.top = '14px';
  }
  function hide(){ hair.setAttribute('opacity', 0); tip.style.opacity = 0; }
  hit.addEventListener('pointermove', move);
  hit.addEventListener('pointerleave', hide);

  const lg = $(legId); lg.textContent = '';
  series.forEach(s => {
    const sp = document.createElement('span'), i = document.createElement('i');
    i.style.background = s.dash ? 'repeating-linear-gradient(90deg,' + cv(s.color) + ' 0 5px,transparent 5px 9px)' : cv(s.color);
    sp.appendChild(i);
    const t = document.createElement('span'); t.innerHTML = s.legend || s.label;
    sp.appendChild(t); lg.appendChild(sp);
  });
}

function table(id, cols, rows){
  const t = $(id); t.textContent = '';
  const th = document.createElement('thead'), tr = document.createElement('tr');
  cols.forEach(c => { const e = document.createElement('th'); e.innerHTML = c; tr.appendChild(e); });
  th.appendChild(tr); t.appendChild(th);
  const tb = document.createElement('tbody');
  rows.forEach(r => { const e = document.createElement('tr');
    r.forEach(c => { const d = document.createElement('td'); d.innerHTML = c; e.appendChild(d); });
    tb.appendChild(e); });
  t.appendChild(tb);
}

function compute(){
  const feed = $('feed').value, P = +$('P').value || 1, V = +$('vol').value || 25;
  const nTi = (+$('mass').value || 100) / 1000 / 79.87;
  const pO2imp = (+$('imp').value || 0) * 1e-6 * P;
  const phases = phaseList(), rows = [];
  for (let TC = 500; TC <= 1500.001; TC += 12.5){
    const T = TC + 273.15, g = P * V / (D.R_ATM * T);
    let rec = { TC };
    if (feed === 'inert'){
      const r = E.inert(T, P, V, nTi, pO2imp, phases, SRC);
      const mu = pO2imp > 0 ? E.pO2ToMuO(pO2imp, T) : null;
      rec = Object.assign(rec, { muO: mu, phase: r.phase, crit: r.crit, ti3: r.ti3,
        margin: mu == null ? null : mu - r.crit, pO2: pO2imp, y: null, reduces: r.reduces && r.ti3 > 1e-6 });
    } else if (feed === 'CO2H2'){
      const n0 = g / 2, gas = E.gasEquilibrium(n0, 2 * n0, 2 * n0, T, P);
      if (!gas) continue;
      const mu = E.muOFromRatio(T, gas.y.H2O / gas.y.H2);
      const f = E.firstPhase(mu, T, phases, SRC);
      let ti3 = 0, ph = f.phase;
      if (f.margin < 0){ const b = E.buffered('CO2H2', T, P, V, nTi, phases, SRC); ti3 = b.ti3; ph = b.phase; }
      rec = Object.assign(rec, { muO: mu, phase: ph, crit: f.crit, margin: f.margin, ti3,
        y: gas.y, pO2: E.muOToPO2(mu, T), reduces: f.margin < 0 });
    } else {
      const b = E.buffered(feed, T, P, V, nTi, phases, SRC);
      const mu = isFinite(b.muO) ? b.muO : b.crit;
      rec = Object.assign(rec, { muO: mu, phase: b.phase, crit: b.crit, ti3: b.ti3,
        margin: b.reduces ? 0 : (b.margin || 0), y: b.gas ? b.gas.y : null,
        pO2: E.muOToPO2(mu, T), reduces: b.reduces });
    }
    rows.push(rec);
  }
  return { rows, feed, phases };
}

function render(){
  const { rows, feed, phases } = compute();
  if (!rows.length) return;
  LAST = rows;
  const T = rows.map(r => r.TC), Tsel = +$('Tpick').value;
  const at = rows.reduce((a,b) => Math.abs(b.TC - Tsel) < Math.abs(a.TC - Tsel) ? b : a);
  const anyRed = rows.some(r => r.reduces);

  const bd = $('badge');
  /* danger is the reduced state - that is the outcome being tested for, and
     the brief reserves the colour for it. Rutile surviving is the safe read. */
  bd.className = 'badge' + (anyRed ? '' : ' safe');
  bd.innerHTML = anyRed
    ? 'Reduces. <em>At ' + at.TC.toFixed(0) + ' °C the solid is ' + sub(at.phase) + '</em>'
    : 'No reduction. <em>Rutile survives the whole range</em>';

  /* A buffered solid sits exactly on a boundary, so its margin is zero by
     construction. Printing 0.0 reads as a near miss; it is not a measurement. */
  const marginVals = rows.filter(r => !r.reduces).map(r => r.margin).filter(v => v != null);
  const worst = marginVals.length ? Math.min(...marginVals) : null;
  const dash = '\u2014';
  const kpis = [
    ['Margin',
     at.reduces ? 'on boundary'
                : (at.margin == null ? dash : (at.margin > 0 ? '+' : '') + at.margin.toFixed(1)),
     at.reduces ? 'the solid holds the gas there'
                : 'kJ per mol O \u00b7 closest ' + (worst == null ? dash : worst.toFixed(1)),
     at.reduces ? 'is-danger is-word' : (at.margin == null ? '' : 'is-safe')],
    ['\u03bc<sub>O</sub> (gas)', at.muO == null ? '\u2212\u221e' : at.muO.toFixed(1), 'kJ per mol O', ''],
    ['pO<sub>2</sub>', at.pO2 > 0 ? sci(at.pO2) : '0', 'atm', ''],
    ['Ti(3+)', at.ti3.toFixed(at.ti3 < 1 ? 3 : 2) + '%', 'of total titanium',
     at.ti3 > 0.01 ? 'is-danger' : ''],
  ];
  $('kpis').innerHTML = kpis.map(c =>
    '<div class="kpi ' + c[3] + '"><div class="k">' + c[0] + '</div><div class="v">' + c[1] + '</div><div class="u">' + c[2] + '</div></div>').join('');

  /* boundaries are ordinal - one hue, stepped by depth in the series */
  const ords = ['--ord-1','--ord-2','--ord-3','--ord-4','--ord-5'];
  const bounds = phases.filter(p => p !== 'TiO2').slice(0, 5).map((p, i) => ({
    x: T, y: rows.map(r => E.muOCritical(p, r.TC + 273.15, SRC)),
    color: ords[i], label: 'TiO2 -> ' + p, legend: sub('TiO2') + ' → ' + sub(p),
    short: p, dash: '5 4' }));
  const gasLine = { x: T, y: rows.map(r => r.muO), color: '--accent', label: 'gas',
                    legend: 'gas', short: 'gas', fill: .1 };
  chart('chMu','lgMu','tipMu',[gasLine, ...bounds],
        { x0:500, x1:1500, ylabel:'μO (kJ per mol O)', yfix:0, endLabels:true, labelRoom:96,
          fmt:v => v.toFixed(1) });
  $('capMu').textContent = anyRed
    ? 'The gas line sits on a boundary: the solid holds it there while both phases coexist.'
    : 'The gas line stays above every boundary. The vertical distance is the margin.';

  const gs = ['--series-1','--series-2','--series-3','--series-4','--series-5'];
  const names = ['CO2','H2','CO','H2O','CH4'];
  if (rows.some(r => r.y)){
    chart('chGas','lgGas','tipGas', names.map((n,i) => ({
      x:T, y: rows.map(r => r.y ? r.y[n] * 100 : null), color: gs[i],
      label: n, legend: sub(n), short: n })),
      { x0:500, x1:1500, ylabel:'gas mol %', y0:0, y1:100, endLabels:true,
        fmt: v => v.toFixed(2) + '%' });
  } else {
    $('chGas').textContent=''; $('lgGas').textContent='';
    const t = el('text',{ x:450, y:150, 'text-anchor':'middle', fill: cv('--muted'), 'font-size':13 });
    t.textContent = 'inert carrier — no C–H–O gas to equilibrate';
    $('chGas').appendChild(t);
  }

  chart('chTi','lgTi','tipTi',[{ x:T, y: rows.map(r => r.ti3), color:'--danger',
    label:'Ti(3+)', legend:'Ti(3+) % of total Ti', fill:.13 }],
    { x0:500, x1:1500, ylabel:'Ti(3+) %', y0:0, fmt: v => v.toFixed(3) + '%' });

  const pick = rows.filter(r => Math.abs(r.TC % 200) < 7);
  table('tbl',['T (°C)','\u03bc<sub>O</sub>','boundary','margin','pO<sub>2</sub> (atm)','Ti(3+) %','phase'],
    pick.map(r => [r.TC.toFixed(0), r.muO == null ? '−∞' : r.muO.toFixed(1), r.crit.toFixed(1),
      r.reduces ? 'buffered'
        : (r.margin == null ? '\u2014' : (r.margin > 0 ? '+' : '') + r.margin.toFixed(1)),
      r.pO2 > 0 ? r.pO2.toExponential(2) : '0', r.ti3.toFixed(3), sub(r.phase)]));

  const notes = [];
  if (feed === 'CH4') notes.push('<strong>CH4 is provisional.</strong> Without graphite in the species list all carbon must leave as CO, and the H2:CO = 2:1 that comes out may be an artefact of that.');
  if ($('deep').value === 'full') notes.push('<strong>Ti20O39 is the model edge, not a phase you would find.</strong> Nineteen of twenty titaniums are still Ti(4+); there the series stands in for oxygen-deficient rutile.');
  if (SRC === 'nist') notes.push('<strong>NIST stops at Ti4O7.</strong> The shallow Magnéli members are missing, so margins here are measured against too hard a step and read high.');
  if (feed === 'inert') notes.push('<strong>The impurity level decides this, not the temperature.</strong> Set O2 to zero and the oxide still barely moves: a sealed headspace cannot hold enough O2 for the reduction to go anywhere.');
  $('caveat').innerHTML = notes.join(' ') || '<strong>Nothing to flag.</strong> These conditions sit inside every fit range and use a single consistent dataset.';
  $('srcNote').textContent = SRC === 'waldner' ? 'full Magnéli series' : 'stops at Ti4O7';
  $('railnote').innerHTML = 'Boundaries take an ordinal ramp — they are ordered by O/Ti. Gas species take separate hues; a one-hue ramp fails colour-vision separation for them by measurement.';
  $('prov').textContent = 'Computed live in this page: ' + rows.length + ' temperatures × ' + phases.length + ' phases. ';
}

/* -- reported values ------------------------------------------------------
   Rendered straight out of the embedded results.json. This table never asks
   the live solver below for anything: the solver runs the oxygen-potential
   route, and these numbers are the Gibbs minimisation. Keeping them separate
   in code is what keeps them from quietly becoming the same thing. */
const RES = JSON.parse(document.getElementById('results').textContent);

function fmtTi(v){
  if (v === 0) return '0';
  if (v < 1e-3) return sci(v);
  return v.toFixed(v < 1 ? 4 : 2);
}

function gibbsTable(){
  const rows = RES.rows.map(r => {
    const gas = Object.entries(r.gas_pct).filter(([, v]) => v > 0.05).slice(0, 3)
      .map(([k, v]) => sub(k) + ' ' + v.toFixed(1)).join(' \u00b7 ');
    const solid = Object.entries(r.phase_split_pct).filter(([, v]) => v > 1e-6)
      .map(([k, v]) => sub(k) + (v > 99.999 ? '' : ' ' + v.toFixed(1) + '%')).join(' + ');
    const red = r.ti3_pct > 0.01;
    return ['<span class="' + (red ? 'is-danger' : '') + '">' + r.label + '</span>',
            r.T_C, fmtTi(r.ti3_pct) + '%', solid || sub('TiO2'), gas || '\u2014',
            r.pO2_atm > 0 ? sci(r.pO2_atm) : '\u2014',
            '<code>' + r.method + '</code>'];
  });
  table('gibbsTbl', ['Atmosphere', 'T (&deg;C)', 'Ti(3+)', 'Solid', 'Gas (mol%)',
                     'pO<sub>2</sub> (atm)', 'method'], rows);

  /* The worst disagreement with the cross-check route, stated rather than
     buried: a check nobody can see the result of is not a check. */
  const d = RES.rows.map(r => r.crosscheck && r.crosscheck.d_mu_O_kJ)
                    .filter(v => v != null).map(Math.abs);
  const worst = d.length ? Math.max(...d) : null;
  $('gibbsProv').innerHTML =
    RES.rows.length + ' points, all <code>' + RES.method + '</code>. Solid data: '
    + RES.conditions.solid_data + '; gases: ' + RES.conditions.gas_data + '. '
    + RES.conditions.m_TiO2_g * 1000 + ' mg TiO2 in ' + RES.conditions.V_cm3
    + ' mL at ' + RES.conditions.P_atm + ' atm. '
    + (worst == null ? ''
       : 'Largest disagreement with the oxygen-potential cross-check where two phases fix the gas: '
         + worst.toFixed(3) + ' kJ per mol O.');
}

gibbsTable();

['feed','deep','P','mass','vol','imp','Tpick'].forEach(id =>
  $(id).addEventListener('input', () => {
    $('impWrap').hidden = $('feed').value !== 'inert';
    $('TpickVal').innerHTML = $('Tpick').value + ' &deg;C';
    render();
  }));
$('src').addEventListener('click', e => {
  const b = e.target.closest('button'); if (!b) return;
  SRC = b.dataset.v;
  [...$('src').children].forEach(x => x.setAttribute('aria-pressed', String(x === b)));
  render();
});
matchMedia('(prefers-color-scheme:dark)').addEventListener('change', render);
$('impWrap').hidden = true;
render();
