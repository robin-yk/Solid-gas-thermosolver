"""Load the built page in a real browser and read back what it drew.

The test suite cannot assume a browser, and a gate that skips itself is
worse than no gate, so this lives here and is run by hand or by CI.

    PW_DIR=/path/with/node_modules/playwright python3 scripts/check_page.py

Exits non-zero on any page error, any figure that did not draw, or any
reading that disagrees with a fresh calculation by the Python package.
"""

import json
import os
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
PAGE = os.path.join(ROOT, 'docs', 'index.html')

from solidgas import activeset as A            # noqa: E402
from solidgas import vacancy_population as VP  # noqa: E402

FEED = {'CO2': 10, 'H2': 10, 'N2': 80}

SCRIPT = r'''
const P = require(process.env.PW + '/node_modules/playwright');
(async () => {
  const b = await P.chromium.launch();
  const pg = await b.newPage();

  const errs = [];
  pg.on('pageerror', e => errs.push(String(e)));
  pg.on('console', m => { if (m.type() === 'error') errs.push('console: ' + m.text()); });
  await pg.goto('file://' + process.argv[2], { waitUntil: 'load' });
  await pg.waitForTimeout(400);
  const svg = id => pg.evaluate(i => {
    const el = document.getElementById(i);
    return el && el.querySelector('svg') ? 1 : 0;
  }, id);
  const txt = id => pg.evaluate(i => {
    const el = document.getElementById(i);
    return el ? el.textContent.replace(/\s+/g, ' ').trim() : null;
  }, id);

  const out = { errors: errs, drew: {}, eq: {} };
  /* the page opens on the home screen: three cards, no workspace */
  out.home = await pg.evaluate(() => ({
    cards: document.querySelectorAll('.card').length,
    thermo: document.getElementById('ws-thermo').style.display,
    topbar: document.getElementById('topbar').style.display }));
  await pg.click('.card[data-ws="ws-thermo"]');
  await pg.waitForTimeout(200);
  out.hash = await pg.evaluate(() => location.hash);
  await pg.click('#run');
  await pg.waitForTimeout(900);
  for (const id of ['figOptimality', 'figBoundary'])
    out.drew[id] = await svg(id);
  /* the two swept figures cost sixty solves, so the page keeps them behind
     a button; the gate presses it rather than pretending they draw free */
  out.beforeSweep = { margin: await svg('figMargin'),
                      conversion: await svg('figConversion') };
  await pg.click('#sweepBtn');
  await pg.waitForTimeout(2500);
  for (const id of ['figMargin', 'figConversion']) out.drew[id] = await svg(id);

  const gasIds = await pg.evaluate(() =>
    Array.from(document.querySelectorAll('#gasInputs input')).map(i => i.id));
  const solve = async (T, gas) => {
    await pg.fill('#T', String(T));
    for (const id of gasIds)
      await pg.fill('#' + id, String(gas[id.replace('gas_', '')] || 0));
    await pg.click('#run');
    await pg.waitForTimeout(700);
    return txt('kpis');
  };
  for (const T of ['600', '800', '1000'])
    out.eq[T] = await solve(T, { CO2: 10, H2: 10, N2: 80 });
  out.pure_h2 = await solve(1200, { H2: 1 });

  await pg.click('.wstab[data-ws="ws-population"]');
  await pg.waitForTimeout(200);
  out.tabSwitch = await pg.evaluate(() => ({
    population: document.getElementById('ws-population').style.display,
    equilibrium: document.getElementById('ws-thermo').style.display }));
  await pg.waitForTimeout(250);
  for (const id of ['figPopulation', 'figRates']) out.drew[id] = await svg(id);
  /* the scatter is the claim, so count the marks the browser actually put
     on it rather than trusting that an <svg> appeared */
  out.popMarks = await pg.evaluate(() => {
    const s = document.querySelector('#figPopulation svg');
    if (!s) return null;
    const n = sel => s.querySelectorAll(sel).length;
    return { circle: n('circle[stroke="#0072B2"]'),
             square: n('rect[stroke="#D55E00"]'),
             triangle: n('path.marker[stroke="#777777"]') };
  });
  out.popKpis = await txt('pnKpis');
  out.popFit = await txt('pnIdent');

  /* the distribution workspace: the TOF figure, the rendered panels, and the
     numbers it prints for the starting scenario and after a change */
  await pg.click('.wstab[data-ws="ws-distribution"]');
  await pg.waitForTimeout(6000);
  out.drew.figTofRange = await pg.evaluate(() =>
    document.querySelector('#figTofRange svg') ? 1 : 0);
  out.tofMarks = await pg.evaluate(() => {
    const s = document.querySelector('#figTofRange svg');
    return { fixed: s.querySelectorAll('rect[stroke="#777777"]').length,
             chosen: s.querySelectorAll('circle[stroke="#0072B2"]').length };
  });
  /* the rendered panels: both images decoded, and a hover over the cut
     face finds a layer in the mask and shows its amount */
  out.renders = await pg.evaluate(() => ['vdParticleImg', 'vdSlabImg'].map(id => {
    const im = document.getElementById(id); return im && im.complete ? im.naturalWidth : 0; }));
  out.hover = await pg.evaluate(() => {
    const st = document.getElementById('vdParticleStage').getBoundingClientRect();
    const found = {};
    for (let fy = 0.3; fy < 0.8; fy += 0.02) for (let fx = 0.3; fx < 0.8; fx += 0.02) {
      const ev = { clientX: st.left + fx * st.width, clientY: st.top + fy * st.height };
      const L = window.VacancyDistribution.layerAt(ev);
      if (L && !found[L]) found[L] = [fx, fy];
    }
    return found;
  });
  if (out.hover.bulk) {
    const box = await pg.locator('#vdParticleStage').boundingBox();
    await pg.mouse.move(box.x + out.hover.bulk[0] * box.width, box.y + out.hover.bulk[1] * box.height);
    await pg.waitForTimeout(200);
    out.tip = await txt('vdTip');
  }
  out.vdKpis = await txt('vdKpis');
  out.vdState = await pg.evaluate(() => window.VacancyDistribution.state());
  await pg.selectOption('#vdSample', 'R1000');
  await pg.selectOption('#vdReactive', 'ISO_z4');
  await pg.waitForTimeout(3000);
  out.vdKpis2 = await txt('vdKpis');
  out.vdState2 = await pg.evaluate(() => window.VacancyDistribution.state());
  out.slabShown = await pg.evaluate(() => window.VacancyDistribution.slabShown());

  out.scope = await pg.evaluate(() =>
    document.body.textContent.replace(/\s+/g, ' ').indexOf('is inferred from the measured rates') >= 0);
  console.log(JSON.stringify(out));
  await b.close();
})();
'''


def _pct(v, d=1):
    return '%.*f%%' % (d, v)


def check(out):
    fails = []

    def ok(name, cond, got=''):
        print('%-46s %s%s' % (name, 'ok' if cond else 'FAIL',
                              '' if cond else '   ' + str(got)))
        if not cond:
            fails.append(name)

    ok('no page errors', not out['errors'], out['errors'])
    for k, v in sorted(out['drew'].items()):
        ok('%s drew' % k, v == 1)
    ok('the swept figures wait for the button',
       out['beforeSweep'] == {'margin': 0, 'conversion': 0}, out['beforeSweep'])
    ok('the page opens on the home screen',
       out['home'] == {'cards': 3, 'thermo': 'none', 'topbar': 'none'}, out['home'])
    ok('a card opens its workspace under its hash', out['hash'] == '#equilibrium', out['hash'])
    ok('the workspace switch works',
       out['tabSwitch']['population'] == '' and out['tabSwitch']['equilibrium'] == 'none',
       out['tabSwitch'])
    ok('the page states its inputs and limits', out['scope'])

    n = len(VP.load_series())
    ok('the population scatter plots all three assignments',
       out['popMarks'] == {'circle': n + 1, 'square': n + 1,
                           'triangle': n + 1}, out['popMarks'])
    rows = VP.partition(VP.load_series(), VP.REPORTED['ns'],
                        VP.REPORTED['E_loss'],
                        10.0 ** VP.REPORTED['log10_nu'])
    hi = max(q['n_iso_umol_g'] for q in rows)
    ok('the population workspace opens at the reported partition',
       ('%.2f' % hi) in out['popKpis'], out['popKpis'])
    ok('mass balance closes exactly',
       max(abs(q['closure_umol_g']) for q in rows) == 0.0
       and 'closes exactly' in out['popFit'], out['popFit'])

    for T, text in sorted(out['eq'].items()):
        r = A.solve(FEED, float(T) + 273.15)
        ok('%s C: the page solved it and got the package answer' % T,
           ('%.2f' % r['conversion_CO2_pct']) in text, text)
        ok('%s C: the page names the same assemblage' % T,
           '+'.join(r['active_condensed_phases']).replace('TiO2', 'TiO')
           in text.replace('₂', '2').replace('TiO2', 'TiO'), text)

    r = A.solve({'H2': 1}, 1473.15)
    ok('pure H2 at 1200 C reduces, and the page says so',
       'Ti10O19' in out['pure_h2'].replace('₁', '1').replace('₀', '0')
       .replace('₉', '9').replace('Ti10O19', 'Ti10O19'), out['pure_h2'])

    tof = json.load(open(os.path.join(ROOT, 'paper_outputs', 'tof_range.json')))
    n = len(tof['samples']) + 1                  # one of each in the legend
    ok('the TOF figure marks every sample twice',
       out['tofMarks'] == {'fixed': n, 'chosen': n},
       out['tofMarks'])
    ok('both rendered panels decoded', min(out['renders']) >= 1000, out['renders'])
    ok('hover finds all three layers in the mask',
       set(out['hover']) == {'surface', 'subsurface', 'bulk'}, out['hover'])
    ok('hover on the bulk shows its calculated amount',
       'Bulk, calculated' in (out.get('tip') or ''), out.get('tip'))
    for key, state in (('vdKpis', out['vdState']), ('vdKpis2', out['vdState2'])):
        ax = tof['axes']
        i = ax['energy_map'].index(state['map'])
        i = i * len(ax['cutoff_nm']) + ax['cutoff_nm'].index(state['cutoff'])
        i = i * len(ax['dG_eV']) + ax['dG_eV'].index(state['dG'])
        i = i * len(ax['f110']) + ax['f110'].index(state['f110'])
        i = i * len(ax['eps']) + [e['key'] for e in ax['eps']].index(state['eps'])
        smp = [q for q in tof['samples'] if q['sample'] == state['sample']][0]
        n = smp['cases']['sites'][i][tof['reactive'].index(state['reactive'])]
        want = float('%.3g' % (smp['rate_co_umol_g_s'] / n))
        ok('%s %s: the page prints the stored TOF' % (state['sample'], state['reactive']),
           ('%g' % want) in out[key], (want, out[key]))
    ok('the surface panel follows the sample', out['slabShown'] == 'R1000', out['slabShown'])

    return fails


def main():
    pw = os.environ.get('PW_DIR')
    if not pw or not os.path.isdir(os.path.join(pw, 'node_modules', 'playwright')):
        print('set PW_DIR to a directory holding node_modules/playwright')
        return 2
    with tempfile.NamedTemporaryFile('w', suffix='.js', delete=False) as fh:
        fh.write(SCRIPT)
        path = fh.name
    try:
        r = subprocess.run(['node', path, PAGE], capture_output=True,
                           text=True, env=dict(os.environ, PW=pw), timeout=300)
    finally:
        os.unlink(path)
    if r.returncode != 0:
        print(r.stdout[-2000:])
        print(r.stderr[-2000:])
        return 1
    fails = check(json.loads(r.stdout))
    print()
    print('%d checks failed' % len(fails) if fails else 'all checks passed')
    return 1 if fails else 0


if __name__ == '__main__':
    sys.exit(main())
