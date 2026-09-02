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
from solidgas import vacancy_inventory as V    # noqa: E402
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

  const out = { errors: errs, drew: {}, eq: {}, sf: {} };
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

  await pg.click('.wstab[data-ws="ws-surface"]');
  await pg.waitForTimeout(200);
  for (const id of ['figCapacity', 'figBound']) out.drew[id] = await svg(id);
  out.tabSwitch = await pg.evaluate(() => ({
    surface: document.getElementById('ws-surface').style.display,
    equilibrium: document.getElementById('ws-thermo').style.display }));

  for (const d of ['0.9', '1.6']) {
    await pg.fill('#sfD', d);
    for (const v of ['95', '190']) {
      await pg.fill('#sfV', v);
      await pg.waitForTimeout(150);
      out.sf[d + '/' + v] = await txt('sfKpis');
    }
  }
  // a measured BET area must override the sphere estimate
  await pg.fill('#sfBET', '3.25');
  await pg.fill('#sfV', '95');
  await pg.waitForTimeout(150);
  out.bet = await txt('sfKpis');
  await pg.click('.wstab[data-ws="ws-population"]');
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

  out.scope = await pg.evaluate(() =>
    document.body.textContent.indexOf('not assigned here') >= 0);
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
    ok('the workspace switch works',
       out['tabSwitch']['surface'] == '' and out['tabSwitch']['equilibrium'] == 'none',
       out['tabSwitch'])
    ok('the page states what it does not claim', out['scope'])

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
       and 'exact by construction' in out['popKpis'], out['popKpis'])

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

    g = V.load()
    for key, text in sorted(out['sf'].items()):
        d, inv = (float(x) for x in key.split('/'))
        bd = V.surface_fraction_bound(g, inv, d_um=d)
        ok('d=%s, %s umol: capacity matches' % (d, inv),
           '%.2f µmol-O g⁻¹' % bd['bridging_capacity_umol_o_g'] in text, text)
        ok('d=%s, %s umol: the bound matches' % (d, inv),
           _pct(bd['surface_fraction_max'] * 100) in text, text)

    bd = V.surface_fraction_bound(g, 95.0, bet_m2_g=3.25)
    ok('a measured BET area overrides the sphere estimate',
       '%.2f µmol-O g⁻¹' % bd['bridging_capacity_umol_o_g'] in out['bet'],
       out['bet'])
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
