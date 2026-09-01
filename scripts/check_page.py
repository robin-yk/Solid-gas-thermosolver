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
  for (const id of ['figConversion', 'figMargin', 'figBoundary'])
    out.drew[id] = await svg(id);

  for (const T of ['600', '800', '1000']) {
    await pg.selectOption('#eqT', T);
    await pg.waitForTimeout(150);
    out.eq[T] = await txt('eqKpis');
  }

  await pg.click('.wstab[data-ws="ws-surface"]');
  await pg.waitForTimeout(200);
  for (const id of ['figCapacity', 'figBound']) out.drew[id] = await svg(id);
  out.tabSwitch = await pg.evaluate(() => ({
    surface: document.getElementById('ws-surface').style.display,
    equilibrium: document.getElementById('ws-equilibrium').style.display }));

  for (const d of ['0.9', '1.6']) {
    await pg.selectOption('#sfD', d);
    for (const v of ['95', '190']) {
      await pg.selectOption('#sfV', v);
      await pg.waitForTimeout(120);
      out.sf[d + '/' + v] = await txt('sfKpis');
    }
  }
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
    ok('the workspace switch works',
       out['tabSwitch']['surface'] == '' and out['tabSwitch']['equilibrium'] == 'none',
       out['tabSwitch'])
    ok('the page states what it does not claim', out['scope'])

    for T, text in sorted(out['eq'].items()):
        r = A.solve(FEED, float(T) + 273.15)
        ok('%s C conversion matches a fresh solve' % T,
           _pct(r['conversion_CO2_pct']) in text, text)
        b = A.reduction_boundary(float(T) + 273.15)
        ok('%s C rutile boundary matches the closed form' % T,
           '%.4f mol%%' % (b['y_CO2'] * 100) in text, text)

    g = V.load()
    for key, text in sorted(out['sf'].items()):
        d, inv = (float(x) for x in key.split('/'))
        bd = V.surface_fraction_bound(g, inv, d_um=d)
        ok('d=%s, %s umol: capacity matches' % (d, inv),
           '%.2f µmol-O g⁻¹' % bd['bridging_capacity_umol_o_g'] in text, text)
        ok('d=%s, %s umol: the bound matches' % (d, inv),
           _pct(bd['surface_fraction_max'] * 100) in text, text)
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
