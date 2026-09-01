"""Load the built page in a real browser and read what it drew back.

The test suite cannot assume a browser, and a gate that skips itself is
worse than a script: it reports green while checking nothing.  So this
lives here and is run by hand (or by CI where a browser exists).

    PW_DIR=/path/with/node_modules/playwright python3 scripts/check_page.py

Exits non-zero on any page error, any figure that did not draw, or any
reading that disagrees with a fresh solve of the Python engine.
"""

import json
import os
import re
import subprocess
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from solidgas import activeset as A
from solidgas import particle as P

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PAGE = os.path.join(ROOT, 'docs', 'ti_solver.html')

SCRIPT = r'''
const {chromium} = require(process.env.PW + '/node_modules/playwright');
(async () => {
  const b = await chromium.launch({executablePath: '/opt/pw-browsers/chromium'});
  const pg = await b.newPage();
  const errs = [];
  pg.on('pageerror', e => errs.push(String(e)));
  pg.on('console', m => { if (m.type() === 'error') errs.push('console: ' + m.text()); });
  await pg.goto('file://' + process.argv[2] + '#defect', {waitUntil: 'load'});
  await pg.waitForTimeout(2500);
  const out = await pg.evaluate(() => ({
    depth: !!document.querySelector('#figDepth svg'),
    where: !!document.querySelector('#figWhere svg'),
    state: (document.getElementById('pxState') || {}).textContent || '',
    reach: (document.getElementById('pxReach') || {}).textContent || '',
    charge: (document.getElementById('pxCharge') || {}).textContent || ''
  }));

  /* the equilibrium workspace: the line draws with no sweep, the
     conversion curve only after one, and both headers have to carry the
     numbers a fresh Python solve gives at the panel's own condition */
  await pg.goto('file://' + process.argv[2] + '#equilibrium', {waitUntil: 'load'});
  await pg.waitForTimeout(2500);
  out.boundary_before_sweep = await pg.evaluate(
    () => !!document.querySelector('#figBoundary svg'));
  out.conversion_before_sweep = await pg.evaluate(
    () => !!document.querySelector('#figConversion svg'));
  await pg.locator('text=Sweep temperature').first().click();
  await pg.waitForTimeout(12000);
  Object.assign(out, await pg.evaluate(() => {
    const txt = id => {
      const e = document.querySelector('#' + id + ' svg');
      return e ? e.textContent : '';
    };
    return { boundary: txt('figBoundary'), conversion: txt('figConversion'),
             margin: !!document.querySelector('#figMargin svg'),
             T_C: document.getElementById('T').value };
  }));
  out.errors = errs;
  process.stdout.write(JSON.stringify(out));
  await b.close();
})();
'''


_SUB = str.maketrans('0123456789', '\u2080\u2081\u2082\u2083\u2084'
                                   '\u2085\u2086\u2087\u2088\u2089')


def _chem(name):
    """Same subscripting FigKit.chem does, so the check reads the glyphs."""
    return re.sub(r'([A-Za-z])(\d+)',
                  lambda m: m.group(1) + m.group(2).translate(_SUB), name)


def _pct(v):
    """Same three branches web/figures_thermo.js writes the header with."""
    if v >= 1:
        return '%.1f%%' % v
    if v >= 0.01:
        return '%.3f%%' % v
    return '%s%%' % ('%.2g' % v)


def _times(r):
    if r >= 1e4:
        return ('%.2g' % r).replace('e+', '\u00d710^').replace('e', '\u00d710^') + '\u00d7'
    if r >= 100:
        return '%d\u00d7' % round(r)
    if r >= 10:
        return '%.0f\u00d7' % r
    return '%.1f\u00d7' % r


def main():
    pw = os.environ.get('PW_DIR')
    if not pw or not os.path.isdir(os.path.join(pw, 'node_modules',
                                                'playwright')):
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
        print(r.stderr)
        return 1
    got = json.loads(r.stdout)

    part = P.particle(P.load())
    res = P.equilibrate(95.0, 600.0, part)

    T_C = float(got['T_C'])
    feed = {'CO2': 1.0, 'H2': 1.0}
    line = A.reduction_boundary(T_C + 273.15)
    eq = A.solve(feed, T_C + 273.15)
    mu = A.gas_oxygen_potential(feed, T_C + 273.15)
    y = A.equivalent_co2_fraction(mu, T_C + 273.15)
    ratio = y / line['y_CO2']

    checks = [
        ('no page errors', got['errors'] == [], got['errors']),
        ('depth figure drew', got['depth'], got['depth']),
        ('where figure drew', got['where'], got['where']),
        ('mu matches', '%.6f eV' % res['mu_eV'] in got['state'],
         '%.6f eV' % res['mu_eV']),
        ('first layer matches',
         '%.5f' % res['theta_first_layer'] in got['state'],
         '%.5f' % res['theta_first_layer']),
        ('ceiling matches',
         '%.1f' % res['point_defect_ceiling_umol_g'] in got['state'],
         '%.1f' % res['point_defect_ceiling_umol_g']),
        ('Debye length reported', '0.569 nm' in got['charge'], None),
        ('breaking barrier reported', '1.645 eV' in got['reach'], None),
        ('electroneutrality verdict shown',
         'cannot be separated' in got['charge'], None),

        ('boundary drew without a sweep', got['boundary_before_sweep'], None),
        ('conversion waited for the sweep',
         not got['conversion_before_sweep'], None),
        ('margin figure drew', got['margin'], None),
        ('boundary names the first phase',
         _chem(line['phase']) in got['boundary'], _chem(line['phase'])),
        ('boundary quotes the closed form',
         _pct(line['y_CO2'] * 100.0) in got['boundary'],
         _pct(line['y_CO2'] * 100.0)),
        ('boundary quotes the distance to it',
         _times(ratio) in got['boundary'], _times(ratio)),
        ('conversion quotes the panel solve',
         '%.1f%%' % eq['conversion_CO2_pct'] in got['conversion'],
         '%.1f%%' % eq['conversion_CO2_pct']),
    ]
    bad = 0
    for name, ok, detail in checks:
        print('%-32s %s%s' % (name, 'ok' if ok else 'FAILED',
                              '' if ok or detail is None else '  ' + str(detail)))
        bad += 0 if ok else 1
    return 1 if bad else 0


if __name__ == '__main__':
    sys.exit(main())
