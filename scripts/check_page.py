"""Load the built page in a real browser and read the particle card back.

The test suite cannot assume a browser, and a gate that skips itself is
worse than a script: it reports green while checking nothing.  So this
lives here and is run by hand (or by CI where a browser exists).

    PW_DIR=/path/with/node_modules/playwright python3 scripts/check_page.py

Exits non-zero on any page error, any figure that did not draw, or any
reading that disagrees with a fresh solve of the Python engine.
"""

import json
import os
import subprocess
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

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
  out.errors = errs;
  process.stdout.write(JSON.stringify(out));
  await b.close();
})();
'''


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
    ]
    bad = 0
    for name, ok, detail in checks:
        print('%-32s %s%s' % (name, 'ok' if ok else 'FAILED',
                              '' if ok or detail is None else '  ' + str(detail)))
        bad += 0 if ok else 1
    return 1 if bad else 0


if __name__ == '__main__':
    sys.exit(main())
