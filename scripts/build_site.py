"""Build tiox.html, the oxygen-potential dashboard, from the solver's JSON output.

The site root is the portal (`portal.html`, copied to `index.html` by
build_portal.py); this page sits under it.

No backend and no build toolchain: the numbers are computed here, embedded in
the page as JSON, and drawn client-side as inline SVG. GitHub Pages serves the
single file as-is.

    python3 scripts/run_all.py
    python3 scripts/build_site.py
"""

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
WEB = ROOT / 'web'
DATA = ROOT / 'data'
DOCS = ROOT / 'docs'
sys.path.insert(0, str(ROOT))

import numpy as np

from solidgas import (ACTIVE_TI_PHASES, FORMULAS, SHOMATE, VALID_RANGE,
                      R_KJ, mu0)

def main():
    """Assemble index.html from the template, the solver, and the data."""
    with (WEB / 'site_template.html').open() as fh:
        html = fh.read()
    with (WEB / 'solver.js').open() as fh:
        solver = fh.read()
    with (WEB / 'page.js').open() as fh:
        page = fh.read()
    with (DATA / 'thermo_data.json').open() as fh:
        data = fh.read()
    # The page must not carry its own copy of the answer. It gets the same file
    # the table is read from, verbatim, and a test compares the two byte for byte.
    with (DATA / 'results.json').open() as fh:
        results = fh.read()
    html = (html.replace('__SOLVER__', solver).replace('__PAGE__', page)
                .replace('__THERMO__', data).replace('__RESULTS__', results))
    out = DOCS / 'tiox.html'
    with out.open('w') as fh:
        fh.write(html)
    print(f'{out.relative_to(ROOT)} written ({out.stat().st_size / 1024:.0f} kB)')


if __name__ == '__main__':
    main()
