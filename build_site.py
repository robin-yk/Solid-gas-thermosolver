"""Build tiox.html, the oxygen-potential dashboard, from the solver's JSON output.

The site root is the portal (`portal.html`, copied to `index.html` by
build_portal.py); this page sits under it.

No backend and no build toolchain: the numbers are computed here, embedded in
the page as JSON, and drawn client-side as inline SVG. GitHub Pages serves the
single file as-is.

    python3 run_all.py              # writes results.json - every reported number
    python3 build_site.py           # writes tiox.html, embedding that same file
"""

import json
import os

import numpy as np

from solidgas import (ACTIVE_TI_PHASES, FORMULAS, SHOMATE, VALID_RANGE,
                      R_KJ, mu0)

HERE = os.path.dirname(os.path.abspath(__file__))


def main():
    """Assemble index.html from the template, the solver, and the data."""
    with open(os.path.join(HERE, 'site_template.html')) as fh:
        html = fh.read()
    with open(os.path.join(HERE, 'solver.js')) as fh:
        solver = fh.read()
    with open(os.path.join(HERE, 'page.js')) as fh:
        page = fh.read()
    with open(os.path.join(HERE, 'thermo_data.json')) as fh:
        data = fh.read()
    # The page must not carry its own copy of the answer. It gets the same file
    # the table is read from, verbatim, and a test compares the two byte for byte.
    with open(os.path.join(HERE, 'results.json')) as fh:
        results = fh.read()
    html = (html.replace('__SOLVER__', solver).replace('__PAGE__', page)
                .replace('__THERMO__', data).replace('__RESULTS__', results))
    out = os.path.join(HERE, 'tiox.html')
    with open(out, 'w') as fh:
        fh.write(html)
    print(f'tiox.html written ({os.path.getsize(out) / 1024:.0f} kB)')


if __name__ == '__main__':
    main()
