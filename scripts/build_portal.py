"""Publish the workspace selector as the site root.

The solver remains self-contained in ``ti_solver.html``. The root page offers
three entries into its equilibrium, defect-statistical-mechanics and
steady-state-redox workspaces.
"""

import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
WEB = ROOT / 'web'
DOCS = ROOT / 'docs'


def main():
    src = WEB / 'portal.html'
    if not src.exists():
        raise SystemExit('web/portal.html is missing')
    dst = DOCS / 'index.html'
    shutil.copyfile(src, dst)
    with dst.open() as fh:
        html = fh.read()
    for target in ('ti_solver.html#thermo', 'ti_solver.html#defect',
                   'ti_solver.html#redox'):
        assert f'href="{target}"' in html, f'portal does not link {target}'
    if not (DOCS / 'ti_solver.html').exists():
        raise SystemExit('docs/ti_solver.html is not built; run scripts/build_ti_solver.py')
    print(f'{dst.relative_to(ROOT)} written ({dst.stat().st_size / 1024:.0f} kB)')


if __name__ == '__main__':
    main()
