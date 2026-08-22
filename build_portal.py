"""Publish the workspace selector as the site root.

The solver remains self-contained in ``ti_solver.html``. The root page offers
two entries into its equilibrium and defect-statistical-mechanics workspaces.
"""

import os
import shutil

HERE = os.path.dirname(os.path.abspath(__file__))


def main():
    src = os.path.join(HERE, 'portal.html')
    if not os.path.exists(src):
        raise SystemExit('portal.html is missing')
    dst = os.path.join(HERE, 'index.html')
    shutil.copyfile(src, dst)
    with open(dst) as fh:
        html = fh.read()
    for target in ('ti_solver.html#thermo', 'ti_solver.html#defect'):
        assert f'href="{target}"' in html, f'portal does not link {target}'
    if not os.path.exists(os.path.join(HERE, 'ti_solver.html')):
        raise SystemExit('ti_solver.html is not built - run build_ti_solver.py')
    print(f'index.html written ({os.path.getsize(dst) / 1024:.0f} kB) = portal.html')


if __name__ == '__main__':
    main()
