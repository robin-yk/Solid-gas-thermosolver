"""Publish the solver as the site root: index.html = ti_solver.html.

The root URL used to be a portal over two tools; it now opens the active-set
solver directly, and the legacy pages (tiox.html, oxide_tool.html) are
reachable from the solver's Legacy tab, labelled with why they are retired.
portal.html is kept in the tree for history but is no longer built into the
root.
"""

import os
import shutil

HERE = os.path.dirname(os.path.abspath(__file__))


def main():
    src = os.path.join(HERE, 'ti_solver.html')
    if not os.path.exists(src):
        raise SystemExit('ti_solver.html is not built - run build_ti_solver.py')
    dst = os.path.join(HERE, 'index.html')
    shutil.copyfile(src, dst)
    with open(dst) as fh:
        html = fh.read()
    # The Legacy tab links these; the root must not link into a 404.
    for target in ('tiox.html', 'oxide_tool.html'):
        assert f'href="{target}"' in html, f'Legacy tab does not link {target}'
        if not os.path.exists(os.path.join(HERE, target)):
            raise SystemExit(f'index links {target} but it is not built')
    print(f'index.html written ({os.path.getsize(dst) / 1024:.0f} kB) '
          f'= ti_solver.html; legacy pages linked from the Legacy tab')


if __name__ == '__main__':
    main()
