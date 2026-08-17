"""Copy the portal to index.html, the site root.

Kept as a copy rather than a redirect so the root is a real page: GitHub Pages
serves index.html directly, and a redirect would put a hop in front of every
visit for no gain.
"""

import os
import shutil

HERE = os.path.dirname(os.path.abspath(__file__))


def main():
    src = os.path.join(HERE, 'portal.html')
    dst = os.path.join(HERE, 'index.html')
    shutil.copyfile(src, dst)
    with open(dst) as fh:
        html = fh.read()
    for target in ('ti_solver.html', 'tiox.html'):
        assert f'href="{target}"' in html, f'portal does not link {target}'
        if not os.path.exists(os.path.join(HERE, target)):
            raise SystemExit(f'portal links {target} but it is not built')
    print(f'index.html written ({os.path.getsize(dst) / 1024:.0f} kB), links 2 tools')


if __name__ == '__main__':
    main()
