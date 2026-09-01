"""Assemble the site: one self-contained page, no external requests.

    python3 scripts/reproduce_paper.py    # produces the numbers
    python3 scripts/build_site.py         # inlines them into docs/index.html

The page displays and does not calculate, so the build is a substitution:
the JSON that `reproduce_paper.py` wrote, the drawing kit, the two figure
modules and the page controller all go inside one file. Nothing is
fetched at load time and no solver ships to the browser.
"""

import json
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
WEB = os.path.join(ROOT, 'web')
OUT = os.path.join(ROOT, 'docs', 'index.html')
DATA = os.path.join(ROOT, 'paper_outputs', 'site_data.json')

PARTS = [
    ('/*CSS*/', os.path.join(WEB, 'site.css')),
    ('/*FIGKIT*/', os.path.join(WEB, 'figkit.js')),
    ('/*THERMO*/', os.path.join(WEB, 'figures_thermo.js')),
    ('/*SURFACE*/', os.path.join(WEB, 'figures_surface.js')),
    ('/*UI*/', os.path.join(WEB, 'site_ui.js')),
]


def read(path):
    with open(path) as fh:
        return fh.read()


def build():
    if not os.path.exists(DATA):
        raise SystemExit('run scripts/reproduce_paper.py first')
    html = read(os.path.join(WEB, 'template.html'))
    with open(DATA) as fh:
        data = json.load(fh)
    html = html.replace('/*DATA*/', json.dumps(data, separators=(',', ':')))
    for token, path in PARTS:
        html = html.replace(token, read(path))
    for token, _ in PARTS + [('/*DATA*/', None)]:
        if token in html:
            raise SystemExit('unsubstituted token left in the page: ' + token)
    if '<script src=' in html or 'href="http' in html.replace(
            'href="https://github.com', ''):
        raise SystemExit('the page must not reach outside itself')
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, 'w') as fh:
        fh.write(html)
    return OUT, len(html)


if __name__ == '__main__':
    path, size = build()
    print('%s written (%d kB)' % (os.path.relpath(path, ROOT), size // 1024))
    sys.exit(0)
