"""Assemble the single self-contained oxide solver page.

No external requests: the data, the engine, the UI and the prose all end up
inside one file. That is a hard requirement rather than a preference - the
template this borrows its structure from pulls two libraries off a CDN, and a
page that stops working when a CDN does is not a self-contained tool.

    python3 scripts/export_oxide.py
    python3 scripts/oxide_reference.py
    python3 scripts/build_oxide.py
"""

from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
WEB = ROOT / 'web'
DATA = ROOT / 'data'
DOCS = ROOT / 'docs'
PARTS = {'__DATA__': DATA / 'oxide_data.json',
         '__ENGINE__': WEB / 'oxide.js',
         '__PAGE__': WEB / 'oxide_page.js',
         '__HOW__': WEB / 'oxide_how.html',
         '__DERIV__': WEB / 'oxide_deriv.html'}


def main():
    with (WEB / 'oxide_template.html').open() as fh:
        html = fh.read()
    for token, path in PARTS.items():
        with path.open() as fh:
            body = fh.read()
        assert token in html, f'{token} missing from the template'
        html = html.replace(token, body)
    out = DOCS / 'oxide_tool.html'
    with out.open('w') as fh:
        fh.write(html)
    for bad in ('http://', 'https://cdn', 'src="http'):
        n = html.count(bad)
        if bad == 'http://' and n:
            raise SystemExit(f'external reference found: {bad}')
    print(f'{out.relative_to(ROOT)} written ({out.stat().st_size / 1024:.0f} kB)')


if __name__ == '__main__':
    main()
