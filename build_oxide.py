"""Assemble the single self-contained oxide solver page.

No external requests: the data, the engine, the UI and the prose all end up
inside one file. That is a hard requirement rather than a preference - the
template this borrows its structure from pulls two libraries off a CDN, and a
page that stops working when a CDN does is not a self-contained tool.

    python3 export_oxide.py    # oxide_data.json, straight from the package
    python3 oxide_reference.py # oxide_reference.json, the Python answer
    python3 build_oxide.py     # oxide_tool.html
"""

import os

HERE = os.path.dirname(os.path.abspath(__file__))
PARTS = {'__DATA__': 'oxide_data.json', '__ENGINE__': 'oxide.js',
         '__PAGE__': 'oxide_page.js', '__HOW__': 'oxide_how.html',
         '__DERIV__': 'oxide_deriv.html'}


def main():
    with open(os.path.join(HERE, 'oxide_template.html')) as fh:
        html = fh.read()
    for token, fname in PARTS.items():
        with open(os.path.join(HERE, fname)) as fh:
            body = fh.read()
        assert token in html, f'{token} missing from the template'
        html = html.replace(token, body)
    out = os.path.join(HERE, 'oxide_tool.html')
    with open(out, 'w') as fh:
        fh.write(html)
    for bad in ('http://', 'https://cdn', 'src="http'):
        n = html.count(bad)
        if bad == 'http://' and n:
            raise SystemExit(f'external reference found: {bad}')
    print(f'oxide_tool.html written ({os.path.getsize(out) / 1024:.0f} kB)')


if __name__ == '__main__':
    main()
