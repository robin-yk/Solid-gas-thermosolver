"""Assemble the self-contained manuscript figure proof sheet.

Everything ends up inside one file: the frozen figure data, the resolved
registry, and the plate kit. No external request, which the gate checks.
Caption slots are resolved here, so a caption that names a number the
extractor did not produce fails the build rather than the reader.

    python3 scripts/export_plates.py
    python3 scripts/build_plates.py
"""

import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
WEB = ROOT / 'web'
DATA = ROOT / 'data'
DOCS = ROOT / 'docs'
sys.path.insert(0, str(Path(__file__).resolve().parent))

import figure_registry as REG                 # noqa: E402

SLOT = re.compile(r'\{\{([^}]+)\}\}')
BOLD = re.compile(r'\*\*(.+?)\*\*', re.S)


def esc(s):
    return (s.replace('&', '&amp;').replace('<', '&lt;')
            .replace('>', '&gt;'))


def caption_html(text, slots, fig_id):
    missing = [k for k in SLOT.findall(text) if k not in slots]
    if missing:
        raise SystemExit(f'{fig_id}: caption names slots the extractor did '
                         f'not produce: {missing}')
    out = esc(text)
    out = SLOT.sub(lambda m: esc(slots[m.group(1)]), out)
    out = BOLD.sub(r'<b>\1</b>', out)
    return out


def main():
    doc = json.loads((DATA / 'figure_data.json').read_text())
    slots = doc['slots']

    figures = []
    used_symbols, used_methods = set(), set()
    for f in REG.FIGURES:
        figures.append({
            'id': f['id'], 'track': f['track'], 'label': f['label'],
            'width': f['width'], 'section': esc(f['section']),
            'claim': esc(f['claim']),
            'caption_html': caption_html(f['caption'], slots, f['id']),
        })
        used_symbols.update(f['symbols'])
        used_methods.update(f['methods'])

    unknown = ([s for s in used_symbols if s not in REG.SYMBOLS]
               + [m for m in used_methods if m not in REG.METHODS])
    if unknown:
        raise SystemExit(f'undefined in the registry vocabulary: {unknown}')

    registry = {
        'spec': REG.SPEC,
        'roles': [{'key': k, 'hue': h, 'where': esc(w)}
                  for k, h, w in REG.ROLES],
        'figures': figures,
        'symbols': sorted((esc(REG.SYMBOLS[s][0]), esc(REG.SYMBOLS[s][1]))
                          for s in used_symbols),
        'methods': sorted((esc(m), esc(REG.METHODS[m]))
                          for m in used_methods),
    }

    html = (WEB / 'figures_template.html').read_text()
    for token, payload in (
            ('__FIG_DATA__', json.dumps(doc, separators=(',', ':'))),
            ('__FIG_REGISTRY__', json.dumps(registry,
                                            separators=(',', ':'))),
            ('__PLATES_JS__', (WEB / 'plates.js').read_text())):
        if token not in html:
            raise SystemExit(f'template lost its {token} slot')
        html = html.replace(token, payload)

    # XML namespace URIs are identifiers, not requests
    scan = html.replace('http://www.w3.org/2000/svg', '') \
               .replace('http://www.w3.org/1999/xlink', '')
    for bad in ('http://', 'https://'):
        if bad in scan:
            ctx = scan[max(0, scan.find(bad) - 60):scan.find(bad) + 60]
            raise SystemExit(f'page must be self-contained, found {bad}: '
                             f'...{ctx}...')

    out = DOCS / 'figures.html'
    out.write_text(html)
    print(f'figures.html written ({len(html) // 1024} kB, '
          f'{len(figures)} figure entries, {len(slots)} slots resolved)')

    # the same sheet without the outer document tags, for publishing as a
    # hosted page: the host supplies <html>/<head>/<body>
    body = html.split('<body>', 1)[1].rsplit('</body>', 1)[0]
    head = html.split('<head>', 1)[1].split('</head>', 1)[0]
    title = re.search(r'<title>.*?</title>', head, re.S).group(0)
    style = re.search(r'<style>.*?</style>', head, re.S).group(0)
    frag = DOCS / 'figures_fragment.html'
    frag.write_text(title + '\n' + style + '\n' + body)
    print(f'figures_fragment.html written ({frag.stat().st_size // 1024} kB)')


if __name__ == '__main__':
    main()
