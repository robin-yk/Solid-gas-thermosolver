"""Render the equation registry as an inline-SVG HTML fragment for the page.

Each equation is typeset by matplotlib mathtext (Computer Modern) and saved
as a standalone SVG with the ink colour replaced by currentColor, so the
glyphs follow the page theme. Element ids are prefixed per equation so the
inlined SVGs cannot collide. Output: equations_method.html, inlined into
ti_solver.html by build_ti_solver.py.
"""

import io
import re

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from equations_registry import R

INK = '#0b0b0b'

plt.rcParams.update({
    'mathtext.fontset': 'cm',
    'font.family': 'DejaVu Sans',
    'text.color': INK,
    'svg.fonttype': 'path',
})


GLYPHS = {}          # glyph id -> <path .../> definition, shared page-wide


def render_svg(tex, fontsize, prefix):
    """One equation as SVG. Glyph outlines are hoisted into the shared GLYPHS
    pool (matplotlib's glyph ids are deterministic per glyph, so identical
    ids carry identical outlines); only the structural ids are prefixed."""
    fig = plt.figure(figsize=(0.1, 0.1))
    fig.patch.set_alpha(0.0)
    t = fig.text(0, 0, tex, fontsize=fontsize)
    fig.canvas.draw()
    bbox = t.get_window_extent()
    pad = 2.0
    w = (bbox.width + 2 * pad) / fig.dpi
    h = (bbox.height + 2 * pad) / fig.dpi
    fig.set_size_inches(w, h)
    t.set_position((pad / (w * fig.dpi), pad / (h * fig.dpi)))
    buf = io.StringIO()
    fig.savefig(buf, format='svg', transparent=True, bbox_inches=None)
    plt.close(fig)
    svg = buf.getvalue()
    svg = svg[svg.index('<svg'):]
    svg = re.sub(r'<!--.*?-->', '', svg, flags=re.S)
    svg = re.sub(r'<metadata>.*?</metadata>', '', svg, flags=re.S)
    svg = svg.replace(INK, 'currentColor')

    def hoist(mdefs):
        for pm in re.finditer(r'<path id="([^"]+)"[^>]*/>', mdefs.group(1)):
            GLYPHS.setdefault(pm.group(1), pm.group(0))
        return ''
    svg = re.sub(r'<defs>(.*?)</defs>', hoist, svg, flags=re.S)

    svg = re.sub(r'id="((?:figure|axes|patch|text)[^"]*)"',
                 f'id="{prefix}-\\1"', svg)
    # Intrinsic size in px (SVG pt -> CSS px at the same number keeps the
    # 13 pt math close to the page's 15 px body text).
    m = re.search(r'width="([\d.]+)pt" height="([\d.]+)pt"', svg)
    if m:
        wpt, hpt = float(m.group(1)), float(m.group(2))
        svg = svg.replace(m.group(0),
                          f'width="{wpt:.0f}" height="{hpt:.0f}" '
                          f'style="max-width:100%;height:auto"', 1)
    return svg


def build_fragment():
    out = []
    n = 0
    for row in R:
        kind = row[0]
        if kind in ('title', 'sub', 'gap'):
            continue                     # the page supplies its own heading
        n += 1
        if kind == 'head':
            out.append(f'<h3>{row[1]}</h3>')
        elif kind == 'note':
            out.append('<div class="eqnote">'
                       + render_svg(row[1], 10.5, f'eqn{n}') + '</div>')
        elif kind == 'eq':
            out.append('<div class="eqrow">'
                       + render_svg(row[1], 13, f'eqn{n}')
                       + f'<span class="eqno">({row[2]})</span></div>')
        elif kind == 'boxeq':
            out.append('<div class="eqbox"><div class="eqrow">'
                       + render_svg(row[1], 13, f'eqn{n}')
                       + f'<span class="eqno">({row[2]})</span></div></div>')
    pool = ('<svg width="0" height="0" style="position:absolute" '
            'aria-hidden="true"><defs>' + ''.join(GLYPHS.values())
            + '</defs></svg>')
    return pool + '\n' + '\n'.join(out)


def main():
    frag = build_fragment()
    with open('equations_method.html', 'w') as fh:
        fh.write(frag)
    print(f'equations_method.html written ({len(frag) / 1024:.0f} kB)')


if __name__ == '__main__':
    main()
