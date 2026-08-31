"""Re-render the frozen manuscript plates from the current extract.

docs/figures/*.svg are the manuscript figures. They are frozen on purpose
- a gate holds them against a fresh render, so a change to the drawing
kit or to the numbers behind them shows up as a diff in a figure rather
than in a proof - and this is the one way they are allowed to move.

Run it after scripts/export_plates.py has refreshed data/figure_data.json,
and read the diff before committing it: every line that moved is a number
a plate is asserting.

    python3 scripts/export_plates.py
    python3 scripts/freeze_plates.py
"""

import json
import os
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
HARNESS = os.path.join(ROOT, 'tests', 'js', 'plates_harness.js')
OUT = os.path.join(ROOT, 'docs', 'figures')


def main():
    r = subprocess.run(['node', HARNESS, ROOT], capture_output=True,
                       text=True, timeout=600)
    if r.returncode != 0:
        sys.stderr.write(r.stderr)
        return 1
    doc = json.loads(r.stdout)
    if doc['errors']:
        for k, v in doc['errors'].items():
            sys.stderr.write('%s: %s\n' % (k, v))
        return 1
    moved = []
    for fid, svg in sorted(doc['svg'].items()):
        path = os.path.join(OUT, fid + '.svg')
        old = None
        if os.path.exists(path):
            with open(path) as fh:
                old = fh.read()
        body = svg.rstrip('\n') + '\n'
        if old != body:
            moved.append(fid)
            with open(path, 'w') as fh:
                fh.write(body)
    print('%d plates rendered, %d moved%s'
          % (len(doc['svg']), len(moved),
             (': ' + ', '.join(moved)) if moved else ''))
    return 0


if __name__ == '__main__':
    sys.exit(main())
