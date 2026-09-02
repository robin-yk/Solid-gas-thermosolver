"""web/vacancy.js against solidgas/vacancy_inventory.py, to the last place.

The bound is a dozen lines of arithmetic, which is exactly why it is worth
mirroring rather than precomputing: the page can then answer a BET area
the user actually measured. A mirror is only safe while something holds it
to the original, so this compares them on a grid and to the ulp.

There is no exponential anywhere in it, so with the multiplications
grouped the same way on both sides the two agree bit for bit and the
budget is zero, not a tolerance."""

import json
import pathlib
import shutil
import struct
import subprocess

import pytest

from solidgas import vacancy_inventory as V

ROOT = pathlib.Path(__file__).resolve().parents[1]
HARNESS = ROOT / 'tests' / 'js' / 'vacancy_harness.js'

DIAMETERS = [0.2, 0.9, 1.6, 5.0, 25.0]
BETS = [0.5, 1.569, 12.0]
INVENTORIES = [1.0, 6.0, 30.0, 95.0, 220.0, 4000.0]


def _ulps(a, b):
    if a == b:
        return 0
    ia, ib = (struct.unpack('<q', struct.pack('<d', v))[0] for v in (a, b))
    if (ia < 0) != (ib < 0):
        return 1 << 62
    return abs(ia - ib)


@pytest.fixture(scope='module')
def js():
    assert shutil.which('node') is not None, 'node is required for this gate'
    spec = {'diameters': DIAMETERS, 'bets': BETS, 'inventories': INVENTORIES}
    r = subprocess.run(['node', str(HARNESS), str(ROOT), json.dumps(spec)],
                       capture_output=True, text=True, timeout=120)
    assert r.returncode == 0, r.stderr[-2000:]
    return json.loads(r.stdout)


def test_the_geometry_is_bit_identical(js):
    g = V.load()
    assert _ulps(js['cell_area_cm2'], V.surface_cell_area_cm2(g)) == 0
    assert _ulps(js['areal_cm2'], V.bridging_areal_density_cm2(g)) == 0
    assert _ulps(js['density_g_cm3'], V.crystal_density_g_cm3(g)) == 0
    assert _ulps(js['site_density'], V.oxygen_site_density_umol_g(g)) == 0


def test_the_area_matches_for_every_diameter(js):
    g = V.load()
    for row in js['areas']:
        want = V.specific_surface_area(g, d_um=row['d_um'])
        assert _ulps(row['area'], want) == 0, row


def test_the_capacity_matches_for_every_size(js):
    g = V.load()
    for row in js['capacities']:
        kw = ({'bet_m2_g': row['bet']} if row.get('bet') is not None
              else {'d_um': row['d_um']})
        assert _ulps(row['full'], V.bridging_oxygen_capacity(g, **kw)) == 0, row
        assert _ulps(row['recon'], V.reconstruction_deficiency(g, **kw)) == 0, row


def test_every_field_of_the_bound_matches(js):
    g = V.load()
    worst = 0
    for row in js['bounds']:
        kw = ({'bet_m2_g': row['bet']} if row.get('bet') is not None
              else {'d_um': row['d_um']})
        want = V.surface_fraction_bound(g, row['inv'], **kw)
        for k, v in want.items():
            if not isinstance(v, float):
                continue
            worst = max(worst, _ulps(row['out'][k], v))
            assert _ulps(row['out'][k], v) == 0, (k, row)
    assert worst == 0, worst


def test_the_mirror_refuses_what_the_module_refuses(js):
    assert js['rejects_zero_inventory'] is True
    assert js['rejects_missing_area'] is True
