"""Python-to-JavaScript parity for the particle engine.

The browser runs a copy of the engine.  A copy that is allowed to drift
is worse than no copy, so every branch of the Python side is put through
the JavaScript side on cases generated here - one spec file, both
engines - and held to a few units in the last place.

Why not bit-identical.  The bracket algebra, the bisection counts and
the layer stack are deterministic and run in the same order on both
sides, so everything that depends on arithmetic alone does agree
exactly.  What does not is the logarithm: neither glibc nor V8 promises
a correctly rounded log, and they disagree in the last bit for some
arguments.  That difference propagates through a bisection and comes out
as one or two units in the last place.  Bounding it in ulps rather than
hiding it behind a relative tolerance keeps the gate able to see a real
divergence, which would be many orders larger.

The budget is not one ulp because one place amplifies it.  The surface
transition potential is kT times a difference of two logarithms of
similar size - phi(0.17) is -1.586 + 1.508 - so a last-bit disagreement
in either log lands on a much smaller number.  The worst case observed
across every branch here is nineteen units in the last place of a
quantity near 0.08 eV, which is under a femto-electronvolt.  Sixty-four
leaves room and is still twelve orders tighter than anything a genuine
divergence in the physics would produce.
"""

ULP_BUDGET = 64

import json
import math
import os
import shutil
import subprocess
import tempfile

import pytest

from solidgas import particle as P
from solidgas.particle import freeenergy as F
from solidgas.particle import material as M
from solidgas.particle import profile as PR
from solidgas.particle import spacecharge as SC
from solidgas.particle import transport as T

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
HARNESS = os.path.join(ROOT, 'tests', 'js', 'particle_harness.js')

SHELL_DEPTHS = [1.0, 2.0, 5.0, 10.0]
EXPOSURE_S = 300.0

# Every branch: dilute and saturated inversion, both compensation forms,
# three polaron ratios, both presets, the two pinned windows, the
# unpinned stretches between them, and a second particle size.
PHI_CASES = [{'theta': t, 'ratio': r}
             for r in (1.0, 2.0, 4.0, 6.0)
             for t in (1e-14, 1e-7, 1e-3, 0.02)
             if t < 1.0 / r] + [
    {'theta': 0.2, 'ratio': 4.0}, {'theta': 0.2499, 'ratio': 4.0},
    {'theta': 0.124, 'ratio': 4.0}, {'theta': 0.126, 'ratio': 4.0}]

INVERSE_CASES = [{'x': x, 'compensated': True, 'ratio': r}
                 for r in (1.0, 4.0, 6.0)
                 for x in (-200.0, -60.0, -1.9459101490553132, -1.0,
                           0.0, 5.0, 30.0)] + [
    {'x': x, 'compensated': False, 'ratio': 4.0}
    for x in (-40.0, 0.0, 40.0)]

EQ_CASES = [
    {'vo': v, 'T_C': t, 'd_um': d, 'preset': pr, 'ratio': ra,
     'compensated': co, 'reconstruction': rc, 'shear_limit': sh}
    for (v, t, d, pr, ra, co, rc, sh) in [
        (95.0, 600.0, 0.9, 'sx', 4.0, True, True, True),
        (1.0, 600.0, 0.9, 'sx', 4.0, True, True, True),
        (35.0, 600.0, 0.9, 'sx', 4.0, True, True, True),      # pinned by 1x2
        (400.0, 600.0, 0.9, 'sx', 4.0, True, True, True),     # pinned by shear
        (95.0, 400.0, 0.9, 'sx', 4.0, True, True, True),
        (95.0, 800.0, 0.9, 'sx', 4.0, True, True, True),
        (95.0, 600.0, 0.9, 'gga', 4.0, True, True, True),
        (95.0, 600.0, 0.15, 'sx', 4.0, True, True, True),     # more surface
        (95.0, 600.0, 0.9, 'sx', 1.0, True, True, True),
        (95.0, 600.0, 0.9, 'sx', 4.0, False, True, True),     # uncompensated
        (95.0, 600.0, 0.9, 'sx', 4.0, True, False, True),     # no 1x2
        (400.0, 600.0, 0.9, 'sx', 4.0, True, True, False),    # no shear cap
    ]]

CRANK = [0.2, 0.5, 0.9, 0.99]
DEPTH_LAYERS = [0, 1, 2, 3, 8]
PRESETS = ['sx', 'gga']


@pytest.fixture(scope='module')
def js():
    assert shutil.which('node') is not None, \
        'node is required to run the parity gate'
    spec = {'phi': PHI_CASES, 'inverse': INVERSE_CASES,
            'presets': PRESETS, 'depth_layers': DEPTH_LAYERS,
            'equilibrium': EQ_CASES, 'crank': CRANK,
            'shell_depths_nm': SHELL_DEPTHS, 'exposure_s': EXPOSURE_S}
    with tempfile.NamedTemporaryFile('w', suffix='.json', delete=False) as fh:
        json.dump(spec, fh)
        path = fh.name
    try:
        r = subprocess.run(['node', HARNESS, ROOT, path],
                           capture_output=True, text=True)
        assert r.returncode == 0, r.stderr
        return json.loads(r.stdout)
    finally:
        os.unlink(path)


def _ulps(a, b):
    """Distance in representable doubles.  Zero means bit-identical."""
    if a == b:
        return 0
    if math.isnan(a) or math.isnan(b):
        return float('inf')
    if math.isinf(a) or math.isinf(b):
        return 0 if a == b else float('inf')
    n = 0
    lo, hi = (a, b) if a < b else (b, a)
    while lo < hi and n <= ULP_BUDGET:
        lo = math.nextafter(lo, math.inf)
        n += 1
    return n


def _same(a, b, what):
    n = _ulps(a, b)
    assert n <= ULP_BUDGET, (
        '%s: python %.17g, javascript %.17g (%s ulp apart)'
        % (what, a, b, n))


def test_the_free_energy_agrees_to_the_last_bits(js):
    for c, got in zip(PHI_CASES, js['phi']):
        _same(F.phi(c['theta'], c['ratio']), got,
              'phi(%g, %g)' % (c['theta'], c['ratio']))


def test_the_inversion_agrees_on_both_branches(js):
    for c, got in zip(INVERSE_CASES, js['inverse']):
        _same(F.theta_of_phi(c['x'], c['compensated'], c['ratio']), got,
              'theta_of_phi(%g, %s, %g)' % (c['x'], c['compensated'],
                                            c['ratio']))


def test_the_decay_length_agrees_to_the_last_bits(js):
    p = P.load()
    for key in PRESETS:
        s = M.segregation(p, key)
        want = js['segregation'][key]
        _same(s['xi_nm'], want['xi_nm'], '%s xi' % key)
        _same(s['dE_seg_eV'], want['dE_seg_eV'], '%s dE' % key)
        for n, e in zip(DEPTH_LAYERS, want['E']):
            _same(M.formation_energy(n * s['d110_nm'], s), e,
                  '%s E(%d layers)' % (key, n))


def test_every_equilibrium_branch_agrees_to_the_last_bits(js):
    p = P.load()
    for c, got in zip(EQ_CASES, js['equilibrium']):
        part = P.particle(p, d_um=c['d_um'], preset=c['preset'])
        r = P.equilibrate(c['vo'], c['T_C'], part, ratio=c['ratio'],
                          compensated=c['compensated'],
                          reconstruction=c['reconstruction'],
                          shear_limit=c['shear_limit'])
        tag = '%g umol %g C %g um %s ratio %g' % (
            c['vo'], c['T_C'], c['d_um'], c['preset'], c['ratio'])
        for k in ('mu_eV', 'theta_interior_mean', 'theta_bulk_plateau',
                  'theta_first_layer', 'umol_surface', 'umol_interior',
                  'umol_extended'):
            _same(r[k], got[k], '%s: %s' % (tag, k))
        _same(r['surface']['theta'], got['theta_surface'],
              '%s: theta_surface' % tag)
        _same(r['surface']['reconstructed_fraction'],
              got['reconstructed_fraction'], '%s: reconstructed' % tag)
        assert r['surface']['phase'] == got['surface_phase'], tag
        assert r['limited_by'] == got['limited_by'], tag
        _same(part['N_surface_umol_g'], got['N_surface_umol_g'],
              '%s: N_surface' % tag)
        _same(part['radius_nm'], got['radius_nm'], '%s: radius' % tag)
        if r['point_defect_ceiling_umol_g'] is not None:
            _same(r['point_defect_ceiling_umol_g'],
                  got['point_defect_ceiling_umol_g'], '%s: ceiling' % tag)


def test_the_downstream_readings_agree_to_the_last_bits(js):
    p = P.load()
    for c, shells, sc, tr in zip(EQ_CASES, js['shells'], js['spacecharge'],
                                 js['transport']):
        part = P.particle(p, d_um=c['d_um'], preset=c['preset'])
        r = P.equilibrate(c['vo'], c['T_C'], part, ratio=c['ratio'],
                          compensated=c['compensated'],
                          reconstruction=c['reconstruction'],
                          shear_limit=c['shear_limit'])
        for row, want in zip(PR.shell_split(r, part, SHELL_DEPTHS), shells):
            _same(row['inventory_fraction'], want,
                  'shell %g nm' % row['depth_nm'])
        _same(SC.assess(r, part, p)['ratio_debye_over_xi'], sc, 'debye/xi')
        _same(T.verdict(c['T_C'], part, EXPOSURE_S, p=p)[
            'E_m_that_would_break_it_eV'], tr, 'breaking barrier')


def test_the_crank_series_agrees_to_the_last_bits(js):
    for f, got in zip(CRANK, js['crank']):
        _same(T.tau_for_uptake(f), got, 'tau for %g' % f)
