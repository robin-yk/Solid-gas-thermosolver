"""The particle engine has to actually be on the page.

An engine nobody can reach is a library, not a workspace.  These gates
hold three things: that the built page carries the engine, its figures
and the card that mounts them; that the typeset formulation states the
physics the engine implements rather than the one it replaced; and that
every number the formulation quotes as prose is one the engine returns.

What is deliberately NOT here is loading the page in a browser and
reading the card back.  That check exists - scripts/check_page.py runs
it - but it needs a browser this suite cannot assume, and a test that
quietly skips itself on a machine without one is worse than a script
somebody has to run: the suite would report green while checking
nothing.  The static gates below cover the wiring; the script covers the
running page.
"""

import os
import re

import pytest

from solidgas import particle as P
from solidgas.particle import material as M

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PAGE = os.path.join(ROOT, 'docs', 'ti_solver.html')
TEMPLATE = os.path.join(ROOT, 'web', 'ti_solver_template.html')
EQS = os.path.join(ROOT, 'web', 'generated', 'equations_statmech.html')
# The rendered sheet is glyph outlines, so its prose is not greppable.
# The authored source is, and it is the thing an edit would change; a
# freshness gate in test_statmech_port.py already holds the two together.
REG = os.path.join(ROOT, 'scripts', 'equations_registry_statmech.py')


def _read(path):
    with open(path) as fh:
        return fh.read()


@pytest.fixture(scope='module')
def page():
    return _read(PAGE)


@pytest.fixture(scope='module')
def template():
    return _read(TEMPLATE)


def test_the_page_carries_the_engine_the_card_needs(page):
    """All three files, inlined, and no leftover build token."""
    for marker in ('root.Particle = api',
                   'root.ParticleFigures = api',
                   'The depth-resolved workspace card'):
        assert marker in page, marker
    for token in ('__PX_ENGINE__', '__PX_FIGS__', '__PX_UI__'):
        assert token not in page, '%s was not substituted' % token


def test_the_card_mounts_where_the_page_says_it_does(template):
    """Every host the UI writes into exists in the defect workspace, and
    every input it listens to is an input that workspace has."""
    ui = _read(os.path.join(ROOT, 'web', 'particle_ui.js'))
    for host in re.findall(r"\$\('([a-zA-Z]+)'\)", ui):
        if host.startswith('px') or host.startswith('fig'):
            assert 'id="%s"' % host in template, 'no host %s on the page' % host
    for inp in ('smVO', 'smT', 'smGeoVal', 'smGeoMode', 'smPreset', 'smExp'):
        assert 'id="%s"' % inp in template, 'no input %s to listen to' % inp

    # both figures live inside the defect workspace, not another one
    defect = template[template.index('id="ws-defect"'):]
    defect = defect[:defect.index('id="ws-redox"')]
    assert 'id="figDepth"' in defect and 'id="figWhere"' in defect


def test_the_sheet_is_actually_rendered():
    eqs = _read(EQS)
    assert len(eqs) > 100_000
    assert eqs.count('<svg') >= 21


def test_the_formulation_states_the_compensated_free_energy():
    """The sheet used to print the plain Fermi-Dirac isotherm, which is a
    neutral vacancy.  If it goes back to that, this fails."""
    eqs = _read(REG)
    for phrase in ('The free energy of a compensated vacancy',
                   'two electrons behind',
                   'runs out of',
                   'NEUTRAL vacancy',
                   'the cube of the right one'):
        assert phrase in eqs, phrase
    for gone in ('analytic site-exclusion',):
        assert gone not in eqs, 'the pre-compensation wording is back: ' + gone


def test_the_formulation_states_the_derived_decay_length():
    eqs = _read(REG)
    assert 'no profile shape is assumed' in eqs
    p = P.load()
    for key, want in (('sx', 0.543), ('gga', 0.725)):
        assert M.segregation(p, key)['xi_nm'] == pytest.approx(want, abs=5e-4)
        assert ('%.3f' % want) in eqs, 'the sheet does not quote xi for ' + key


def test_the_formulation_discloses_the_open_items():
    """Two things the sheet must not leave out: that the full space-charge
    solve is missing, and that the kinetic layer's equilibrium measure no
    longer agrees with the thermodynamic one."""
    eqs = _read(REG)
    assert 'Poisson-Boltzmann solve is not implemented' in eqs
    assert 'OPEN ITEM' in eqs
    assert 'no longer agree' in eqs


def test_every_number_the_card_and_sheet_quote_is_the_engine_s():
    """Numbers printed as prose in the formulation, checked against a
    fresh solve rather than trusted."""
    p = P.load()
    part = P.particle(p)
    r = P.equilibrate(95.0, 600.0, part)
    eqs = _read(REG)
    claims = [
        ('112.4', r['point_defect_ceiling_umol_g'], '%.1f'),
        ('13.6', part['N_surface_umol_g'], '%.1f'),
        ('0.57$ nm', 0.0, None),      # against xi, checked in test_particle
        ('0.54$ nm', 0.0, None),
        ('0.5 %', 0.0, None),         # the cost of smearing the layer stack
        ('85 %', 100.0 * 95.0 / r['point_defect_ceiling_umol_g'], '%.0f %%'),
        ('0.995', 0.0, None),
    ]
    for text, value, fmt in claims:
        assert text in eqs, 'the sheet no longer quotes ' + text
        if fmt:
            assert (fmt % value) == text, (
                'the sheet says %r; the engine says %r' % (text, fmt % value))
