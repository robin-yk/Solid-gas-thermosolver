"""Gates for the fixed-inventory equilibrium TOF range."""
import csv
import filecmp
import itertools
import math
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

HERE = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(HERE))

from tofrange import particle as pt                       # noqa: E402
from tofrange import scenarios as sc                      # noqa: E402
from tofrange.engine import Model                         # noqa: E402

T = 873.15
KT = pt.KB * T


def registry():
    with open(HERE / 'specification' / 'parameter_registry.csv') as f:
        return {r['id']: r['value'] for r in csv.DictReader(f)}


def test_constants_match_the_registry():
    r = registry()
    assert (pt.A_NM, pt.C_NM, pt.U) == (float(r['rutile_a_nm']), float(r['rutile_c_nm']), float(r['rutile_u']))
    assert (pt.RHO, pt.M) == (float(r['rutile_density_g_cm3']), float(r['rutile_molar_mass_g_mol']))
    assert (pt.NA, pt.KB) == (float(r['avogadro']), float(r['boltzmann_eV_K']))
    assert T == float(r['equilibrium_T_K'])
    assert sc.CONT_XI == (0.25, 0.5, 1.0) and sc.A_LI == float(r['li_continuum_A_eV'])
    assert sc.DIAMETERS == (900.0,) + tuple(float(d) for d in r['diameter_sensitivity_nm'].split(';'))[::-1]


def test_capacities_partition_every_site_once():
    assert pt.bridging_capacity(900) == pytest.approx(13.558, abs=5e-4)     # manuscript Note 2a: ~13.6
    for nl, emap in sc.discrete_maps().values():
        m = sc.build_discrete(nl, emap, 'LOCAL', 900.0)
        o = sum(m.C[i] for i, t in enumerate(m.tags) if t != 'Ti')
        ti = sum(m.C[i] for i, t in enumerate(m.tags) if t == 'Ti')
        assert o == pytest.approx(pt.O_TOTAL, rel=1e-14) and ti == pytest.approx(pt.TI_TOTAL, rel=1e-14)
    m = sc.build_continuum(0.5, 'LOCAL', 900.0)
    assert m.C[0::2].sum() == pytest.approx(pt.O_TOTAL, rel=1e-13)


def test_neutral_solution_matches_an_independent_bisection():
    nl, emap = sc.discrete_maps()['PAB']
    m = sc.build_discrete(nl, emap, 'NEUTRAL', 900.0)
    E, C = m.E[:, 1], m.C
    lo, hi = -5.0, 5.0
    for _ in range(200):
        mu = 0.5 * (lo + hi)
        n = float(np.sum(C / (1 + np.exp((E - mu) / KT))))
        lo, hi = (lo, mu) if n > 94.0 else (mu, hi)
    sol = m.solve(94.0, T)
    assert sol.mu_eV == pytest.approx(mu, abs=1e-12)
    assert sol.x[:, 1] == pytest.approx(C / (1 + np.exp((E - mu) / KT)), rel=1e-9)


@pytest.mark.parametrize('name', ['PAB', 'HAM', 'LI_SBR1', 'LI_L2'])
def test_local_solution_satisfies_the_stationarity_condition(name):
    """eps + kT ln theta/(1-theta) + 2 kT ln y/(1-y) = mu at every O site."""
    nl, emap = sc.discrete_maps()[name]
    m = sc.build_discrete(nl, emap, 'LOCAL', 900.0)
    sol = m.solve(193.1, T)
    frac = sol.x[:, 1] / m.C
    y = {m.region[i]: frac[i] for i, t in enumerate(m.tags) if t == 'Ti'}
    for i, t in enumerate(m.tags):
        if t == 'Ti':
            continue
        th, yy = frac[i], y[m.region[i]]
        lhs = m.E[i, 1] + KT * math.log(th / (1 - th)) + 2 * KT * math.log(yy / (1 - yy))
        assert lhs == pytest.approx(sol.mu_eV, abs=1e-9)
    assert abs(sol.inventory_residual) < 1e-9 and sol.charge_residual < 1e-8


def test_note_2b_is_reproduced():
    """Manuscript SI Note 2b: x(R), R600 interior and top-2-nm share, xi scan."""
    m = sc.build_continuum(0.5, 'LOCAL', 900.0)
    x = [sc.continuum_surface_x(m.solve(n, T).mu_eV, T, True) for n in (13.45, 36.85, 94.0, 193.1, 781.2)]
    assert [round(100 * v, 2) for v in x[:3]] == [10.13, 17.59, 22.49]
    assert [round(100 * v, 1) for v in x[3:]] == [24.1, 24.9]
    sol = m.solve(94.0, T)
    assert round(sol.x[-2, 1] / m.C[-2], 5) == 0.00338
    z = sc.continuum_edges(450.0)
    top = lambda s: s.x[0::2, 1][z[1:] <= 2.0].sum() / 94.0   # noqa: E731
    assert round(top(sol), 3) == 0.111
    for xi, ms in ((0.25, 0.064), (1.0, 0.193)):
        assert round(top(sc.build_continuum(xi, 'LOCAL', 900.0).solve(94.0, T)), 3) == ms
    ratio = [x[4] * (1 - x[4]) ** k / (x[2] * (1 - x[2]) ** k) for k in (2, 4, 8)]
    assert (round(min(ratio), 2), round(max(ratio), 2)) == (0.86, 1.04)


def _chain(n, periodic):
    """All 2^n patterns of an n-site bridging row with additive energy -1."""
    E, v, iso = [], [], []
    for s in itertools.product((0, 1), repeat=n):
        E.append(-1.0 * sum(s)); v.append(sum(s))
        nb = lambda i: [(i - 1) % n, (i + 1) % n] if periodic else [j for j in (i - 1, i + 1) if 0 <= j < n]  # noqa: E731
        iso.append(sum(s[i] and not any(s[j] for j in nb(i)) for i in range(n)))
    return dict(C=1.0, E=E, v=v), np.array(iso, float)


def test_isolated_count_is_exact_for_independent_sites():
    """A periodic row of 12 sites, enumerated as 4096 states, gives exactly
    theta(1-theta)^2 isolated vacancies per site; a 4-site open chain (the
    coupled_v1 surface cell) counts end sites with one neighbour and gives
    (2-theta)/(2(1-theta)) times more."""
    for n, periodic in ((12, True), (4, False)):
        dom, iso = _chain(n, periodic)
        sol = Model([dom]).solve(0.3 * n, T)
        th = 0.3
        per_site = float(sol.x[0] @ iso) / n
        expect = th * (1 - th) ** 2 * (1 if periodic else (2 - th) / (2 * (1 - th)))
        assert per_site == pytest.approx(expect, rel=1e-10)


def test_coupled_v1_neutral_r600_is_reproduced():
    """coupled_v1 inputs/project.json, Model.equilibrium(94.0, 873.15)."""
    cb = pt.bridging_capacity(900)
    dom, _ = _chain(4, False)
    dom = dict(dom, C=cb / 4, E=[1.31 * e for e in dom['E']], tag='surface')
    m = Model([dom, dict(C=cb, E=[0, -0.72], v=[0, 1], tag='sub'),
               dict(C=pt.O_TOTAL - 2 * cb, E=[0, 0], v=[0, 1], tag='bulk')])
    sol = m.solve(94.0, T)
    got = [sol.vacancies(t) for t in ('surface', 'sub', 'bulk')]
    assert got == pytest.approx([13.557923761144504, 13.215478361519827, 67.2265978773728], rel=1e-9)


def test_outputs_regenerate_byte_for_byte(tmp_path):
    work = tmp_path / 'p'
    shutil.copytree(HERE, work, ignore=shutil.ignore_patterns('outputs', '__pycache__', '.pytest_cache'))
    (work / 'outputs').mkdir()
    subprocess.run([sys.executable, 'run.py'], cwd=work, check=True)
    for name in ('scenario_tof.csv', 'sample_tof_range.csv'):
        assert filecmp.cmp(work / 'outputs' / name, HERE / 'outputs' / name, shallow=False), name
