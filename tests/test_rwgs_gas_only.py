"""Gas-only RWGS: high-precision oracle, Python, and browser mirror."""

import json
import pathlib
import shutil
import subprocess

import pytest

from solidgas import activeset as A


ROOT = pathlib.Path(__file__).resolve().parents[1]
REF = ROOT / 'data' / 'reference_rwgs_gas_only.json'
HARNESS = ROOT / 'tests' / 'js' / 'rwgs_harness.js'
FIELDS = ('Kp', 'extent_per_mol_feed', 'CO2_conversion_pct',
          'H2_utilization_pct', 'CO_yield_per_mol_feed_pct')


@pytest.fixture(scope='module')
def reference():
    assert REF.exists(), 'run python3 scripts/oracle_rwgs.py'
    return json.loads(REF.read_text())


@pytest.fixture(scope='module')
def browser(reference):
    assert shutil.which('node') is not None, 'node is required for this gate'
    run = subprocess.run(['node', str(HARNESS), str(ROOT)],
                         capture_output=True, text=True, timeout=120)
    assert run.returncode == 0, run.stderr[-2000:]
    rows = json.loads(run.stdout)
    assert len(rows) == len(reference['rows'])
    return rows


def test_python_and_browser_match_the_80_digit_reference(reference, browser):
    for oracle, js in zip(reference['rows'], browser):
        py = A.rwgs_gas_only(float(oracle['T_C']) + 273.15,
                             float(oracle['H2_CO2_ratio']))
        key = (oracle['T_C'], oracle['H2_CO2_ratio'])
        for field in FIELDS:
            target = float(oracle[field])
            assert py[field] == pytest.approx(target, rel=2e-13, abs=2e-13), \
                (key, field, 'python')
            assert js[field] == pytest.approx(target, rel=2e-13, abs=2e-13), \
                (key, field, 'browser')


@pytest.mark.parametrize('temperature_c,ratio', [
    (400, 0.1), (600, 1), (900, 3), (1200, 100), (1500, 1e7)
])
def test_fresh_feed_mass_balance_and_stoichiometry(temperature_c, ratio):
    result = A.rwgs_gas_only(temperature_c + 273.15, ratio)
    gas = result['gas_fractions']
    assert sum(gas.values()) == pytest.approx(1.0, abs=2e-15)
    assert gas['CO'] == pytest.approx(gas['H2O'], abs=2e-15)
    assert gas['CO2'] + gas['CO'] == pytest.approx(1 / (1 + ratio),
                                                   abs=2e-15)
    assert gas['H2'] + gas['H2O'] == pytest.approx(ratio / (1 + ratio),
                                                   abs=2e-15)


def test_the_one_to_one_value_is_the_workspace_value():
    result = A.rwgs_gas_only(1173.15, 1)
    assert result['CO2_conversion_pct'] == pytest.approx(53.0095951542,
                                                         abs=1e-10)


def test_invalid_ratio_is_rejected():
    with pytest.raises(ValueError, match='positive'):
        A.rwgs_gas_only(1173.15, 0)
