"""The operating atlas, built from the browser's live solvers."""

import json
import pathlib
import shutil
import subprocess

import pytest

from solidgas import activeset as A


ROOT = pathlib.Path(__file__).resolve().parents[1]
HARNESS = ROOT / 'tests' / 'js' / 'atlas_figures_harness.js'


@pytest.fixture(scope='module')
def atlas():
    assert shutil.which('node') is not None, 'node is required for this gate'
    run = subprocess.run(['node', str(HARNESS), str(ROOT), '900'],
                         capture_output=True, text=True, timeout=180)
    assert run.returncode == 0, run.stderr[-2000:]
    return json.loads(run.stdout)


def test_all_three_atlas_figures_draw(atlas):
    assert atlas['bytes']['map'] > 100_000
    assert atlas['bytes']['trade'] > 4_000
    assert atlas['bytes']['coupling'] > 4_000
    assert atlas['viewBoxes'] == {
        'map': '0 0 720 360', 'trade': '0 0 360 360',
        'coupling': '0 0 360 360'}


def test_map_contains_every_computed_cell_and_boundary_point(atlas):
    payload = atlas['payload']['map']
    assert len(payload['cells']) == payload['nT'] * payload['nR'] == 1584
    assert len(payload['boundary']) == 181
    assert all(0 <= q['conv_pct'] <= 100 for q in payload['cells'])


def test_map_marker_is_the_python_gas_only_result(atlas):
    marker = atlas['payload']['map']['current']
    expected = A.rwgs_gas_only(1173.15, 1)
    assert marker['conv_pct'] == pytest.approx(
        expected['CO2_conversion_pct'], rel=2e-13)


def test_coupling_boundary_and_both_models_are_live(atlas):
    payload = atlas['payload']['coupling']
    boundary = A.reduction_boundary(1173.15)
    expected_ratio = (1 - boundary['y_CO2']) / boundary['y_CO2']
    assert payload['boundary_ratio'] == pytest.approx(expected_ratio,
                                                       rel=1e-10)
    assert any(q['phases'] != ['TiO2'] for q in payload['rows'])
    assert max(q['gas_only'] - q['gas_solid'] for q in payload['rows']) > 0.1


def test_figure_language_names_the_two_calculations(atlas):
    all_text = ' '.join(v for rows in atlas['text'].values() for v in rows)
    assert 'gas phase only' in all_text
    assert 'gas + finite Ti–O charge' in all_text
    assert 'reduction begins' in all_text
