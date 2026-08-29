"""Redox cycling with memory, and the one experiment that measures it.

The model adds exactly one parameter to the memoryless picture, so the
gates are mostly about what that parameter does at its limits: at zero the
model must be the old one to the last digit, and at one it must give the
least recovery the construction can produce. In between, mass has to be
conserved cycle by cycle, because a bookkeeping model that loses oxygen is
worth nothing.

The two reported cycles are carried as a case with their own numbers, and
the discriminating prediction is gated on the property that makes it
discriminating: a memoryless inventory predicts no gap at all.
"""

import math

import pytest

from solidgas import cycling as CY

# the reported integrals, in umol-O per g TiO2
REPORTED = [{'created': 94.8, 'recovered': 81.5},
            {'created': 61.0, 'recovered': 51.6}]


# --------------------------------------------------------- the limits

def test_no_locking_is_the_memoryless_model(p=None):
    """At lock = 0 the total inventory is again the whole state.

    Two schedules that reach the same total by different paths have to
    give the same recovery, which is exactly the property the memoryless
    partition has and the reason it cannot see cycle history."""
    f = 0.86
    a = CY.run([90.0], f, lock=0.0)
    b = CY.run([30.0, 30.0, 30.0], f, lock=0.0)
    assert CY.total({'accessible': a[-1]['accessible_umol_g'],
                     'locked': a[-1]['locked_umol_g']}) > 0
    for row in b:
        assert row['locked_umol_g'] == 0.0
    q = CY.two_sample_prediction(95.0, f, lock=0.0)
    assert q['recovery_gap_umol_g'] == pytest.approx(0.0, abs=1e-12)
    assert q['recovery_ratio'] == pytest.approx(1.0, abs=1e-12)


def test_complete_locking_is_the_least_recovery_the_model_can_give(p=None):
    """At lock = 1 only the freshly created inventory is reachable."""
    f = 0.86
    rows = CY.run([94.8, 61.0], f, lock=1.0)
    assert rows[1]['pre_accessible_umol_g'] == pytest.approx(61.0)
    assert rows[1]['recovered_umol_g'] == pytest.approx(f * 61.0)
    for lock in (0.0, 0.25, 0.5, 0.75, 1.0):
        r = CY.run([94.8, 61.0], f, lock=lock)[1]['recovered_umol_g']
        assert r <= f * (94.8 - f * 94.8 + 61.0) + 1e-12
        assert r >= f * 61.0 - 1e-12


def test_recovery_falls_monotonically_with_locking(p=None):
    prev = None
    for lock in (0.0, 0.2, 0.4, 0.6, 0.8, 1.0):
        r = CY.run([94.8, 61.0], 0.86, lock=lock)[1]['recovered_umol_g']
        if prev is not None:
            assert r < prev, lock
        prev = r


# --------------------------------------------------------- bookkeeping

def test_every_cycle_closes_on_oxygen(p=None):
    """created = recovered + the change in residual, exactly."""
    for lock in (0.0, 0.37, 1.0):
        for f in (0.1, 0.5, 0.86, 1.0):
            rows = CY.run([94.8, 61.0, 40.0, 12.0], f, lock=lock)
            prev = 0.0
            for r in rows:
                closes = (r['created_umol_g'] - r['recovered_umol_g']
                          - (r['residual_umol_g'] - prev))
                assert abs(closes) < 1e-12, (lock, f, r['cycle'])
                prev = r['residual_umol_g']
            assert rows[-1]['accessible_umol_g'] >= 0.0
            assert rows[-1]['locked_umol_g'] >= 0.0


def test_the_residual_never_falls(p=None):
    rows = CY.run([94.8, 61.0, 40.0, 12.0], 0.86, lock=1.0)
    res = [r['residual_umol_g'] for r in rows]
    assert all(b > a for a, b in zip(res, res[1:])), res


def test_the_two_normalisations_diverge_once_anything_locks(p=None):
    """Why the definition of 'recovery fraction' has to be stated.

    Against the newly created inventory the fraction stays at f forever;
    against the total present it falls as the locked pool grows. The two
    agree only on the first cycle, and only the first."""
    rows = CY.run([94.8, 61.0, 61.0], 0.86, lock=1.0)
    assert rows[0]['f_vs_created'] == pytest.approx(
        rows[0]['f_vs_pre_total'], abs=1e-12)
    for r in rows[1:]:
        assert r['f_vs_created'] == pytest.approx(0.86, abs=1e-12)
        assert r['f_vs_pre_total'] < r['f_vs_created']
    tot = [r['f_vs_pre_total'] for r in rows]
    assert all(b < a for a, b in zip(tot, tot[1:])), tot


# ------------------------------------------------------------- fitting

def test_the_fit_recovers_parameters_it_was_given(p=None):
    """Synthetic cycles in, the same (f, lock) out."""
    for f in (0.4, 0.7, 0.86):
        for lock in (0.0, 0.3, 0.65, 1.0):
            rows = CY.run([94.8, 61.0, 61.0], f, lock=lock)
            obs = [{'created': r['created_umol_g'],
                    'recovered': r['recovered_umol_g']} for r in rows]
            got = CY.fit(obs)
            assert got['f_rec'] == pytest.approx(f, abs=1e-12), (f, lock)
            assert got['lock'] == pytest.approx(lock, abs=1e-9), (f, lock)
            assert got['rms_umol_g'] < 1e-9, (f, lock)


def test_a_single_cycle_cannot_constrain_the_memory(p=None):
    with pytest.raises(ValueError, match='two cycles'):
        CY.fit(REPORTED[:1])


def test_the_reported_cycles_pin_the_locking_at_its_bound():
    """What the measurement actually says, with the caveat it comes with.

    The first cycle fixes f at 0.860. The second then needs an accessible
    pool of 60.0 umol-O/g, which is less than the 61.0 it was just given -
    so not only did none of the 13.3 residual come back, the fresh
    inventory itself recovered slightly less well. Complete locking is the
    least recovery this construction can predict and it still over-shoots
    by 1.6%, which is the model telling you where it ends."""
    got = CY.fit(REPORTED)
    assert got['f_rec'] == pytest.approx(0.8597, abs=5e-4)
    assert got['lock'] == 1.0
    assert got['lock_at_bound'] is True
    assert got['lock_unclamped'][0] > 1.0
    assert got['notes'] and 'over-predicts' in got['notes'][0]
    assert got['residual_umol_g'][0] == pytest.approx(0.0, abs=1e-12)
    over = got['residual_umol_g'][1]
    assert 0.0 < over < 1.0
    assert over / REPORTED[1]['recovered'] < 0.02
    # and the per-cycle recovery against fresh inventory barely moves,
    # which is the reading the locking makes correct rather than assumed
    f1, f2 = got['f_per_cycle']
    assert abs(f1 - f2) / f1 < 0.02


def test_the_memoryless_model_gets_the_trend_backwards(REPORTED_=REPORTED):
    """Why the memory had to be added at all.

    Against the total present, the measurement falls from 0.860 to 0.694.
    A memoryless model re-equilibrates the whole inventory every time, so
    with a smaller total it can only predict a recovery fraction that is
    the same or larger. The sign of the change is the evidence."""
    obs_ratio = [o['recovered'] for o in REPORTED_]
    pre = [94.8, 94.8 - 81.5 + 61.0]
    measured = [obs_ratio[i] / pre[i] for i in range(2)]
    assert measured[1] < measured[0]
    memoryless = CY.run([o['created'] for o in REPORTED_], 0.8597, lock=0.0)
    modelled = [r['f_vs_pre_total'] for r in memoryless]
    assert modelled[1] == pytest.approx(modelled[0], abs=1e-12), \
        'without memory the fraction cannot move at all'


# ------------------------------------------------- the discriminating test

def test_the_two_sample_gap_is_the_locking_read_off_the_experiment():
    """Same total inventory, different history: the gap is the memory.

    Zero at lock = 0 for any f, because total V_O is then the whole state.
    Monotone in lock, so a measured gap reads the parameter directly."""
    prev = -1.0
    for lock in (0.0, 0.25, 0.5, 0.75, 1.0):
        q = CY.two_sample_prediction(95.0, 0.8597, lock=lock)
        gap = q['recovery_gap_umol_g']
        assert gap > prev, lock
        prev = gap
        assert q['cycled']['locked'] >= 0.0
        assert (q['cycled']['accessible'] + q['cycled']['locked']
                == pytest.approx(95.0, rel=1e-12)), \
            'both samples must be compared at the same total inventory'
    q = CY.two_sample_prediction(95.0, 0.8597, lock=1.0)
    assert q['recovery_gap_umol_g'] == pytest.approx(7.64, abs=0.05)
    assert q['recovery_ratio'] == pytest.approx(0.906, abs=5e-3)
    for f in (0.3, 0.6, 0.9):
        assert CY.two_sample_prediction(
            95.0, f, lock=0.0)['recovery_gap_umol_g'] == pytest.approx(
                0.0, abs=1e-12)


def test_building_the_cycled_sample_more_slowly_widens_the_gap():
    """More cycles to the same total means more chances to lock."""
    gaps = [CY.two_sample_prediction(95.0, 0.8597, lock=1.0,
                                     split=n)['recovery_gap_umol_g']
            for n in (1, 2, 4, 8)]
    assert all(b > a for a, b in zip(gaps, gaps[1:])), gaps


# ---------------------------------------------------------- input guards

def test_the_inputs_are_checked():
    with pytest.raises(ValueError):
        CY.reduce_step(CY.initial_state(), -1.0)
    with pytest.raises(ValueError):
        CY.reoxidise_step(CY.initial_state(10.0), 1.5)
    with pytest.raises(ValueError):
        CY.reoxidise_step(CY.initial_state(10.0), 0.5, lock=-0.1)
    with pytest.raises(ValueError):
        CY.two_sample_prediction(95.0, 0.86, split=0)
    with pytest.raises(ValueError):
        CY.fit([{'created': 0.0, 'recovered': 0.0},
                {'created': 1.0, 'recovered': 0.5}])
