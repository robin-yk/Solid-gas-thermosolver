"""Redox cycling with memory: two reservoirs instead of one number.

The defect partition in statmech.py has a single scalar state, the total
vacancy inventory, and it re-equilibrates that total on every call. That is
the right model for one measurement at one condition and the wrong one for
a sequence of them, because it cannot represent the difference between a
vacancy that has just been made and one that has already refused to come
out once.

So the state here is a pair:

    A   accessible   inventory a CO2 exposure can still reach
    L   locked       inventory that survived an exposure and stopped
                     responding to them

with one parameter joining them. A reduction step adds its oxygen deficit
to A. A reoxidation step removes a fraction f of A, and a fraction `lock`
of what is left of A becomes L. At lock = 0 nothing is ever locked, L stays
empty, and the model is the memoryless one exactly - so the memory is a
declared addition that can be switched off, not a new assumption baked in.

What the parameters mean, and do not mean. f is the fraction of the
accessible pool one exposure recovers; it is an empirical per-exposure
number here, and the placeholder kinetics or the CGMC layer can be asked
to predict it instead. `lock` is a phenomenological survival probability
and not a mechanism: it is equally consistent with vacancies migrating
beyond reach, with aggregation into an ordered defect, and with a local
structural relaxation that this model does not resolve. Which of those it
is cannot be decided from cycle integrals alone, and this module does not
pretend otherwise.

Two cycles constrain two parameters and no more. An unlocking rate - locked
inventory slowly returning to A - is the obvious third, and it is left out
rather than fitted to data that cannot see it.
"""

import math


def initial_state(accessible=0.0, locked=0.0):
    return {'accessible': float(accessible), 'locked': float(locked)}


def total(state):
    return state['accessible'] + state['locked']


def reduce_step(state, created):
    """A reduction adds its oxygen deficit to the accessible pool.

    Newly formed vacancies are made where the reductant reaches, which is
    the same place the oxidant reaches, so they start accessible. If a
    reduction were to create locked inventory directly, that would be a
    different mechanism and would need saying."""
    if created < 0.0:
        raise ValueError('a reduction step cannot create negative inventory')
    return {'accessible': state['accessible'] + float(created),
            'locked': state['locked']}


def reoxidise_step(state, f_rec, lock=1.0):
    """One exposure: recover f of the accessible pool, lock some of the rest.

    Returns the new state and how much oxygen went in, which is the
    quantity a CO2 integral measures."""
    if not 0.0 <= f_rec <= 1.0:
        raise ValueError('f_rec is a fraction of the accessible pool')
    if not 0.0 <= lock <= 1.0:
        raise ValueError('lock is a probability')
    a = state['accessible']
    recovered = f_rec * a
    left = a - recovered
    return ({'accessible': left * (1.0 - lock),
             'locked': state['locked'] + left * lock},
            recovered)


def run(schedule, f_rec, lock=1.0, state=None):
    """A sequence of reduce / reoxidise pairs, carrying the state forward.

    schedule is a list of created inventories, one per cycle. f_rec may be
    a number or a callable taking (cycle_index, state) so a kinetic layer
    can supply it."""
    st = state or initial_state()
    rows = []
    for i, created in enumerate(schedule):
        before_total = total(st)
        st = reduce_step(st, created)
        pre = dict(st)
        f = f_rec(i, st) if callable(f_rec) else f_rec
        st, recovered = reoxidise_step(st, f, lock)
        rows.append({
            'cycle': i + 1, 'created_umol_g': float(created),
            'carried_in_umol_g': before_total,
            'pre_total_umol_g': total(pre),
            'pre_accessible_umol_g': pre['accessible'],
            'pre_locked_umol_g': pre['locked'],
            'f_rec_used': f,
            'recovered_umol_g': recovered,
            'residual_umol_g': total(st),
            'accessible_umol_g': st['accessible'],
            'locked_umol_g': st['locked'],
            'f_vs_created': recovered / created if created else 0.0,
            'f_vs_pre_total': (recovered / total(pre)
                               if total(pre) else 0.0)})
    return rows


# ------------------------------------------------------------------ fitting

def fit(observations, lock_max=1.0):
    """Back out (f, lock) from measured per-cycle integrals.

    observations is a list of {'created', 'recovered'} in umol-O/g, in
    order. The first cycle starts from a fresh sample, so it fixes f
    directly; the remaining cycles then determine how much of each
    residual survived, because that is the only thing left free.

    lock is a probability, so it is capped at 1. When the cap binds, the
    honest reading is that complete locking still over-predicts the later
    cycles and f itself must be falling - which the return says explicitly
    rather than absorbing into an out-of-range parameter."""
    obs = list(observations)
    if len(obs) < 2:
        raise ValueError('two cycles are the minimum that constrains lock')
    first = obs[0]
    if not first['created'] > 0:
        raise ValueError('the first cycle must create some inventory')
    f = first['recovered'] / first['created']
    if not 0.0 < f <= 1.0:
        raise ValueError('the first cycle recovers more than it created')

    # with f fixed, each later cycle says what its accessible pool was,
    # and the difference from the fresh inventory is what did not survive
    # Only the part of the ACCESSIBLE pool that an exposure failed to
    # recover is available to lock; inventory locked in an earlier cycle
    # is already gone from A and must not be counted again. So the
    # denominator is A_n - R_n and not the cumulative residual - those
    # coincide for two cycles and diverge from the third onwards, which is
    # exactly where a fit against the wrong one goes quietly wrong.
    locks, raw, notes = [], [], []
    left = first['created'] - first['recovered']
    for n, row in enumerate(obs[1:], start=2):
        need = row['recovered'] / f              # accessible pool implied
        carried = need - row['created']          # ... that came from before
        lk = lock_max if left <= 0.0 else 1.0 - carried / left
        raw.append(lk)
        if lk > lock_max:
            # at the bound only the fresh inventory is accessible, so the
            # smallest recovery the model can predict is f * created
            over = f * row['created'] - row['recovered']
            notes.append(
                'cycle %d: complete locking is the least recovery this '
                'model can predict and it still over-predicts by %.2f '
                'umol-O/g (%.1f%%), so f is falling as well as the '
                'residual surviving'
                % (n, over, 100.0 * over / row['recovered']))
            lk = lock_max
        if lk < 0.0:
            notes.append('cycle %d: implies negative locking, i.e. the '
                         'residual is more accessible than fresh inventory'
                         % n)
        locks.append(lk)
        left = need - row['recovered']           # A_n - R_n for the next

    lock = sum(locks) / len(locks)
    lock = min(lock_max, max(0.0, lock))
    pred = run([o['created'] for o in obs], f, lock)
    resid = [pred[i]['recovered_umol_g'] - obs[i]['recovered']
             for i in range(len(obs))]
    return {'f_rec': f, 'lock': lock, 'lock_per_cycle': locks,
            'lock_unclamped': raw,
            'lock_at_bound': any(x >= lock_max for x in locks),
            'predicted': pred,
            'residual_umol_g': resid,
            'rms_umol_g': math.sqrt(sum(r * r for r in resid) / len(resid)),
            'notes': notes,
            'f_per_cycle': [o['recovered'] / o['created'] if o['created']
                            else 0.0 for o in obs]}


# -------------------------------------------------- the discriminating test

def two_sample_prediction(total_vo, f_rec, lock=1.0, split=None):
    """Two samples at the same total inventory, reached by different paths.

    The experiment that separates a memoryless inventory from one with a
    locked reservoir: make the same total V_O in one long reduction and in
    several short reduce / reoxidise cycles, then measure the CO2 recovery
    of each. A memoryless model says they must agree, because total V_O is
    its whole state. A locked reservoir says the cycled sample recovers
    less, by exactly the amount that has stopped responding.

    split, when given, is the number of cycles used to build the cycled
    sample; the default builds it in two."""
    fresh = initial_state(accessible=total_vo)
    _, rec_fresh = reoxidise_step(fresh, f_rec, lock)

    n = 2 if split is None else int(split)
    if n < 1:
        raise ValueError('the cycled sample needs at least one cycle')
    st = initial_state()
    # build up to the same total by repeated reduce / reoxidise, then a
    # final reduction that brings the total to the target
    per = total_vo / float(n + 1)
    for _ in range(n):
        st = reduce_step(st, per)
        st, _ = reoxidise_step(st, f_rec, lock)
    st = reduce_step(st, total_vo - total(st))
    cycled = dict(st)
    _, rec_cycled = reoxidise_step(cycled, f_rec, lock)

    return {'total_umol_g': total_vo, 'f_rec': f_rec, 'lock': lock,
            'cycles_to_build': n,
            'fresh': {'accessible': total_vo, 'locked': 0.0,
                      'recovered_umol_g': rec_fresh},
            'cycled': {'accessible': cycled['accessible'],
                       'locked': cycled['locked'],
                       'recovered_umol_g': rec_cycled},
            'recovery_gap_umol_g': rec_fresh - rec_cycled,
            'recovery_ratio': (rec_cycled / rec_fresh if rec_fresh
                               else float('nan')),
            'note': ('A memoryless inventory predicts a gap of zero at any '
                     'f, because total V_O is its entire state. Any '
                     'measured gap is the locked reservoir, and its size '
                     'is the lock parameter read directly off the '
                     'experiment.')}
