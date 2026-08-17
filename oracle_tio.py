"""High-precision numerical oracle for the Ti-O gas/solid Gibbs equilibrium.

This is the independent reference the production active-set solver is tested
against. It is not a port of the production code and it must never become one:

  * it does NOT import the production GibbsSystem class, the legacy
    minimisers, or any solver code from `solidgas` - only the raw
    thermodynamic coefficient tables (SHOMATE, WALDNER), which both sides
    are required to share;
  * every arithmetic step runs in mpmath at 80 significant digits, so a
    solid present at 1e-41 mol and a Gibbs difference of 1e-34 kJ are
    ordinary numbers here, not casualties of double precision;
  * it solves each candidate active set of condensed phases as a square
    root-finding problem in element-potential coordinates - a different
    formulation from the production minimiser - and classifies every
    candidate by its KKT conditions.

Formulation, per candidate active set A of condensed phases:

    gas species      mu_i = mu0_i + RT ln(y_i P)  = sum_e a_ie lambda_e
    active solids    mu0_j                        = t_j lambda_Ti + o_j lambda_O
    balances         sum_i a_ie n_i + sum_j a_je m_j = b_e
    normalisation    sum_i y_i = 1

The solid mole numbers are recovered from the two balances that contain them
(Ti, and the combination O' = O - 2*Ti) *linearly*, so no Newton variable ever
sits at 1e-41; and the O' combination cancels the 2*n_Ti term in b_O
symbolically, so the oxygen actually available to the gas is formed without
subtracting two nearly equal numbers.

|A| = 2 pins (lambda_Ti, lambda_O) exactly (two line compounds, two unknowns:
the two-phase buffer). |A| = 1 leaves lambda_O to the gas. |A| = 3 is
over-determined at generic temperature and is solved only as a degeneracy
probe. Gibbs' phase rule for these M-O condensed phases allows two generically
- C and H never enter a solid here, so they do not raise the count.

KKT reduced costs are reported in two bases, and they are not the same number
once the registry leaves the Magneli series:

    r_j (per formula) = mu0_j - t_j lambda_Ti - o_j lambda_O
    r_j (per mol O)   = r_j / def_j,   def_j = 2 t_j - o_j

For every TiO2 -> Magneli step def_j = 1, which is the accident that made the
two bases agree in earlier work on this system.

Zero-inventory rule: a gas element whose available inventory is exactly zero
has its species removed (their equilibrium amount is exactly zero) and its
element potential is undefined. If the candidate then leaves a reducible
phase inactive, the directional derivative of G toward creating a trace of
O-bearing gas plus that phase is -infinity (the xi*ln(xi) edge), so the
candidate fails KKT with r = -inf. This is why a bone-dry reducing feed always
sits on a two-phase buffer rather than on pure rutile.

Run:  python3 oracle_tio.py          -> reference_results_high_precision.json
"""

import json
import sys

from mpmath import mp, mpf, matrix, lu_solve, log, exp, sqrt, floor

# Data-only imports: coefficient tables. No solver code crosses this line.
from solidgas.shomate import SHOMATE, BREAKPOINT, DEFAULT_BREAKPOINT
from solidgas.waldner import WALDNER

mp.dps = 90          # working precision; the contract is 80 significant digits
REPORT_DPS = 80

R_KJ = mpf('8.314472') / 1000       # kJ / mol K
R_ATM = mpf('82.06')                # cm3 atm / mol K

GASES = ['CO2', 'H2', 'CO', 'H2O', 'N2', 'O2']
SOLIDS = ['TiO2', 'Ti10O19', 'Ti8O15', 'Ti6O11', 'Ti5O9',
          'Ti4O7', 'Ti3O5', 'Ti2O3']

GAS_FORMULA = {
    'CO2': {'C': 1, 'O': 2}, 'H2': {'H': 2}, 'CO': {'C': 1, 'O': 1},
    'H2O': {'H': 2, 'O': 1}, 'N2': {'N': 2}, 'O2': {'O': 2},
}
STOICH = {'TiO2': (1, 2), 'Ti10O19': (10, 19), 'Ti8O15': (8, 15),
          'Ti6O11': (6, 11), 'Ti5O9': (5, 9), 'Ti4O7': (4, 7),
          'Ti3O5': (3, 5), 'Ti2O3': (2, 3)}
# Ti(3+) per formula unit: two for every Magneli member and for Ti2O3.
TI3 = {p: (0 if p == 'TiO2' else 2) for p in SOLIDS}

ATOMIC_WEIGHT = {'Ti': mpf('47.867'), 'O': mpf('15.999')}

P_ATM = mpf(1)
V_CM3 = mpf(25)
MASS_G = mpf('0.100')
TEMPS_C = [500, 700, 900, 1100, 1300, 1500]

# KKT / classification tolerances, in kJ. At 80 digits the algebra closes to
# ~1e-70, so these separate physics from arithmetic with room to spare.
TOL_KKT = mpf('1e-30')       # an inactive phase is stable if r_j > -TOL_KKT
TOL_DEG = mpf('1e-6')        # |r_j| below this flags a degenerate boundary
TOL_RES = mpf('1e-60')       # Newton convergence, scaled residual
TRACE_MOL = mpf('1e-12')     # below this, log10(n) is reported alongside

# --------------------------------------------------------------- mu0, mpmath


def _mp(x):
    return mpf(x) if not isinstance(x, mpf) else x


def mu0_shomate(name, T):
    """NIST Shomate G(T) in kJ/mol, evaluated at working precision."""
    lo, hi, dHf = SHOMATE[name]
    brk = BREAKPOINT.get(name, DEFAULT_BREAKPOINT)
    A, B, C, D, E, F, G, H = [_mp(c) for c in (lo if T <= brk else hi)]
    dHf = _mp(dHf)
    t = T / 1000
    Ht = dHf + (A * t + B * t**2 / 2 + C * t**3 / 3 + D * t**4 / 4
                - E / t + F - H)
    St = (A * log(t) + B * t + C * t**2 / 2 + D * t**3 / 3
          - E / (2 * t**2) + G)
    return Ht - T * St / 1000


def _cp_int(form, c, T0, T1):
    a, b, cc, d = [_mp(x) for x in c]
    if form == 'poly':
        F = lambda T: a * T + b * T**2 / 2 + cc * T**3 / 3 - d / T
    else:  # magneli: Cp = k0 + k1 T^-1/2 + k2 T^-2 + k3 T^-3
        F = lambda T: a * T + 2 * b * sqrt(T) - cc / T - d / (2 * T**2)
    return F(T1) - F(T0)


def _cp_over_t_int(form, c, T0, T1):
    a, b, cc, d = [_mp(x) for x in c]
    if form == 'poly':
        F = lambda T: a * log(T) + b * T + cc * T**2 / 2 - d / (2 * T**2)
    else:
        F = lambda T: a * log(T) - 2 * b / sqrt(T) - cc / (2 * T**2) - d / (3 * T**3)
    return F(T1) - F(T0)


def mu0_waldner(name, T):
    """Waldner & Eriksson G(T) in kJ/mol, piecewise Cp integrated from 298.15 K.

    Transition enthalpies enter H as dHt and S as dHt/T_transition at the
    range join, exactly as the assessment prescribes.
    """
    d = WALDNER[name]
    H = _mp(d['dHf298'])
    S = _mp(d['S298'])
    lo = mpf('298.15')
    segs = []
    for (r0, r1, form, c, dHt) in d['ranges']:
        r0m, r1m = _mp(r0), _mp(r1)
        if T <= r0m:
            break
        hi = T if T < r1m else r1m
        if hi > (lo if lo > r0m else r0m):
            segs.append((max(lo, r0m), hi, form, c, _mp(dHt), r0m))
            lo = hi
        if T <= r1m:
            break
    if segs and segs[-1][1] < T:            # extrapolate the last range
        r0, r1, form, c, dHt = d['ranges'][-1]
        segs.append((segs[-1][1], T, form, c, mpf(0), _mp(r0)))
    for (t0, t1, form, c, dHt, T_tr) in segs:
        H += dHt + _cp_int(form, c, t0, t1)
        S += (dHt / T_tr if dHt else mpf(0)) + _cp_over_t_int(form, c, t0, t1)
    return (H - T * S) / 1000


def cp_waldner(name, T):
    for (r0, r1, form, c, _dHt) in WALDNER[name]['ranges']:
        if T <= _mp(r1):
            a, b, cc, d = [_mp(x) for x in c]
            if form == 'poly':
                return a + b * T + cc * T**2 + d / T**2
            return a + b / sqrt(T) + cc / T**2 + d / T**3
    r0, r1, form, c, _dHt = WALDNER[name]['ranges'][-1]
    a, b, cc, d = [_mp(x) for x in c]
    if form == 'poly':
        return a + b * T + cc * T**2 + d / T**2
    return a + b / sqrt(T) + cc / T**2 + d / T**3


def cp_shomate(name, T):
    lo, hi, _ = SHOMATE[name]
    brk = BREAKPOINT.get(name, DEFAULT_BREAKPOINT)
    A, B, C, D, E, _F, _G, _H = [_mp(c) for c in (lo if T <= brk else hi)]
    t = T / 1000
    return A + B * t + C * t**2 + D * t**3 + E / t**2


def mu_table(T):
    """Every chemical potential needed at one temperature, computed once."""
    mu = {g: mu0_shomate(g, T) for g in GASES}
    for s in SOLIDS:
        mu[s] = mu0_waldner(s, T)
    return mu


# ------------------------------------------------------------ feed and cases

CASES = {
    'rwgs_1_1':     dict(label='RWGS feed CO2:H2 = 1:1',
                         feed={'CO2': 1, 'H2': 1}),
    'h2_rich':      dict(label='H2-rich RWGS feed CO2:H2 = 1:3',
                         feed={'CO2': 1, 'H2': 3}),
    'product_mix':  dict(label='RWGS product/recycle CO2:H2:CO:H2O = 3:3:2:2',
                         feed={'CO2': 3, 'H2': 3, 'CO': 2, 'H2O': 2}),
    'pure_h2':      dict(label='pure H2',
                         feed={'H2': 1}),
    'pure_n2':      dict(label='pure N2',
                         feed={'N2': 1}),
    'n2_10ppm_o2':  dict(label='N2 with 10 ppm O2',
                         feed={'N2': 999990, 'O2': 10}),
}


def build_point(feed_rel, T, P=P_ATM, V=V_CM3, mass=MASS_G, T_charge=None):
    """Feed moles, element inventory, and the cancellation-free split.

    n_g0 = P V / (R T_charge); T_charge defaults to the reaction temperature,
    which is the sealed-vessel convention this repository has used throughout.
    """
    Tc = T if T_charge is None else T_charge
    n_g0 = P * V / (R_ATM * Tc)
    tot = sum(_mp(v) for v in feed_rel.values())
    y0 = {g: _mp(feed_rel.get(g, 0)) / tot for g in GASES}
    n0 = {g: y0[g] * n_g0 for g in GASES}
    mw = ATOMIC_WEIGHT['Ti'] + 2 * ATOMIC_WEIGHT['O']
    n_Ti = mass / mw                      # charge loaded as rutile
    b = {'C': mpf(0), 'H': mpf(0), 'N': mpf(0), 'O': mpf(0), 'Ti': n_Ti}
    bO_gasfeed = mpf(0)
    for g in GASES:
        for e, k in GAS_FORMULA[g].items():
            b[e] += k * n0[g]
    bO_gasfeed = sum(GAS_FORMULA[g].get('O', 0) * n0[g] for g in GASES)
    b['O'] += 2 * n_Ti                    # host solid oxygen
    return dict(n_g0=n_g0, y0=y0, n0=n0, n_Ti=n_Ti, b=b,
                bO_gasfeed=bO_gasfeed, P=P, V=V, T=T, T_charge=Tc, mass=mass)


# ---------------------------------------------------------- candidate solver


def solve_candidate(active, pt, mu, seed=None):
    """One active set: square Newton for the gas, linear recovery of solids.

    Returns a record with feasibility, KKT reduced costs of every inactive
    phase, moles, potentials, and the reason for any rejection.
    """
    T, P, b = pt['T'], pt['P'], pt['b']
    RT = R_KJ * T
    tj = {s: mpf(STOICH[s][0]) for s in active}
    oj = {s: mpf(STOICH[s][1]) for s in active}
    dj = {s: 2 * tj[s] - oj[s] for s in active}   # O freed per formula vs TiO2

    lam = {}
    triple_residual = None
    if len(active) >= 2:
        a1, a2 = active[0], active[1]
        det = tj[a1] * oj[a2] - tj[a2] * oj[a1]
        lam['Ti'] = (mu[a1] * oj[a2] - mu[a2] * oj[a1]) / det
        lam['O'] = (tj[a1] * mu[a2] - tj[a2] * mu[a1]) / det
        if len(active) == 3:
            a3 = active[2]
            triple_residual = mu[a3] - tj[a3] * lam['Ti'] - oj[a3] * lam['O']

    # Which gas elements exist at all, and how much oxygen the gas can hold.
    have = {e: b[e] > 0 for e in ('C', 'H', 'N')}
    if len(active) == 1:
        s0 = active[0]
        m_fixed = b['Ti'] / tj[s0]
        avail_O = pt['bO_gasfeed'] + dj[s0] * m_fixed
        o_free = avail_O > 0
    else:
        m_fixed = None
        avail_O = None
        o_free = False                    # lambda_O pinned by the buffer

    def elem_ok(e):
        if e in ('C', 'H', 'N'):
            return have[e]
        if e == 'O':
            return o_free or len(active) >= 2
        return False

    species = [g for g in GASES
               if all(elem_ok(e) for e in GAS_FORMULA[g])]
    if not species:
        return dict(active=list(active), feasible=False, kkt_ok=False,
                    reason='no gas species with inventory')

    newton_elems = [e for e in ('C', 'H', 'N') if have[e]]
    if o_free:
        newton_elems.append('O')

    # ln y_i = sum_e a_ie theta_e - g_i, with pinned lambda_O folded into g_i.
    gi = {}
    for s in species:
        v = mu[s] / RT + log(P)
        if len(active) >= 2 and 'O' in GAS_FORMULA[s]:
            v -= GAS_FORMULA[s]['O'] * lam['O'] / RT
        gi[s] = v

    rhs = {e: b[e] for e in ('C', 'H', 'N')}
    if o_free:
        rhs['O'] = avail_O

    nE = len(newton_elems)
    if seed is not None and len(seed) == nE + 1:
        x = matrix(seed)
    else:
        # Seed: feed-shaped gas, floored; theta from least squares on ln y.
        yg = {}
        for s in species:
            yg[s] = pt['y0'].get(s, mpf(0))
        floor_ = mpf('1e-8')
        tot = sum(yg.values())
        for s in species:
            yg[s] = (yg[s] / tot if tot > 0 else mpf(1) / len(species))
            if yg[s] < floor_:
                yg[s] = floor_
        tot = sum(yg.values())
        for s in species:
            yg[s] /= tot
        # least squares A theta = ln y + g
        if nE:
            A = matrix(len(species), nE)
            rv = matrix(len(species), 1)
            for i, s in enumerate(species):
                for k, e in enumerate(newton_elems):
                    A[i, k] = GAS_FORMULA[s].get(e, 0)
                rv[i] = log(yg[s]) + gi[s]
            AtA = A.T * A
            Atr = A.T * rv
            try:
                th = lu_solve(AtA, Atr)
            except Exception:
                th = matrix([[mpf(0)] for _ in range(nE)])
        else:
            th = matrix(0, 1)
        x = matrix(nE + 1, 1)
        for k in range(nE):
            x[k] = th[k]
        x[nE] = log(pt['n_g0'])

    scale = {e: max(rhs[e], mpf('1e-300')) for e in newton_elems}

    def gas_state(xv):
        th = {e: xv[k] for k, e in enumerate(newton_elems)}
        lnNg = xv[nE]
        lny = {}
        for s in species:
            v = -gi[s]
            for e, k in GAS_FORMULA[s].items():
                if e in th:
                    v += k * th[e]
            lny[s] = v
        y = {s: exp(lny[s]) for s in species}
        return th, lnNg, y

    def residual(xv):
        th, lnNg, y = gas_state(xv)
        Ng = exp(lnNg)
        res = matrix(nE + 1, 1)
        for k, e in enumerate(newton_elems):
            tot = sum(GAS_FORMULA[s].get(e, 0) * y[s] for s in species)
            res[k] = (tot * Ng - rhs[e]) / scale[e]
        res[nE] = sum(y.values()) - 1
        return res, y, Ng

    def jacobian(xv, y, Ng):
        J = matrix(nE + 1, nE + 1)
        for k, e in enumerate(newton_elems):
            for k2, f in enumerate(newton_elems):
                J[k, k2] = sum(GAS_FORMULA[s].get(e, 0) * GAS_FORMULA[s].get(f, 0)
                               * y[s] for s in species) * Ng / scale[e]
            J[k, nE] = sum(GAS_FORMULA[s].get(e, 0) * y[s]
                           for s in species) * Ng / scale[e]
        for k2, f in enumerate(newton_elems):
            J[nE, k2] = sum(GAS_FORMULA[s].get(f, 0) * y[s] for s in species)
        J[nE, nE] = mpf(0)
        return J

    res, y, Ng = residual(x)
    rn = max(abs(v) for v in res)
    it = 0
    while rn > TOL_RES and it < 300:
        J = jacobian(x, y, Ng)
        try:
            step = lu_solve(J, -res)
        except Exception:
            return dict(active=list(active), feasible=False, kkt_ok=False,
                        reason='singular Jacobian')
        lam_step = mpf(1)
        for _ in range(80):
            xn = x + lam_step * step
            resn, yn, Ngn = residual(xn)
            rnn = max(abs(v) for v in resn)
            if rnn < rn:
                break
            lam_step /= 2
        else:
            return dict(active=list(active), feasible=False, kkt_ok=False,
                        reason=f'stalled at residual {mp.nstr(rn, 3)}')
        x, res, y, Ng, rn = xn, resn, yn, Ngn, rnn
        it += 1
    if rn > TOL_RES:
        return dict(active=list(active), feasible=False, kkt_ok=False,
                    reason=f'no convergence, residual {mp.nstr(rn, 3)}')

    th, lnNg, y = gas_state(x)
    n_gas = {s: y[s] * Ng for s in species}
    for g in GASES:
        n_gas.setdefault(g, mpf(0))

    # Solid moles, linearly. O' = O - 2*Ti balance keeps this cancellation-free.
    gasO = sum(GAS_FORMULA[s].get('O', 0) * n_gas[s] for s in species)
    if len(active) == 1:
        m = {active[0]: m_fixed}
    elif len(active) == 2:
        a1, a2 = active
        det = tj[a1] * dj[a2] - tj[a2] * dj[a1]
        if det == 0:
            return dict(active=list(active), feasible=False, kkt_ok=False,
                        reason='degenerate stoichiometry pair')
        rO = gasO - pt['bO_gasfeed']
        m = {a1: (b['Ti'] * dj[a2] - tj[a2] * rO) / det,
             a2: (tj[a1] * rO - dj[a1] * b['Ti']) / det}
    else:
        # Degeneracy probe only: the m-family is 1-dimensional when consistent.
        m = {s: None for s in active}

    # Element potentials still missing: lambda_O (single, gas-defined) and
    # lambda_Ti; lambda_C/H/N from theta.
    if o_free:
        lam['O'] = RT * th['O']
    if len(active) == 1 and 'O' in lam:
        s0 = active[0]
        lam['Ti'] = (mu[s0] - oj[s0] * lam['O']) / tj[s0]
    for k, e in enumerate(newton_elems):
        lam[e] = RT * th[e]

    # KKT reduced costs for every phase left out of the active set.
    r_costs = {}
    lam_O_defined = 'O' in lam
    for s in SOLIDS:
        if s in active:
            continue
        t_s, o_s = mpf(STOICH[s][0]), mpf(STOICH[s][1])
        d_s = 2 * t_s - o_s
        if not lam_O_defined:
            r_costs[s] = dict(per_formula='-inf' if d_s > 0 else None,
                              per_mol_O='-inf' if d_s > 0 else None,
                              oxygen_removed_per_extent=float(d_s))
            continue
        rf = mu[s] - t_s * lam['Ti'] - o_s * lam['O']
        r_costs[s] = dict(per_formula=rf,
                          per_mol_O=(rf / d_s if d_s > 0 else None),
                          oxygen_removed_per_extent=float(d_s))

    feasible = True
    reason = 'ok'
    if len(active) == 2:
        if not (m[active[0]] > 0 and m[active[1]] > 0):
            feasible = False
            reason = 'negative condensed amount: ' + ', '.join(
                f'{s}={mp.nstr(m[s], 6)}' for s in active)
    if len(active) == 3:
        feasible = False
        reason = ('degeneracy probe; solid-equation residual '
                  f'{mp.nstr(triple_residual, 6)} kJ')

    kkt_ok = feasible
    worst = None
    for s, rc in r_costs.items():
        v = rc['per_formula']
        if v == '-inf':
            kkt_ok = False
            worst = '-inf'
            break
        if v is not None:
            if worst is None or v < worst:
                worst = v
            if v < -TOL_KKT:
                kkt_ok = False

    return dict(active=list(active), feasible=feasible, kkt_ok=kkt_ok,
                reason=reason, y=y, n_gas=n_gas, m=m, lam=lam, Ng=Ng,
                r_costs=r_costs, worst_inactive_r=worst,
                triple_residual=triple_residual, iterations=it,
                newton_residual=rn, seed_out=[x[k] for k in range(nE + 1)])


# ------------------------------------------------------------ point pipeline


def candidates():
    singles = [[s] for s in SOLIDS]
    pairs = [[SOLIDS[i], SOLIDS[j]]
             for i in range(len(SOLIDS)) for j in range(i + 1, len(SOLIDS))]
    return singles + pairs


def triples():
    n = len(SOLIDS)
    return [[SOLIDS[i], SOLIDS[j], SOLIDS[k]]
            for i in range(n) for j in range(i + 1, n) for k in range(j + 1, n)]


def g_total(n_gas, m, mu, T, P):
    RT = R_KJ * T
    Ng = sum(v for v in n_gas.values() if v > 0)
    G = mpf(0)
    for s, n in n_gas.items():
        if n > 0:
            G += n * (mu[s] + RT * log(n / Ng * P))
    for s, v in m.items():
        if v:
            G += v * mu[s]
    return G


def g_feed(pt, mu):
    RT = R_KJ * pt['T']
    G = mpf(0)
    for s, n in pt['n0'].items():
        if n > 0:
            G += n * (mu[s] + RT * log(n / pt['n_g0'] * pt['P']))
    G += pt['n_Ti'] * mu['TiO2']
    return G


def jnum(x):
    """mp number -> JSON: float when a double can carry it, string otherwise."""
    if x is None or isinstance(x, str):
        return x
    if x == 0:
        return 0.0
    a = abs(x)
    if mpf('1e-290') < a < mpf('1e290'):
        return float(x)
    return mp.nstr(x, 17)


def solve_point(case, cfg, T_C, seeds):
    T = mpf(str(T_C)) + mpf('273.15')
    pt = build_point(cfg['feed'], T)
    mu = mu_table(T)

    results = []
    passing = []
    for cand in candidates():
        key = tuple(cand)
        r = solve_candidate(cand, pt, mu, seed=seeds.get(key))
        if 'seed_out' in r:
            seeds[key] = r['seed_out']
        results.append(r)
        if r['feasible'] and r['kkt_ok']:
            passing.append(r)

    if not passing:
        raise RuntimeError(f'{case} {T_C}C: no candidate satisfies KKT')
    passing.sort(key=lambda r: len(r['active']))
    win = passing[0]

    degenerate = []
    for s, rc in win['r_costs'].items():
        v = rc['per_mol_O']
        if isinstance(v, mpf) and abs(v) < TOL_DEG:
            degenerate.append(s)
    if len(passing) > 1 and not degenerate:
        raise RuntimeError(
            f'{case} {T_C}C: {len(passing)} candidates pass KKT without '
            f'a degenerate reduced cost - convexity says this cannot happen')

    # Degeneracy probes: consistent triples only ever appear at special T.
    probe_hits = []
    for cand in triples():
        tj = {s: mpf(STOICH[s][0]) for s in cand}
        oj = {s: mpf(STOICH[s][1]) for s in cand}
        a1, a2, a3 = cand
        det = tj[a1] * oj[a2] - tj[a2] * oj[a1]
        lam_Ti = (mu[a1] * oj[a2] - mu[a2] * oj[a1]) / det
        lam_O = (tj[a1] * mu[a2] - tj[a2] * mu[a1]) / det
        res3 = mu[a3] - tj[a3] * lam_Ti - oj[a3] * lam_O
        if abs(res3) < TOL_DEG:
            probe_hits.append({'phases': cand, 'residual_kJ': jnum(res3)})

    m = {s: (v if v else mpf(0)) for s, v in win['m'].items()}
    n_gas = win['n_gas']
    G_tot = g_total(n_gas, m, mu, T, pt['P'])
    G_ref = g_feed(pt, mu)

    # Element residuals, reconstructed from scratch.
    resid = {}
    for e in ('C', 'H', 'N', 'O', 'Ti'):
        tot = sum(GAS_FORMULA[g].get(e, 0) * n_gas[g] for g in GASES)
        for s, v in m.items():
            t_s, o_s = STOICH[s]
            tot += (t_s if e == 'Ti' else o_s if e == 'O' else 0) * v
        resid[e] = mp.nstr(tot - pt['b'][e], 5) if pt['b'][e] or tot else '0.0'

    solid_O = sum(mpf(STOICH[s][1]) * v for s, v in m.items())
    solid_Ti = sum(mpf(STOICH[s][0]) * v for s, v in m.items())
    ti3 = sum(TI3[s] * v for s, v in m.items())
    gasO = sum(GAS_FORMULA[g].get('O', 0) * n_gas[g] for g in GASES)

    traces = {}
    for s, v in list(n_gas.items()) + list(m.items()):
        if 0 < v < TRACE_MOL:
            traces[s] = float(log(v) / log(mpf(10)))

    Ngtot = sum(n_gas.values())
    row = {
        'system': 'Ti-O', 'case': case, 'label': cfg['label'],
        'T_C': T_C, 'T_K': float(T), 'P_atm': float(pt['P']),
        'V_cm3': float(pt['V']), 'mass_g': float(pt['mass']),
        'T_charge': 'reaction T',
        'initial_solid': 'TiO2',
        'feed_mole_fractions': {g: jnum(pt['y0'][g]) for g in GASES
                                if pt['y0'][g] > 0},
        'n_gas_charge_mol': jnum(pt['n_g0']),
        'n_Ti_mol': jnum(pt['n_Ti']),
        'method': 'gibbs_min',
        'solver': f'element-potential Newton, mpmath dps={mp.dps}',
        'active_condensed_phases': win['active'],
        'species_moles': {s: jnum(n_gas[s]) for s in GASES},
        'solid_moles': {s: jnum(m.get(s, mpf(0))) for s in SOLIDS},
        'gas_fractions': {s: jnum(n_gas[s] / Ngtot) for s in GASES},
        'lambda_kJ_per_mol': {e: jnum(v) for e, v in win['lam'].items()},
        'inactive_phase_reduced_costs': {
            s: {'per_formula_kJ': jnum(rc['per_formula']),
                'per_mol_O_kJ': jnum(rc['per_mol_O']),
                'oxygen_removed_per_extent': rc['oxygen_removed_per_extent']}
            for s, rc in win['r_costs'].items()},
        'reduced_cost_basis': ('per_formula: mu0_j - t_j*lambda_Ti - o_j*lambda_O; '
                               'per_mol_O: divided by def_j = 2*t_j - o_j'),
        'G_total_kJ': jnum(G_tot),
        'common_reference_G_kJ': jnum(G_ref),
        'G_rel_kJ': jnum(G_tot - G_ref),
        'reduced_pct': jnum(100 * ti3 / pt['n_Ti']),
        'x_in_TiO2x': jnum(2 - solid_O / solid_Ti) if solid_Ti > 0 else None,
        'oxygen_to_gas_mol_O': jnum(gasO - pt['bO_gasfeed']),
        'delta_n_gas': {g: jnum(n_gas[g] - pt['n0'][g]) for g in GASES},
        'element_residuals_mol': resid,
        'log10_trace_amounts': traces,
        'boundary_or_degenerate': bool(degenerate) or bool(probe_hits),
        'degenerate_phases': degenerate,
        'triple_probe_hits': probe_hits,
        'solver_tolerance': {'newton_scaled_residual': mp.nstr(TOL_RES, 2),
                             'kkt_kJ': mp.nstr(TOL_KKT, 2),
                             'degeneracy_kJ_per_mol_O': mp.nstr(TOL_DEG, 2)},
        'reporting_detection_limit': ('none at oracle precision; production '
                                      'floors documented in its own record'),
        'converged': True,
        'candidates': [
            {'active': r['active'], 'feasible': r['feasible'],
             'kkt_ok': r['kkt_ok'], 'reason': r['reason'],
             'worst_inactive_r_kJ': (jnum(r['worst_inactive_r'])
                                     if r.get('worst_inactive_r') is not None
                                     else None),
             'solid_moles': ({s: jnum(v) for s, v in r['m'].items()}
                             if r.get('m') else None)}
            for r in results],
    }
    return row


# ------------------------------------------------- boundary data validation


def boundary_block():
    """Waldner ladder, NIST cross-check, and the composition boundary anchor.

    The independent-data leg: NIST-JANAF carries TiO2, Ti4O7, Ti3O5, Ti2O3
    from a different assessment lineage. The two are compared as data against
    data - critical potentials, formation-enthalpy offsets, heat capacities -
    with no solver anywhere in the loop.
    """
    out = {'temperatures_C': TEMPS_C, 'waldner_mu_O_critical_kJ_per_mol_O': {},
           'nist_mu_O_critical_kJ_per_mol_O': {},
           'equilibrium_pO2_atm': {}, 'dHf_offsets_waldner_minus_nist_kJ': {},
           'cp_relative_difference': {}, 'h2_h2o_boundary': {}}
    shared = ['Ti4O7', 'Ti3O5', 'Ti2O3']
    for name in shared + ['TiO2']:
        w = mpf(WALDNER[name]['dHf298']) / 1000
        n = mpf(SHOMATE[name][2])
        out['dHf_offsets_waldner_minus_nist_kJ'][name] = float(w - n)
    for T_C in TEMPS_C:
        T = mpf(str(T_C)) + mpf('273.15')
        muW = {s: mu0_waldner(s, T) for s in SOLIDS}
        muN = {s: mu0_shomate(s, T) for s in shared + ['TiO2']}
        wl, nl, po = {}, {}, {}
        for s in SOLIDS:
            if s == 'TiO2':
                continue
            t_s, o_s = STOICH[s]
            dO = 2 * t_s - o_s
            crit = (t_s * muW['TiO2'] - muW[s]) / dO
            wl[s] = float(crit)
            po[s] = jnum(exp((2 * crit - mu0_shomate('O2', T)) / (R_KJ * T)))
        for s in shared:
            t_s, o_s = STOICH[s]
            dO = 2 * t_s - o_s
            nl[s] = float((t_s * muN['TiO2'] - muN[s]) / dO)
        out['waldner_mu_O_critical_kJ_per_mol_O'][str(T_C)] = wl
        out['nist_mu_O_critical_kJ_per_mol_O'][str(T_C)] = nl
        out['equilibrium_pO2_atm'][str(T_C)] = po
        cp = {}
        for s in shared + ['TiO2']:
            w, n = cp_waldner(s, T), cp_shomate(s, T)
            cp[s] = float(abs(w - n) / n)
        out['cp_relative_difference'][str(T_C)] = cp

    # The composition boundary used by the above/below gate test: at 900 C,
    # the H2O/H2 ratio that sits exactly on the TiO2/Ti10O19 buffer.
    T = mpf('1173.15')
    crit = (10 * mu0_waldner('TiO2', T) - mu0_waldner('Ti10O19', T))  # per O
    ratio = exp((crit - (mu0_shomate('H2O', T) - mu0_shomate('H2', T)))
                / (R_KJ * T))
    out['h2_h2o_boundary'] = {
        'T_C': 900,
        'mu_O_critical_kJ_per_mol_O': float(crit),
        'H2O_over_H2_critical_ratio': jnum(ratio),
        'note': ('feeds H2 + r*H2O at r = ratio*(1 -/+ 1e-3) must land on '
                 'the two-phase buffer / single-phase rutile respectively'),
    }
    return out


def main():
    rows = []
    for case, cfg in CASES.items():
        seeds = {}
        for T_C in TEMPS_C:
            row = solve_point(case, cfg, T_C, seeds)
            rows.append(row)
            print(f"{case:14s} {T_C:5d}C  active={'+'.join(row['active_condensed_phases']):20s} "
                  f"red={row['reduced_pct'] if isinstance(row['reduced_pct'], float) else row['reduced_pct']}"
                  , flush=True)

    # Boundary flip points for the gate test, solved like any other case.
    bb = boundary_block()
    ratio = mpf(bb['h2_h2o_boundary']['H2O_over_H2_critical_ratio']) \
        if isinstance(bb['h2_h2o_boundary']['H2O_over_H2_critical_ratio'], str) \
        else mpf(repr(bb['h2_h2o_boundary']['H2O_over_H2_critical_ratio']))
    for tag, fac, expect in (('h2_h2o_below_boundary', mpf(1) - mpf('1e-3'),
                              'TiO2+Ti10O19'),
                             ('h2_h2o_above_boundary', mpf(1) + mpf('1e-3'),
                              'TiO2')):
        cfg = dict(label=f'H2 with H2O/H2 = crit*{mp.nstr(fac, 6)}',
                   feed={'H2': 1, 'H2O': ratio * fac})
        row = solve_point(tag, cfg, 900, {})
        row['expected_assemblage'] = expect
        rows.append(row)
        print(f"{tag:24s} 900C  active={'+'.join(row['active_condensed_phases'])}",
              flush=True)

    doc = {
        'title': 'High-precision reference: Ti-O gas/solid Gibbs equilibrium',
        'generator': 'oracle_tio.py',
        'precision': f'mpmath, {mp.dps} working digits, {REPORT_DPS} certified',
        'method': 'gibbs_min',
        'formulation': ('element-potential root solve per condensed active '
                        'set; KKT classification over all candidates; '
                        'independent of the production solver code'),
        'scope': ('species-limited C-H-O-N / Ti-O equilibrium. Gases: '
                  'CO2 H2 CO H2O N2 O2. CH4, graphite, nitrides and the '
                  'rutile TiO2-x solution phase are excluded. Finite-'
                  'inventory closed contact, not a flow reactor.'),
        'conditions': {'P_atm': 1.0, 'V_cm3': 25.0, 'mass_g': 0.1,
                       'initial_solid': 'TiO2',
                       'gas_charge': 'n_g0 = P*V/(R_atm*T) at reaction T',
                       'atomic_weights': {k: float(v) for k, v
                                          in ATOMIC_WEIGHT.items()}},
        'gas_data': 'NIST-JANAF Shomate (solidgas/shomate.py tables)',
        'solid_data': ('Waldner & Eriksson, Calphad 23 (1999) 189-218 '
                       '(solidgas/waldner.py tables)'),
        'cases': {k: v['label'] for k, v in CASES.items()},
        'temperatures_C': TEMPS_C,
        'rows': rows,
        'boundary_validation': bb,
    }
    with open('reference_results_high_precision.json', 'w') as fh:
        json.dump(doc, fh, indent=1)
    print(f'\n{len(rows)} points -> reference_results_high_precision.json')
    return 0


if __name__ == '__main__':
    sys.exit(main())
