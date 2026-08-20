"""Production gas/solid Gibbs equilibrium: condensed phases by active set.

This is the solver the page runs. The rules it exists to enforce:

  * gas species live in logarithmic coordinates and are never floored - a
    species with zero elemental inventory is removed outright (its amount is
    exactly zero by balance), everything else stays strictly positive;
  * condensed phases are enumerated as active sets: a phase left out of the
    candidate is not a variable at all, so its mole number is exactly 0.0 in
    the result, not a residue at some floor;
  * every inactive phase gets a KKT reduced cost at the solution, in two
    bases (per formula unit formed, and per mole of oxygen moved), and the
    stability decision is made by those signs - never by the size of a mole
    number, and never by comparing total-G values that double precision
    cannot resolve;
  * all candidates of one charge share one Gibbs reference: the loaded,
    unreacted feed. G_rel below double resolution is flagged, not printed as
    if it were a measurement.

The algorithm is the element-potential form of the same convex program the
legacy minimisers attacked head-on. For a candidate active set A:

    gases          ln y_i = sum_e a_ie theta_e - mu0_i/RT - ln P
    active solids  t_j lambda_Ti + o_j lambda_O = mu0_j
    balances       per element; the solids enter only Ti and O, and O is used
                   as the combination (O - rho_host * Ti), which removes the
                   host oxygen symbolically instead of by subtraction.

|A| = 2 pins (lambda_Ti, lambda_O) as a 2x2 linear solve - the two-phase
buffer. |A| = 1 leaves lambda_O to the gas. Newton runs only over the gas
unknowns (at most five), and the solid amounts are recovered linearly, so a
condensed phase present at 1e-41 mol costs the optimiser nothing.

Because the program is convex, a candidate whose solution is feasible and
whose inactive reduced costs are all nonnegative IS the global minimum; the
enumeration cannot pick a plausible-but-wrong basin. Degeneracies (a reduced
cost at zero within tolerance) are reported as boundaries, not resolved by
fiat.

The legacy solvers (`gibbs.GibbsSystem.minimise`, `equilibrium.solve`) are
kept for regression only and are not imported here.
"""

import math
from itertools import combinations

from .shomate import SHOMATE, BREAKPOINT, DEFAULT_BREAKPOINT, R_KJ, R_ATM
from . import waldner

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
TI3 = {p: (0 if p == 'TiO2' else 2) for p in SOLIDS}

ATOMIC_WEIGHT = {'Ti': 47.867, 'O': 15.999}

# Classification tolerances. KKT signs decide stability; sizes only display.
TOL_NEWTON = 1e-13          # scaled balance residual
TOL_KKT_KJ = 1e-9           # inactive r_j below -this fails the candidate
TOL_DEG_KJ = 1e-6           # |r_j| per mol O below this flags a boundary
TRACE_MOL = 1e-12           # below this, log10(n) is reported alongside
EXP_CAP = 700.0             # exp() guard; beyond it the step is damped back

MODE_NOTE = ('C-H-O-N / Ti-O Gibbs equilibrium; gas registry in this '
             'version: CO2 H2 CO H2O N2 O2; CH4, graphite, nitrides and the '
             'rutile TiO2-x solution phase are not included in this version; '
             'closed finite-inventory contact, not a flow reactor')


def mu0_gas(name, T):
    """NIST Shomate G(T), kJ/mol."""
    lo, hi, dHf = SHOMATE[name]
    A, B, C, D, E, F, G, H = (lo if T <= BREAKPOINT.get(name, DEFAULT_BREAKPOINT)
                              else hi)
    t = T / 1000.0
    Ht = dHf + (A * t + B * t * t / 2 + C * t**3 / 3 + D * t**4 / 4
                - E / t + F - H)
    St = (A * math.log(t) + B * t + C * t * t / 2 + D * t**3 / 3
          - E / (2 * t * t) + G)
    return Ht - T * St / 1000.0


def mu0_solid(name, T):
    """Waldner & Eriksson G(T), kJ/mol."""
    return waldner.mu0(name, T)


def mole_weight(solid):
    t, o = STOICH[solid]
    return t * ATOMIC_WEIGHT['Ti'] + o * ATOMIC_WEIGHT['O']


# ----------------------------------------------------------------- the charge


def build_charge(feed_rel, T_K, P_atm=1.0, V_cm3=25.0, mass_g=0.100,
                 initial_solid='TiO2', T_charge_K=None):
    """Feed moles and element inventory for one condition.

    `feed_rel` maps gas species to nonnegative relative amounts; it is
    normalised here and the normalised fractions are part of the record, so
    nothing is rescaled silently. n_g0 = P*V/(R*T_charge); T_charge defaults
    to the reaction temperature (sealed-vessel convention).
    """
    tot = sum(float(v) for v in feed_rel.values())
    if tot <= 0:
        raise ValueError('gas feed sums to zero')
    if mass_g <= 0:
        raise ValueError('solid mass must be positive')
    if initial_solid not in STOICH:
        raise ValueError(f'unknown initial solid {initial_solid}')
    for k, v in feed_rel.items():
        if k not in GAS_FORMULA:
            raise ValueError(f'unknown gas species {k}')
        if float(v) < 0:
            raise ValueError(f'negative feed amount for {k}')
    Tc = T_K if T_charge_K is None else T_charge_K
    n_g0 = P_atm * V_cm3 / (R_ATM * Tc)
    y0 = {g: float(feed_rel.get(g, 0.0)) / tot for g in GASES}
    n0 = {g: y0[g] * n_g0 for g in GASES}
    t_h, o_h = STOICH[initial_solid]
    n_solid0 = mass_g / mole_weight(initial_solid)
    n_Ti = n_solid0 * t_h
    rho = o_h / t_h                      # host O per Ti; the O' combination
    b = {'C': 0.0, 'H': 0.0, 'N': 0.0, 'O': 0.0, 'Ti': n_Ti}
    bO_gas = 0.0
    for g in GASES:
        for e, k in GAS_FORMULA[g].items():
            b[e] += k * n0[g]
    bO_gas = sum(GAS_FORMULA[g].get('O', 0) * n0[g] for g in GASES)
    b['O'] += o_h * n_solid0
    return dict(T=T_K, T_charge=Tc, P=P_atm, V=V_cm3, mass=mass_g,
                initial_solid=initial_solid, n_solid0=n_solid0,
                n_g0=n_g0, y0=y0, n0=n0, n_Ti=n_Ti, rho=rho, b=b,
                bO_gasfeed=bO_gas)


# ------------------------------------------------------------ one active set


def _solve_candidate(active, ch, mu, registry, seed=None):
    """Square Newton over the gas; solids recovered linearly. Pure doubles.

    `registry` is the full list of condensed phases allowed in this solve;
    reduced costs are reported for every registry member left inactive.
    """
    T, P, b = ch['T'], ch['P'], ch['b']
    RT = R_KJ * T
    rho = ch['rho']
    tj = {s: float(STOICH[s][0]) for s in active}
    oj = {s: float(STOICH[s][1]) for s in active}
    dj = {s: rho * tj[s] - oj[s] for s in active}   # O freed per formula

    lam = {}
    if len(active) >= 2:
        a1, a2 = active[0], active[1]
        det = tj[a1] * oj[a2] - tj[a2] * oj[a1]
        lam['Ti'] = (mu[a1] * oj[a2] - mu[a2] * oj[a1]) / det
        lam['O'] = (tj[a1] * mu[a2] - tj[a2] * mu[a1]) / det

    have = {e: b[e] > 0 for e in ('C', 'H', 'N')}
    if len(active) == 1:
        s0 = active[0]
        m_fixed = b['Ti'] / tj[s0]
        avail_O = ch['bO_gasfeed'] + dj[s0] * m_fixed
        o_free = avail_O > 0
    else:
        m_fixed, avail_O, o_free = None, None, False

    def elem_ok(e):
        if e in ('C', 'H', 'N'):
            return have[e]
        return o_free or len(active) >= 2      # e == 'O'

    species = [g for g in GASES if all(elem_ok(e) for e in GAS_FORMULA[g])]
    if not species:
        return dict(active=list(active), feasible=False, kkt_ok=False,
                    reason='no gas species with inventory')

    elems = [e for e in ('C', 'H', 'N') if have[e]]
    if o_free:
        elems.append('O')
    nE = len(elems)

    gi = {}
    for s in species:
        v = mu[s] / RT + math.log(P)
        if len(active) >= 2 and 'O' in GAS_FORMULA[s]:
            v -= GAS_FORMULA[s]['O'] * lam['O'] / RT
        gi[s] = v

    rhs = {e: b[e] for e in ('C', 'H', 'N')}
    if o_free:
        rhs['O'] = avail_O
    scale = {e: max(rhs[e], 1e-300) for e in elems}

    if seed is not None and len(seed) == nE + 1:
        x = list(seed)
    else:
        yg = {s: max(ch['y0'].get(s, 0.0), 0.0) for s in species}
        tot = sum(yg.values())
        for s in species:
            yg[s] = (yg[s] / tot) if tot > 0 else 1.0 / len(species)
            yg[s] = max(yg[s], 1e-8)
        tot = sum(yg.values())
        for s in species:
            yg[s] /= tot
        # theta seed: normal equations of ln y + g = A theta
        AtA = [[0.0] * nE for _ in range(nE)]
        Atr = [0.0] * nE
        for s in species:
            r = math.log(yg[s]) + gi[s]
            for i, e in enumerate(elems):
                ae = GAS_FORMULA[s].get(e, 0)
                if not ae:
                    continue
                Atr[i] += ae * r
                for k, f in enumerate(elems):
                    AtA[i][k] += ae * GAS_FORMULA[s].get(f, 0)
        th = _lin_solve(AtA, Atr) or [0.0] * nE
        x = th + [math.log(ch['n_g0'])]

    def gas_y(xv):
        y = {}
        for s in species:
            v = -gi[s]
            for i, e in enumerate(elems):
                v += GAS_FORMULA[s].get(e, 0) * xv[i]
            y[s] = math.exp(min(v, EXP_CAP))
        return y

    def residual(xv):
        y = gas_y(xv)
        Ng = math.exp(min(xv[nE], EXP_CAP))
        res = []
        for e in elems:
            tot = sum(GAS_FORMULA[s].get(e, 0) * y[s] for s in species)
            res.append((tot * Ng - rhs[e]) / scale[e])
        res.append(sum(y.values()) - 1.0)
        return res, y, Ng

    res, y, Ng = residual(x)
    rn = max(abs(v) for v in res)
    it = 0
    while rn > TOL_NEWTON and it < 200:
        J = [[0.0] * (nE + 1) for _ in range(nE + 1)]
        for i, e in enumerate(elems):
            for k, f in enumerate(elems):
                J[i][k] = sum(GAS_FORMULA[s].get(e, 0) * GAS_FORMULA[s].get(f, 0)
                              * y[s] for s in species) * Ng / scale[e]
            J[i][nE] = sum(GAS_FORMULA[s].get(e, 0) * y[s]
                           for s in species) * Ng / scale[e]
        for k, f in enumerate(elems):
            J[nE][k] = sum(GAS_FORMULA[s].get(f, 0) * y[s] for s in species)
        J[nE][nE] = 0.0
        step = _lin_solve(J, [-v for v in res])
        if step is None:
            return dict(active=list(active), feasible=False, kkt_ok=False,
                        reason='singular Jacobian')
        damp = 1.0
        for _ in range(60):
            xn = [x[i] + damp * step[i] for i in range(nE + 1)]
            resn, yn, Ngn = residual(xn)
            rnn = max(abs(v) for v in resn)
            if rnn < rn:
                break
            damp /= 2
        else:
            return dict(active=list(active), feasible=False, kkt_ok=False,
                        reason=f'stalled at residual {rn:.2e}')
        x, res, y, Ng, rn = xn, resn, yn, Ngn, rnn
        it += 1
    if rn > TOL_NEWTON:
        return dict(active=list(active), feasible=False, kkt_ok=False,
                    reason=f'no convergence, residual {rn:.2e}')

    n_gas = {g: 0.0 for g in GASES}
    for s in species:
        n_gas[s] = y[s] * Ng

    gasO = sum(GAS_FORMULA[s].get('O', 0) * n_gas[s] for s in species)
    if len(active) == 1:
        m = {active[0]: m_fixed}
    else:
        a1, a2 = active[0], active[1]
        det = tj[a1] * dj[a2] - tj[a2] * dj[a1]
        if det == 0:
            return dict(active=list(active), feasible=False, kkt_ok=False,
                        reason='degenerate stoichiometry pair')
        rO = gasO - ch['bO_gasfeed']
        m = {a1: (b['Ti'] * dj[a2] - tj[a2] * rO) / det,
             a2: (tj[a1] * rO - dj[a1] * b['Ti']) / det}

    if o_free:
        lam['O'] = RT * x[elems.index('O')]
    if len(active) == 1 and 'O' in lam:
        s0 = active[0]
        lam['Ti'] = (mu[s0] - oj[s0] * lam['O']) / tj[s0]
    for i, e in enumerate(elems):
        lam[e] = RT * x[i]

    r_costs = {}
    lam_O_def = 'O' in lam
    for s in registry:
        if s in active:
            continue
        t_s, o_s = STOICH[s]
        d_s = ch['rho'] * t_s - o_s
        if not lam_O_def:
            r_costs[s] = dict(per_formula=(-math.inf if d_s > 0 else None),
                              per_mol_O=(-math.inf if d_s > 0 else None),
                              dO=d_s)
            continue
        rf = mu[s] - t_s * lam['Ti'] - o_s * lam['O']
        r_costs[s] = dict(per_formula=rf,
                          per_mol_O=(rf / d_s if d_s != 0 else None), dO=d_s)

    feasible = all(v is not None and v > 0 for v in m.values())
    reason = 'ok' if feasible else (
        'negative condensed amount: '
        + ', '.join(f'{s}={m[s]:.3e}' for s in active if m[s] is not None))

    kkt_ok = feasible
    worst = None
    for s, rc in r_costs.items():
        v = rc['per_formula']
        if v is None:
            continue
        if worst is None or v < worst:
            worst = v
        if v < -TOL_KKT_KJ:
            kkt_ok = False

    return dict(active=list(active), feasible=feasible, kkt_ok=kkt_ok,
                reason=reason, y={s: y[s] for s in species}, n_gas=n_gas,
                m=m, lam=lam, r_costs=r_costs, worst_inactive_r=worst,
                iterations=it, newton_residual=rn,
                seed_out=list(x))


def _lin_solve(A, rhs):
    """Gaussian elimination with partial pivoting; None on singularity."""
    n = len(rhs)
    M = [row[:] + [rhs[i]] for i, row in enumerate(A)]
    for c in range(n):
        p = max(range(c, n), key=lambda r: abs(M[r][c]))
        if abs(M[p][c]) < 1e-300:
            return None
        M[c], M[p] = M[p], M[c]
        piv = M[c][c]
        for r in range(n):
            if r == c:
                continue
            f = M[r][c] / piv
            if f:
                for k in range(c, n + 1):
                    M[r][k] -= f * M[c][k]
    return [M[r][n] / M[r][r] for r in range(n)]


# --------------------------------------------------------------- full solve


def enumerate_candidates(solids):
    """Singles and pairs, in a canonical order independent of input order."""
    ss = sorted(solids, key=lambda s: (-STOICH[s][1] / STOICH[s][0], s))
    return [[s] for s in ss] + [list(c) for c in combinations(ss, 2)]


def triple_probes(solids, mu, tol=TOL_DEG_KJ):
    """Three line compounds share (lambda_Ti, lambda_O) only at special T."""
    hits = []
    ss = sorted(solids, key=lambda s: (-STOICH[s][1] / STOICH[s][0], s))
    for cand in combinations(ss, 3):
        a1, a2, a3 = cand
        t1, o1 = STOICH[a1]
        t2, o2 = STOICH[a2]
        t3, o3 = STOICH[a3]
        det = t1 * o2 - t2 * o1
        lam_Ti = (mu[a1] * o2 - mu[a2] * o1) / det
        lam_O = (t1 * mu[a2] - t2 * mu[a1]) / det
        res = mu[a3] - t3 * lam_Ti - o3 * lam_O
        if abs(res) < tol:
            hits.append({'phases': list(cand), 'residual_kJ': res})
    return hits


def _g_absolute(n_gas, m, mu, T, P):
    RT = R_KJ * T
    Ng = sum(v for v in n_gas.values() if v > 0)
    G = 0.0
    for s, n in n_gas.items():
        if n > 0:
            G += n * (mu[s] + RT * math.log(n / Ng * P))
    for s, v in m.items():
        G += v * mu[s]
    return G


def solve(feed_rel, T_K, P_atm=1.0, V_cm3=25.0, mass_g=0.100,
          initial_solid='TiO2', T_charge_K=None, solids=None):
    """The full answer for one condition. Every value is gibbs_min.

    Returns a record carrying the winning active set, exact zeros for every
    excluded phase, KKT reduced costs in both bases for every inactive phase,
    the shared Gibbs reference, per-element balance residuals, and the
    boundary/degeneracy flags. Raises if no candidate satisfies KKT.
    """
    ch = build_charge(feed_rel, T_K, P_atm, V_cm3, mass_g,
                      initial_solid, T_charge_K)
    # Canonical phase order at the door: every output dict, candidate list
    # and mole vector is then independent of how the caller ordered the input.
    solids = sorted(solids or SOLIDS,
                    key=lambda s: (-STOICH[s][1] / STOICH[s][0], s))
    mu = {g: mu0_gas(g, T_K) for g in GASES}
    for s in solids:
        mu[s] = mu0_solid(s, T_K)

    results = []
    passing = []
    for cand in enumerate_candidates(solids):
        r = _solve_candidate(cand, ch, mu, solids)
        results.append(r)
        if r['feasible'] and r['kkt_ok']:
            passing.append(r)

    if not passing:
        raise RuntimeError('no candidate active set satisfies KKT: '
                           + '; '.join(f"{'+'.join(r['active'])}: {r['reason']}"
                                       for r in results))
    passing.sort(key=lambda r: len(r['active']))
    win = passing[0]

    degenerate = [s for s, rc in win['r_costs'].items()
                  if rc['per_mol_O'] is not None
                  and not math.isinf(rc['per_mol_O'])
                  and abs(rc['per_mol_O']) < TOL_DEG_KJ]
    probes = triple_probes(solids, mu)

    m = {s: 0.0 for s in solids}             # excluded phases: exactly zero
    for s, v in win['m'].items():
        m[s] = v
    n_gas = win['n_gas']

    G_tot = _g_absolute(n_gas, m, mu, T_K, P_atm)
    RT = R_KJ * T_K
    G_ref = ch['n_solid0'] * mu[initial_solid] if initial_solid in mu else \
        ch['n_solid0'] * mu0_solid(initial_solid, T_K)
    for s, n in ch['n0'].items():
        if n > 0:
            G_ref += n * (mu[s] + RT * math.log(n / ch['n_g0'] * P_atm))
    G_rel = G_tot - G_ref
    g_floor = 1e-12 * max(abs(G_tot), abs(G_ref))

    balance = {}
    for e in ('C', 'H', 'N', 'O', 'Ti'):
        tot = sum(GAS_FORMULA[g].get(e, 0) * n_gas[g] for g in GASES)
        for s, v in m.items():
            t_s, o_s = STOICH[s]
            tot += (t_s if e == 'Ti' else o_s if e == 'O' else 0) * v
        balance[e] = tot - ch['b'][e]

    solid_O = sum(STOICH[s][1] * v for s, v in m.items())
    solid_Ti = sum(STOICH[s][0] * v for s, v in m.items())
    ti3 = sum(TI3[s] * v for s, v in m.items())
    gasO = sum(GAS_FORMULA[g].get('O', 0) * n_gas[g] for g in GASES)
    Ngtot = sum(n_gas.values())

    traces = {}
    for s, v in list(n_gas.items()) + list(m.items()):
        if 0 < v < TRACE_MOL:
            traces[s] = math.log10(v)

    return {
        'system': 'Ti-O',
        'mode': MODE_NOTE,
        'method': 'gibbs_min',
        'solver': 'active-set element-potential Newton, double precision',
        'T_C': T_K - 273.15, 'T_K': T_K, 'P_atm': P_atm, 'V_cm3': V_cm3,
        'T_charge_K': ch['T_charge'], 'mass_g': mass_g,
        'initial_solid': initial_solid,
        'feed_mole_fractions': {g: ch['y0'][g] for g in GASES if ch['y0'][g] > 0},
        'n_gas_charge_mol': ch['n_g0'], 'n_Ti_mol': ch['n_Ti'],
        'active_condensed_phases': win['active'],
        'species_moles': dict(n_gas),
        'solid_moles': m,
        'gas_fractions': {g: (n_gas[g] / Ngtot if Ngtot > 0 else 0.0)
                          for g in GASES},
        'phase_split_pct': {s: (STOICH[s][0] * v / ch['n_Ti'] * 100.0)
                            for s, v in m.items()},
        'lambda_kJ_per_mol': dict(win['lam']),
        'inactive_phase_reduced_costs': {
            s: {'per_formula_kJ': rc['per_formula'],
                'per_mol_O_kJ': rc['per_mol_O'],
                'oxygen_removed_per_extent': rc['dO']}
            for s, rc in win['r_costs'].items()},
        'reduced_cost_basis': ('per_formula: mu0_j - t_j*lambda_Ti - '
                               'o_j*lambda_O; per_mol_O: divided by '
                               'def_j = rho_host*t_j - o_j'),
        'G_total_kJ': G_tot,
        'common_reference_G_kJ': G_ref,
        'G_rel_kJ': G_rel,
        'G_rel_below_double_resolution': abs(G_rel) < g_floor,
        'reduced_pct': 100.0 * ti3 / ch['n_Ti'],
        'x_in_TiO2x': (2.0 - solid_O / solid_Ti) if solid_Ti > 0 else None,
        'oxygen_to_gas_mol_O': gasO - ch['bO_gasfeed'],
        'conversion_CO2_pct': ((ch['n0']['CO2'] - n_gas['CO2'])
                               / ch['n0']['CO2'] * 100.0
                               if ch['n0']['CO2'] > 0 else None),
        'delta_n_gas': {g: n_gas[g] - ch['n0'][g] for g in GASES},
        'balance_by_element': balance,
        'log10_trace_amounts': traces,
        'boundary_or_degenerate': bool(degenerate) or bool(probes),
        'degenerate_phases': degenerate,
        'triple_probe_hits': [{'phases': p['phases'],
                               'residual_kJ': p['residual_kJ']}
                              for p in probes],
        'solver_tolerance': {'newton_scaled_residual': TOL_NEWTON,
                             'kkt_kJ': TOL_KKT_KJ,
                             'degeneracy_kJ_per_mol_O': TOL_DEG_KJ},
        'reporting_detection_limit': {
            'moles': 'exact zeros for excluded phases; doubles elsewhere',
            'G_rel_kJ_floor': g_floor},
        'converged': True,
        'iterations': win['iterations'],
        'newton_residual': win['newton_residual'],
        'candidates': [{'active': r['active'], 'feasible': r['feasible'],
                        'kkt_ok': r['kkt_ok'], 'reason': r['reason'],
                        'worst_inactive_r_kJ': r.get('worst_inactive_r'),
                        'solid_moles': r.get('m')}
                       for r in results],
    }
