"""Isolated surface vacancies, as an effective population.

The titration gives the total oxygen removed. It does not give the number
of isolated vacancies at the outermost surface, which is what activates
CO2, and the two do not move together: across the reduction series the
measured inventory rises eightfold while the initial CO rate peaks after
600 C. This is the smallest model that produces that shape.

    n_V,tot = n_iso + n_assoc + n_below

Removal starts at the surface and is generated uniformly, g = n_V,tot/t_red.
A vacancy joins an active region of capacity n_s with probability
f_free = 1 - (n_iso + n_assoc)/n_s and otherwise goes below it; two isolated
vacancies associate at second order.

Only one equation has to be integrated. Writing q = n_iso + n_assoc, the
loss term cancels out of the sum,

    dq/dt = g (1 - q/n_s)      ->   q(t)       = n_s [1 - exp(-g t / n_s)]
                                    n_below(t) = g t - q(t)

so f_free is exp(-g t / n_s) in closed form and what is left is scalar:

    dn_iso/dt = g exp(-g t / n_s) - (2 k_loss / n_s) n_iso^2,   n_iso(0) = 0
    k_loss    = nu_eff exp(-E_loss / kB T)
    n_assoc   = q - n_iso

That matters for more than speed. Integrating all three explicitly is
unstable wherever a search tries a small E_loss with a large nu_eff -
k_loss dt exceeds one and the populations go negative or overflow - and the
closure n_iso + n_assoc + n_below = n_V,tot becomes something to hope for
rather than an identity. Here the closure is exact by construction and the
scalar equation goes to an implicit solver, so the same code answers at the
fitted point and out where the optimiser wanders.

The initial CO rate is taken proportional to n_iso with one per-site
coefficient for every sample, so rates are predicted as ratios normalised
to R600 and no absolute turnover number is claimed.

What this is not. The measured rates are used to infer n_iso, so the model
does not predict them independently, and n_iso must not be reused as a
turnover-frequency denominator. The rates constrain n_iso but say nothing
that separates n_assoc from n_below. E_loss and nu_eff are an effective
temperature dependence for the net loss of isolated vacancies, not an
elementary barrier.
"""

import csv
import math
import os

from scipy.integrate import solve_ivp

KB_EV = 8.617333262e-5
_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SERIES_PATH = os.path.join(_ROOT, 'data', 'reduction_series.csv')

# The active region cannot be smaller than the bridging row nor larger than
# the bridging row plus the first subsurface oxygen layer.
NS_BOUNDS = (13.6, 67.8)
NORMALISE_TO = 'R600'
CLOSURE_TOL = 1e-6

# What the supplement reports, at the digits it was fitted with.
REPORTED = {'ns': 15.33381126, 'E_loss': 2.18450850, 'log10_nu': 6.74333249}

RTOL = 1e-9
ATOL = 1e-11
METHOD = 'Radau'


def load_series(path=SERIES_PATH):
    with open(path) as fh:
        rows = list(csv.DictReader(fh))
    for r in rows:
        for k in ('T_red_C', 't_red_s', 'nV_total_umol_g', 'r_CO_umol_g_s'):
            r[k] = float(r[k])
        r['T_red_K'] = r['T_red_C'] + 273.15
    return rows


def rate_constant(T_K, E_loss, nu_eff):
    """Net loss constant for isolated vacancies, per second."""
    x = -E_loss / (KB_EV * T_K)
    return nu_eff * math.exp(x) if x > -700.0 else 0.0


def integrate(nV_total, t_red_s, T_K, ns, E_loss, nu_eff,
              rtol=RTOL, atol=ATOL, method=METHOD):
    """The three populations at the end of the reduction.

    Only n_iso is integrated; q and n_below are closed form, so the closure
    n_iso + n_assoc + n_below = n_V,tot holds identically.
    """
    if not (nV_total > 0 and t_red_s > 0 and ns > 0):
        raise ValueError('inventory, reduction time and capacity must be positive')
    g = nV_total / t_red_s
    k = rate_constant(T_K, E_loss, nu_eff)
    lam = g / ns
    a = 2.0 * k / ns

    def rhs(t, y):
        return [g * math.exp(-lam * t) - a * y[0] * y[0]]

    def jac(t, y):
        return [[-2.0 * a * y[0]]]

    sol = solve_ivp(rhs, (0.0, t_red_s), [0.0], method=method,
                    jac=jac if method in ('Radau', 'BDF') else None,
                    rtol=rtol, atol=atol)
    if not sol.success:
        raise RuntimeError('the isolated-population solve failed: %s' % sol.message)
    iso = float(sol.y[0][-1])
    q = ns * (1.0 - math.exp(-lam * t_red_s))
    return iso, q - iso, g * t_red_s - q


def partition(series, ns, E_loss, nu_eff, normalise_to=NORMALISE_TO, **kw):
    """Every sample's populations and its predicted CO rate."""
    out = []
    for r in series:
        iso, ass, bel = integrate(r['nV_total_umol_g'], r['t_red_s'],
                                  r['T_red_K'], ns, E_loss, nu_eff, **kw)
        out.append({
            'sample': r['sample'], 'treatment': r['treatment'],
            'T_red_C': r['T_red_C'], 'T_red_K': r['T_red_K'],
            't_red_s': r['t_red_s'],
            'nV_total_umol_g': r['nV_total_umol_g'],
            'r_CO_measured': r['r_CO_umol_g_s'],
            'n_iso_umol_g': iso, 'n_assoc_umol_g': ass, 'n_below_umol_g': bel,
            'closure_umol_g': iso + ass + bel - r['nV_total_umol_g'],
            'below_fraction': bel / r['nV_total_umol_g'],
        })
    ref = [q for q in out if q['sample'] == normalise_to]
    if not ref:
        raise ValueError('no sample named %r to normalise against' % normalise_to)
    ref = ref[0]
    scale = ref['r_CO_measured'] / ref['n_iso_umol_g'] if ref['n_iso_umol_g'] > 0 else 0.0
    for q in out:
        q['r_CO_fitted'] = scale * q['n_iso_umol_g']
        q['residual'] = q['r_CO_fitted'] - q['r_CO_measured']
        q['residual_relative'] = q['residual'] / q['r_CO_measured']
        q['residual_log'] = (math.log(q['r_CO_fitted'] / q['r_CO_measured'])
                             if q['r_CO_fitted'] > 0 else float('inf'))
    return out


# --------------------------------------------------------------- residuals

RESIDUALS = ('log', 'relative', 'raw')


def residuals(rows, kind='log'):
    """The residual vector under one of three definitions.

    The measured rates span 86-fold, so this choice is not cosmetic: raw
    least squares is dominated by R600 and barely sees R1000, while the log
    residual weights every sample by its own size. `log` is the default
    because the model is a statement about ratios. The supplement does not
    record which was used, so every fitted number here is reported under
    all three.
    """
    if kind == 'raw':
        return [q['residual'] for q in rows]
    if kind == 'relative':
        return [q['residual_relative'] for q in rows]
    if kind == 'log':
        return [q['residual_log'] for q in rows]
    raise ValueError('residual kind must be one of %r' % (RESIDUALS,))


def objective(series, ns, E_loss, nu_eff, kind='log', **kw):
    """Sum of squared residuals, penalised where the model is not valid."""
    try:
        rows = partition(series, ns, E_loss, nu_eff, **kw)
    except (ValueError, RuntimeError, OverflowError, FloatingPointError):
        return 1e30
    for q in rows:
        if (q['n_iso_umol_g'] < -1e-9 or q['n_assoc_umol_g'] < -1e-9
                or q['n_below_umol_g'] < -1e-9
                or not math.isfinite(q['n_iso_umol_g'])
                or abs(q['closure_umol_g']) > CLOSURE_TOL):
            return 1e30
    r = residuals(rows, kind)
    if any(not math.isfinite(v) for v in r):
        return 1e30
    return float(sum(v * v for v in r))


# ------------------------------------------------------------------- the fit

# Search tolerances. The optimiser only has to rank candidates, and Radau at
# 1e-7 ranks them identically to 1e-9 while costing a fifth as much; every
# number that gets reported is recomputed at RTOL/ATOL afterwards.
SEARCH_RTOL = 1e-7
SEARCH_ATOL = 1e-9


def _nelder_mead(cost, start, steps, iters=600, tol=1e-12):
    n = len(start)
    pts = [list(start)]
    for i in range(n):
        p = list(start)
        p[i] += steps[i]
        pts.append(p)
    vals = [cost(p) for p in pts]
    for _ in range(iters):
        order = sorted(range(n + 1), key=lambda i: vals[i])
        pts = [pts[i] for i in order]
        vals = [vals[i] for i in order]
        if abs(vals[-1] - vals[0]) <= tol * (abs(vals[0]) + tol):
            break
        cen = [sum(p[j] for p in pts[:-1]) / n for j in range(n)]
        refl = [cen[j] + (cen[j] - pts[-1][j]) for j in range(n)]
        fr = cost(refl)
        if fr < vals[0]:
            exp_ = [cen[j] + 2.0 * (cen[j] - pts[-1][j]) for j in range(n)]
            fe = cost(exp_)
            pts[-1], vals[-1] = (exp_, fe) if fe < fr else (refl, fr)
        elif fr < vals[-2]:
            pts[-1], vals[-1] = refl, fr
        else:
            con = [cen[j] + 0.5 * (pts[-1][j] - cen[j]) for j in range(n)]
            fc = cost(con)
            if fc < vals[-1]:
                pts[-1], vals[-1] = con, fc
            else:
                for i in range(1, n + 1):
                    pts[i] = [pts[0][j] + 0.5 * (pts[i][j] - pts[0][j])
                              for j in range(n)]
                    vals[i] = cost(pts[i])
    best = min(range(n + 1), key=lambda i: vals[i])
    return pts[best], vals[best]


def _unpack(p, ns_bounds):
    """(log ns, E_loss, log10 nu) -> the physical triple, ns held in bounds."""
    lo, hi = ns_bounds
    ns = min(max(math.exp(p[0]), lo), hi)
    return ns, max(p[1], 0.0), 10.0 ** p[2]


def fit(series, kind='log', ns_bounds=NS_BOUNDS, start=None, iters=600):
    """Fit (n_s, E_loss, nu_eff) under one residual definition.

    The search coordinates are log n_s, E_loss and log10 nu_eff, as the
    supplement's reviewer asked: n_s is a positive scale, and it is the
    logarithm of the prefactor that trades against the barrier, not the
    prefactor itself.
    """
    lo, hi = ns_bounds
    if start is None:
        start = [math.log(20.0), 2.2, 6.7]

    def cost(p):
        ns, E, nu = _unpack(p, ns_bounds)
        return objective(series, ns, E, nu, kind=kind,
                         rtol=SEARCH_RTOL, atol=SEARCH_ATOL)

    best, val = _nelder_mead(cost, start, [0.35, 0.25, 0.6], iters=iters)
    ns, E, nu = _unpack(best, ns_bounds)
    return {'residual_kind': kind,
            'ns_umol_g': ns, 'E_loss_eV': E, 'nu_eff_s1': nu,
            'log10_nu_eff': math.log10(nu),
            'objective': objective(series, ns, E, nu, kind=kind),
            'ns_at_bound': abs(ns - lo) < 1e-6 or abs(ns - hi) < 1e-6,
            'ns_bounds': list(ns_bounds)}


def fit_at_fixed_nu(series, nu_eff, kind='log', ns_bounds=NS_BOUNDS, iters=400):
    """The same fit with the prefactor pinned - the identifiability probe."""
    lo, hi = ns_bounds

    def cost(p):
        ns = min(max(math.exp(p[0]), lo), hi)
        return objective(series, ns, max(p[1], 0.0), nu_eff, kind=kind,
                         rtol=SEARCH_RTOL, atol=SEARCH_ATOL)

    best, _v = None, float('inf')
    for ns0 in (14.0, 20.0, 35.0, 60.0):
        for E0 in (2.0, 3.5, 5.0):
            p, v = _nelder_mead(cost, [math.log(ns0), E0], [0.3, 0.4], iters=iters)
            if v < _v:
                best, _v = p, v
    ns = min(max(math.exp(best[0]), lo), hi)
    E = max(best[1], 0.0)
    return {'residual_kind': kind, 'ns_umol_g': ns, 'E_loss_eV': E,
            'nu_eff_s1': nu_eff,
            'objective': objective(series, ns, E, nu_eff, kind=kind)}


def rate_ratio(series, ns, E_loss, nu_eff, hi='R800', lo='R600', **kw):
    """One sample's predicted rate over another's."""
    part = {q['sample']: q for q in partition(series, ns, E_loss, nu_eff, **kw)}
    return part[hi]['r_CO_fitted'] / part[lo]['r_CO_fitted']


def measured_ratio(series, hi='R800', lo='R600'):
    by = {r['sample']: r['r_CO_umol_g_s'] for r in series}
    return by[hi] / by[lo]


# ------------------------------------------------------------ identifiability

def correlation(series, ns, E_loss, nu_eff, kind='log',
                rel_step=(0.02, 0.01, 0.01)):
    """Correlations from a finite-difference Gauss-Newton covariance.

    This is the Hessian estimate and it is named as such: a profile or a
    resampling estimate of the same correlations is a different object and
    can give a different number. The supplement does not say which it
    quotes, so this reports the method alongside the value.
    """
    p0 = [ns, E_loss, math.log10(nu_eff)]
    h = [abs(p0[i]) * rel_step[i] or rel_step[i] for i in range(3)]

    def resid(p):
        return residuals(partition(series, p[0], p[1], 10.0 ** p[2]), kind)

    J = []
    for i in range(3):
        up, dn = list(p0), list(p0)
        up[i] += h[i]
        dn[i] -= h[i]
        ru, rd = resid(up), resid(dn)
        J.append([(a - b) / (2 * h[i]) for a, b in zip(ru, rd)])

    JTJ = [[sum(J[a][k] * J[b][k] for k in range(len(J[0]))) for b in range(3)]
           for a in range(3)]
    cov = _invert3(JTJ)
    if cov is None:
        return None
    rho = [[cov[a][b] / math.sqrt(cov[a][a] * cov[b][b]) for b in range(3)]
           for a in range(3)]
    return {'method': 'Gauss-Newton Hessian, central differences',
            'residual_kind': kind,
            'order': ['ns', 'E_loss', 'log10_nu_eff'], 'matrix': rho,
            'rho_E_lognu': rho[1][2], 'rho_ns_E': rho[0][1],
            'rho_ns_lognu': rho[0][2]}


def _invert3(m):
    a, b, c = m[0]
    d, e, f = m[1]
    g, h, i = m[2]
    det = a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g)
    if abs(det) < 1e-300:
        return None
    return [[(e * i - f * h) / det, (c * h - b * i) / det, (b * f - c * e) / det],
            [(f * g - d * i) / det, (a * i - c * g) / det, (c * d - a * f) / det],
            [(d * h - e * g) / det, (b * g - a * h) / det, (a * e - b * d) / det]]


# -------------------------------------------------------------- the step rule

def partition_step(series, ns, E_loss, nu_eff, normalise_to=NORMALISE_TO, **kw):
    """The alternative the supplement rejects, written out so it can be tested.

    f_free is 1 until the active region is full and 0 afterwards, instead of
    falling linearly. q then fills at rate g up to n_s and stops:

        q(t) = min(g t, n_s),   n_below(t) = g t - q(t)

    and the isolated population obeys the same scalar equation with that
    f_free in place of the exponential.
    """
    out = []
    for r in series:
        nV, t_red, T_K = r['nV_total_umol_g'], r['t_red_s'], r['T_red_K']
        g = nV / t_red
        k = rate_constant(T_K, E_loss, nu_eff)
        a = 2.0 * k / ns
        t_full = ns / g

        def rhs(t, y):
            return [(g if t < t_full else 0.0) - a * y[0] * y[0]]

        sol = solve_ivp(rhs, (0.0, t_red), [0.0], method='Radau',
                        rtol=kw.get('rtol', RTOL), atol=kw.get('atol', ATOL),
                        max_step=max(t_full / 8.0, t_red / 2000.0))
        if not sol.success:
            raise RuntimeError(sol.message)
        iso = float(sol.y[0][-1])
        q = min(g * t_red, ns)
        out.append({'sample': r['sample'], 'n_iso_umol_g': iso,
                    'n_assoc_umol_g': q - iso,
                    'n_below_umol_g': g * t_red - q,
                    'nV_total_umol_g': nV,
                    'r_CO_measured': r['r_CO_umol_g_s']})
    ref = [q for q in out if q['sample'] == normalise_to][0]
    scale = ref['r_CO_measured'] / ref['n_iso_umol_g'] if ref['n_iso_umol_g'] > 0 else 0.0
    for q in out:
        q['r_CO_fitted'] = scale * q['n_iso_umol_g']
        q['residual'] = q['r_CO_fitted'] - q['r_CO_measured']
        q['residual_relative'] = q['residual'] / q['r_CO_measured']
        q['residual_log'] = (math.log(q['r_CO_fitted'] / q['r_CO_measured'])
                             if q['r_CO_fitted'] > 0 else float('inf'))
    return out
