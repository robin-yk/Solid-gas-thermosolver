"""Coarse-grained KMC of reoxidation transport, after Katsoulakis & Vlachos.

Radial coarse cells over one spherical particle: cell 0 is the surface
monolayer (capacity capped at the reconstruction threshold), cell 1 the
first subsurface O layer, cells 2..m-1 equal-width bulk shells from the
outside inward. A vacancy exchange between adjacent cells k -> l fires at

    R(k->l) = Lambda_kl * eta_k * (q_l - eta_l) * exp(beta * eps_k),

with the symmetric kernel

    Lambda_kl = nu0 * a^2 * A_kl * exp(-beta * E_TS) / (d_kl * v_O * q_k * q_l),

and a common saddle E_TS = E_m + eps_bulk, so the activation energy to
leave cell k is E_TS - eps_k. Two properties pin the construction:

1. Detailed balance holds exactly against the product-binomial Gibbs
   measure with H = sum_k eps_k eta_k (the generalised-exclusion identity
   eta_k (q_l - eta_l) rho(eta_k) rho(eta_l) =
   (eta_l + 1)(q_k - eta_k + 1) rho(eta_k - 1) rho(eta_l + 1)), so the
   closed system relaxes to the same Fermi-Dirac split over site classes
   that the equilibrium workspace computes - a cross-check the tests run.
2. In the dilute limit the mean equation is the finite-volume
   discretisation of spherical diffusion with D = a^2 nu0 exp(-beta E_m),
   which pins the kernel's magnitude (and carries the framework's q^2
   time coarse-graining: per-vacancy escape scales as (a/W)^2).

Surface refill by CO2 is the one irreversible channel (open system):
eta_0 -> eta_0 - 1 at k_fill = A_s p^m exp(-beta Ea_s) per vacancy.

Deterministic mode integrates the mean master equation with an adaptive
exponential-Euler stepper whose fixed point is exact balance; stochastic
mode tau-leaps a scaled-down column with the same rates and a seeded RNG
shared with the lattice engine. cgmc.js is a line-for-line port.
"""

import math

from . import statmech as base

KB_EV = base.KB_EV
N_AVO = base.N_AVO


def _kfill(p, kin):
    d = base.co2_kinetics_params(p, kin)
    kt = KB_EV * (d['T_reox_C'] + 273.15)
    pf = d['p_CO2_atm'] ** d['p_order_m']
    return (d['A_eff_per_s_atm']['surface'] * pf
            * math.exp(-d['Ea_eV']['surface'] / kt), d)


def build(p, vo_total=None, d_um=None, dist_fractions=None, eps=None,
          kin=None, cg=None):
    """Radial cell system for one particle; all counts are absolute."""
    cgp = dict(p['kinetics_cgmc'])
    if cg:
        cgp.update(cg)
    if vo_total is None:
        vo_total = p['defaults']['VO_total_umol_g']
    if d_um is None:
        d_um = p['defaults']['particle_diameter_um']
    if eps is None:
        eps = base.energetics(p)['eps']
    a_ang = p['lattice_constants_A']['a']
    c_ang = p['lattice_constants_A']['c']
    a_nm = a_ang / math.sqrt(2.0) / 10.0            # d110 layer spacing
    n_o = 4.0 / (a_ang * a_ang * c_ang) * 1000.0    # O sites per nm^3
    v_o = 1.0 / n_o
    sig_b = 1.0 / (c_ang * a_ang * math.sqrt(2.0)) * 100.0   # per nm^2
    sig_l = 4.0 * sig_b
    r_nm = d_um * 1000.0 / 2.0
    kf, kin_used = _kfill(p, kin)
    kt = KB_EV * (kin_used['T_reox_C'] + 273.15)
    beta = 1.0 / kt

    m = max(6, int(cgp['n_cells']))
    n_bulk_cells = m - 2
    r_core = r_nm - 2.0 * a_nm
    w = r_core / n_bulk_cells

    q = [4.0 * math.pi * r_nm * r_nm * sig_b,
         4.0 * math.pi * (r_nm - a_nm) * (r_nm - a_nm) * sig_l]
    centers = [r_nm - 0.5 * a_nm, r_nm - 1.5 * a_nm]
    eps_c = [eps[0], eps[1]]
    for i in range(n_bulk_cells):                   # outermost first
        r_out = r_core - i * w
        r_in = r_out - w
        q.append(n_o * 4.0 * math.pi / 3.0
                 * (r_out * r_out * r_out - r_in * r_in * r_in))
        centers.append(0.5 * (r_out + r_in))
        eps_c.append(eps[2])
    ifaces = [r_nm - a_nm, r_core]
    for i in range(1, n_bulk_cells):
        ifaces.append(r_core - i * w)

    e_ts = cgp['E_m_bulk_eV'] + eps[2]
    nu0 = cgp['nu0_per_s']
    lam = []
    for j in range(m - 1):
        area = 4.0 * math.pi * ifaces[j] * ifaces[j]
        dist = centers[j] - centers[j + 1]
        lam.append(nu0 * a_nm * a_nm * area * math.exp(-beta * e_ts)
                   / (dist * v_o * q[j] * q[j + 1]))

    n_o_particle = n_o * 4.0 * math.pi / 3.0 * r_nm * r_nm * r_nm
    n_o_per_g = 2.0 / p['molar_mass_g_mol'] * 1e6
    v_total = vo_total / n_o_per_g * n_o_particle

    if dist_fractions is None:
        geom = base.geometry(p, d_um=d_um)
        d0 = base.distribute(p, kin_used['T_reox_C'], vo_total, geom,
                             eps=eps)
        dist_fractions = d0['fractions']
    eta = [0.0] * m
    want_s = dist_fractions['surface'] * v_total
    eta[0] = min(want_s, 0.999999 * q[0])
    want_ss = dist_fractions['subsurface'] * v_total + (want_s - eta[0])
    eta[1] = min(want_ss, 0.999999 * q[1])
    rest = v_total - eta[0] - eta[1]
    q_bulk = sum(q[2:])
    for j in range(2, m):
        eta[j] = rest * q[j] / q_bulk

    return {'m': m, 'q': q, 'eps_c': eps_c, 'lam': lam, 'beta': beta,
            'k_fill': kf, 'eta0': eta, 'V0': v_total, 'R_nm': r_nm,
            'a_nm': a_nm, 'centers': centers, 'T_reox_C':
            kin_used['T_reox_C'], 'E_m_eV': cgp['E_m_bulk_eV'],
            'eps_ctrl': cgp['eps_ctrl'],
            'column_counts': cgp['stochastic_column_counts']}


def _expm2(a, b, c, d, dt, x0, y0):
    """exp(M dt) (x0,y0) for M = [[a,b],[c,d]], real non-positive spectrum.

    Written in eigenvalue form: both exp(lambda dt) factors are <= 1 for
    this dissipative block, so there is no overflow at any dt."""
    mu = 0.5 * (a + d)
    disc = 0.25 * (a - d) * (a - d) + b * c
    root = math.sqrt(disc if disc > 0.0 else 0.0)
    lp = (mu + root) * dt
    lm = (mu - root) * dt
    ep = math.exp(lp) if lp < 700.0 else math.exp(700.0)
    em = math.exp(lm) if lm > -700.0 else 0.0
    half_sum = 0.5 * (ep + em)
    if root * dt < 1e-8:
        sh_over = dt * ep
    else:
        sh_over = (ep - em) / (2.0 * root)
    x = half_sum * x0 + sh_over * ((a - mu) * x0 + b * y0)
    y = half_sum * y0 + sh_over * (c * x0 + (d - mu) * y0)
    return (x if x > 0.0 else 0.0), (y if y > 0.0 else 0.0)


def _sweep(sys, q, eta, dt, forward):
    """One conservative sweep of exact pair relaxations.

    For each adjacent pair with blocking factors frozen at the pair's
    entry values, the exchange is a linear two-state system whose exact
    solution relaxes eta_k exponentially toward the pair equilibrium
    x* = (eta_k + eta_l) r_l / (r_k + r_l). Mass moves only within the
    pair, so conservation is exact at any dt; a pair already in detailed
    balance sits at x* and is left untouched, so the closed-system
    equilibrium is a fixed point of the stepper, not a tolerance."""
    m = sys['m']
    lam = sys['lam']
    eps_c = sys['eps_c']
    beta = sys['beta']
    kf = sys['k_fill']
    removed = 0.0
    new = list(eta)
    order = range(m - 1) if forward else range(m - 2, -1, -1)
    for j in order:
        k = j
        l = j + 1
        a1 = math.exp(beta * eps_c[k])
        a2 = math.exp(beta * eps_c[l])
        if j == 0 and kf > 0.0:
            # surface pair with the refill sink, solved exactly at frozen
            # blocking: the slaved surface cell then follows its smooth
            # quasi-steady trajectory instead of a drain/rebalance sawtooth
            r01 = lam[0] * a1 * (q[1] - new[1])
            r10 = lam[0] * a2 * (q[0] - new[0])
            tot0 = new[0] + new[1]
            x, y = _expm2(-(r01 + kf), r10, r01, -r10, dt, new[0], new[1])
            if x > q[0]:
                x = q[0]
            if y > q[1]:
                y = q[1]
            if x + y > tot0:
                s2 = tot0 / (x + y)
                x *= s2
                y *= s2
            removed += tot0 - x - y
            new[0] = x
            new[1] = y
            continue
        rk = lam[j] * a1 * (q[l] - new[l])
        rl = lam[j] * a2 * (q[k] - new[k])
        s = rk + rl
        if s <= 0.0:
            continue
        tot = new[k] + new[l]
        # exact pair balance with the true blocking factors:
        # a1 x (q_l - tot + x) = a2 (tot - x)(q_k - x), quadratic in x.
        # Relaxing toward the exact balance (not a frozen-blocking one)
        # makes every pair update a no-op at global equilibrium, so the
        # stepper cannot orbit it; the entry rates still set the pace.
        aa = a1 - a2
        bb = a1 * (q[l] - tot) + a2 * (q[k] + tot)
        cc = -a2 * tot * q[k]
        lo = tot - q[l]
        if lo < 0.0:
            lo = 0.0
        hi = tot if tot < q[k] else q[k]
        if abs(aa) * (abs(lo) + abs(hi)) < 1e-14 * abs(bb):
            xstar = -cc / bb
        else:
            disc = bb * bb - 4.0 * aa * cc
            sq = math.sqrt(disc if disc > 0.0 else 0.0)
            qq = -0.5 * (bb + (sq if bb >= 0.0 else -sq))
            xstar = cc / qq
            if not (lo <= xstar <= hi):
                xstar = qq / aa
        if xstar < lo:
            xstar = lo
        if xstar > hi:
            xstar = hi
        x = s * dt
        w = 1.0 - (math.exp(-x) if x < 700.0 else 0.0)
        xk = new[k] + (xstar - new[k]) * w
        if xk < lo:
            xk = lo
        if xk > hi:
            xk = hi
        new[k] = xk
        new[l] = tot - xk
    return new, removed


def run(sys, t_end, mode='det', seed=902, n_out=90, profile_times=None,
        eps_ctrl=None):
    """Integrate to t_end; log-spaced outputs of recovery and profiles.

    Deterministic: adaptive alternating-direction pair-relaxation sweeps
    (exact per pair at frozen blocking) with exact exponential surface
    refill. Stochastic: tau-leap of the same frozen channel rates on a
    column scaled to ~column_counts vacancies, Poisson below 30 expected
    events per channel, Gaussian above, seeded with the shared RNG."""
    m = sys['m']
    eta = list(sys['eta0'])
    scale = 1.0
    rng = None
    if mode == 'stoch':
        scale = max(1.0, sys['V0'] / sys['column_counts'])
        eta = [x / scale for x in eta]
        rng = base.make_rng(seed)
    q = [x / scale for x in sys['q']]
    v0 = sys['V0'] / scale
    kf = sys['k_fill']
    if eps_ctrl is None:
        eps_ctrl = sys['eps_ctrl']

    t_first = 1e-11
    outs = [0.0]
    for i in range(n_out):
        outs.append(t_first * (t_end / t_first) ** (i / (n_out - 1.0)))
    if profile_times is None:
        profile_times = [t_end * f for f in (1e-6, 1e-3, 3e-2, 1.0)]
    prof_left = sorted(profile_times)
    profiles = []

    times = [0.0]
    rec = [0.0]
    t = 0.0
    dt = 1e-12
    oi = 1
    removed = 0.0
    guard = 0
    step_i = 0
    while t < t_end and guard < 400000:
        guard += 1
        if mode == 'stoch':
            # tau-leap validity: expected events per channel stay a small
            # fraction of the source population (stiff channels pin dt)
            cap = _leap_cap(sys, q, eta, kf)
            if dt > cap:
                dt = cap
        while True:
            if mode == 'stoch':
                fill = _poisson(kf * eta[0] * dt, rng)
                if fill > eta[0]:
                    fill = eta[0]
                work = list(eta)
                work[0] -= fill
                new = _leap(sys, q, work, dt, rng)
            else:
                new, fill = _sweep(sys, q, eta, dt, step_i % 2 == 0)
            worst = 0.0
            floor = 1e-6 * v0
            for k in range(m):
                # resolve both the occupancy and the blocking factor:
                # whichever of eta and (q - eta) is smaller is the
                # sensitive quantity; cells holding less than ~1e-6 of
                # the inventory do not constrain the step
                ref = eta[k] if eta[k] < q[k] - eta[k] else q[k] - eta[k]
                d = abs(new[k] - eta[k]) / (ref + floor)
                if d > worst:
                    worst = d
            if worst <= eps_ctrl or dt <= 1e-13:
                break
            dt *= 0.5
        eta = new
        removed += fill
        t += dt
        step_i += 1
        if worst < 0.5 * eps_ctrl:
            dt *= 1.3
        if t + dt > t_end:
            dt = max(t_end - t, 1e-13)
        while oi < len(outs) and t >= outs[oi]:
            times.append(t)
            rec.append(removed / v0)
            oi += 1
        while prof_left and t >= prof_left[0]:
            profiles.append({'t_s': t,
                             'theta': [eta[k] / q[k] for k in range(m)]})
            prof_left.pop(0)
    while prof_left:
        profiles.append({'t_s': t,
                         'theta': [eta[k] / q[k] for k in range(m)]})
        prof_left.pop(0)
    times.append(t)
    rec.append(removed / v0)
    return {'t_s': times, 'recovery': rec, 'final_eta': eta,
            'final_theta': [eta[k] / q[k] for k in range(m)],
            'profiles': profiles, 'steps': guard, 'mode': mode,
            'recovered_fraction': rec[-1], 'V0_column': v0,
            'mass_check': (sum(eta) + removed) / v0 - 1.0,
            'completed': t >= 0.999 * t_end,
            'centers_nm': sys['centers'], 'depth_nm':
            [sys['R_nm'] - c for c in sys['centers']]}


def _poisson(lam_, rng):
    if lam_ <= 0.0:
        return 0.0
    if lam_ < 30.0:
        limit = math.exp(-lam_)
        k = 0
        prod = rng()
        while prod > limit:
            k += 1
            prod *= rng()
        return float(k)
    u1 = rng()
    u2 = rng()
    z = math.sqrt(-2.0 * math.log(1.0 - u1)) * math.cos(2.0 * math.pi * u2)
    n = round(lam_ + math.sqrt(lam_) * z)
    return float(n) if n > 0 else 0.0


def _leap_cap(sys, q, eta, kf):
    m = sys['m']
    lam = sys['lam']
    eps_c = sys['eps_c']
    beta = sys['beta']
    cap = 1e30
    for j in range(m - 1):
        k = j
        l = j + 1
        for (rate, pop) in (
                (lam[j] * math.exp(beta * eps_c[k]) * eta[k]
                 * (q[l] - eta[l]), eta[k]),
                (lam[j] * math.exp(beta * eps_c[l]) * eta[l]
                 * (q[k] - eta[k]), eta[l])):
            if rate > 0.0:
                c = 0.2 * max(pop, 30.0) / rate
                if c < cap:
                    cap = c
    if kf > 0.0 and eta[0] > 0.0:
        c = 0.2 * max(eta[0], 30.0) / (kf * eta[0])
        if c < cap:
            cap = c
    return cap


def _leap(sys, q, eta, dt, rng):
    """One tau-leap of the directed channel rates, blocking frozen."""
    m = sys['m']
    lam = sys['lam']
    eps_c = sys['eps_c']
    beta = sys['beta']
    new = list(eta)
    for j in range(m - 1):
        k = j
        l = j + 1
        fkl = lam[j] * math.exp(beta * eps_c[k]) * eta[k] * (q[l] - eta[l])
        flk = lam[j] * math.exp(beta * eps_c[l]) * eta[l] * (q[k] - eta[k])
        for pair in ((fkl, k, l), (flk, l, k)):
            x = pair[0] * dt
            n = x if x > 1e5 else _poisson(x, rng)
            if n > new[pair[1]]:
                n = new[pair[1]]
            room = q[pair[2]] - new[pair[2]]
            if n > room:
                n = room if room > 0.0 else 0.0
            new[pair[1]] -= n
            new[pair[2]] += n
    return new


def effective_barrier(p, target_recovery, t_exp, vo_total=None, d_um=None,
                      kin=None, cg=None, dist_fractions=None, eps=None,
                      lo=0.3, hi=3.0):
    """The E_m that reproduces a measured recovery at t_exp (bisection).

    Recovery falls monotonically with E_m once transport is limiting; the
    readout is the panel's bridge from the measured transient to an
    effective migration barrier."""
    for _ in range(36):
        mid = 0.5 * (lo + hi)
        sys = build(p, vo_total=vo_total, d_um=d_um, kin=kin,
                    cg=dict(cg or {}, E_m_bulk_eV=mid),
                    dist_fractions=dist_fractions, eps=eps)
        r = run(sys, t_exp, n_out=12)['recovered_fraction']
        if r > target_recovery:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)
