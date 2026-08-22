"""Independent high-precision reference for the defect stat-mech engine.

Shares only rutile_dft.json with the production code - no import of
solidgas.statmech, no shared solver code. Three kinds of certified values:

1. Analytic matching grid: the uncapped three-class site-exclusion
   distribution, solved by mpmath bisection at 50 digits. The reported
   reconstruction interval is retained as a validity warning.
2. Exact ideal-lattice canonical ensemble (generating polynomial): finite-
   lattice class occupancies and the exact chemical potential
   mu = -kT ln(Z(N+1)/Z(N)) that Widom insertion estimates.
3. Exact enumeration of small interacting lattices (a 10-site ring and a
   24-site three-class slab): canonical mean energy, gap statistics, class
   occupancies, chemical potential.

Output: reference_statmech.json, consumed by tests/test_statmech.py and by
the JS parity harness.
"""

import itertools
import json
import math

from mpmath import mp, mpf, exp as mexp, log as mlog

mp.dps = 50

KB_EV = mpf('8.617333262e-5')
N_AVO = mpf('6.02214076e23')

with open('rutile_dft.json') as fh:
    P = json.load(fh)


# ------------------------------------------------------------ geometry

def geometry(d_um=None, bet_m2_g=None, ss_layers=1):
    a = mpf(str(P['lattice_constants_A']['a']))
    c = mpf(str(P['lattice_constants_A']['c']))
    m = mpf(str(P['molar_mass_g_mol']))
    cell = c * a * mp.sqrt(2)
    sig_b = mpf('1e16') / cell
    sig_l = mpf('4e16') / cell
    rho = 2 * m / (N_AVO * a * a * c * mpf('1e-24'))
    if bet_m2_g is not None:
        area = mpf(str(bet_m2_g)) * mpf('1e4')
    else:
        area = 6 / (rho * mpf(str(d_um)) * mpf('1e-4'))
    n_o = 2 / m * mpf('1e6')
    n_s = area * sig_b / N_AVO * mpf('1e6')
    n_ss = area * sig_l * ss_layers / N_AVO * mpf('1e6')
    return {'area_m2_g': area * mpf('1e-4'), 'N_s': n_s, 'N_ss': n_ss,
            'N_b': n_o - n_s - n_ss, 'N_O_total': n_o}


def fd(mu, eps, kt):
    return 1 / (1 + mexp((eps - mu) / kt))


def eps_of(preset):
    q = P['vacancy_energetics']['presets'][preset]
    return [mpf(str(q['surface_eV'])), mpf(str(q['subsurface_eV'])),
            mpf(str(q['bulk_eV']))]


def boundaries(t_c, geom, preset):
    """Closed-form phase boundaries: (1x2) onset, reconstruction complete,
    CS-precipitation ceiling. Independent re-derivation at 50 digits."""
    kt = KB_EV * (mpf(str(t_c)) + mpf('273.15'))
    eps = eps_of(preset)
    sp = P['surface_phases']
    th_t = mpf(str(sp['theta_transition']))
    th_r = mpf(str(sp['theta_reconstructed_eff']))
    th_sol = mpf(str(P['saturation']['x_max_shear'])) / 2
    mu_t = eps[0] + kt * mlog(th_t / (1 - th_t))
    mu_cs = eps[2] + kt * mlog(th_sol / (1 - th_sol))
    deep_t = (geom['N_ss'] * fd(mu_t, eps[1], kt)
              + geom['N_b'] * fd(mu_t, eps[2], kt))
    deep_cs = (geom['N_ss'] * fd(mu_cs, eps[1], kt)
               + geom['N_b'] * fd(mu_cs, eps[2], kt))
    return {'mu_t': mu_t, 'mu_cs': mu_cs, 'th_t': th_t, 'th_r': th_r,
            'th_sol': th_sol, 'kt': kt, 'eps': eps,
            'vo_onset': geom['N_s'] * th_t + deep_t,
            'vo_recon': geom['N_s'] * th_r + deep_t,
            'vo_cs': geom['N_s'] * th_r + deep_cs}


def match(t_c, vo, geom, preset):
    """Phase-aware matching: (1x1) lattice gas, first-order (1x2)
    coexistence riser, reconstructed line phase, mu pinned at the CS
    coexistence estimate with lever-rule precipitation."""
    b = boundaries(t_c, geom, preset)
    kt, eps = b['kt'], b['eps']
    vo_mp = mpf(str(vo))

    def deep(mu):
        return (geom['N_ss'] * fd(mu, eps[1], kt)
                + geom['N_b'] * fd(mu, eps[2], kt))
    ext = mpf(0)
    if vo_mp <= b['vo_onset']:
        regime, phase = 'dilute', '1x1'
        lo, hi = mpf(-40), b['mu_t']
        for _ in range(240):
            mid = (lo + hi) / 2
            if geom['N_s'] * fd(mid, eps[0], kt) + deep(mid) < vo_mp:
                lo = mid
            else:
                hi = mid
        mu = (lo + hi) / 2
        ts = fd(mu, eps[0], kt)
        us = geom['N_s'] * ts
        phi = mpf(0)
    elif vo_mp < b['vo_recon']:
        regime, phase = 'surface_coexistence', 'coexistence'
        mu = b['mu_t']
        us = vo_mp - deep(mu)
        ts = us / geom['N_s']
        phi = (ts - b['th_t']) / (b['th_r'] - b['th_t'])
    else:
        if vo_mp <= b['vo_cs']:
            regime, phase = 'reconstructed', '1x2'
            lo, hi = b['mu_t'], b['mu_cs']
            target = vo_mp - geom['N_s'] * b['th_r']
            for _ in range(240):
                mid = (lo + hi) / 2
                if deep(mid) < target:
                    lo = mid
                else:
                    hi = mid
            mu = (lo + hi) / 2
        else:
            regime, phase = 'cs_coexistence', '1x2'
            mu = b['mu_cs']
            ext = vo_mp - geom['N_s'] * b['th_r'] - deep(mu)
        ts = b['th_r']
        us = geom['N_s'] * b['th_r']
        phi = mpf(1)
    tss = fd(mu, eps[1], kt)
    tb = fd(mu, eps[2], kt)
    uss, ub = geom['N_ss'] * tss, geom['N_b'] * tb
    tot = us + uss + ub + ext
    x = vo_mp / (mpf('1e6') / mpf(str(P['molar_mass_g_mol'])))
    sat = P['saturation']
    warn = []
    if phase != '1x1':
        warn.append('surface_phase')
    if tss > mpf(str(sat['subsurface_dilute_warn'])):
        warn.append('subsurface_dense')
    if regime == 'cs_coexistence':
        warn.append('cs_precipitation')
    elif x >= mpf('0.8') * mpf(str(sat['x_max_shear'])):
        warn.append('x_near')
    return {'T_C': t_c, 'VO_total_umol_g': vo, 'preset': preset,
            'mu_V_eV': float(mu),
            'theta': {'surface': float(ts), 'subsurface': float(tss),
                      'bulk': float(tb)},
            'umol_g': {'surface': float(us), 'subsurface': float(uss),
                       'bulk': float(ub)},
            'extended_defects_umol_g': float(ext),
            'fractions': {'surface': float(us / tot),
                          'subsurface': float(uss / tot),
                          'bulk': float(ub / tot),
                          'extended': float(ext / tot)},
            'x_TiO2mx': float(x), 'Ti3_frac': float(2 * x),
            'regime': regime,
            'surface_phase': {'phase': phase,
                              'phi_reconstructed': float(phi)},
            'boundaries': {'mu_t_eV': float(b['mu_t']),
                           'mu_cs_eV': float(b['mu_cs']),
                           'VO_onset_umol_g': float(b['vo_onset']),
                           'VO_recon_complete_umol_g': float(b['vo_recon']),
                           'VO_cs_umol_g': float(b['vo_cs'])},
            'surface_reconstruction_regime': phase != '1x1',
            'warn_kinds': warn}


# ------------------------------------------------ exact ideal ensemble

def ideal_ensemble(m_by_class, eps, t_c, n_vac):
    """Coefficients of prod_c (1 + z w_c)^{M_c}: exact Z(N), <n_c>, mu."""
    kt = KB_EV * (mpf(str(t_c)) + mpf('273.15'))
    w = [mexp(-mpf(str(e)) / kt) for e in eps]

    def convolve(pa, pb, limit):
        out = [mpf(0)] * min(len(pa) + len(pb) - 1, limit + 1)
        for i, ai in enumerate(pa):
            if i > limit:
                break
            for j, bj in enumerate(pb):
                if i + j > limit:
                    break
                out[i + j] += ai * bj
        return out

    def binom_poly(m_c, w_c, limit):
        poly = [mpf(1)]
        base = [mpf(1), w_c]
        k = m_c
        while k:
            if k & 1:
                poly = convolve(poly, base, limit)
            k >>= 1
            if k:
                base = convolve(base, base, limit)
        return poly

    limit = n_vac + 1

    def z_poly(ms):
        poly = [mpf(1)]
        for m_c, w_c in zip(ms, w):
            poly = convolve(poly, binom_poly(m_c, w_c, limit), limit)
        return poly

    zp = z_poly(m_by_class)
    z_n = zp[n_vac]
    z_n1 = zp[n_vac + 1] if n_vac + 1 < len(zp) else mpf(0)
    theta = []
    for c, m_c in enumerate(m_by_class):
        ms = list(m_by_class)
        ms[c] = m_c - 1
        zc = z_poly(ms)
        mean_n = m_c * w[c] * zc[n_vac - 1] / z_n
        theta.append(float(mean_n / m_c))
    mu = -kt * mlog(z_n1 / z_n)
    return {'theta_exact': theta, 'mu_exact_eV': float(mu)}


# --------------------------------------------- exact small enumerations

def ring_energy(sites, ls, pair):
    e = mpf(0)
    ss = sorted(sites)
    for i in range(len(ss)):
        for j in range(i + 1, len(ss)):
            d = ss[j] - ss[i]
            d = min(d, ls - d)
            if d <= len(pair):
                e += mpf(str(pair[d - 1]))
    return e


def ring_case(ls, n_vac, pair, t_c):
    kt = KB_EV * (mpf(str(t_c)) + mpf('273.15'))
    z = mpf(0)
    e_sum = mpf(0)
    gap_counts = [mpf(0)] * (ls + 1)
    for combo in itertools.combinations(range(ls), n_vac):
        e = ring_energy(combo, ls, pair)
        w = mexp(-e / kt)
        z += w
        e_sum += e * w
        ss = sorted(combo)
        for t in range(len(ss)):
            gap = (ss[t + 1] - ss[t]) if t + 1 < len(ss) else ss[0] + ls - ss[t]
            gap_counts[gap] += w
    z1 = mpf(0)
    for combo in itertools.combinations(range(ls), n_vac + 1):
        z1 += mexp(-ring_energy(combo, ls, pair) / kt)
    tot_gaps = sum(gap_counts)
    return {'row_sites': ls, 'n_vac': n_vac, 'pair_eV': pair, 'T_C': t_c,
            'E_mean_eV': float(e_sum / z),
            'gap_probs': [float(g / tot_gaps) for g in gap_counts],
            'mu_exact_eV': float(-kt * mlog(z1 / z))}


def micro3_case(ls, n_vac, pair, eps, t_c):
    """One surface ring of ls sites plus ls subsurface and ls bulk sites."""
    kt = KB_EV * (mpf(str(t_c)) + mpf('273.15'))
    n = 3 * ls
    epsm = [mpf(str(e)) for e in eps]

    def config_energy(combo):
        e = mpf(0)
        surf = []
        for s in combo:
            c = s // ls
            e += epsm[c]
            if c == 0:
                surf.append(s)
        return e + ring_energy(surf, ls, pair)

    z = mpf(0)
    e_sum = mpf(0)
    n_cls = [mpf(0)] * 3
    for combo in itertools.combinations(range(n), n_vac):
        e = config_energy(combo)
        w = mexp(-e / kt)
        z += w
        e_sum += e * w
        for s in combo:
            n_cls[s // ls] += w
    z1 = mpf(0)
    for combo in itertools.combinations(range(n), n_vac + 1):
        z1 += mexp(-config_energy(combo) / kt)
    return {'row_sites': ls, 'n_vac': n_vac, 'pair_eV': pair,
            'eps_eV': eps, 'T_C': t_c,
            'theta_exact': [float(n_cls[c] / z / ls) for c in range(3)],
            'E_mean_eV': float(e_sum / z),
            'mu_exact_eV': float(-kt * mlog(z1 / z))}


# ------------------------------------------------- CO2 accessibility

def co2_accessibility(umol_by_class, t_reox_c, exposure_s, p_co2_atm):
    """Same closed form as the engine, evaluated at 50 digits."""
    kin = P['co2_kinetics']
    kt = KB_EV * (mpf(str(t_reox_c)) + mpf('273.15'))
    pf = mpf(str(p_co2_atm)) ** mpf(str(kin['p_order_m']))
    probs = {}
    acc = mpf(0)
    tot = mpf(0)
    for c in ('surface', 'subsurface', 'bulk'):
        k = (mpf(str(kin['A_eff_per_s_atm'][c])) * pf
             * mexp(-mpf(str(kin['Ea_eV'][c])) / kt))
        prob = 1 - mexp(-k * mpf(str(exposure_s)))
        probs[c] = float(prob)
        acc += mpf(str(umol_by_class[c])) * prob
        tot += mpf(str(umol_by_class[c]))
    return {'T_reox_C': t_reox_c, 'exposure_s': exposure_s,
            'p_CO2_atm': p_co2_atm, 'P': probs,
            'accessible_umol_g': float(acc),
            'f_recoverable': float(acc / tot)}


# ------------------------------------------------------------ sfc32 spec

def sfc32_first(seed, count):
    state = seed & 0xFFFFFFFF
    words = []
    for _ in range(4):
        state = (state + 0x9E3779B9) & 0xFFFFFFFF
        zz = state
        zz = ((zz ^ (zz >> 16)) * 0x21F0AAAD) & 0xFFFFFFFF
        zz = ((zz ^ (zz >> 15)) * 0x735A2D97) & 0xFFFFFFFF
        words.append(zz ^ (zz >> 15))
    a, b, c, d = words
    out = []
    for i in range(12 + count):
        t = (a + b) & 0xFFFFFFFF
        a = b ^ (b >> 9)
        b = (c + ((c << 3) & 0xFFFFFFFF)) & 0xFFFFFFFF
        c = ((c << 21) | (c >> 11)) & 0xFFFFFFFF
        d = (d + 1) & 0xFFFFFFFF
        t = (t + d) & 0xFFFFFFFF
        c = (c + t) & 0xFFFFFFFF
        if i >= 12:
            out.append(t / 4294967296.0)
    return out


def main():
    geom = geometry(d_um=P['defaults']['particle_diameter_um'])
    grid = []
    for t_c in (400, 600, 800, 1000):
        for vo in (1.0, 10.0, 38.0, 95.0, 150.0):
            grid.append(match(t_c, vo, geom, 'sx'))
    flagship = match(P['defaults']['T_C'], P['defaults']['VO_total_umol_g'],
                     geom, 'sx')
    flagship_gga = match(P['defaults']['T_C'],
                         P['defaults']['VO_total_umol_g'], geom, 'gga')
    bet_geom = geometry(bet_m2_g=1.60)
    bet_case = match(600, 95.0, bet_geom, 'sx')

    eps_sx = [float(e) for e in eps_of('sx')]
    pair = P['surface_ordering']['pair_eV']['along_row']
    kin_d = P['co2_kinetics']['defaults']
    co2_rows = [co2_accessibility(flagship['umol_g'], kin_d['T_reox_C'],
                                  t, kin_d['p_CO2_atm'])
                for t in (30, 300, 3000)]
    co2_rows.append(co2_accessibility(flagship['umol_g'], 700, 300,
                                      kin_d['p_CO2_atm']))

    doc = {
        'generated_by': 'oracle_statmech.py',
        'precision_dps': mp.dps,
        'geometry_0p9um': {k: float(v) for k, v in geom.items()},
        'geometry_bet_1p60': {k: float(v) for k, v in bet_geom.items()},
        'rng_sfc32_seed12345_first8': sfc32_first(12345, 8),
        'analytic_grid': grid,
        'flagship': flagship,
        'flagship_gga': flagship_gga,
        'bet_case': bet_case,
        'ideal_ensemble': dict(
            {'m_by_class': [48, 48, 96], 'eps_eV': eps_sx, 'T_C': 600,
             'n_vac': 14},
            **ideal_ensemble([48, 48, 96], eps_sx, 600, 14)),
        'ring_case': ring_case(10, 3, pair, 600),
        'micro3_case': micro3_case(8, 4, pair, eps_sx, 600),
        'co2_flagship': co2_rows,
    }
    with open('reference_statmech.json', 'w') as fh:
        json.dump(doc, fh, indent=1)
    print('reference_statmech.json written')
    print('flagship: mu=%.9f f=(%.6f, %.6f, %.6f) warn=%s' % (
        flagship['mu_V_eV'], flagship['fractions']['surface'],
        flagship['fractions']['subsurface'], flagship['fractions']['bulk'],
        flagship['warn_kinds']))


if __name__ == '__main__':
    main()
