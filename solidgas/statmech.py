"""Canonical defect statistical mechanics for oxygen vacancies in rutile.

The measured total vacancy inventory is distributed over three site classes -
rutile (110) bridging O, first subsurface O layer, bulk - by equalising the
vacancy chemical potential. Surface vacancies interact through a short-range
pair table, so the surface isotherm theta_s(mu) is sampled by canonical
Metropolis Monte Carlo on a representative slab (Widom insertion supplies mu);
subsurface and bulk are non-interacting in this version and follow the exact
site-exclusion isotherm. The particle-scale balance then uses the real
geometric site counts, not the slab's class ratio: a 0.9 um particle has
~0.05% of its O sites in the surface layer, the slab has 10%.

The reported 17-20% onset interval for the (1x2) reconstruction is a model
boundary, not a saturation capacity. Surface coverage is never clipped. A
solution inside or above that interval is marked as a conditional extrapolation
of the unreconstructed (1x1) Hamiltonian because reconstructed-phase
energetics are not implemented.

statmech.js is a line-for-line port of this file. Loops, operation order and
the RNG are kept identical so that a same-seed run is bit-reproducible across
the two implementations (the only library call inside the sampling loop is
exp(), whose last-ulp differences are compared against, never accumulated,
so a decision flip has probability ~1e-16 per move).
"""

import json
import math
import os

KB_EV = 8.617333262e-5          # eV / K
N_AVO = 6.02214076e23

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_PATH = os.path.join(_ROOT, 'data', 'rutile_dft.json')


def load_params(path=DATA_PATH):
    with open(path) as fh:
        return json.load(fh)


# ---------------------------------------------------------------- RNG
# sfc32 seeded by splitmix32: fast, 128-bit state, and expressible in
# identical uint32 arithmetic in Python and JavaScript.

def _splitmix32(seed):
    state = seed & 0xFFFFFFFF

    def nxt():
        nonlocal state
        state = (state + 0x9E3779B9) & 0xFFFFFFFF
        z = state
        z = ((z ^ (z >> 16)) * 0x21F0AAAD) & 0xFFFFFFFF
        z = ((z ^ (z >> 15)) * 0x735A2D97) & 0xFFFFFFFF
        return z ^ (z >> 15)
    return nxt


def make_rng(seed):
    sm = _splitmix32(seed)
    st = [sm(), sm(), sm(), sm()]

    def rng():
        a, b, c, d = st
        t = (a + b) & 0xFFFFFFFF
        a = b ^ (b >> 9)
        b = (c + ((c << 3) & 0xFFFFFFFF)) & 0xFFFFFFFF
        c = ((c << 21) | (c >> 11)) & 0xFFFFFFFF
        d = (d + 1) & 0xFFFFFFFF
        t = (t + d) & 0xFFFFFFFF
        c = (c + t) & 0xFFFFFFFF
        st[0], st[1], st[2], st[3] = a, b, c, d
        return t / 4294967296.0
    for _ in range(12):
        rng()
    return rng


# ---------------------------------------------------------------- geometry

def site_densities(p):
    """(bridging-O, full-O-layer) areal densities in cm^-2 for the (110) face.

    The 1x1 (110) cell is c x a*sqrt(2) and holds 1 bridging O and 4 O in a
    full layer of thickness d110 = a/sqrt(2)."""
    a = p['lattice_constants_A']['a']
    c = p['lattice_constants_A']['c']
    cell_A2 = c * a * math.sqrt(2.0)
    return 1e16 / cell_A2, 4e16 / cell_A2


def density_g_cm3(p):
    a = p['lattice_constants_A']['a']
    c = p['lattice_constants_A']['c']
    vcell_cm3 = a * a * c * 1e-24
    return 2.0 * p['molar_mass_g_mol'] / (N_AVO * vcell_cm3)


def geometry(p, d_um=None, bet_m2_g=None, ss_layers=None):
    """Real site counts in umol-O/g for a spherical particle or a BET area."""
    if ss_layers is None:
        ss_layers = p['defaults']['subsurface_layers']
    if bet_m2_g is not None:
        area_cm2_g = bet_m2_g * 1e4
    else:
        if d_um is None:
            d_um = p['defaults']['particle_diameter_um']
        area_cm2_g = 6.0 / (density_g_cm3(p) * d_um * 1e-4)
    sig_b, sig_l = site_densities(p)
    n_o = 2.0 / p['molar_mass_g_mol'] * 1e6
    n_s = area_cm2_g * sig_b / N_AVO * 1e6
    n_ss = area_cm2_g * sig_l * ss_layers / N_AVO * 1e6
    return {'area_m2_g': area_cm2_g * 1e-4, 'N_s': n_s, 'N_ss': n_ss,
            'N_b': n_o - n_s - n_ss, 'N_O_total': n_o,
            'ss_layers': ss_layers}


# ---------------------------------------------------------------- lattice

def build_lattice(rows, row_sites, ss_layers, bulk_layers, ordering):
    """Slab site classes and the surface pair-interaction neighbour lists.

    along_row[k] is the pair energy at separation k+1 bridging sites; the
    cross_row list covers the two inter-row offsets inside the 10 A cutoff.
    Lists are bidirectional and built in a fixed order for float parity."""
    along = ordering['pair_eV']['along_row']
    cross = ordering['pair_eV']['cross_row']
    if along and row_sites < 2 * len(along) + 1:
        raise ValueError('row too short for the along-row interaction range')
    if any(v != 0.0 for v in cross) and rows < 3:
        raise ValueError('cross-row interactions need at least 3 rows')
    n_surf = rows * row_sites
    n_ss = ss_layers * rows * row_sites
    n_bulk = bulk_layers * rows * row_sites
    n = n_surf + n_ss + n_bulk
    cls = [0] * n_surf + [1] * n_ss + [2] * n_bulk
    nbrs = [None] * n_surf

    def sid(r, k):
        return (r % rows) * row_sites + (k % row_sites)
    for r in range(rows):
        for k in range(row_sites):
            lst = []
            for dz in range(1, len(along) + 1):
                v = along[dz - 1]
                if v != 0.0:
                    lst.append((sid(r, k + dz), v))
                    lst.append((sid(r, k - dz), v))
            for off in range(len(cross)):
                v = cross[off]
                if v == 0.0:
                    continue
                if off == 0:
                    lst.append((sid(r + 1, k), v))
                    lst.append((sid(r - 1, k), v))
                else:
                    lst.append((sid(r + 1, k + off), v))
                    lst.append((sid(r + 1, k - off), v))
                    lst.append((sid(r - 1, k + off), v))
                    lst.append((sid(r - 1, k - off), v))
            nbrs[sid(r, k)] = lst
    # flat neighbour arrays: same iteration order, no per-pair objects in
    # the sampling loop (and the same layout the JS port uses)
    nbr_off = [0] * (n_surf + 1)
    nbr_idx = []
    nbr_v = []
    for i in range(n_surf):
        for j, v in nbrs[i]:
            nbr_idx.append(j)
            nbr_v.append(v)
        nbr_off[i + 1] = len(nbr_idx)
    return {'n': n, 'n_surf': n_surf, 'n_ss': n_ss, 'n_bulk': n_bulk,
            'rows': rows, 'row_sites': row_sites, 'cls': cls,
            'nbr_off': nbr_off, 'nbr_idx': nbr_idx, 'nbr_v': nbr_v}


def _inter(i, occ, lat):
    """Pair energy of a vacancy at i with the occupied vacancies around it."""
    if i >= lat['n_surf']:
        return 0.0
    s = 0.0
    off = lat['nbr_off']
    idx = lat['nbr_idx']
    vv = lat['nbr_v']
    for t in range(off[i], off[i + 1]):
        if occ[idx[t]]:
            s += vv[t]
    return s


def total_energy(lat, occ, eps_cls):
    e = 0.0
    for i in range(lat['n']):
        if occ[i]:
            e += eps_cls[lat['cls'][i]]
    p = 0.0
    off = lat['nbr_off']
    idx = lat['nbr_idx']
    vv = lat['nbr_v']
    for i in range(lat['n_surf']):
        if occ[i]:
            for t in range(off[i], off[i + 1]):
                if occ[idx[t]]:
                    p += vv[t]
    return e + 0.5 * p


# ---------------------------------------------------------------- MC

def _gap_hist(lat, occ, hist):
    rows, ls = lat['rows'], lat['row_sites']
    for r in range(rows):
        base = r * ls
        ks = [k for k in range(ls) if occ[base + k]]
        m = len(ks)
        if m < 2:
            continue
        for t in range(m):
            if t + 1 < m:
                gap = ks[t + 1] - ks[t]
            else:
                gap = ks[0] + ls - ks[t]
            hist[gap] += 1


def adjust_filling(lat, occ, vac, target, rng):
    n = lat['n']
    while len(vac) > target:
        iv = int(rng() * len(vac))
        occ[vac[iv]] = 0
        vac.pop(iv)
    while len(vac) < target:
        s = int(rng() * n)
        if not occ[s]:
            occ[s] = 1
            vac.append(s)


def mc_run(lat, eps_cls, kt, occ, vac, sweeps_eq, sweeps_sample, rng,
           widom_every=10, hist_every=25, blocks=20):
    """Canonical swap Metropolis on the slab; occ/vac are mutated in place.

    A move draws a vacancy and a uniform site (occupied targets are counted
    null moves, keeping the proposal symmetric); swaps are non-local because
    the target is the equilibrium ensemble, not diffusion dynamics. Widom
    insertion over the empty lattice gives mu = -kT ln(<W>/(N_V+1))."""
    n = lat['n']
    cls = lat['cls']
    n_surf = lat['n_surf']
    off = lat['nbr_off']
    idx = lat['nbr_idx']
    vv = lat['nbr_v']
    exp = math.exp
    nv = len(vac)
    beta = 1.0 / kt
    cnt = [0, 0, 0]
    for v in vac:
        cnt[cls[v]] += 1
    e_state = [total_energy(lat, occ, eps_cls)]

    def attempt():
        iv = int(rng() * nv)
        v = vac[iv]
        s = int(rng() * n)
        if occ[s]:
            return
        occ[v] = 0
        i_s = 0.0
        if s < n_surf:
            for t in range(off[s], off[s + 1]):
                if occ[idx[t]]:
                    i_s += vv[t]
        i_v = 0.0
        if v < n_surf:
            for t in range(off[v], off[v + 1]):
                if occ[idx[t]]:
                    i_v += vv[t]
        de = (eps_cls[cls[s]] - eps_cls[cls[v]]) + (i_s - i_v)
        if de <= 0.0 or rng() < exp(-de * beta):
            occ[s] = 1
            vac[iv] = s
            cnt[cls[v]] -= 1
            cnt[cls[s]] += 1
            e_state[0] += de
        else:
            occ[v] = 1

    for _ in range(sweeps_eq * n):
        attempt()

    bsize = max(1, sweeps_sample // blocks)
    th_blocks = []
    wid_sum = 0.0
    wid_n = 0
    hist = [0] * (lat['row_sites'] + 1)
    e_sum = 0.0
    for b in range(blocks):
        acc0 = 0
        acc1 = 0
        acc2 = 0
        for sw in range(bsize):
            for _ in range(n):
                attempt()
            acc0 += cnt[0]
            acc1 += cnt[1]
            acc2 += cnt[2]
            e_sum += e_state[0]
            gsw = b * bsize + sw
            if widom_every and gsw % widom_every == 0:
                w = 0.0
                for i in range(n):
                    if not occ[i]:
                        w += math.exp(-beta * (eps_cls[cls[i]]
                                               + _inter(i, occ, lat)))
                wid_sum += w
                wid_n += 1
            if hist_every and gsw % hist_every == 0:
                _gap_hist(lat, occ, hist)
        th_blocks.append([acc0 / (bsize * lat['n_surf']),
                          acc1 / (bsize * lat['n_ss']) if lat['n_ss'] else 0.0,
                          acc2 / (bsize * lat['n_bulk']) if lat['n_bulk'] else 0.0])

    theta = [0.0, 0.0, 0.0]
    for tb in th_blocks:
        for k in range(3):
            theta[k] += tb[k]
    for k in range(3):
        theta[k] /= blocks
    err = [0.0, 0.0, 0.0]
    if blocks > 1:
        for k in range(3):
            s2 = 0.0
            for tb in th_blocks:
                s2 += (tb[k] - theta[k]) ** 2
            err[k] = math.sqrt(s2 / (blocks - 1) / blocks)
    mu_w = -kt * math.log((wid_sum / wid_n) / (nv + 1)) if wid_n else None
    e_check = total_energy(lat, occ, eps_cls)
    return {'theta_s': theta[0], 'theta_ss': theta[1], 'theta_b': theta[2],
            'err_s': err[0], 'err_ss': err[1], 'err_b': err[2],
            'mu_eV': mu_w, 'E_mean_eV': e_sum / (blocks * bsize),
            'E_drift_eV': abs(e_state[0] - e_check),
            'hist': hist, 'n_vac': nv,
            'vac_sorted': sorted(vac),
            'surf_vac': sorted(v for v in vac if v < lat['n_surf'])}


def isotherm_scan(p, t_c, seed=None, quality=1.0, preset=None, eps=None,
                  zero_pairs=False, ordering=None, mc=None, progress=None):
    """theta_s(mu) by warm-started canonical runs at increasing filling.

    The isotherm depends on temperature and energetics only, so a single scan
    serves every total inventory and particle size."""
    d = dict(p['defaults']['mc'])
    if mc:
        d.update(mc)
    if ordering is None:
        ordering = p['surface_ordering']
    if zero_pairs:
        ordering = {'pair_eV': {'along_row': [], 'cross_row': []}}
    lat = build_lattice(d['rows'], d['row_sites'], d['ss_layers'],
                        d['bulk_layers'], ordering)
    eps_cls = eps if eps is not None else energetics(p, preset)['eps']
    kt = KB_EV * (t_c + 273.15)
    rng = make_rng(d['seed'] if seed is None else seed)
    occ = bytearray(lat['n'])
    vac = []
    sw_eq = max(1, int(round(d['equil_sweeps'] * quality)))
    sw_warm = max(1, int(round(d['warm_equil_sweeps'] * quality)))
    sw_sm = max(d['blocks'], int(round(d['sample_sweeps'] * quality)))
    pts = []
    for idx, target in enumerate(d['fillings']):
        adjust_filling(lat, occ, vac, target, rng)
        res = mc_run(lat, eps_cls, kt, occ, vac,
                     sw_eq if idx == 0 else sw_warm, sw_sm, rng,
                     d['widom_every'], d['hist_every'], d['blocks'])
        res['filling'] = target
        del res['vac_sorted']
        # the final surface configuration travels with the point: the
        # particle view draws the sampled arrangement, not synthetic noise
        res['rows'] = d['rows']
        res['row_sites'] = d['row_sites']
        pts.append(res)
        if progress:
            progress(idx + 1, len(d['fillings']))
    return pts


# ---------------------------------------------------------------- matching

def fd(mu, eps_i, kt):
    x = (eps_i - mu) / kt
    if x > 700.0:
        return 0.0
    if x < -700.0:
        return 1.0
    return 1.0 / (1.0 + math.exp(x))


def energetics(p, preset=None):
    ve = p['vacancy_energetics']
    key = preset or ve['default_preset']
    q = ve['presets'][key]
    return {'preset': key,
            'eps': [q['surface_eV'], q['subsurface_eV'], q['bulk_eV']]}


def theta_s_interp(points, eps_s, kt):
    """Monotone interpolant of the sampled surface isotherm.

    Below the sampled range the surface is dilute and near-ideal, so the
    exact non-interacting isotherm continues the curve. Above the sampled
    range, the last two points continue linearly in logit(theta). This avoids
    inventing a saturation plateau at the numerical scan boundary. Results
    past the reconstruction onset are separately marked as extrapolations of
    the unreconstructed surface Hamiltonian."""
    pts = [q for q in points if q['mu_eV'] is not None]
    pts = sorted(pts, key=lambda q: q['mu_eV'])
    mus = []
    ths = []
    run = 0.0
    for q in pts:
        if q['theta_s'] > run:
            run = q['theta_s']
        mus.append(q['mu_eV'])
        ths.append(run)

    def f(mu):
        if not mus:
            return fd(mu, eps_s, kt)
        if mu <= mus[0]:
            t = fd(mu, eps_s, kt)
            return t if t < ths[0] else ths[0]
        if mu >= mus[-1]:
            if len(mus) < 2 or mus[-1] <= mus[-2]:
                return ths[-1]
            tiny = 1e-12
            t0 = min(1.0 - tiny, max(tiny, ths[-2]))
            t1 = min(1.0 - tiny, max(tiny, ths[-1]))
            l0 = math.log(t0 / (1.0 - t0))
            l1 = math.log(t1 / (1.0 - t1))
            slope = max(0.0, (l1 - l0) / (mus[-1] - mus[-2]))
            z = l1 + slope * (mu - mus[-1])
            if z > 700.0:
                return 1.0
            return 1.0 / (1.0 + math.exp(-z))
        lo, hi = 0, len(mus) - 1
        while hi - lo > 1:
            mid = (lo + hi) // 2
            if mus[mid] <= mu:
                lo = mid
            else:
                hi = mid
        w = (mu - mus[lo]) / (mus[hi] - mus[lo])
        return ths[lo] + w * (ths[hi] - ths[lo])
    return f


def distribute(p, t_c, vo_total, geom, theta_s_fn=None, preset=None, eps=None):
    """Split the total inventory over the real site counts at common mu."""
    kt = KB_EV * (t_c + 273.15)
    if eps is None:
        eps = energetics(p, preset)['eps']
    if theta_s_fn is None:
        theta_s_fn = lambda mu: fd(mu, eps[0], kt)  # noqa: E731

    def parts(mu):
        return theta_s_fn(mu), fd(mu, eps[1], kt), fd(mu, eps[2], kt)

    def total(mu):
        ts, tss, tb = parts(mu)
        return geom['N_s'] * ts + geom['N_ss'] * tss + geom['N_b'] * tb

    lo, hi = -3.0, 3.0
    for _ in range(40):                 # bracket even extreme temperatures
        if total(lo) <= vo_total:
            break
        lo *= 2.0
    for _ in range(40):
        if total(hi) >= vo_total:
            break
        hi *= 2.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if total(mid) < vo_total:
            lo = mid
        else:
            hi = mid
    mu = 0.5 * (lo + hi)
    ts, tss, tb = parts(mu)
    us = geom['N_s'] * ts
    uss = geom['N_ss'] * tss
    ub = geom['N_b'] * tb
    tot = us + uss + ub
    mol_ti = 1e6 / p['molar_mass_g_mol']
    x = vo_total / mol_ti
    sat = p['saturation']
    onset_lo = p['surface_ordering']['critical_coverage']
    onset_hi = p['surface_ordering']['reconstruction_coverage']
    warnings = []
    if ts >= onset_lo:
        relation = ('is inside' if ts < onset_hi else 'exceeds')
        warnings.append({
            'kind': 'surface_reconstruction',
            'text': ('Surface bridging-O vacancy coverage %.1f%% %s the '
                     'reported %.0f-%.0f%% onset interval for the rutile '
                     'TiO2(110)-(1x2) reconstruction. The value is not '
                     'clipped. This distribution is a conditional '
                     'extrapolation of the unreconstructed (1x1) surface '
                     'Hamiltonian; reconstructed-phase energetics are not '
                     'implemented.'
                     % (ts * 100, relation, onset_lo * 100,
                        onset_hi * 100))})
    if tss > sat['subsurface_dilute_warn']:
        warnings.append({'kind': 'subsurface_dense',
                         'text': 'Subsurface layer occupancy %.0f%% is far '
                                 'beyond the dilute-defect regime; read this '
                                 'as a heavily reduced near-surface shell '
                                 '(extended defects), not isolated vacancies.'
                                 % (tss * 100)})
    if x >= sat['x_max_shear']:
        warnings.append({'kind': 'x_max',
                         'text': 'x = %.4f in TiO2-x exceeds the rutile '
                                 'point-defect range (~%.3f); crystallographic '
                                 'shear planes are expected.'
                                 % (x, sat['x_max_shear'])})
    elif x >= 0.8 * sat['x_max_shear']:
        warnings.append({'kind': 'x_near',
                         'text': 'x = %.4f in TiO2-x is near the rutile '
                                 'point-defect range (~%.3f).'
                                 % (x, sat['x_max_shear'])})
    if vo_total > 0 and abs(tot - vo_total) > 1e-6 * vo_total:
        warnings.append({'kind': 'unresolved',
                         'text': 'Inventory not matched within tolerance; '
                                 'check the isotherm range.'})
    return {'method': 'canonical_statmech',
            'T_C': t_c, 'VO_total_umol_g': vo_total,
            'mu_V_eV': mu,
            'theta': {'surface': ts, 'subsurface': tss, 'bulk': tb},
            'umol_g': {'surface': us, 'subsurface': uss, 'bulk': ub},
            'fractions': {'surface': us / tot if tot else 0.0,
                          'subsurface': uss / tot if tot else 0.0,
                          'bulk': ub / tot if tot else 0.0},
            'matched_umol_g': tot,
            'x_TiO2mx': x, 'Ti3_frac': 2.0 * x,
            'surface_reconstruction_regime': ts >= onset_lo,
            'warnings': warnings}


def spacing_summary(points, theta_s):
    """Gap statistics of the scan point nearest the operating coverage.

    The literature constraints (modal spacing 4-5 sites, 1-2 suppressed) are
    statements about the arrangement at a given coverage, so the comparison
    point is chosen by theta_s. Solutions past the reconstruction onset remain
    conditional on the unreconstructed surface Hamiltonian."""
    best = None
    for q in points:
        if q['mu_eV'] is None:
            continue
        d = abs(q['theta_s'] - theta_s)
        if best is None or d < best[0]:
            best = (d, q)
    if best is None:
        return None
    q = best[1]
    tot = 0
    mode = 0
    mode_n = -1
    near = 0
    for gap in range(1, len(q['hist'])):
        c = q['hist'][gap]
        tot += c
        if c > mode_n:
            mode_n = c
            mode = gap
        if gap <= 2:
            near += c
    if tot == 0:
        return None
    return {'at_theta_s': q['theta_s'], 'modal_gap_sites': mode,
            'P_gap_le_2': near / tot, 'pairs_sampled': tot}


# ------------------------------------------------- thermodynamic bridge

def phase_validity(active, reduced_costs, is_reduced):
    """Tier the gas-solid equilibrium verdict for the rutile defect model.

    invalid  - the equilibrium assemblage is a reduced phase with no
               rutile: the host lattice itself is not the stable phase,
               so a fixed-inventory rutile description is metastable at
               best.
    boundary - rutile coexists with a reduced phase: the charge sits on
               the reduction boundary (this includes the trace-buffer
               states under near-inert gas).
    stable   - rutile alone; the margin is the smallest reduced-phase
               KKT reduced cost per mole of oxygen.

    Pure function of the equilibrium result, ported verbatim to the
    browser; the caller supplies which phases count as reduced."""
    red_active = [s for s in active if is_reduced.get(s)]
    if red_active:
        tier = 'boundary' if 'TiO2' in active else 'invalid'
        return {'tier': tier, 'reduced_active': red_active,
                'nearest_phase': red_active[0],
                'margin_per_mol_O_kJ': 0.0 if tier == 'boundary' else None}
    best_s = None
    best_r = None
    for s in sorted(reduced_costs):
        if not is_reduced.get(s):
            continue
        rc = reduced_costs[s]
        r = rc['per_mol_O_kJ'] if isinstance(rc, dict) else rc
        if r is None or math.isinf(r):
            continue
        if best_r is None or r < best_r:
            best_r = float(r)
            best_s = s
    return {'tier': 'stable', 'reduced_active': [],
            'nearest_phase': best_s, 'margin_per_mol_O_kJ': best_r}


# ------------------------------------------------------------ CO2 layer

def co2_kinetics_params(p, kin=None):
    """Kinetics parameters with per-call overrides; all placeholders until
    fitted - the JSON note says so and the UI repeats it."""
    base = p['co2_kinetics']
    d = {'Ea_eV': dict(base['Ea_eV']),
         'A_eff_per_s_atm': dict(base['A_eff_per_s_atm']),
         'p_order_m': base['p_order_m'],
         'exposure_s': base['defaults']['exposure_s'],
         'p_CO2_atm': base['defaults']['p_CO2_atm'],
         'T_reox_C': base['defaults']['T_reox_C']}
    if kin:
        for k, v in kin.items():
            if k in ('Ea_eV', 'A_eff_per_s_atm'):
                d[k].update(v)
            else:
                d[k] = v
    return d


def accessibility(p, umol_by_class, kin=None):
    """First-order refill of each site class during one CO2 exposure:
    k_c = A_c p^m exp(-Ea_c/kT), P_c(t) = 1 - exp(-k_c t)."""
    d = co2_kinetics_params(p, kin)
    kt = KB_EV * (d['T_reox_C'] + 273.15)
    pf = d['p_CO2_atm'] ** d['p_order_m']
    probs = {}
    acc = 0.0
    tot = 0.0
    for c in ('surface', 'subsurface', 'bulk'):
        k = d['A_eff_per_s_atm'][c] * pf * math.exp(-d['Ea_eV'][c] / kt)
        x = k * d['exposure_s']
        prob = 1.0 if x > 700.0 else 1.0 - math.exp(-x)
        probs[c] = prob
        acc += umol_by_class[c] * prob
        tot += umol_by_class[c]
    return {'P': probs, 'accessible_umol_g': acc,
            'f_recoverable': acc / tot if tot else 0.0,
            'params': d}


def dielectric(p, vo_total):
    """Empirical loss factor measured in this work.

    The fit is a linear regression of log10(eps'') on log10(VO), so the
    correlation is the power law eps'' = 10**a * VO**b; outside the fitted
    inventory range the value is flagged as an extrapolation."""
    c = p['dielectric_correlation']
    if c.get('log10_intercept') is None or c.get('exponent') is None \
            or vo_total <= 0:
        return None
    e2 = 10.0 ** (c['log10_intercept']
                  + c['exponent'] * math.log10(vo_total))
    rng = c.get('fit_range_umol_g')
    out = {'eps2': e2,
           'extrapolated': bool(rng and (vo_total < rng[0]
                                         or vo_total > rng[1]))}
    if c.get('eps_prime'):
        out['tan_delta'] = e2 / c['eps_prime']
    return out


def sweep(p, t_c, geom, vo_list, theta_s_fn=None, preset=None, eps=None,
          kin=None):
    """Distribution and accessibility along a total-inventory axis.

    Instant once the isotherm exists - each point is one scalar matching
    solve. The surface follows the sampled isotherm without a hard coverage
    cap; subsurface and bulk take the remaining inventory, and the accessible
    curve falls away from the total."""
    rows = []
    for vo in vo_list:
        d = distribute(p, t_c, vo, geom, theta_s_fn=theta_s_fn,
                       preset=preset, eps=eps)
        a = accessibility(p, d['umol_g'], kin)
        rows.append({'VO_total_umol_g': vo,
                     'surface': d['umol_g']['surface'],
                     'subsurface': d['umol_g']['subsurface'],
                     'bulk': d['umol_g']['bulk'],
                     'accessible': a['accessible_umol_g'],
                     'f_recoverable': a['f_recoverable'],
                     'mu_V_eV': d['mu_V_eV']})
    return rows


def run_full(p, t_c=None, vo_total=None, d_um=None, bet_m2_g=None,
             ss_layers=None, preset=None, eps=None, seed=None, quality=1.0,
             zero_pairs=False, ordering=None, mc=None, kin=None,
             progress=None):
    """Isotherm scan plus particle-scale matching in one call."""
    if t_c is None:
        t_c = p['defaults']['T_C']
    if vo_total is None:
        vo_total = p['defaults']['VO_total_umol_g']
    geom = geometry(p, d_um=d_um, bet_m2_g=bet_m2_g, ss_layers=ss_layers)
    kt = KB_EV * (t_c + 273.15)
    if eps is None:
        eps = energetics(p, preset)['eps']
    pts = isotherm_scan(p, t_c, seed=seed, quality=quality, eps=eps,
                        zero_pairs=zero_pairs, ordering=ordering, mc=mc,
                        progress=progress)
    if zero_pairs:
        fn = lambda mu: fd(mu, eps[0], kt)  # exact ideal surface isotherm
    else:
        fn = theta_s_interp(pts, eps[0], kt)
    out = distribute(p, t_c, vo_total, geom, theta_s_fn=fn, eps=eps)
    out['geometry'] = geom
    out['isotherm'] = pts
    out['spacing'] = spacing_summary(pts, out['theta']['surface'])
    out['accessibility'] = accessibility(p, out['umol_g'], kin)
    return out
