"""Every case, its reactive-site count, its apparent TOF and its class.

PRIMARY cases use sourced inputs only:
    900 nm sphere; energy maps PAB, HAM, LI_SBR1, LI_L2, LI_CONT (xi 0.5 nm);
    NEUTRAL and LOCAL closures; no bulk pairs; 600 C equilibrium;
    reactive sites BRI, ISO_z2, ISO_z4, ISO_z8.
SENSITIVITY families change one thing at a time on every PRIMARY model:
    diameter      600, 300 nm (v12 Q1 scan; no measured size distribution)
    decay_length  LI_CONT xi 0.25, 1 nm (manuscript Note 2b scan)
    global        GLOBAL closure, eps 64 (a axis, fields normal to (110)) or 107
                  (c axis), Parker 1961 Fig. 1 at 873 K; S0 penalty 0.2 eV
    pairs_eq      ZHA2017 bulk pairs, equilibrated at 600 C
    pairs_frozen  ZHA2017 bulk pairs formed at the treatment temperature and
                  frozen; everything else re-equilibrates at 600 C
    basal         alpha = 1: layer-1 in-plane vacancies also count (discrete maps)
    inventory     R600 at the Origin alternative, 94.6 umol/g
Surface limit (every case)
    An unreconstructed (110) surface holds at most 17 % bridging vacancies
    (BIR2024; manuscript Note 2a, 0.17 ML). Above it the coverage used for
    reactive sites is capped at 0.17 and the case is flagged RECON_CAP: the
    TOF is then a lower bound, because reconstruction can only remove sites.
Boundary flags (class BOUNDARY, value is a lower bound on TOF)
    YUAN_STRONG_REDUCTION  treatment at or above 900 C (Yuan 2024: Ti2O3-(1x2),
                           kinetically trapped)
    PAIRS_BEYOND_DILUTE    more than half the eligible bulk vacancies paired:
                           larger clusters expected, which only lower the
                           surface population further
INVALID: fewer than 0.01 umol/g or 1 % of the bridging capacity reactive.
"""
import csv
from pathlib import Path

from .build import build, surface_coverage, pair_fraction

HERE = Path(__file__).resolve().parent.parent
T_EQ = 873.15
THETA_C = 0.17
MIN_SITES, MIN_FRACTION = 0.01, 0.01
STRONG_REDUCTION_C = 900.0
PAIR_BOUNDARY = 0.5
MAPS = ('PAB', 'HAM', 'LI_SBR1', 'LI_L2', 'LI_CONT')
CLOSURES = ('NEUTRAL', 'LOCAL')
REACTIVE = {'BRI': None, 'ISO_z2': 2, 'ISO_z4': 4, 'ISO_z8': 8}
EPS = {'a_axis': 64.0, 'c_axis': 107.0}


def samples():
    with open(HERE / 'experimental-data' / 'samples.csv') as f:
        rows = {r['sample']: r for r in csv.DictReader(f)}
    with open(HERE / 'experimental-data' / 'treatment_histories.csv') as f:
        for r in csv.DictReader(f):
            rows[r['sample']]['treatment_T_C'] = r['treatment_T_C']
    return list(rows.values())


def model_specs():
    """(family, label, kwargs for build) for every equilibrium model."""
    base = [dict(name=m, closure=c, d_nm=900.0) for m in MAPS for c in CLOSURES]
    out = [('PRIMARY', '', b) for b in base]
    for d in (600.0, 300.0):
        out += [('diameter', f'{d:g} nm', dict(b, d_nm=d)) for b in base]
    for xi in (0.25, 1.0):
        out += [('decay_length', f'xi {xi:g} nm', dict(b, xi=xi)) for b in base if b['name'] == 'LI_CONT']
    for axis, eps in EPS.items():
        out += [('global', f'eps {eps:g} ({axis})', dict(name=m, closure='GLOBAL', d_nm=900.0, eps_r=eps))
                for m in MAPS]
    out += [('pairs_eq', 'ZHA2017 at 600 C', dict(b, pairs=True)) for b in base]
    return out


def reactive_sites(theta, c_bri, z):
    th = min(theta, THETA_C)
    return c_bri * th * (1.0 if z is None else (1.0 - th) ** z)


def valid_count(n, c_bri):
    return n >= MIN_SITES and n >= MIN_FRACTION * c_bri


def evaluate(sol, lay, sample, family, label, spec, pf=0.0, theta=None):
    """Cases (one per reactive definition) for one solved sample.

    theta overrides the equilibrium coverage (transport cases pass the
    coverage reached at the observation time).
    """
    if theta is None:
        theta = surface_coverage(sol, lay)
    rate = float(sample['rate_co_umol_g_s'])
    T_treat = float(sample['treatment_T_C'])
    rows = []
    defs = dict(REACTIVE)
    if family == 'basal':
        defs = {'BRI+BASAL': None}
    for rname, z in defs.items():
        n = reactive_sites(theta, lay.c_bri, z)
        n_raw = lay.c_bri * theta * (1.0 if z is None else (1.0 - theta) ** z)
        if family == 'basal':
            basal = sol.occupancy(lay.basal) * sol.model.C[lay.basal]
            n, n_raw = n + basal, n_raw + basal
        flags = []
        if theta > THETA_C:
            flags.append('RECON_CAP')
        if T_treat >= STRONG_REDUCTION_C:
            flags.append('YUAN_STRONG_REDUCTION')
        if pf > PAIR_BOUNDARY:
            flags.append('PAIRS_BEYOND_DILUTE')
        if not valid_count(n, lay.c_bri):
            cls = 'INVALID'
        elif 'YUAN_STRONG_REDUCTION' in flags or 'PAIRS_BEYOND_DILUTE' in flags:
            cls = 'BOUNDARY'
        else:
            cls = 'PRIMARY' if family == 'PRIMARY' else 'SENSITIVITY'
        rows.append(dict(
            sample=sample['sample'], family=family, variant=label,
            diameter_nm=f"{spec['d_nm']:g}", energy_map=spec['name']
            + (f" xi {spec.get('xi', 0.5):g}" if spec['name'] == 'LI_CONT' else ''),
            closure=spec['closure'] + (f" eps {spec['eps_r']:g}" if spec.get('eps_r') else ''),
            reactive=rname, cls=cls, flags=' '.join(flags),
            theta_bri=f'{theta:.6g}', bridging_capacity_umol_g=f'{lay.c_bri:.6g}',
            N_react_umol_g=f'{n:.6g}', TOF_s_1=f'{rate / n:.6g}',
            N_react_uncapped_umol_g=f'{n_raw:.6g}',
            TOF_uncapped_s_1=f'{rate / n_raw:.6g}' if valid_count(n_raw, lay.c_bri) else '',
            pair_fraction=f'{pf:.4g}', iterations=sol.iterations if sol is not None else '',
            inventory_residual=f'{sol.inventory_residual:.1e}' if sol is not None else ''))
    return rows


def run_all(progress=None):
    S = samples()
    measured = [s for s in S if s['inventory_status'] == 'SELECTED']
    rows = []
    for family, label, spec in model_specs():
        model, lay = build(**spec)
        sols = {s['sample']: model.solve(float(s['inventory_umol_g']), T_EQ) for s in measured}
        for s in measured:
            rows += evaluate(sols[s['sample']], lay, s, family, label, spec,
                             pair_fraction(sols[s['sample']], lay))
        if spec.get('pairs'):
            rows += frozen_pairs(spec, measured)
        if family == 'PRIMARY':
            r600 = next(s for s in measured if s['sample'] == 'R600')
            alt = dict(r600, inventory_umol_g=r600['inventory_alt_umol_g'])
            rows += evaluate(model.solve(float(alt['inventory_umol_g']), T_EQ), lay, alt,
                             'inventory', 'R600 94.6', spec)
            if lay.kind == 'discrete':
                for s in measured:
                    rows += evaluate(sols[s['sample']], lay, s, 'basal', 'alpha 1', spec)
        if progress:
            progress(family, label, spec)
    return S, rows


def frozen_pairs(spec, measured):
    """Pairs formed at the treatment temperature stay; the rest re-equilibrates."""
    rows = []
    plain, lay0 = build(**dict(spec, pairs=False))
    eligible = [o for o in lay0.o if o.get('eligible')]
    for s in measured:
        T_treat = float(s['treatment_T_C']) + 273.15
        if abs(T_treat - T_EQ) < 1e-9:
            continue
        m_t, lay_t = build(**dict(spec, T=T_treat))
        sol_t = m_t.solve(float(s['inventory_umol_g']), T_treat)
        v_fix, q_reg = 0.0, {}
        for i in lay_t.pair_idx:
            n_pairs = float(sol_t.x[i][1:].sum())
            v_fix += 2 * n_pairs
            if m_t.region[i] >= 0:
                r = m_t.region_ids[m_t.region[i]]
                q_reg[r] = q_reg.get(r, 0.0) + 4 * n_pairs
        sol = plain.solve(float(s['inventory_umol_g']), T_EQ, fixed=dict(v=v_fix, q_region=q_reg))
        mono = sum(sol.occupancy(o['idx']) * o['C'] for o in eligible)
        rows += evaluate(sol, lay0, s, 'pairs_frozen', f"ZHA2017 frozen at {s['treatment_T_C']} C",
                         spec, v_fix / (v_fix + mono))
    return rows
