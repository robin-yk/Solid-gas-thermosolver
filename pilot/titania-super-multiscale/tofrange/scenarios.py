"""Every case, its reactive-site count, its apparent TOF and its class.

PRIMARY cases use sourced inputs only:
    900 nm sphere; energy maps PAB, HAM, LI_SBR1, LI_L2, LI_CONT (xi 0.5 nm);
    NEUTRAL and LOCAL closures; no bulk pairs; 600 C equilibrium;
    reactive sites BRI, ISO_z2, ISO_z4, ISO_z8.
SENSITIVITY families change one thing at a time on every PRIMARY model:
    diameter      1250, 1600 nm (the samples' size range, 900-1600 nm)
    size_mix      equal mass per diameter over 900-1600 nm (8 diameters), same
                  inventory per gram in every size; reactive sites mass-weighted
    recon_fixed   explicit Ti2O3-(1x2) area fraction 0.25, 0.5, 0.75
    recon_state   explicit (1x2) state with relative energy dG = -0.4 ... +0.4 eV;
                  the reconstructed fraction is an output (discrete maps)
    aggregates    bulk aggregates of every size from Monte Carlo, pairwise-
                  additive ZHA2017 energy, cutoff scanned 0.28-1.0 nm
    facet         (110) share of the surface 0.75, 0.5; the rest has no
                  explicit sites (no sourced energies for other facets)
    matsunaga     Li maps with a basal site 0.11 eV above bridging
                  (Matsunaga 2014); basal vacancies also counted
    decay_length  LI_CONT xi 0.25, 1 nm (manuscript Note 2b scan)
    global        GLOBAL closure, eps 64 (a axis, fields normal to (110)) or 107
                  (c axis), Parker 1961 Fig. 1 at 873 K; S0 penalty 0.2 eV
    pairs_eq      ZHA2017 bulk pairs, equilibrated at 600 C
    pairs_frozen  ZHA2017 bulk pairs formed at the treatment temperature and
                  frozen; everything else re-equilibrates at 600 C
    basal         alpha = 1: layer-1 in-plane vacancies also count (discrete maps)
    inventory     R600 at the Origin alternative, 94.6 umol/g
COMBINED family (several changes at once, so their interaction is computed):
    combined      aggregates (cutoff 0.34, 0.40 nm) + (1x2) reconstruction state
                  (dG -0.2, 0, +0.2 eV) + equal-mass 900-1600 nm mixture,
                  discrete maps, NEUTRAL and LOCAL
Surface limit (every family except recon_fixed, recon_state and combined)
    An unreconstructed (110) surface holds at most 17 % bridging vacancies
    (BIR2024; manuscript Note 2a, 0.17 ML). Above it the coverage used for
    reactive sites is capped at 0.17 and the case is flagged RECON_CAP. This is
    a rule; it does not compute a reconstructed population. The recon families
    carry the reconstruction explicitly, uncapped, flagged ABOVE_17PCT where the
    unreconstructed part still exceeds 17 %.
Boundary flags (class BOUNDARY, value is a lower bound on TOF)
    YUAN_STRONG_REDUCTION  treatment at or above 900 C (Yuan 2024: Ti2O3-(1x2),
                           kinetically trapped)
    PAIRS_BEYOND_DILUTE    more than half the eligible bulk vacancies paired:
                           larger clusters expected, which only lower the
                           surface population further
BELOW_SITE_THRESHOLD: fewer than 0.01 umol/g or 1 % of the bridging capacity
reactive. This is a reporting threshold of this pilot. The TOF (rate / sites)
and Q are still computed and kept; such cases only stay out of the summary
range, which they would otherwise set by near-zero denominators.
"""
import csv
from pathlib import Path

import numpy as np

from .build import build, pair_fraction, surface_state

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
SIZES = (1250.0, 1600.0)
SIZE_MIX = tuple(np.linspace(900.0, 1600.0, 8))
RECON_F = (0.25, 0.5, 0.75)
RECON_DG = (-0.4, -0.2, 0.0, 0.2, 0.4)
MC_CUTOFFS = (0.28, 0.34, 0.40, 0.45, 0.6, 1.0)
FACET_110 = (0.75, 0.5)
MATSUNAGA = 0.11
EXPLICIT_RECON = ('recon_fixed', 'recon_state', 'combined')
COMBINED_CUTOFFS = (0.34, 0.40)
COMBINED_DG = (-0.2, 0.0, 0.2)


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
    for d in SIZES:
        out += [('diameter', f'{d:g} nm', dict(b, d_nm=d)) for b in base]
    for f in RECON_F:
        out += [('recon_fixed', f'f_rec {f:g}', dict(b, recon=('fixed', f))) for b in base]
    for dG in RECON_DG:
        out += [('recon_state', f'dG {dG:+g} eV', dict(b, recon=('state', dG)))
                for b in base if b['name'] != 'LI_CONT']
    for f in FACET_110:
        out += [('facet', f'f110 {f:g}', dict(b, f110=f)) for b in base]
    out += [('matsunaga', 'basal +0.11 eV', dict(b, basal_shift=MATSUNAGA))
            for b in base if b['name'] in ('LI_SBR1', 'LI_L2')]
    for cut in MC_CUTOFFS:
        out += [('aggregates', f'MC cutoff {cut:g} nm', dict(b, aggregates=cut)) for b in base]
    for xi in (0.25, 1.0):
        out += [('decay_length', f'xi {xi:g} nm', dict(b, xi=xi)) for b in base if b['name'] == 'LI_CONT']
    for axis, eps in EPS.items():
        out += [('global', f'eps {eps:g} ({axis})', dict(name=m, closure='GLOBAL', d_nm=900.0, eps_r=eps))
                for m in MAPS]
    out += [('pairs_eq', 'ZHA2017 at 600 C', dict(b, pairs=True)) for b in base]
    return out


def reactive_sites(theta, c_bri, z, cap=True):
    th = min(theta, THETA_C) if cap else theta
    return c_bri * th * (1.0 if z is None else (1.0 - th) ** z)


def valid_count(n, c_bri):
    return n >= MIN_SITES and n >= MIN_FRACTION * c_bri


def site_counts(sol, lay, family, theta=None):
    """Reactive sites per definition: {name: (capped or explicit, uncapped)}.

    Also returns (theta, reconstructed fraction). theta overrides the
    equilibrium coverage (transport cases pass the coverage at t_obs)."""
    if sol is None:
        th, c_react, f_rec = theta, lay.c_bri * (1 - lay.f_fixed), lay.f_fixed
    else:
        th, c_react, f_rec = surface_state(sol, lay)
    if theta is not None:
        th = theta
    cap = family not in EXPLICIT_RECON
    defs = {'BRI+BASAL': None} if family == 'basal' else dict(REACTIVE)
    if family == 'matsunaga':
        defs['BRI+BASAL'] = None
    out = {}
    for rname, z in defs.items():
        n = reactive_sites(th, c_react, z, cap)
        n_raw = reactive_sites(th, c_react, z, False)
        if rname == 'BRI+BASAL':
            basal = sol.occupancy(lay.basal) * sol.model.C[lay.basal]
            n, n_raw = n + basal, n_raw + basal
        out[rname] = (n, n_raw)
    return out, th, f_rec


def evaluate(sol, lay, sample, family, label, spec, pf=0.0, theta=None, counts=None):
    """Cases (one per reactive definition) for one solved sample."""
    if counts is None:
        counts = site_counts(sol, lay, family, theta)
    sites, th, f_rec = counts
    rate = float(sample['rate_co_umol_g_s'])
    T_treat = float(sample['treatment_T_C'])
    rows = []
    for rname, (n, n_raw) in sites.items():
        flags = []
        if th > THETA_C:
            flags.append('ABOVE_17PCT' if family in EXPLICIT_RECON else 'RECON_CAP')
        if T_treat >= STRONG_REDUCTION_C:
            flags.append('YUAN_STRONG_REDUCTION')
        if pf > PAIR_BOUNDARY:
            flags.append('PAIRS_BEYOND_DILUTE')
        if not valid_count(n, lay.c_bri):
            cls = 'BELOW_SITE_THRESHOLD'
        elif 'YUAN_STRONG_REDUCTION' in flags or 'PAIRS_BEYOND_DILUTE' in flags:
            cls = 'BOUNDARY'
        else:
            cls = 'PRIMARY' if family == 'PRIMARY' else 'SENSITIVITY'
        rows.append(dict(
            sample=sample['sample'], family=family, variant=label,
            diameter_nm=f"{spec['d_nm']:g}" if spec.get('d_nm') else 'mix 900-1600',
            energy_map=spec['name'] + (f" xi {spec.get('xi', 0.5):g}" if spec['name'] == 'LI_CONT' else ''),
            closure=spec['closure'] + (f" eps {spec['eps_r']:g}" if spec.get('eps_r') else ''),
            reactive=rname, cls=cls, flags=' '.join(flags),
            theta_bri=f'{th:.6g}', reconstructed_fraction=f'{f_rec:.4g}',
            bridging_capacity_umol_g=f'{lay.c_bri:.6g}',
            N_react_umol_g=f'{n:.6g}', TOF_s_1=f'{rate / n:.6g}' if n > 0 else 'inf',
            N_react_uncapped_umol_g=f'{n_raw:.6g}',
            TOF_uncapped_s_1=f'{rate / n_raw:.6g}' if valid_count(n_raw, lay.c_bri) else '',
            pair_fraction=f'{pf:.4g}', iterations=sol.iterations if sol is not None else '',
            inventory_residual=f'{sol.inventory_residual:.1e}' if sol is not None else ''))
    return rows


def size_mixture(measured, family='size_mix', label='equal mass, 900-1600 nm', maps=MAPS, **extra):
    """Equal mass at each diameter in SIZE_MIX; reactive sites add by mass.
    extra: further build options applied at every diameter (combined family)."""
    rows = []
    for m in maps:
        for c in CLOSURES:
            per = []
            for d in SIZE_MIX:
                model, lay = build(m, c, d, **extra)
                per.append([(site_counts(model.solve(float(s['inventory_umol_g']), T_EQ, fixed=lay.fixed),
                                         lay, family), lay) for s in measured])
            lay900 = per[0][0][1]
            for k, s in enumerate(measured):
                sites = {r: tuple(np.mean([p[k][0][0][r][i] for p in per]) for i in (0, 1)) for r in REACTIVE}
                th = float(np.mean([p[k][0][1] for p in per]))
                f_rec = float(np.mean([p[k][0][2] for p in per]))
                spec = dict(name=m, closure=c, d_nm=None)
                rows += evaluate(None, lay900, s, family, label, spec, counts=(sites, th, f_rec))
    return rows


def combined(measured):
    """Aggregates, the (1x2) reconstruction state and the 900-1600 nm mass
    mixture switched on together (discrete maps: the state needs explicit
    bridging sites)."""
    rows = []
    for cut in COMBINED_CUTOFFS:
        for dG in COMBINED_DG:
            rows += size_mixture(measured, 'combined', f'MC cutoff {cut:g} nm, dG {dG:+g} eV, 900-1600 nm',
                                 maps=[m for m in MAPS if m != 'LI_CONT'],
                                 aggregates=cut, recon=('state', dG))
    return rows


def run_all(progress=None):
    S = samples()
    measured = [s for s in S if s['inventory_status'] == 'SELECTED']
    rows = []
    for family, label, spec in model_specs():
        model, lay = build(**spec)
        sols = {s['sample']: model.solve(float(s['inventory_umol_g']), T_EQ, fixed=lay.fixed)
                for s in measured}
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
    rows += size_mixture(measured)
    rows += combined(measured)
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
