"""Apparent TOF range for every sample: python3 run.py (writes outputs/).

outputs/cases.csv               every case: inputs, coverage, reactive sites, TOF, class, flags
outputs/sample_tof_range.csv    per sample: PRIMARY range, with sensitivities, boundary bound
outputs/sensitivity_effects.csv how far each one-at-a-time change moves the TOF
outputs/transport_gate.csv      time for the surface to reach its 600 C state
"""
import csv
import math
import statistics
from pathlib import Path

from tofrange.build import build
from tofrange.kinetics import BULK_BARRIERS, SURFACE_SETS, relax
from tofrange.scenarios import CLOSURES, MAPS, T_EQ, evaluate, run_all

HERE = Path(__file__).resolve().parent
FIXED_SITES = 2.31                 # SI Note 2a fixed denominator, umol/g
T_OBS = (1.0, 10.0, 60.0, 600.0)   # initial-rate observation time grid (v12, assumed)


def fmt(x):
    return f'{x:.6g}'


def write(name, rows, keys=None):
    keys = keys or list(dict.fromkeys(k for r in rows for k in r))
    with open(HERE / 'outputs' / name, 'w', newline='') as f:
        w = csv.DictWriter(f, keys, lineterminator='\n', restval='')
        w.writeheader(); w.writerows(rows)


def label(r):
    return '/'.join((r['diameter_nm'], r['energy_map'], r['closure'], r['reactive'], r['variant'] or '-'))


def base_key(r):
    m = r['energy_map'].split(' xi')[0]
    c = 'LOCAL' if r['closure'].startswith('GLOBAL') else r['closure']
    return r['sample'], m, c, r['reactive'].replace('+BASAL', '')


def summary(S, rows):
    """Per sample. core = the 40 PRIMARY-family cases (900 nm, sourced maps,
    NEUTRAL/LOCAL, four reactive definitions); all = every family except
    BELOW_SITE_THRESHOLD cases, whose TOF and Q range is reported apart.
    Flags say which values are lower bounds."""
    out = []
    for s in S:
        rate = float(s['rate_co_umol_g_s'])
        rec = dict(sample=s['sample'], rate_co_umol_g_s=s['rate_co_umol_g_s'],
                   fixed_TOF_SI2a_s_1=fmt(rate / FIXED_SITES))
        mine = [r for r in rows if r['sample'] == s['sample'] and r['variant'] != 'R600 94.6']
        if not mine:
            out.append(dict(rec, status='MISSING_INPUT: inventory'))
            continue
        rec['status'] = 'CALCULATED'
        valid = [r for r in mine if r['cls'] != 'BELOW_SITE_THRESHOLD']
        for tag, sel in (('core', [r for r in valid if r['family'] == 'PRIMARY']), ('all', valid)):
            sel = sorted(sel, key=lambda r: float(r['TOF_s_1']))
            v = [float(r['TOF_s_1']) for r in sel]
            rec[f'{tag}_n'] = len(sel)
            rec[f'{tag}_TOF_min_s_1'] = fmt(v[0])
            rec[f'{tag}_TOF_median_s_1'] = fmt(statistics.median(v))
            rec[f'{tag}_TOF_max_s_1'] = fmt(v[-1])
            q = sorted(float(r['Q_vs_R600']) for r in sel if r['Q_vs_R600'] and not r['Q_ref_below_threshold'])
            if q:
                rec[f'{tag}_Q_min'] = fmt(q[0])
                rec[f'{tag}_Q_max'] = fmt(q[-1])
            rec[f'{tag}_at_min'] = label(sel[0])
            rec[f'{tag}_at_max'] = label(sel[-1])
            flags = {}
            for r in sel:
                for f in r['flags'].split():
                    flags[f] = flags.get(f, 0) + 1
            rec[f'{tag}_flags'] = ' '.join(f'{f}:{n}/{len(sel)}' for f, n in sorted(flags.items()))
        core = [r for r in mine if r['family'] == 'PRIMARY']
        raw = sorted(float(r['TOF_uncapped_s_1']) for r in core if r['TOF_uncapped_s_1'])
        rec['core_uncapped_TOF_min_s_1'] = fmt(raw[0])
        rec['core_uncapped_TOF_max_s_1'] = fmt(raw[-1])
        below = [r for r in mine if r['cls'] == 'BELOW_SITE_THRESHOLD']
        rec['below_threshold_n'] = len(below)
        v = sorted(float(r['TOF_s_1']) for r in below if r['TOF_s_1'] != 'inf')
        q = sorted(float(r['Q_vs_R600']) for r in below if r['Q_vs_R600'])
        if v:
            rec['below_threshold_TOF_min_s_1'], rec['below_threshold_TOF_max_s_1'] = fmt(v[0]), fmt(v[-1])
        if q:
            rec['below_threshold_Q_min'], rec['below_threshold_Q_max'] = fmt(q[0]), fmt(q[-1])
        inv = {}
        for r in below:
            inv[r['family']] = inv.get(r['family'], 0) + 1
        rec['below_threshold_by_family'] = ' '.join(f'{f}:{n}' for f, n in sorted(inv.items()))
        out.append(rec)
    return out


def effects(rows):
    """Change against the matched core case, in decades. BELOW_SITE_THRESHOLD
    cases are counted instead (their near-zero denominators would set the range)."""
    prim = {base_key(r): float(r['TOF_s_1']) for r in rows if r['family'] == 'PRIMARY'}
    groups = {}
    for r in rows:
        if r['family'] == 'PRIMARY':
            continue
        g = groups.setdefault((r['sample'], r['family'], r['variant']), dict(d=[], below=0))
        k = base_key(r)
        if r['cls'] == 'BELOW_SITE_THRESHOLD':
            g['below'] += 1
        elif k in prim:
            g['d'].append(math.log10(float(r['TOF_s_1']) / prim[k]))
    out = []
    for (s, f, v), g in sorted(groups.items()):
        d = g['d']
        out.append(dict(sample=s, family=f, variant=v, n=len(d), below_threshold=g['below'],
                        median_log10_change=f'{statistics.median(d):+.3f}' if d else '',
                        min_log10_change=f'{min(d):+.3f}' if d else '',
                        max_log10_change=f'{max(d):+.3f}' if d else ''))
    return out


def transport_gate(S):
    """Relax from equilibrium at the treatment temperature, hold at 600 C.

    Where the surface has not reached its 600 C coverage (within 1 %) by an
    observation time, that coverage becomes a transport case of its own.
    """
    gate, rows = [], []
    moved = [s for s in S if s['inventory_status'] == 'SELECTED'
             and abs(float(s['treatment_T_C']) + 273.15 - T_EQ) > 1e-9]
    for name in MAPS:
        for closure in CLOSURES:
            kw = dict(n_cont=600, z0_cont=1e-2) if name == 'LI_CONT' else {}
            model, lay = build(name, closure, 900.0, **kw)
            spec = dict(name=name, closure=closure, d_nm=900.0)
            sets = list(SURFACE_SETS) if lay.kind == 'discrete' else ['-']
            for sset in sets:
                for bname, B in BULK_BARRIERS.items():
                    for s in moved:
                        T0 = float(s['treatment_T_C']) + 273.15
                        r = relax(model, lay, float(s['inventory_umol_g']), T0, T_EQ, B,
                                  surface_set=sset if sset != '-' else 'WU', t_obs=T_OBS)
                        g = dict(sample=s['sample'], energy_map=name, closure=closure,
                                 surface_edges=sset, bulk_barrier=f'{bname} {B:g} eV',
                                 coverage_at_treatment_T=fmt(r['coverage_start']),
                                 coverage_eq_600C=fmt(r['coverage_eq']),
                                 t_within_1pct_s=fmt(r['t_1pct']),
                                 mass_drift=f"{r['mass_drift']:.1e}")
                        for t in T_OBS:
                            g[f'coverage_at_{t:g}s'] = fmt(r['coverage'][t])
                            if abs(r['coverage'][t] / r['coverage_eq'] - 1) > 0.01:
                                rows += evaluate(None, lay, s, 'transport',
                                                 f'{sset} {bname} from {s["treatment_T_C"]} C, t {t:g} s',
                                                 spec, theta=r['coverage'][t])
                        gate.append(g)
    return gate, rows


def add_q(rows):
    """Q = TOF / TOF(R600) within the same scenario (Note 11.2), R600 at 94.0."""
    key = lambda r: (r['family'], r['variant'], r['diameter_nm'], r['energy_map'], r['closure'], r['reactive'])  # noqa: E731
    ok = lambda r: r['TOF_s_1'] != 'inf'                                   # noqa: E731
    R600 = [r for r in rows if r['sample'] == 'R600' and ok(r) and r['variant'] != 'R600 94.6']
    ref = {key(r): float(r['TOF_s_1']) for r in R600}
    ref_below = {key(r) for r in R600 if r['cls'] == 'BELOW_SITE_THRESHOLD'}
    for r in rows:
        k = key(r)
        r['Q_vs_R600'] = fmt(float(r['TOF_s_1']) / ref[k]) if ok(r) and k in ref else ''
        r['Q_ref_below_threshold'] = 'yes' if r['Q_vs_R600'] and k in ref_below else ''


def main():
    S, rows = run_all()
    gate, lagging = transport_gate(S)
    rows += lagging
    add_q(rows)
    write('cases.csv', rows)
    write('sample_tof_range.csv', summary(S, rows))
    write('sensitivity_effects.csv', effects(rows))
    write('transport_gate.csv', gate)


if __name__ == '__main__':
    main()
