"""Solve the model over its parameters for every sample: python3 run.py.

outputs/model_cases.csv       one row per sample, parameter point and reactive
                              definition: vacancies in every pool, reactive
                              sites, TOF, Q = TOF / TOF(R600) at the same point
outputs/sample_tof_range.csv  per sample: TOF and Q range over the parameters,
                              and the range of every pool
"""
import csv
import itertools
import statistics
import sys
from multiprocessing import Pool
from pathlib import Path

import numpy as np

from tofrange import model as md

HERE = Path(__file__).resolve().parent
FIXED_SITES = 2.31                 # SI Note 2a fixed denominator, umol/g
MIN_SITES, MIN_FRACTION = 0.01, 0.01
STRONG_REDUCTION_C = 900.0


def fmt(x):
    return f'{x:.6g}'


def samples():
    with open(HERE / 'experimental-data' / 'samples.csv') as f:
        rows = {r['sample']: r for r in csv.DictReader(f)}
    with open(HERE / 'experimental-data' / 'treatment_histories.csv') as f:
        for r in csv.DictReader(f):
            rows[r['sample']]['treatment_T_C'] = r['treatment_T_C']
    return list(rows.values())


def grid():
    return list(itertools.product(md.MAPS, md.CUTOFFS, md.DG, md.F110, md.EPS))


def point(args):
    """Every sample at one parameter point, mass-averaged over the diameters."""
    (name, cut, dG, f110, axis), S = args
    per = []
    for d in md.DIAMETERS:
        m, lay = md.particle(name, d, cutoff=cut, dG=dG, f110=f110, eps_r=md.EPS[axis])
        per.append([(md.reactive_sites(s, lay), md.populations(s, lay))
                    for s in (m.solve(float(x['inventory_umol_g']), md.T_EQ) for x in S)])
    c_bri900 = md.pt.bridging_capacity(md.DIAMETERS[0]) * f110
    rows = []
    for k, s in enumerate(S):
        sites = {r: float(np.mean([p[k][0][0][r] for p in per])) for r in md.REACTIVE}
        th = float(np.mean([p[k][0][1] for p in per]))
        f_rec = float(np.mean([p[k][0][2] for p in per]))
        pops = {q: float(np.mean([p[k][1][q] for p in per])) for q in md.POOLS}
        rate = float(s['rate_co_umol_g_s'])
        for rname, n in sites.items():
            below = n < MIN_SITES or n < MIN_FRACTION * c_bri900
            rows.append(dict(
                sample=s['sample'], energy_map=name, cutoff_nm=f'{cut:g}', dG_eV=f'{dG:+g}',
                f110=f'{f110:g}', eps=f'{md.EPS[axis]:g} ({axis})', reactive=rname,
                theta_bri=fmt(th), reconstructed_fraction=fmt(f_rec),
                **{f'N_{q}_umol_g': fmt(v) for q, v in pops.items()},
                N_react_umol_g=fmt(n), TOF_s_1=fmt(rate / n) if n > 0 else 'inf',
                site_class='BELOW_SITE_THRESHOLD' if below else 'OK',
                flags='YUAN_STRONG_REDUCTION' if float(s['treatment_T_C']) >= STRONG_REDUCTION_C else ''))
    return rows


def key(r):
    return tuple(r[k] for k in ('energy_map', 'cutoff_nm', 'dG_eV', 'f110', 'eps', 'reactive'))


def label(r):
    return f"{r['energy_map']} / cutoff {r['cutoff_nm']} nm / dG {r['dG_eV']} eV / f110 {r['f110']} / eps {r['eps']} / {r['reactive']}"


def add_q(rows):
    """Q = TOF / TOF(R600) at the same parameter point and reactive definition."""
    ref = {key(r): r for r in rows if r['sample'] == 'R600' and r['TOF_s_1'] != 'inf'}
    for r in rows:
        k = key(r)
        ok = r['TOF_s_1'] != 'inf' and k in ref
        r['Q_vs_R600'] = fmt(float(r['TOF_s_1']) / float(ref[k]['TOF_s_1'])) if ok else ''
        r['Q_ref_below_threshold'] = 'yes' if ok and ref[k]['site_class'] != 'OK' else ''


def summary(S, rows):
    out = []
    for s in S:
        mine = [r for r in rows if r['sample'] == s['sample']]
        rec = dict(sample=s['sample'], rate_co_umol_g_s=s['rate_co_umol_g_s'],
                   inventory_umol_g=s['inventory_umol_g'],
                   fixed_TOF_SI2a_s_1=fmt(float(s['rate_co_umol_g_s']) / FIXED_SITES))
        ok = sorted((r for r in mine if r['site_class'] == 'OK'), key=lambda r: float(r['TOF_s_1']))
        v = [float(r['TOF_s_1']) for r in ok]
        rec.update(n_cases=len(mine), n_ok=len(ok), TOF_min_s_1=fmt(v[0]), TOF_median_s_1=fmt(statistics.median(v)),
                   TOF_max_s_1=fmt(v[-1]), TOF_at_min=label(ok[0]), TOF_at_max=label(ok[-1]))
        qs = sorted((r for r in ok if r['Q_vs_R600'] and not r['Q_ref_below_threshold']),
                    key=lambda r: float(r['Q_vs_R600']))
        if qs:
            rec.update(Q_min=qs[0]['Q_vs_R600'], Q_max=qs[-1]['Q_vs_R600'],
                       Q_at_min=label(qs[0]), Q_at_max=label(qs[-1]))
        below = [r for r in mine if r['site_class'] != 'OK']
        vb = sorted(float(r['TOF_s_1']) for r in below if r['TOF_s_1'] != 'inf')
        qb = sorted(float(r['Q_vs_R600']) for r in below if r['Q_vs_R600'])
        rec['n_below_threshold'] = len(below)
        if vb:
            rec.update(below_TOF_min_s_1=fmt(vb[0]), below_TOF_max_s_1=fmt(vb[-1]))
        if qb:
            rec.update(below_Q_min=fmt(qb[0]), below_Q_max=fmt(qb[-1]))
        pts = [r for r in mine if r['reactive'] == 'BRI']
        for q in md.POOLS + ('react_BRI',):
            col = 'N_react_umol_g' if q == 'react_BRI' else f'N_{q}_umol_g'
            vals = sorted(float(r[col]) for r in pts)
            rec[f'{q}_min_umol_g'], rec[f'{q}_median_umol_g'], rec[f'{q}_max_umol_g'] = (
                fmt(vals[0]), fmt(statistics.median(vals)), fmt(vals[-1]))
        out.append(rec)
    return out


def write(name, rows):
    keys = list(dict.fromkeys(k for r in rows for k in r))
    with open(HERE / 'outputs' / name, 'w', newline='') as f:
        w = csv.DictWriter(f, keys, lineterminator='\n', restval='')
        w.writeheader(); w.writerows(rows)


def solve_all(points=None, procs=4):
    S = samples()
    points = grid() if points is None else points
    with Pool(procs) as pool:
        rows = [r for block in pool.imap(point, [(p, S) for p in points], chunksize=1) for r in block]
    return S, rows


def main():
    S, rows = solve_all()
    add_q(rows)
    write('model_cases.csv', rows)
    write('sample_tof_range.csv', summary(S, rows))


if __name__ == '__main__':
    sys.exit(main())
