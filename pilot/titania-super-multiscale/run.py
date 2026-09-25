"""Apparent TOF range for every sample: python3 run.py (writes outputs/)."""
import csv
import math
import statistics
from pathlib import Path

from tofrange import particle as pt
from tofrange.scenarios import scenarios

HERE = Path(__file__).resolve().parent
T = 873.15                         # common equilibrium temperature (Note 2b)
FIXED_SITES = 2.31                 # SI Note 2a fixed denominator, umol/g
REACTIVE = {'BRI': None, 'ISO_z2': 2, 'ISO_z4': 4, 'ISO_z8': 8}
# v12 boundary rule 2: a denominator below 0.01 umol/g or below 1 % of the
# bridging capacity gives a TOF set by rounding, so the case is INVALID.
MIN_SITES, MIN_FRACTION = 0.01, 0.01


def samples():
    with open(HERE / 'experimental-data' / 'samples.csv') as f:
        return list(csv.DictReader(f))


def fmt(x):
    return f'{x:.6g}'


def main():
    rows, S = [], samples()
    measured = [s for s in S if s['inventory_status'] == 'SELECTED']
    for emap, closure, d_nm, build, theta_of in scenarios():
        model = build()
        cb = pt.bridging_capacity(d_nm)
        for s in measured:
            sol = model.solve(float(s['inventory_umol_g']), T)
            th = theta_of(sol, T)
            for rname, z in REACTIVE.items():
                n_react = cb * th * (1 if z is None else (1 - th) ** z)
                ok = n_react >= MIN_SITES and n_react >= MIN_FRACTION * cb
                rows.append(dict(sample=s['sample'], status='VALID' if ok else 'INVALID', diameter_nm=fmt(d_nm), energy_map=emap, closure=closure,
                                 reactive=rname, theta_bri=fmt(th), N_react_umol_g=fmt(n_react),
                                 TOF_s_1=fmt(float(s['rate_co_umol_g_s']) / n_react),
                                 iterations=sol.iterations, inventory_residual=f'{sol.inventory_residual:.1e}',
                                 charge_residual=f'{sol.charge_residual:.1e}'))
    out = HERE / 'outputs'
    with open(out / 'scenario_tof.csv', 'w', newline='') as f:
        w = csv.DictWriter(f, rows[0].keys(), lineterminator='\n'); w.writeheader(); w.writerows(rows)

    dims = ('diameter_nm', 'energy_map', 'closure', 'reactive')
    summary = []
    for s in S:
        r = float(s['rate_co_umol_g_s'])
        base = dict(sample=s['sample'], rate_co_umol_g_s=s['rate_co_umol_g_s'],
                    fixed_TOF_SI2a_s_1=fmt(r / FIXED_SITES))
        mine = [x for x in rows if x['sample'] == s['sample'] and x['status'] == 'VALID']
        n_invalid = sum(x['sample'] == s['sample'] and x['status'] == 'INVALID' for x in rows)
        if not mine:
            summary.append(dict(base, status='MISSING_INPUT: inventory'))
            continue
        tof = sorted(mine, key=lambda x: float(x['TOF_s_1']))
        vals = [float(x['TOF_s_1']) for x in tof]
        rec = dict(base, status='CALCULATED', n_valid=len(vals), n_invalid=n_invalid, TOF_min_s_1=fmt(vals[0]),
                   TOF_median_s_1=fmt(statistics.median(vals)), TOF_max_s_1=fmt(vals[-1]),
                   at_min='/'.join(tof[0][d] for d in dims), at_max='/'.join(tof[-1][d] for d in dims))
        # Span each assumption controls: median, over all other settings, of
        # log10(max/min) when only that assumption varies.
        for d in dims:
            groups = {}
            for x in mine:
                groups.setdefault(tuple(x[o] for o in dims if o != d), []).append(float(x['TOF_s_1']))
            spans = [math.log10(max(g) / min(g)) for g in groups.values() if len(g) > 1]
            rec['decades_' + d] = f'{statistics.median(spans):.3f}'
        summary.append(rec)
    keys = list(dict.fromkeys(k for x in summary for k in x))
    with open(out / 'sample_tof_range.csv', 'w', newline='') as f:
        w = csv.DictWriter(f, keys, lineterminator='\n'); w.writeheader(); w.writerows(summary)


if __name__ == '__main__':
    main()
