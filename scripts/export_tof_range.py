"""Pack the one-model results into paper_outputs/tof_range.json.

    python3 scripts/export_tof_range.py

Reads pilot/titania-super-multiscale/outputs/{model_cases,sample_tof_range}.csv
(written by that directory's run.py, about 40 min on four cores) and the
explicit-layer site capacities of the same model. The page reads this file
and computes nothing of its own: per-case numbers are rounded to four
significant figures, the per-sample ranges are copied as written.
"""
import csv
import json
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PILOT = os.path.join(ROOT, 'pilot', 'titania-super-multiscale')
OUT = os.path.join(ROOT, 'paper_outputs', 'tof_range.json')
sys.path.insert(0, PILOT)

from tofrange import model as md  # noqa: E402

RUN = __import__('run')


def g4(x):
    return float(f'{float(x):.4g}')


def capacities(f110):
    """Explicit-layer O sites (umol/g) of each pool, mass-averaged over the
    diameters the way run.point averages the pools."""
    acc = dict(bridging=0.0, basal_L1=0.0, L1_subbridging=0.0, subsurface_L2_4=0.0)
    for d in md.DIAMETERS:
        o_sites, _ = md.pt.layer_sites(d, md.K_EXPL)
        for k, site, c, _ in o_sites:
            pool = ('subsurface_L2_4' if k > 1 else
                    {'BRI': 'bridging', 'IPL': 'basal_L1', 'SBR': 'L1_subbridging'}[site])
            acc[pool] += c * f110 / len(md.DIAMETERS)
    return {k: g4(v) for k, v in acc.items()}


def main():
    with open(os.path.join(PILOT, 'outputs', 'model_cases.csv')) as f:
        rows = list(csv.DictReader(f))
    with open(os.path.join(PILOT, 'outputs', 'sample_tof_range.csv')) as f:
        summary = {r['sample']: r for r in csv.DictReader(f)}
    samples = RUN.samples()
    grid = RUN.grid()
    at = {(m, f'{c:g}', f'{g:+g}', f'{f:g}', f'{md.EPS[e]:g} ({e})'): i
          for i, (m, c, g, f, e) in enumerate(grid)}
    react = list(md.REACTIVE)
    cases = {s['sample']: dict(theta=[None] * len(grid), f_rec=[None] * len(grid),
                               pools=[None] * len(grid), sites=[[None] * len(react) for _ in grid],
                               below=[[0] * len(react) for _ in grid])
             for s in samples}
    for r in rows:
        i = at[(r['energy_map'], r['cutoff_nm'], r['dG_eV'], r['f110'], r['eps'])]
        c = cases[r['sample']]
        j = react.index(r['reactive'])
        c['theta'][i] = g4(r['theta_bri'])
        c['f_rec'][i] = g4(r['reconstructed_fraction'])
        c['pools'][i] = [g4(r[f'N_{q}_umol_g']) for q in md.POOLS]
        c['sites'][i][j] = g4(r['N_react_umol_g'])
        c['below'][i][j] = 0 if r['site_class'] == 'OK' else 1
    for name, c in cases.items():
        missing = [i for i, v in enumerate(c['theta']) if v is None]
        if missing:
            raise SystemExit(f'{name}: {len(missing)} parameter points missing from model_cases.csv')
    doc = dict(
        about='Apparent CO turnover frequency of each sample over the parameters of one model. '
              'Written by scripts/export_tof_range.py from pilot/titania-super-multiscale/outputs.',
        fixed_sites_umol_g=RUN.FIXED_SITES,
        site_threshold=dict(min_sites_umol_g=RUN.MIN_SITES, min_fraction_of_bridging=RUN.MIN_FRACTION),
        strong_reduction_C=RUN.STRONG_REDUCTION_C,
        axes=dict(energy_map=list(md.MAPS), cutoff_nm=list(md.CUTOFFS), dG_eV=list(md.DG),
                  f110=list(md.F110), eps=[dict(key=k, value=v) for k, v in md.EPS.items()]),
        reactive=react, pools=list(md.POOLS),
        diameters_nm=[g4(d) for d in md.DIAMETERS], explicit_depth_nm=g4(md.K_EXPL * md.pt.D110),
        capacity_umol_g={f'{f:g}': capacities(f) for f in md.F110},
        samples=[dict(sample=s['sample'], inventory_umol_g=float(s['inventory_umol_g']),
                      rate_co_umol_g_s=float(s['rate_co_umol_g_s']),
                      treatment_T_C=float(s['treatment_T_C']),
                      summary=summary[s['sample']], cases=cases[s['sample']]) for s in samples])
    with open(OUT, 'w') as f:
        json.dump(doc, f, separators=(',', ':'))
        f.write('\n')
    print(f'{os.path.relpath(OUT, ROOT)} written ({os.path.getsize(OUT) // 1024} kB)')


if __name__ == '__main__':
    main()
