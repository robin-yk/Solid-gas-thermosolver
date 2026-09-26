"""Resolve the five displayed scenarios at each sample inventory, D = 900 nm."""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'pilot/titania-super-multiscale'))
from tofrange import model as md
import run

SCENARIOS = [('S1', 'LI_SBR1', .28, 0), ('S2', 'LI_SBR1', .28, .4),
             ('S3', 'HAM', .28, 0), ('S4', 'HAM', .34, 0), ('S5', 'PAB', .34, 0)]

def main():
    result = {s['sample']: [] for s in run.samples()}
    for name, energy, cutoff, dg in SCENARIOS:
        m, layout = md.particle(energy, 900, cutoff=cutoff, dG=dg, f110=.75, eps_r=64)
        for sample in run.samples():
            inventory = float(sample['inventory_umol_g'])
            sol = m.solve(inventory, md.T_EQ)
            rows = []
            for i, tag in enumerate(m.tags):
                if tag == 'Ti':
                    continue
                depth = 450 - m.es['r'][m.shell[i]]
                capacity = m.C[i] * (700 if tag == 'bulk' else 2 if i == layout.idx['cell'] else 1)
                vacancies = float(sol.x[i] @ m.v[i])
                rows.append(dict(depth_nm=float(depth), vacancies=vacancies, capacity=float(capacity)))
            # Combine domains sharing an atomic plane (including the residual BRI sites).
            grouped = {}
            for row in rows:
                key = round(row['depth_nm'], 10)
                q = grouped.setdefault(key, dict(depth_nm=row['depth_nm'], vacancies=0., capacity=0.))
                q['vacancies'] += row['vacancies']; q['capacity'] += row['capacity']
            rows = sorted(grouped.values(), key=lambda r: r['depth_nm'])
            assert abs(sum(r['vacancies'] for r in rows) - inventory) < 1e-5
            for row in rows:
                row['percent'] = 100 * row['vacancies'] / row['capacity']
                assert 0 <= row['percent'] <= 100
            result[sample['sample']].append(dict(id=name, points=rows))
            print(name, sample['sample'], len(rows), flush=True)
    output = dict(diameter_nm=900, temperature_K=md.T_EQ, f110=.75, eps_r=64, samples=result)
    (ROOT / 'web/depth_profiles.json').write_text(json.dumps(output, separators=(',', ':')) + '\n')

if __name__ == '__main__':
    main()
