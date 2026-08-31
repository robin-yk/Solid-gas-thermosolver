"""Print what the particle engine says, at one condition or a ladder.

    python scripts/run_particle.py                    the shipped condition
    python scripts/run_particle.py --vo 40 --T 700
    python scripts/run_particle.py --ladder
"""

import argparse
import sys

sys.path.insert(0, __file__.rsplit('/', 2)[0])

from solidgas import particle as P
from solidgas.particle import spacecharge as SC
from solidgas.particle import transport as T


def one(p, part, vo, t_c, exposure_s):
    r = P.equilibrate(vo, t_c, part)
    seg = part['segregation']
    print('rutile TiO2, %s energies, %.2f um particle, %.0f C, %.4g umol-O/g'
          % (seg['label'], part['diameter_um'], t_c, vo))
    print('  x in TiO(2-x)              %.6f   (shear limit %.4f)'
          % (2 * r['theta_particle'], part['x_max_shear']))
    print('  depth profile              E(z) = %.3f - %.3f exp(-z/%.4f nm)'
          % (seg['e_bulk_eV'], seg['dE_seg_eV'], seg['xi_nm']))
    print()
    print('  chemical potential         %+.6f eV' % r['mu_eV'])
    print('  bridging row               theta = %.4f   (%s)'
          % (r['surface']['theta'], r['surface']['phase']))
    print('  first subsurface layer     theta = %.6f' % r['theta_first_layer'])
    print('  bulk plateau               theta = %.6e' % r['theta_bulk_plateau'])
    print('  enrichment, layer 1 / bulk %.1f x'
          % (r['theta_first_layer'] / r['theta_bulk_plateau']))
    print()
    print('  on the row                 %9.3f umol-O/g' % r['umol_surface'])
    print('  below it                   %9.3f umol-O/g' % r['umol_interior'])
    print('  extended defects           %9.3f umol-O/g' % r['umol_extended'])
    print('  pinned by                  %s' % (r['limited_by'] or 'nothing'))
    print('  point-defect ceiling       %9.3f umol-O/g  (this loading is %.0f%%)'
          % (r['point_defect_ceiling_umol_g'],
             100 * vo / r['point_defect_ceiling_umol_g']))
    print()
    print('  where the oxygen is:')
    for row in P.shell_split(r, part):
        print('    within %5.1f nm          %6.3f%% of the volume, '
              '%6.2f%% of the oxygen  (%4.0f x)'
              % (row['depth_nm'], 100 * row['volume_fraction'],
                 100 * row['inventory_fraction'],
                 row['inventory_fraction'] / row['volume_fraction']))
    print()
    v = T.verdict(t_c, part, exposure_s=exposure_s, p=p)
    print('  is it reachable:')
    print('    D                        %.3e cm2/s' % v['D_cm2_per_s'])
    print('    99%% homogenised in       %.3e s   against %g s of exposure'
          % (v['t_homogenise_s'], exposure_s))
    print('    barrier that would break %.3f eV   (literature %.2f, margin %.3f)'
          % (v['E_m_that_would_break_it_eV'], v['E_m_eV'],
             v['barrier_margin_eV']))
    print()
    a = SC.assess(r, part, p)
    c = SC.ceria_comparison(r, part, p)
    print('  is local electroneutrality allowed:')
    print('    Debye length             %.3f nm   against xi = %.3f nm'
          % (a['debye_nm_at_bulk'], a['xi_nm']))
    print('    %.1fx more dilute than the ceria case it is taken from'
          % c['dilution_factor'])
    print('    %s' % a['verdict'])


def ladder(part):
    print('%10s %10s %8s %10s %10s %10s  %s'
          % ('umol-O/g', 'mu (eV)', 'th_surf', 'th_layer1', 'th_bulk',
             'extended', 'pinned by'))
    for vo in (1, 5, 20, 30, 33.5, 35, 37.5, 50, 95, 112, 150, 400, 2000):
        r = P.equilibrate(float(vo), 600.0, part)
        print('%10.4g %+10.6f %8.4f %10.3e %10.3e %10.2f  %s'
              % (vo, r['mu_eV'], r['surface']['theta'],
                 r['theta_first_layer'], r['theta_bulk_plateau'],
                 r['umol_extended'], r['limited_by'] or '-'))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--vo', type=float, default=None)
    ap.add_argument('--T', type=float, default=None)
    ap.add_argument('--d-um', type=float, default=None)
    ap.add_argument('--preset', default=None)
    ap.add_argument('--exposure', type=float, default=300.0)
    ap.add_argument('--ladder', action='store_true')
    a = ap.parse_args()

    p = P.load()
    part = P.particle(p, d_um=a.d_um, preset=a.preset)
    if a.ladder:
        return ladder(part)
    one(p, part,
        a.vo if a.vo is not None else p['defaults']['VO_total_umol_g'],
        a.T if a.T is not None else float(p['defaults']['T_C']),
        a.exposure)


if __name__ == '__main__':
    main()
