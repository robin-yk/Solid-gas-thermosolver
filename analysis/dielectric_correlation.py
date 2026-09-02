"""Dielectric loss against oxygen removal: an empirical regression.

This is a measurement analysis, not thermodynamics, which is why it lives
outside `solidgas/`. It fits

    log10(eps'') = a + b * log10(VO_total)

by ordinary least squares and reports the slope, the intercept, their
standard errors and R-squared. It is not coupled to the equilibrium
calculation or to the population model, and nothing in the package
reads its output.

The fit reported in the manuscript was produced in OriginLab from four
reduced samples; `stored()` returns it verbatim with its provenance so
the number in the paper and the number here cannot drift apart. `fit()`
reproduces the same regression from a CSV of (VO_umol_g, eps_imag) rows
once those four points are digitised.
"""

import csv
import json
import math
import os

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
STORED_PATH = os.path.join(_ROOT, 'data', 'dielectric_fit.json')


def stored(path=STORED_PATH):
    """The manuscript's own fit, as published."""
    with open(path) as fh:
        return json.load(fh)


def fit(rows):
    """Least squares on log10-log10. `rows` is (VO_umol_g, eps_imag) pairs.

    Returns the same fields `stored()` carries, so the two can be compared
    directly. Standard errors are the usual OLS ones; with four points the
    residual has two degrees of freedom and the errors are large, which is
    the honest reading of the published fit as well.
    """
    pts = [(float(x), float(y)) for x, y in rows if float(x) > 0 and float(y) > 0]
    n = len(pts)
    if n < 3:
        raise ValueError('three points is the minimum that leaves a residual')
    xs = [math.log10(x) for x, _ in pts]
    ys = [math.log10(y) for _, y in pts]
    mx, my = sum(xs) / n, sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    syy = sum((y - my) ** 2 for y in ys)
    if sxx <= 0:
        raise ValueError('every point sits at the same vacancy content')
    b = sxy / sxx
    a = my - b * mx
    rss = sum((y - (a + b * x)) ** 2 for x, y in zip(xs, ys))
    dof = n - 2
    s2 = rss / dof if dof > 0 else float('nan')
    return {
        'form': 'power_law',
        'exponent': b,
        'exponent_err': math.sqrt(s2 / sxx) if dof > 0 else float('nan'),
        'log10_intercept': a,
        'log10_intercept_err': (math.sqrt(s2 * (1.0 / n + mx * mx / sxx))
                                if dof > 0 else float('nan')),
        'fit_stats': {
            'pearson_r': sxy / math.sqrt(sxx * syy) if syy > 0 else float('nan'),
            'r_squared': 1.0 - rss / syy if syy > 0 else float('nan'),
            'residual_ss': rss,
            'n_points': n,
        },
    }


def fit_csv(path):
    """Same, reading a two-column CSV with a header row."""
    with open(path) as fh:
        rdr = csv.DictReader(fh)
        rows = [(r['VO_umol_g'], r['eps_imag']) for r in rdr]
    return fit(rows)


def predict(model, vo_umol_g):
    """eps'' at a given oxygen removal, from a fit or from `stored()`."""
    return 10.0 ** (model['log10_intercept']
                    + model['exponent'] * math.log10(vo_umol_g))
