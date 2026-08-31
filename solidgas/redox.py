"""Steady-state Mars-van Krevelen redox: two half-reactions, one crossing.

The mechanism is Mars-van Krevelen and this is its site balance at steady
state: the reductant takes lattice oxygen and leaves a vacancy, the oxidant
fills a vacancy and puts lattice oxygen back, nothing adsorbs, and there are
no surface intermediates.

theta is the vacancy fraction on the TWO-COORDINATE BRIDGING OXYGEN ROWS of
rutile (110), as a fraction of those sites. That is where the rest of the
model is already referenced: the descriptor is E_vac for (110) bridging O
(4.39 eV, sX, Li/Guo/Robertson 2015) and the face ceiling theta_max = 0.5 is
the (1x2) added-row reconstruction, a bridging-row phase. Three-coordinate
in-plane oxygen is not counted - it is held more strongly, its removal energy
is not the descriptor's, and one theta over both would describe neither site.
The particle-average theta is a different quantity: stoichiometric
bookkeeping, x/2 in MO(2-x), with no claim about which defect carries the
reduction (in rutile at low x that is the titanium interstitial, and past the
solubility the solid shears rather than making more point defects).

The reductant and the oxidant are deliberately unnamed. The rate parameters
in data/redox.json are identical placeholders for both halves, so naming a
gas pair would assert chemistry the numbers do not carry; what the model does
require of the pair is unit order in gas and in site, which is what lets
temperature, pressure and prefactor cancel out of the crossing.

A reductant removes lattice oxygen at a rate that needs oxygen to remove, and
an oxidant refills vacancies at a rate that needs vacancies to fill:

    r_red(theta) = A_r p_r^m_r (1 - theta)^n_r exp(-E_r/kT)     decreasing
    r_ox (theta) = A_o p_o^m_o      theta ^n_o exp(-E_o/kT)     increasing

One is falling and the other rising over theta in [0, 1], so they cross
exactly once and the crossing is the steady state. That is the whole model:
one nonlinear equation in one variable, solved here by bisection and, at unit
site orders, in closed form as theta* = k_r / (k_r + k_o). Nothing in this
file is new - it is 1970s reaction engineering - and the closed form is kept
public precisely so the shipped numeric path can be checked against it.

What makes the pair a volcano rather than a ranking is that both barriers
move with the SAME material descriptor, the oxygen vacancy formation energy,
and they move in opposite directions: an oxide that gives up oxygen easily is
quick to reduce and slow to reoxidise. With linear (Bronsted-Evans-Polanyi)
coupling the steady rate is the harmonic mean of the two rate constants, and
the harmonic mean of two opposed exponentials peaks. The peak sits at
theta* = alpha / (alpha + beta), independent of temperature, pressure and
prefactors; only its position along the descriptor axis moves with those.

The uniform assumption, and what relaxing it actually does. theta is read as
one number for the whole particle, so the reacting surface is assumed to
carry the particle-average vacancy fraction. Supplying a layer-resolved
mapping average -> surface changes the answer in two ways, and only one of
them is the one usually expected:

  The steady RATE does not move at all. Any strictly increasing mapping is a
  reparameterisation of the same crossing in surface coverage, and the
  crossing coverage is fixed by the rate constants alone, so the turnover
  is identical to the last digit. The uniform model gets the rate right.

  The steady STATE moves by orders of magnitude. What changes is the
  particle inventory that holds the surface at that coverage. For rutile at
  600 C the surface runs ~1844x above the particle average in the dilute
  limit - close to its geometric ceiling of N_O_total / N_s = 1847 - so the
  carrier sits at roughly a two-thousandth of the reduction extent the
  uniform model assigns it. Degree of reduction is the measured quantity, so
  this is the half that matters for comparison with experiment.

  And on half the descriptor axis the steady state stops existing. The
  mapping saturates: a reconstructing surface is pinned at a fixed
  occupancy, so above that inventory both half-rates stop depending on how
  reduced the particle is and the imbalance keeps its sign. Every oxide
  easier to reduce than the balance point then has no crossing - it reduces
  until some other process stops it. The correction removes the root rather
  than moving it, which is not something a perturbative reading of the
  uniform model anticipates.

surface_correction() measures all of this rather than asserting it, by
taking the mapping from whatever layer-resolved model is supplied. The
uniform case is the reference, not the answer.

Domain. theta is a point-defect concentration only while the reduced solid
stays a single phase. Past the declared solubility the vacancies order into
shear planes and mass action in theta stops meaning anything; steady_state()
reports that as a flag and a warning, never by clipping.
"""

import json
import math
import os

KB_EV = 8.617333262e-5          # eV / K

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_PATH = os.path.join(_ROOT, 'data', 'redox.json')


def load_params(path=DATA_PATH):
    with open(path) as fh:
        return json.load(fh)


# ------------------------------------------------------------ rate constants

def _arrhenius(log10_a, e_ev, kt):
    """A exp(-E/kT) with the exponent guarded, not the barrier.

    A negative E is left alone here: linear free-energy relations do
    extrapolate below zero, and silently flooring the barrier would hide
    that the descriptor has left the range the relation was drawn for.
    rate_constants() raises the flag instead."""
    x = e_ev / kt
    if x > 700.0:
        return 0.0
    if x < -700.0:
        x = -700.0
    return (10.0 ** log10_a) * math.exp(-x)


def rate_constants(p, t_c=None, dev=0.0, p_red=None, p_ox=None):
    """Both half-reaction rate constants at one point on the descriptor axis.

    dev is the offset of the oxygen vacancy formation energy from the
    reference material, in eV; negative means easier to reduce. The two
    barriers move with it in opposite directions, which is the only place
    the material enters."""
    cond = p['conditions']
    if t_c is None:
        t_c = cond['T_C']
    if p_red is None:
        p_red = cond['p_red_atm']
    if p_ox is None:
        p_ox = cond['p_ox_atm']
    kt = KB_EV * (t_c + 273.15)
    red = p['half_reactions']['reduction']
    ox = p['half_reactions']['oxidation']

    e_red = red['E0_eV'] + red['bep_slope'] * dev
    e_ox = ox['E0_eV'] - ox['bep_slope'] * dev
    k_red = _arrhenius(red['log10_A_s1'], e_red, kt) * p_red ** red['order_gas']
    k_ox = _arrhenius(ox['log10_A_s1'], e_ox, kt) * p_ox ** ox['order_gas']

    warnings = []
    if e_red < 0.0 or e_ox < 0.0:
        warnings.append({
            'kind': 'bep_out_of_range',
            'text': ('The linear free-energy relation extrapolates to a '
                     'negative barrier at this descriptor offset '
                     '(E_red = %.2f eV, E_ox = %.2f eV). The rate constant '
                     'is still evaluated, but the relation is being read '
                     'outside the range it can support.' % (e_red, e_ox))})
    return {'k_red_s1': k_red, 'k_ox_s1': k_ox,
            'E_red_eV': e_red, 'E_ox_eV': e_ox,
            'kT_eV': kt, 'T_C': t_c, 'dev_eV': dev,
            'p_red_atm': p_red, 'p_ox_atm': p_ox,
            'warnings': warnings}


def theta_first_order(k_red, k_ox):
    """Closed-form crossing at unit site orders: theta* = k_r / (k_r + k_o).

    Public because it is the oracle the bisection is checked against, not
    because the solver takes this shortcut - steady_state() always bisects."""
    tot = k_red + k_ox
    if tot <= 0.0:
        return 0.0
    return k_red / tot


# ------------------------------------------------------------- steady state

def half_rates(p, theta, k):
    """(reduction, oxidation) rate per reacting oxygen site at theta."""
    red = p['half_reactions']['reduction']
    ox = p['half_reactions']['oxidation']
    avail = max(0.0, 1.0 - theta)
    vac = max(0.0, theta)
    return (k['k_red_s1'] * avail ** red['order_site'],
            k['k_ox_s1'] * vac ** ox['order_site'])


def _sigmoid_pair(z):
    """(theta, 1 - theta) from the logit, both to full relative precision.

    Bisecting in theta cannot resolve 1 - theta once theta nears 1: the
    subtraction throws away every digit the root has, and the rate residual
    degrades to 1e-5 at the ends of the descriptor axis. Bisecting in the
    logit and returning the pair keeps both tails exact."""
    if z >= 0.0:
        e = math.exp(-z) if z < 700.0 else 0.0
        return 1.0 / (1.0 + e), e / (1.0 + e)
    e = math.exp(z) if z > -700.0 else 0.0
    return e / (1.0 + e), 1.0 / (1.0 + e)


def _mapping_saturates(surf, n):
    """Whether the average -> surface mapping stops responding.

    A reconstructing surface is pinned at a fixed occupancy by its own
    phase, so past a modest inventory it no longer reports how reduced the
    particle is. Both half-rates then go flat in the inventory and the
    balance has nothing left to solve for. Probed over the upper half of the
    axis, where any real mapping has already saturated if it is going to."""
    if n < 2:
        return False
    vals = [surf(0.5 + 0.5 * i / float(n - 1)) for i in range(n)]
    scale = max(abs(v) for v in vals)
    if scale == 0.0:
        return False
    return max(abs(v - vals[-1]) for v in vals) <= 1e-12 * scale


def steady_state(p, t_c=None, dev=0.0, p_red=None, p_ox=None,
                 theta_sol=None, theta_of_average=None, window_n=17,
                 k=None):
    """The crossing of the two half-reaction rates.

    theta_of_average, when given, is the layer-resolved correction: a
    callable mapping the particle-average vacancy fraction to the vacancy
    fraction on the reacting surface. The balance is then solved in the
    average while both rates are evaluated at the surface, which is the
    'one more pass' over the uniform model. Left out, the mapping is the
    identity and the result is the uniform model exactly."""
    if k is None:
        k = rate_constants(p, t_c=t_c, dev=dev, p_red=p_red, p_ox=p_ox)
    warnings = list(k['warnings'])
    st = p['structure']
    if theta_sol is None:
        theta_sol = st['x_solubility'] / st['oxygen_per_formula']
    identity = theta_of_average is None
    surf = theta_of_average if not identity else (lambda t: t)
    red_n = p['half_reactions']['reduction']['order_site']
    ox_n = p['half_reactions']['oxidation']['order_site']

    def imbalance(theta):
        r_red, r_ox = half_rates(p, surf(theta), k)
        return r_red - r_ox

    def imbalance_z(z):
        """The same function of the logit, without the losing subtraction."""
        th, one_minus = _sigmoid_pair(z)
        if not identity:
            th = surf(th)
            one_minus = 1.0 - th
        return (k['k_red_s1'] * max(0.0, one_minus) ** red_n
                - k['k_ox_s1'] * max(0.0, th) ** ox_n)

    # r_red falls and r_ox rises over the interval, so the sign at the ends
    # settles whether the crossing is interior before any bisection.
    lo, hi = -60.0, 60.0
    for _ in range(6):
        if imbalance_z(lo) >= 0.0 >= imbalance_z(hi):
            break
        lo, hi = 2.0 * lo, 2.0 * hi
    if imbalance_z(hi) > 0.0:
        theta, bracketed = 1.0, False
    elif imbalance_z(lo) < 0.0:
        theta, bracketed = 0.0, False
    else:
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            if imbalance_z(mid) > 0.0:
                lo = mid
            else:
                hi = mid
        theta, bracketed = _sigmoid_pair(0.5 * (lo + hi))[0], True

    runaway = None
    if not bracketed:
        runaway = 'reduction' if theta == 1.0 else 'oxidation'
        warnings.append({
            'kind': 'runaway_' + runaway,
            'text': ('The two rates do not cross inside 0 < theta < 1 at '
                     'this descriptor and these pressures: the imbalance '
                     'keeps its sign and the solid runs to theta = %.0f. '
                     'There is no steady state here, which is an answer '
                     'about the carrier and not a solver failure - it '
                     'reduces (or reoxidises) until something other than '
                     'this pair of rate laws stops it.' % theta)})

    # Whether the root LOCATES the state is a local question: the imbalance
    # has to actually move through it. Flatness elsewhere on the axis says
    # nothing about the root, so this is probed one logit unit either side
    # rather than over the whole interval.
    localised = bracketed
    if bracketed:
        zc = 0.5 * (lo + hi)
        f_lo, f_hi = imbalance_z(zc - 1.0), imbalance_z(zc + 1.0)
        span = max(abs(f_lo), abs(f_hi))
        localised = span > 0.0 and min(abs(f_lo), abs(f_hi)) > 1e-9 * span
        if not localised:
            warnings.append({
                'kind': 'root_not_localised',
                'text': ('The rate imbalance does not move through the '
                         'crossing: it reaches zero and stays there, so '
                         'every inventory beyond this point is equally a '
                         'steady state and the root reported below is only '
                         'the near edge of that set. Mass action in a '
                         'single average theta cannot pick one out.')})

    buffered = (not identity) and _mapping_saturates(surf, window_n)
    if buffered:
        warnings.append({
            'kind': 'surface_buffered',
            'text': ('The reacting surface stops responding to the '
                     'inventory over the upper part of the axis: it is held '
                     'at a fixed occupancy by a surface phase. Above that '
                     'point neither half-rate knows how reduced the '
                     'particle is, which is why a carrier on the reducible '
                     'side of the balance point has no steady state to '
                     'settle into rather than merely a different one.')})

    theta_surface = surf(theta)
    if identity and bracketed:
        # read the pair off the converged logit rather than subtracting
        th, one_minus = _sigmoid_pair(0.5 * (lo + hi))
        r_red = k['k_red_s1'] * one_minus ** red_n
        r_ox = k['k_ox_s1'] * th ** ox_n
    else:
        r_red, r_ox = half_rates(p, theta_surface, k)
    rate = 0.5 * (r_red + r_ox)
    scale = max(r_red, r_ox)
    residual = 0.0 if scale == 0.0 else abs(r_red - r_ox) / scale

    band = p['regime_thresholds']
    if theta < band['theta_lo']:
        regime = 'reduction_limited'
    elif theta > band['theta_hi']:
        regime = 'oxidation_limited'
    else:
        regime = 'balanced'

    x_equiv = theta * st['oxygen_per_formula']
    beyond = theta > theta_sol
    if beyond:
        warnings.append({
            'kind': 'beyond_solubility',
            'text': ('The steady state sits at theta = %.4f (x = %.4f in '
                     'MO(2-x)), above the declared point-defect solubility '
                     'theta = %.4f. Mass action in theta does not describe '
                     'this solid: the vacancies order into a second phase '
                     'before the carrier reaches this crossing, so the '
                     'number below is where an unshearing solid WOULD sit, '
                     'and the real carrier changes structure instead.'
                     % (theta, x_equiv, theta_sol))})

    return {'theta': theta, 'theta_surface': theta_surface,
            'rate_s1': rate, 'r_red_s1': r_red, 'r_ox_s1': r_ox,
            'residual': residual, 'interior_crossing': bracketed,
            'runaway': runaway, 'surface_buffered': buffered,
            'localised': localised,
            'regime': regime,
            'theta_solubility': theta_sol, 'x_equiv': x_equiv,
            'beyond_solubility': beyond,
            'k_red_s1': k['k_red_s1'], 'k_ox_s1': k['k_ox_s1'],
            'E_red_eV': k['E_red_eV'], 'E_ox_eV': k['E_ox_eV'],
            'dev_eV': dev, 'T_C': k['T_C'],
            'p_red_atm': k['p_red_atm'], 'p_ox_atm': k['p_ox_atm'],
            'uniform': identity,
            'warnings': warnings}


# ------------------------------------------------- the assumption-free axis

def from_ratio(p, ratio, theta_sol=None, theta_of_average=None,
               window_n=17):
    """The crossing from the ratio of the two rate constants and nothing else.

    K = k_red / k_ox is the whole material dependence of the steady state:
    at unit site orders theta* = K/(1+K) exactly, with no temperature, no
    pressure, no prefactor and no linear free-energy relation anywhere in
    it. Everything else this module offers - the descriptor axis, the
    Bronsted-Evans-Polanyi slopes, the volcano - is a way of predicting K
    from a material property, and each of those steps adds a parameter
    that has to be assumed. This entry point adds none.

    Rates come back in units of k_ox, since only the ratio was given."""
    if not (ratio >= 0.0):
        raise ValueError('K = k_red/k_ox must be non-negative')
    k = {'k_red_s1': float(ratio), 'k_ox_s1': 1.0,
         'E_red_eV': None, 'E_ox_eV': None, 'kT_eV': None, 'T_C': None,
         'dev_eV': None, 'p_red_atm': None, 'p_ox_atm': None,
         'warnings': []}
    out = steady_state(p, theta_sol=theta_sol, k=k,
                       theta_of_average=theta_of_average,
                       window_n=window_n)
    out['K'] = float(ratio)
    out['rate_unit'] = 'k_ox'
    return out


def theta_of_ratio(ratio):
    """theta* = K/(1+K), the closed form the solver is checked against."""
    return ratio / (1.0 + ratio) if ratio >= 0.0 else 0.0


def ratio_sweep(p, k_lo=1e-4, k_hi=1e4, n=161, theta_sol=None,
                theta_of_average=None, window_n=0):
    """One steady state per decade-spaced K: the family, with no assumptions.

    Logarithmic because K is a ratio, so equal factors deserve equal space."""
    rows = []
    lg0, lg1 = math.log10(k_lo), math.log10(k_hi)
    for i in range(n):
        kk = 10.0 ** (lg0 + (lg1 - lg0) * i / float(n - 1))
        s = from_ratio(p, kk, theta_sol=theta_sol,
                       theta_of_average=theta_of_average,
                       window_n=window_n)
        rows.append({'K': kk, 'theta': s['theta'],
                     'theta_surface': s['theta_surface'],
                     'rate_over_k_ox': s['rate_s1'], 'regime': s['regime'],
                     'beyond_solubility': s['beyond_solubility'],
                     'runaway': s['runaway']})
    return rows


def phase_map(p, k_lo=1e-4, k_hi=1e4, e_lo=1.0, e_hi=1e4, nk=61, ne=41,
              theta_sol=None, theta_surface_max=None):
    """Regions of behaviour over (K, surface enrichment).

    The second axis is E = theta_surface / theta_average, which is what a
    layer-resolved partition supplies and what the uniform assumption sets
    to 1. Reading the diagram along E answers the question the uniform
    model cannot ask: how much of the design window was lost by assuming
    the particle is uniformly reduced.

    Two boundaries, both closed forms rather than swept:

      no steady state   theta_s* = K/(1+K) > theta_surface_max
                        - the face is pinned by its own phase, so neither
                          half-rate depends on the inventory any more.
                          A vertical line: it does not involve E.

      second phase      theta_bulk = theta_s*/E > theta_sol
                        - the average composition leaves the point-defect
                          range and the solid orders its vacancies away.
                          A curve: raising E is exactly what buys room
                          here."""
    st = p['structure']
    if theta_sol is None:
        theta_sol = st['x_solubility'] / st['oxygen_per_formula']
    if theta_surface_max is None:
        theta_surface_max = st.get('theta_surface_max', 0.5)
    band = p['regime_thresholds']
    lgk0, lgk1 = math.log10(k_lo), math.log10(k_hi)
    lge0, lge1 = math.log10(e_lo), math.log10(e_hi)
    cells = []
    for j in range(ne):
        e = 10.0 ** (lge0 + (lge1 - lge0) * j / float(ne - 1))
        row = []
        for i in range(nk):
            kk = 10.0 ** (lgk0 + (lgk1 - lgk0) * i / float(nk - 1))
            th = theta_of_ratio(kk)
            if th > theta_surface_max:
                row.append('no_steady_state')
            elif th / e > theta_sol:
                row.append('second_phase')
            elif th < band['theta_lo']:
                row.append('reduction_limited')
            elif th > band['theta_hi']:
                row.append('oxidation_limited')
            else:
                row.append('balanced')
        cells.append({'E': e, 'cells': row})
    k_runaway = (theta_surface_max / (1.0 - theta_surface_max)
                 if theta_surface_max < 1.0 else float('inf'))
    return {'K_lo': k_lo, 'K_hi': k_hi, 'E_lo': e_lo, 'E_hi': e_hi,
            'nk': nk, 'ne': ne, 'theta_sol': theta_sol,
            'theta_surface_max': theta_surface_max,
            'K_no_steady_state': k_runaway,
            'rows': cells}


def usable_window(p, enrichment, theta_sol=None, theta_surface_max=None):
    """The largest K a carrier can run at, at a given surface enrichment.

    Both limits in closed form: the solid must stay one phase, and the face
    must stay below the coverage its own phase pins it at. The smaller of
    the two is the answer, and which one binds is the interesting part."""
    st = p['structure']
    if theta_sol is None:
        theta_sol = st['x_solubility'] / st['oxygen_per_formula']
    if theta_surface_max is None:
        theta_surface_max = st.get('theta_surface_max', 0.5)
    th_phase = min(1.0, enrichment * theta_sol)
    k_phase = (th_phase / (1.0 - th_phase) if th_phase < 1.0
               else float('inf'))
    k_face = (theta_surface_max / (1.0 - theta_surface_max)
              if theta_surface_max < 1.0 else float('inf'))
    return {'enrichment': enrichment,
            'K_max': min(k_phase, k_face),
            'K_second_phase': k_phase, 'K_no_steady_state': k_face,
            'binding': 'second_phase' if k_phase < k_face
                       else 'no_steady_state'}


# ------------------------------------------------------- graphical solution

def crossing_curves(p, t_c=None, dev=0.0, p_red=None, p_ox=None, n=201):
    """The two rate curves over theta, for the plot the crossing is read off.

    This is the figure and the solver looking at the same functions: the
    root returned by steady_state() is where these two columns meet."""
    k = rate_constants(p, t_c=t_c, dev=dev, p_red=p_red, p_ox=p_ox)
    rows = []
    for i in range(n):
        theta = i / float(n - 1)
        r_red, r_ox = half_rates(p, theta, k)
        rows.append({'theta': theta, 'r_red_s1': r_red, 'r_ox_s1': r_ox})
    return rows


# ------------------------------------------------------------ family of runs

def sabatier_optimum(p, t_c=None, p_red=None, p_ox=None):
    """Closed-form peak of the volcano, valid at unit site orders.

    Maximising the harmonic mean k_r k_o / (k_r + k_o) over the descriptor
    gives beta k_r = alpha k_o, hence

        theta*  = alpha / (alpha + beta)
        dev*    = [ (E_ox0 - E_red0) - kT ln(alpha/beta)
                    + kT ln(C_r / C_o) ] / (alpha + beta)

    with C the prefactor times its pressure term. theta* at the peak carries
    no temperature, no pressure and no prefactor: only the two slopes set it.
    At other site orders the peak has to be swept for, and this returns None."""
    red = p['half_reactions']['reduction']
    ox = p['half_reactions']['oxidation']
    if red['order_site'] != 1.0 or ox['order_site'] != 1.0:
        return None
    cond = p['conditions']
    if t_c is None:
        t_c = cond['T_C']
    if p_red is None:
        p_red = cond['p_red_atm']
    if p_ox is None:
        p_ox = cond['p_ox_atm']
    kt = KB_EV * (t_c + 273.15)
    a, b = red['bep_slope'], ox['bep_slope']
    if a + b <= 0.0:
        return None
    c_r = (10.0 ** red['log10_A_s1']) * p_red ** red['order_gas']
    c_o = (10.0 ** ox['log10_A_s1']) * p_ox ** ox['order_gas']
    dev = ((ox['E0_eV'] - red['E0_eV']) - kt * math.log(a / b)
           + kt * math.log(c_r / c_o)) / (a + b)
    return {'dev_eV': dev, 'theta': a / (a + b), 'T_C': t_c}


def descriptor_sweep(p, t_c=None, p_red=None, p_ox=None, n=121,
                     dev_lo=None, dev_hi=None, theta_sol=None,
                     theta_of_average=None, window_n=17):
    """The family of solutions: one steady state per descriptor offset.

    This is the deliverable before any material is fitted. Each row is a
    hypothetical oxide placed by one number, and the rate column is the
    volcano."""
    rng = p['descriptor']['offset_range_eV']
    if dev_lo is None:
        dev_lo = rng[0]
    if dev_hi is None:
        dev_hi = rng[1]
    rows = []
    for i in range(n):
        dev = dev_lo + (dev_hi - dev_lo) * i / float(n - 1)
        s = steady_state(p, t_c=t_c, dev=dev, p_red=p_red, p_ox=p_ox,
                         theta_sol=theta_sol, window_n=window_n,
                         theta_of_average=theta_of_average)
        rows.append({'dev_eV': dev, 'theta': s['theta'],
                     'rate_s1': s['rate_s1'], 'regime': s['regime'],
                     'beyond_solubility': s['beyond_solubility'],
                     'runaway': s['runaway'], 'localised': s['localised'],
                     'k_red_s1': s['k_red_s1'], 'k_ox_s1': s['k_ox_s1']})
    return rows


# ------------------------------------------------- the one-more-pass measure

def surface_correction(p, theta_of_average, t_c=None, dev=0.0,
                       p_red=None, p_ox=None, theta_sol=None):
    """How far the uniform answer is from the layer-resolved one.

    Runs the same crossing twice - once with the reacting surface held at
    the particle average, once with it taken from the supplied layer-
    resolved mapping - and returns both with the ratios. This is the
    quantity that says whether the second pass is a refinement or the
    whole answer."""
    uni = steady_state(p, t_c=t_c, dev=dev, p_red=p_red, p_ox=p_ox,
                       theta_sol=theta_sol)
    cor = steady_state(p, t_c=t_c, dev=dev, p_red=p_red, p_ox=p_ox,
                       theta_sol=theta_sol, theta_of_average=theta_of_average)
    def ratio(a, b):
        return float('inf') if b == 0.0 else a / b
    return {'uniform': uni, 'corrected': cor,
            'rate_ratio': ratio(cor['rate_s1'], uni['rate_s1']),
            'theta_ratio': ratio(cor['theta'], uni['theta']),
            'surface_ratio': ratio(cor['theta_surface'], cor['theta'])}
