"""Free energy of an oxygen vacancy that pays for its own electrons.

Removing a neutral oxygen atom from an oxide leaves two electrons behind.
In rutile they do not delocalise: they localise as two Ti(3+) small
polarons on the cation sublattice.  So a vacancy is not one event on one
sublattice, it is one event on the oxygen sublattice and two on the
cation sublattice, and all three carry configurational entropy.

Rutile has two Ti per four O, so N_Ti = N_O / 2 and

    theta_polaron = 2 * (N_O theta) / N_Ti = 4 * theta.

The free energy per oxygen site is then

    f(theta) = E theta
             + kT [ theta ln theta + (1-theta) ln(1-theta) ]
             + kT (1/2) [ p ln p + (1-p) ln(1-p) ],   p = 4 theta

and its derivative, the vacancy chemical potential, is

    mu(theta) = E + kT PHI(theta),
    PHI(theta) = ln(theta/(1-theta)) + 2 ln(4 theta / (1 - 4 theta)).

Two consequences, and they are the reason this module exists.

The cation sublattice runs out at theta = 1/4.  That composition is
TiO(1.5), which is Ti2O3.  Nothing here imposes a cap: PHI diverges
there because there is no fourth Ti(3+) to make.  A free energy written
on the oxygen sublattice alone has no such divergence and will happily
report theta -> 1, a completely deoxygenated lattice.

In the dilute limit the two logarithms add to 3 kT ln theta, so
theta ~ exp(mu / 3kT) and, with mu_V = -mu_O = -(1/2) kT ln p(O2),
theta ~ p(O2)^(-1/6).  That is the textbook exponent for a doubly
ionised vacancy compensated by electrons.  The oxygen-sublattice-only
form gives p(O2)^(-1/2), the exponent for a NEUTRAL vacancy, and an
enrichment between two sites that is the cube of the right one.

The construction is the one Skorodumova et al. established for ceria
(Phys. Rev. Lett. 89, 166601) and that the four-model ceria study this
engine is built on uses for its equilibrium profile; the mapping to
rutile is exact because fluorite CeO2 also has N_cation = N_O / 2, so
the factor of four and the quarter cap are the same numbers.
"""

import math

# theta_polaron / theta_vacancy.  Two polarons per vacancy, divided by
# the cation-to-oxygen site ratio.  Rutile and fluorite both have one
# cation per two oxygens, so both give four.
POLARON_RATIO = 4.0

# Where the cation sublattice is full.  At ratio 4 that is TiO(1.5) = Ti2O3.
THETA_CAP = 1.0 / POLARON_RATIO

# The ratio is exposed because it is not four everywhere.  Per (1x1)
# (110) cell the bridging row holds one oxygen and the topmost titanium
# layer holds two titaniums, so polarons confined to that layer would
# give a ratio of one, not four.  Ratio four is the statement that they
# are not confined - they localise in the first subsurface layer, which
# is what DFT reports for reduced rutile (110) (Reticcioli et al., Phys.
# Rev. X 7, 031053 and Phys. Rev. B 98, 045306).  Passing a different
# ratio is how that statement gets tested rather than assumed.

_ITERS = 80          # both brackets are under two units wide in log space


def cap(ratio=POLARON_RATIO):
    """Site fraction at which the cation sublattice is full."""
    return 1.0 / ratio


def _brackets(ratio):
    """(phi at the half-cap, 2 ln ratio, dilute slack).

    What the two inversion branches need, derived from the ratio rather
    than tabulated, so a non-default ratio stays exact.

    At theta = half the cap the polaron term vanishes identically, so
    phi there is just the oxygen logit - that value is the branch split.
    The dilute slack is the largest the two (1-x) corrections get on the
    lower branch; both increase with theta, so the half-cap endpoint is
    the maximum, and there 1 - ratio*theta is exactly one half whatever
    the ratio is."""
    half = 0.5 / ratio
    phi_half = math.log(half / (1.0 - half))
    log_ratio2 = 2.0 * math.log(ratio)
    slack = -math.log(1.0 - half) - 2.0 * math.log(0.5)
    return phi_half, log_ratio2, slack


def phi(theta, ratio=POLARON_RATIO):
    """Dimensionless part of the chemical potential, mu = E + kT phi.

    Strictly increasing on (0, 1/ratio): the derivative is
    1/(theta(1-theta)) + 2/(theta(1 - ratio theta)), positive throughout.
    Diverges to -inf at zero and to +inf at the cap."""
    if theta <= 0.0:
        return float('-inf')
    if theta >= 1.0 / ratio:
        return float('inf')
    p = ratio * theta
    return (math.log(theta / (1.0 - theta))
            + 2.0 * math.log(p / (1.0 - p)))


def phi_oxygen_only(theta):
    """The same quantity with the cation sublattice left out.

    Kept because the difference is the argument, not because anything
    should call it: this is the form that has no cap and reports a
    fully deoxygenated surface.  Its dilute limit is ln theta, one
    third of phi's, which is the p(O2)^(-1/2) exponent."""
    if theta <= 0.0:
        return float('-inf')
    if theta >= 1.0:
        return float('inf')
    return math.log(theta / (1.0 - theta))


def mu_of_theta(theta, e_form, kt, compensated=True, ratio=POLARON_RATIO):
    """Vacancy chemical potential, eV, at a site of formation energy e_form."""
    f = phi(theta, ratio) if compensated else phi_oxygen_only(theta)
    return e_form + kt * f


def theta_of_phi(x, compensated=True, ratio=POLARON_RATIO):
    """Invert phi.  Returns theta in (0, 1/ratio).

    Both branches bisect in a logarithmic variable so the tiny
    occupancies reached at deep-bulk formation energies keep their
    relative accuracy, and both brackets are derived rather than
    guessed, which is what keeps the iteration count small.

    Dilute branch, theta <= half the cap.  Write

        phi = 3 ln theta + 2 ln ratio - ln(1-theta) - 2 ln(1 - ratio theta).

    The last two terms are non-negative and increase with theta, so over
    the branch they lie between zero and a slack fixed by the half-cap
    endpoint.  Hence

        (x - 2 ln ratio - slack)/3  <=  ln theta  <=  (x - 2 ln ratio)/3,

    a bracket half a unit wide whatever x is.

    Saturated branch.  With w = cap - theta,

        phi = ln(theta/(1-theta)) + 2 ln(1 - ratio w) - 2 ln ratio - 2 ln w,

    which is monotone decreasing in ln w.  The half-cap is an exact
    upper end for ln w and the lower end is reached by doubling."""
    if not compensated:
        if x > 700.0:
            return 1.0
        if x < -700.0:
            return 0.0
        return 1.0 / (1.0 + math.exp(-x))

    phi_half, log_ratio2, slack = _brackets(ratio)
    top = 1.0 / ratio

    if x <= phi_half:
        hi = (x - log_ratio2) / 3.0
        if hi < -745.0:
            return 0.0
        lo = max((x - log_ratio2 - slack) / 3.0, -745.0)
        for _ in range(_ITERS):
            mid = 0.5 * (lo + hi)
            if phi(math.exp(mid), ratio) < x:
                lo = mid
            else:
                hi = mid
        return math.exp(0.5 * (lo + hi))

    # phi decreases in v = ln(cap - theta): small v means theta at the
    # cap and phi at +inf.  v at the half-cap is an exact upper end,
    # because the branch was entered on x > phi(half-cap).  The lower
    # end is found by doubling the gap, which is a handful of steps and
    # costs nothing next to keeping an analytic bound valid for every
    # ratio - at ratio one the cap collides with full occupancy and the
    # closed-form endpoint diverges.
    hi = math.log(0.5 * top)
    lo = hi - 1.0
    while lo > -745.0 and phi(top - math.exp(lo), ratio) < x:
        lo = hi - 2.0 * (hi - lo)
    if lo <= -745.0:
        if phi(top - math.exp(-745.0), ratio) < x:
            return top
        lo = -745.0
    for _ in range(_ITERS):
        mid = 0.5 * (lo + hi)
        t = top - math.exp(mid)
        if t <= 0.0 or phi(t, ratio) > x:
            lo = mid
        else:
            hi = mid
    return top - math.exp(0.5 * (lo + hi))


def theta_of_mu(mu, e_form, kt, compensated=True, ratio=POLARON_RATIO):
    """Site occupancy at a given vacancy chemical potential."""
    return theta_of_phi((mu - e_form) / kt, compensated, ratio)


def dilute_theta(mu, e_form, kt, ratio=POLARON_RATIO):
    """The theta -> 0 asymptote, exp(((mu-E)/kT - 2 ln ratio)/3).

    Exposed so the power law can be checked against the full inversion
    instead of asserted in a comment.  Three in the denominator is the
    whole point: it is what turns p(O2)^(-1/2) into p(O2)^(-1/6)."""
    return math.exp(((mu - e_form) / kt - 2.0 * math.log(ratio)) / 3.0)


def free_energy_per_site(theta, e_form, kt, compensated=True,
                         ratio=POLARON_RATIO):
    """f(theta) in eV per oxygen site.  mu_of_theta is its derivative.

    The cation term carries the site-count ratio 1/ratio * 2, which is
    the number of cation sites per oxygen site: two polarons per vacancy
    spread over a sublattice that is 'ratio' times coarser."""
    if theta <= 0.0:
        return 0.0
    ideal = theta * math.log(theta) + (1.0 - theta) * math.log(1.0 - theta)
    f = e_form * theta + kt * ideal
    if not compensated:
        return f
    p = ratio * theta
    if p >= 1.0:
        return float('inf')
    cat = p * math.log(p) + (1.0 - p) * math.log(1.0 - p)
    return f + kt * cat * (2.0 / ratio)
