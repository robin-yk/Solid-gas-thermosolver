"""Full-composition Gibbs minimisation in reaction-extent coordinates.

Same objective and same feasible set as the species-space minimisation in
`equilibrium.py` - the change of variables is exact - but the constant part of
G is removed analytically instead of being carried through the arithmetic.

That is not a refinement, it is what makes the inert case computable at all.
Sealed N2 pulls rutile off stoichiometry by a mole fraction of order 1e-38, and
in species space the optimiser never sees it: n(TiO2) = n_Ti - 10*xi with n_Ti
at 1e-3 rounds to n_Ti exactly, because a double carries sixteen digits and the
gap is thirty-nine. Minimising the reaction Gibbs energy directly keeps xi as a
primary variable at its own magnitude, so the same equilibrium condition is
solved without ever forming that difference.

Components are chosen as the unreacted charge: the feed gases plus rutile. Each
remaining species is formed from them by one reaction and carries one extent.
G is reported relative to that unreacted state.
"""

import numpy as np
from scipy.optimize import minimize as _minimize

from .shomate import mu0, R_KJ, R_ATM  # noqa: F401
from .equilibrium import FORMULAS, METAL_PHASES, solid_mu

LOG_FLOOR = -110.0          # exp(-110) = 1.7e-48, below any physical extent
LOG_CEIL = np.log(10.0)     # basis is scaled so no species exceeds ~1


class GibbsSystem:
    """One atmosphere: its species, its elements, and its reference charge.

    `components` must be one species per element and non-singular; the natural
    choice is the charge as loaded, which makes the reference state the
    unreacted feed and every extent a reaction actually running forward.
    """

    def __init__(self, gases, solids, components, source='waldner'):
        self.gases, self.solids = list(gases), list(solids)
        self.species = self.gases + self.solids
        self.source = source

        elements = []
        for s in self.species:
            for e in FORMULAS[s]:
                if e not in elements:
                    elements.append(e)
        self.elements = sorted(elements)

        self.E = np.array([[FORMULAS[s].get(e, 0) for e in self.elements]
                           for s in self.species], dtype=float)
        gasset = set(self.gases)
        self.is_gas = np.array([s in gasset for s in self.species])

        self.comp = [self.species.index(c) for c in components]
        self.free = [i for i in range(len(self.species)) if i not in self.comp]
        Ek = self.E[self.comp].T
        if Ek.shape[0] != Ek.shape[1]:
            raise ValueError(f'need one component per element, got '
                             f'{len(components)} for {self.elements}')
        if abs(np.linalg.det(Ek)) < 1e-12:
            raise ValueError('component set is singular')
        self.Ainv = np.linalg.inv(Ek)
        # dn(component) per unit extent: what forming a free species costs.
        self.dK = -self.Ainv @ self.E[self.free].T

    # -- composition -------------------------------------------------------

    def b_from(self, **moles):
        """Element vector for a charge written as species moles.

        The element order is whatever sorting produced, and it is not the same
        for Ti-O as for Ce-O - sorted puts Ce before H. Building the vector by
        name instead of by position is the only way that cannot go wrong
        silently, which is exactly how it would go wrong.
        """
        b = np.zeros(len(self.elements))
        for name, amount in moles.items():
            for e, k in FORMULAS[name].items():
                b[self.elements.index(e)] += k * amount
        return b

    def n0(self, b):
        """The reference charge: components only, everything else absent.

        A charge of pure CH4 over rutile has to be expressed in four components
        for four elements, so two of them come out at zero; inverting the
        element matrix leaves those at round-off rather than at nothing, and
        the sign of that round-off would decide feasibility. Snap it.
        """
        b = np.asarray(b, dtype=float)
        v = self.Ainv @ b
        v[np.abs(v) < 1e-14 * max(np.abs(b).max(), 1e-300)] = 0.0
        n = np.zeros(len(self.species))
        n[self.comp] = v
        return n

    def delta_n(self, xi):
        d = np.zeros(len(self.species))
        d[self.free] = xi
        d[self.comp] = self.dK @ xi
        return d

    def n_of(self, xi, b):
        return self.n0(b) + self.delta_n(xi)

    def stoichiometry(self, name):
        """How free species `name` is built out of the components."""
        j = self.free.index(self.species.index(name))
        return {self.species[k]: -self.dK[i, j]
                for i, k in enumerate(self.comp) if abs(self.dK[i, j]) > 1e-12}

    # -- energetics --------------------------------------------------------

    def _mu(self, name, T):
        if name in self.solids:
            return solid_mu(name, T, self.source)
        return mu0(name, T)

    def dG_form(self, T):
        """Reaction Gibbs energy for forming each free species, kJ/mol.

        This is the whole standard-state dependence of G on the extents, formed
        directly rather than as a difference of two large sums.
        """
        mu = np.array([self._mu(s, T) for s in self.species])
        return mu[self.free] + self.dK.T @ mu[self.comp]

    def _dmix(self, n0, dn, n, P_tot):
        """Change in the gas mixing term, in units of RT.

        Every piece is written so that a small change against a large reference
        goes through log1p rather than through a subtraction of two logs.
        """
        s = 0.0
        for i in np.where(self.is_gas)[0]:
            if n[i] < 0:
                return np.inf
            if n0[i] > 0 and abs(dn[i]) < 0.5 * n0[i]:
                s += n0[i] * np.log1p(dn[i] / n0[i]) + dn[i] * np.log(n[i])
            else:
                s += (n[i] * np.log(n[i]) if n[i] > 0 else 0.0)
                s -= (n0[i] * np.log(n0[i]) if n0[i] > 0 else 0.0)
        g = np.where(self.is_gas)[0]
        ng0, d = n0[g].sum(), dn[g].sum()
        ng = ng0 + d
        if ng <= 0:
            return np.inf
        if ng0 > 0 and abs(d) < 0.5 * ng0:
            s -= ng0 * np.log1p(d / ng0) + d * np.log(ng)
        else:
            s -= ng * np.log(ng) - (ng0 * np.log(ng0) if ng0 > 0 else 0.0)
        return s + d * np.log(P_tot)

    def G_rel(self, xi, b, T, P_tot=1.0):
        """Gibbs energy relative to the unreacted charge, in kJ."""
        xi = np.asarray(xi, dtype=float)
        n0 = self.n0(b)
        dn = self.delta_n(xi)
        n = n0 + dn
        floor = -1e-12 * max(1.0, float(np.abs(n0).max()))
        if np.any(n < floor):
            return 1e10
        n = np.maximum(n, 0.0)
        m = self._dmix(n0, dn, n, P_tot)
        if not np.isfinite(m):
            return 1e10
        return float(self.dG_form(T) @ xi) + R_KJ * T * m

    def grad_G_rel(self, xi, b, T, P_tot=1.0):
        """dG/dxi, analytic.

        The mixing term differentiates down to a sum over gases of nu * ln(y P):
        the terms in d(ng)/dxi cancel identically, because the mole fractions
        sum to one. Finite differences on this objective are useless - the whole
        of G_rel can sit at 1e-19 kJ - so the gradient is written out.
        """
        xi = np.asarray(xi, dtype=float)
        n = np.maximum(self.n_of(xi, b), 0.0)
        g = np.where(self.is_gas)[0]
        ng = n[g].sum()
        ly = np.zeros(len(self.species))
        if ng > 0:
            for i in g:
                ly[i] = np.log(max(n[i], 1e-300) / ng * P_tot)
        # dn/dxi: identity on the free block, dK on the component block
        dn = np.zeros((len(self.species), len(self.free)))
        dn[self.free, :] = np.eye(len(self.free))
        dn[self.comp, :] = self.dK
        return self.dG_form(T) + R_KJ * T * (dn[g].T @ ly[g])

    # -- solving -----------------------------------------------------------

    def minimise(self, b, T, seeds, P_tot=1.0, ftol=1e-13):
        """Minimise G over the extents from each seed; lowest result wins.

        The basis is scaled to order one first. G is homogeneous of degree one
        in n - mole fractions do not care how much material there is - so this
        is exact, and it keeps the optimiser's tolerances meaningful against a
        charge of a millimole.
        """
        b = np.asarray(b, dtype=float)
        scale = float(np.abs(b).sum())
        if scale <= 0:
            raise ValueError('empty charge')
        bs = b / scale
        n0s = self.n0(bs)
        if np.any(n0s < -1e-14):
            raise ValueError(f'reference charge is infeasible: {n0s}')

        nf = len(self.free)
        cons = [{'type': 'ineq',
                 'fun': lambda u, i=i, k=k: (n0s[k] + self.dK[i] @ np.exp(u)),
                 'jac': lambda u, i=i: self.dK[i] * np.exp(u)}
                for i, k in enumerate(self.comp)]

        best = None
        for seed in seeds:
            u0 = np.clip(np.log(np.maximum(np.asarray(seed, float) / scale,
                                           np.exp(LOG_FLOOR))),
                         LOG_FLOOR, LOG_CEIL)
            f0 = self.G_rel(np.exp(u0), bs, T, P_tot)
            fs = max(abs(f0), 1e-300) if np.isfinite(f0) and f0 < 1e9 else 1.0
            try:
                r = _minimize(lambda u: self.G_rel(np.exp(u), bs, T, P_tot) / fs,
                              u0, method='SLSQP',
                              jac=lambda u: (self.grad_G_rel(np.exp(u), bs, T, P_tot)
                                             * np.exp(u) / fs),
                              bounds=[(LOG_FLOOR, LOG_CEIL)] * nf,
                              constraints=cons,
                              options={'ftol': ftol, 'maxiter': 500})
            except Exception:
                continue
            xi = np.exp(np.clip(r.x, LOG_FLOOR, LOG_CEIL))
            g = self.G_rel(xi, bs, T, P_tot)
            if not np.isfinite(g) or g > 1e9:
                continue
            if best is None or g < best[0]:
                best = (g, xi)
        if best is None:
            return None
        g, xi = best
        n = self.n_of(xi * scale, b)
        return {'n': np.maximum(n, 0.0), 'xi': xi * scale,
                'G_rel_kJ': g * scale,
                'residual': float(np.max(np.abs(self.E.T @ n - b))),
                'method': 'gibbs_min', 'formulation': 'reaction_extent',
                'n_seeds': len(seeds)}

    # -- readout -----------------------------------------------------------

    def gas_fractions(self, n):
        g = np.where(self.is_gas)[0]
        tot = n[g].sum()
        return {self.species[i]: (n[i] / tot if tot > 0 else 0.0) for i in g}

    def reduced_percent(self, n, n_metal_total):
        """Percentage of the metal present in its reduced state, M(3+)."""
        if n_metal_total <= 0:
            return 0.0
        red = sum(n[self.species.index(p)] * METAL_PHASES[p][1]
                  for p in self.solids if p in METAL_PHASES)
        return min(red / n_metal_total, 1.0) * 100.0

    def phase_split(self, n, n_metal_total):
        if n_metal_total <= 0:
            return {}
        return {p: n[self.species.index(p)] * METAL_PHASES[p][0] / n_metal_total * 100.0
                for p in self.solids if p in METAL_PHASES}

    # Kept so the titanium scripts and tests read unchanged.
    ti3_percent = reduced_percent
