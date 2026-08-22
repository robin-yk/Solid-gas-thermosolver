"""The defect stat-mech equation registry, (S1)-(S12).

Rendered by export_equations.py into equations_statmech.html and inlined
into the Model panel of the Defect stat mech workspace. Same row kinds as
equations_registry.py; the S-prefix keeps the numbering apart from the
thermodynamic sheet's (1)-(30).
"""

R = []
R.append(('title', 'The defect model, as implemented'))
R.append(('sub', 'Canonical oxygen-vacancy lattice · solidgas/statmech.py'))
R.append(('gap', 0.10))

R.append(('head', '1 · Canonical lattice model'))
R.append(('eq', r'$E(\mathbf{n}) \;=\; \sum_i \epsilon_{c(i)}\, n_i \;+\; \frac{1}{2} \sum_{i \neq j} V_{ij}\, n_i n_j$', 'S1'))
R.append(('note', r'$n_i \in \{0,1\}$: oxygen vacancy at site $i$;  $N_V = \sum_i n_i$ is fixed by the measured inventory (canonical, not grand-canonical)'))
R.append(('eq', r'$\epsilon_{c} \in \{\epsilon_s,\ \epsilon_{ss},\ \epsilon_b\} \;=\; \{0,\ +0.59,\ +1.31\}\ \mathrm{eV\ (sX)}\,;\qquad \{0,\ +0.39,\ +1.08\}\ \mathrm{eV\ (GGA)}$', 'S2'))
R.append(('note', r'$V_{ij}$: along-row separations $\Delta = 1, 2, 3$ bridging sites and cross-row pairs inside the 10 $\mathrm{\AA}$ cutoff; polaron relaxation is inside $\epsilon_c$'))
R.append(('eq', r'$V_O(i) + O(j) \;\rightarrow\; O(i) + V_O(j)\,,\qquad P_{\mathrm{accept}} \;=\; \min\!\left[ 1,\ e^{-\Delta E / k_B T} \right]$', 'S3'))
R.append(('eq', r'$\mu_V \;=\; -k_B T\, \ln\! \left[ \frac{1}{N_V + 1} \left\langle \sum_{i\ \mathrm{empty}} e^{-\Delta E_i^{\mathrm{ins}} / k_B T} \right\rangle \right]$', 'S4'))
R.append(('note', r'Widom insertion on the sampling slab; verified against exact generating-polynomial ensembles and exact enumerations'))

R.append(('head', '2 · Particle-scale matching'))
R.append(('eq', r'$\theta_c(\mu_V) \;=\; \left[\, 1 + e^{(\epsilon_c - \mu_V)/k_B T} \right]^{-1} \ \ (c = ss,\, b)\,;\qquad \theta_s(\mu_V)\ \mathrm{from\ the\ MC\ isotherm}$', 'S5'))
R.append(('eq', r'$\theta_{\mathrm{rec}} \simeq 0.17\!\text{-}\!0.20 \qquad \mathrm{(reported\ onset\ interval\ for\ the\ 1\!\times\!2\ reconstruction;\ validity\ boundary,\ not\ a\ cap)}$', 'S6'))
R.append(('boxeq', r'$N_s\,\theta_s(\mu_V) \;+\; N_{ss}\,\theta_{ss}(\mu_V) \;+\; N_b\,\theta_b(\mu_V) \;=\; N_V^{\mathrm{tot}}$', 'S7'))
R.append(('eq', r'$N_s = A\,\sigma_{\mathrm{br}}\,,\qquad N_{ss} = A\,\sigma_{\mathrm{layer}}\,,\qquad A = \frac{6}{\rho\, d}\ \ \mathrm{or\ BET}$', 'S8'))
R.append(('note', r'real geometric site counts, not the slab ratio: at $d = 0.9\ \mu\mathrm{m}$, $N_s = 13.6\ \mu\mathrm{mol\text{-}O/g}$ against $25\,042\ \mu\mathrm{mol\text{-}O/g}$ of O sites'))

R.append(('head', '3 · CO<sub>2</sub> accessibility (placeholder kinetics)'))
R.append(('eq', r'$k_c \;=\; A_c\ p_{\mathrm{CO_2}}^{\,m}\ e^{-E_{a,c}/k_B T}$', 'S9'))
R.append(('eq', r'$P_c(t) \;=\; 1 - e^{-k_c\, t}$', 'S10'))
R.append(('boxeq', r'$V_O^{\mathrm{acc}} \;=\; \sum_c N_c\, \theta_c\, P_c\,,\qquad f_{\mathrm{rec}} \;=\; V_O^{\mathrm{acc}} \,/\, V_O^{\mathrm{tot}}$', 'S11'))
R.append(('note', r'first-order refill per class; $A_c$ are effective prefactors (deeper classes are transport-limited), to be replaced by fitted or DFT values'))

R.append(('head', '4 · Dielectric correlation (empirical)'))
R.append(('eq', r"$\varepsilon'' \;=\; a\, V_O^{\mathrm{tot}} + b\,,\qquad \tan\delta \;=\; \varepsilon'' / \varepsilon'$", 'S12'))
R.append(('note', r'coefficients from the experimental fit of this work; blank until set in rutile_dft.json'))
