"""The defect equation registry, (S1)-(S21).

Rendered by export_equations.py into equations_statmech.html and inlined
into the Model panel of the Defect stat mech workspace. Same row kinds as
equations_registry.py; the S-prefix keeps the numbering apart from the
thermodynamic sheet's (1)-(30).

The sheet covers two engines. Sections 1 and 2 are the free energy and
they are shared: solidgas/statmech.py partitions it over three site
classes, solidgas/particle/ over a depth coordinate. Sections 3 onward
say which engine each block belongs to.
"""

R = []
R.append(('title', 'The defect model, as implemented'))
R.append(('sub', 'Oxygen vacancies with electronic compensation · solidgas/statmech.py and solidgas/particle/'))
R.append(('gap', 0.10))

R.append(('head', '1 · The free energy of a compensated vacancy'))
R.append(('note', r'Removing a neutral oxygen leaves two electrons behind and in rutile they localise as two $\mathrm{Ti^{3+}}$ small polarons. A vacancy is therefore one event on the oxygen sublattice and two on the cation sublattice, and all three carry configurational entropy. Rutile has one Ti per two O, so $N_{\mathrm{Ti}} = N_O/2$ and $\theta_p = 4\theta$.'))
R.append(('eq', r'$\frac{f(\theta)}{N_O} \;=\; E\,\theta \;+\; k_B T [\theta \ln \theta + (1-\theta)\ln(1-\theta)] \;+\; \frac{k_B T}{2} [\theta_p \ln \theta_p + (1-\theta_p)\ln(1-\theta_p)]\,,\qquad \theta_p = 4\theta$', 'S1'))
R.append(('boxeq', r'$\mu_V(\theta) \;=\; \frac{\partial f}{\partial \theta} \;=\; E \;+\; k_B T\left[\, \ln\frac{\theta}{1-\theta} \;+\; 2\ln\frac{4\theta}{1-4\theta} \right]$', 'S2'))
R.append(('note', r'The second logarithm diverges at $\theta = 1/4$. That composition is $\mathrm{TiO_{1.5}}$, which is $\mathrm{Ti_2O_3}$: nothing imposes the ceiling, the cation sublattice simply runs out of $\mathrm{Ti^{3+}}$ sites. Halving the gap to the quarter costs $2 k_B T \ln 2$, without bound.'))
R.append(('eq', r'$\theta \;\longrightarrow\; 4^{-2/3}\, e^{(\mu_V - E)/3k_B T}\ \ (\theta \to 0)\,,\qquad \mu_V = -\mu_O = -\frac{1}{2}k_B T \ln p_{\mathrm{O_2}} \;\;\Longrightarrow\;\; \theta \propto p_{\mathrm{O_2}}^{-1/6}$', 'S3'))
R.append(('note', r'The textbook exponent for a doubly ionised vacancy with electronic compensation. Written on the oxygen sublattice alone the two logarithms collapse to one, the exponent becomes $-1/2$ - a NEUTRAL vacancy - and the enrichment between two sites is the cube of the right one.'))

R.append(('head', '2 · Where the energy comes from'))
R.append(('eq', r'$\epsilon_{c} \in \{\epsilon_s,\ \epsilon_{ss},\ \epsilon_b\} \;=\; \{0,\ +0.59,\ +1.31\}\ \mathrm{eV\ (sX)}\,;\qquad \{0,\ +0.39,\ +1.08\}\ \mathrm{eV\ (GGA)}$', 'S4'))
R.append(('note', r'Li, Guo and Robertson, J. Phys. Chem. C 119 (2015): absolute neutral $V_O$ formation energies on rutile (110), surface / first subsurface / bulk. Relative values are used; the oxygen chemical potential cancels in the canonical ensemble. Polaron relaxation is inside $\epsilon_c$, so S1 counts the polaron ENTROPY and not its energy twice.'))
R.append(('boxeq', r'$E(z) \;=\; \epsilon_b \;-\; \Delta E_{\mathrm{seg}}\, e^{-z/\xi}\,,\qquad \Delta E_{\mathrm{seg}} = \epsilon_b - \epsilon_s\,,\qquad \xi \;=\; \frac{-\,d_{110}}{\ln\![(\epsilon_b - \epsilon_{ss})/(\epsilon_b - \epsilon_s)]}$', 'S5'))
R.append(('note', r'Two parameters, and the DFT set has three energies with the bulk as the asymptote, so both are determined and no profile shape is assumed: $\xi = 0.543$ nm on sX, $0.725$ nm on GGA - one and two thirds and two and a quarter $d_{110}$ layers. This replaces a declared span fraction between the subsurface and bulk values, which the dataset labelled assumed; the exponential reproduces it to within a few hundredths without being told to.'))

R.append(('head', '3 · Particle-scale matching · solidgas/particle/'))
R.append(('note', r'One chemical potential is shared by every oxygen site in the particle - that is what equilibrium means - and $\mu_V$ is fixed by the one thing the experiment measures, the total oxygen removed. Two pools: the bridging row, an areal density of one O per $(1{\times}1)$ cell, and the interior, a stack of $d_{110}$ oxygen layers at the bulk density with layer $k$ centred at depth $k\,d_{110}$.'))
R.append(('eq', r'$\langle\theta\rangle_{\mathrm{int}}(\mu_V) \;=\; \frac{\sum_k V_k\ \theta(\mu_V,\, E(k\,d_{110}))}{\sum_k V_k}\,,\qquad V_k \;=\; \frac{1}{3}[(r_{\mathrm{in}} - (k-1)d_{110})^3 - (r_{\mathrm{in}} - k\,d_{110})^3]$', 'S6'))
R.append(('note', r'A sum, not an integral. Smearing the stack into a continuum - which is what the ceria formulation this construction follows does - costs 0.5 % in the interior mean at the rutile decay length, because the energy changes by a large factor from one layer to the next. Both forms are implemented and a gate holds them apart.'))
R.append(('boxeq', r'$N_{\mathrm{row}}\,\theta_{\mathrm{row}}^{\mathrm{phase}}(\mu_V) \;+\; N_{\mathrm{int}}\,\langle\theta\rangle_{\mathrm{int}}(\mu_V) \;+\; N_{\mathrm{ext}} \;=\; N_V^{\mathrm{tot}}$', 'S7'))
R.append(('eq', r'$N_{\mathrm{row}} = A\,\sigma_{\mathrm{br}}\,,\qquad \sigma_{\mathrm{layer}} = 4\,\sigma_{\mathrm{br}}\,,\qquad A = \frac{6}{\rho\, d}\ \ \mathrm{or\ BET}\,,\qquad r = \frac{3}{\rho A}$', 'S8'))
R.append(('note', r'Real geometric site counts. At $d = 0.9\ \mu\mathrm{m}$, $N_{\mathrm{row}} = 13.6\ \mu\mathrm{mol\text{-}O/g}$ against $25\,042\ \mu\mathrm{mol\text{-}O/g}$ of oxygen sites. The row is a quarter of a layer and is kept out of the continuum for that reason: giving it a bulk density would quadruple the cheapest sites in the particle.'))

R.append(('head', '4 · The two phases that cut the profile off'))
R.append(('eq', r'$\theta_{\mathrm{row}}^{\mathrm{phase}} = \theta(\mu_V,\epsilon_s)\ \left[\mu_V\!<\!\mu_t\right];\ \ \ \theta_t \rightarrow 0.5\ \mathrm{lever\ rule}\ \left[\mu_V\!=\!\mu_t\right];\ \ \ 0.5\ \left[\mu_V\!>\!\mu_t\right]\,,\qquad \mu_t = \mu_V(\theta_t,\, \epsilon_s)\,,\ \ \theta_t = 0.17$', 'S9'))
R.append(('note', r'First-order transition to the $1\!\times\!2$ added-row $\mathrm{Ti_2O_3}$ line phase (Onishi-Iwasawa stoichiometry: one O vacancy per two bridging sites). It has no configurational entropy of its own, so $\mu_V$ pins while the row converts; at $0.9\ \mu\mathrm{m}$ and 600 C the window is $33.2 \rightarrow 37.7\ \mu\mathrm{mol\text{-}O/g}$. The high-coverage $1\!\times\!1$ branch is excluded from the phase set.'))
R.append(('eq', r'$\mu_{cs} \;=\; \mu_V(\theta_{\mathrm{sol}},\, \epsilon_b)\,,\qquad \theta_{\mathrm{sol}} = x_{\mathrm{sol}}/2 = 0.004\,,\qquad N_{\mathrm{ext}} \;=\; N_V^{\mathrm{tot}} - N_V(\mu_{cs})$', 'S10'))
R.append(('note', r'Rutile does not dissolve point defects all the way to $\mathrm{Ti_2O_3}$: past $x \approx 0.008$ crystallographic shear planes form, so $\mu_V$ pins again and the excess leaves the point-defect count as extended defects. This is the one part of the construction the ceria system has no analogue for. At 600 C and $0.9\ \mu\mathrm{m}$ the point-defect ceiling is $112.4\ \mu\mathrm{mol\text{-}O/g}$ - the shipped loading is 85 % of it.'))

R.append(('head', '5 · Is the equilibrium reachable? · particle/transport.py'))
R.append(('eq', r'$D \;=\; a^2 \nu_0\, e^{-E_m/k_B T}\,,\qquad \tau = \frac{D t}{R^2}\,,\qquad \frac{M(t)}{M_\infty} \;=\; 1 - \frac{6}{\pi^2}\sum_{n=1}^{\infty} \frac{1}{n^2}\, e^{-n^2 \pi^2 \tau}$', 'S11'))
R.append(('note', r'Crank sphere diffusion. An equilibrium profile is only the answer if the oxygen had time to arrange itself. With the literature bulk barrier $E_m = 0.65$ eV a $0.9\ \mu\mathrm{m}$ particle is 99 % homogenised in 0.54 ms at 600 C, six orders below a 300 s exposure.'))
R.append(('boxeq', r'$E_m^{\ast} \;=\; k_B T\, \ln\!\frac{a^2 \nu_0}{\tau_{0.99} R^2 / t_{\mathrm{exp}}} \;=\; 1.645\ \mathrm{eV}\qquad \mathrm{(the\ barrier\ that\ would\ make\ transport\ the\ limit)}$', 'S12'))
R.append(('note', r'Claiming transport is fast is only worth something with the distance to it not being fast. The margin over the literature barrier is 0.995 eV; the recovery measured in this work implies a higher effective barrier than the literature one, so the margin is finite.'))

R.append(('head', '6 · Is local electroneutrality allowed? · particle/spacecharge.py'))
R.append(('eq', r'$\lambda_D \;=\; \sqrt{\frac{\varepsilon_r \varepsilon_0\, k_B T}{\sum_i z_i^2 e^2 n_i}}\,,\qquad \sum_i z_i^2 n_i \;=\; (4\theta + \frac{1}{2}\cdot 4\theta)\, n_O$', 'S13'))
R.append(('note', r'The profile assumes each vacancy drags its two polarons with it pointwise, which is the ceria study’s Model 3. Its Model 4 lets the two species equilibrate independently in a self-consistent potential and found the difference confined to the outermost nanometre - and said explicitly that the collapse followed from how CONCENTRATED that system was. This one runs 13 times more dilute. At the bulk plateau $\lambda_D = 0.57$ nm against $\xi = 0.54$ nm: the same length, so screening and chemistry cannot be separated, and over $\varepsilon_r = 50 \rightarrow 170$ the ratio only moves from 0.74 to 1.37. The full Poisson-Boltzmann solve is not implemented; the error it would quantify runs towards MORE surface enrichment.'))

R.append(('head', '7 · Site-resolved lattice sampling · solidgas/statmech.py'))
R.append(('eq', r'$E(\mathbf{n}) \;=\; \sum_i \epsilon_{c(i)}\, n_i \;+\; \frac{1}{2} \sum_{i \neq j} V_{ij}\, n_i n_j$', 'S14'))
R.append(('note', r'$n_i \in \{0,1\}$: oxygen vacancy at site $i$;  $N_V = \sum_i n_i$ is fixed by the measured inventory (canonical, not grand-canonical). $V_{ij}$: along-row separations $\Delta = 1, 2, 3$ bridging sites and cross-row pairs inside the 10 $\mathrm{\AA}$ cutoff. This layer produces the ordering exhibit - spacings and the reconstruction onset - not the partition.'))
R.append(('eq', r'$V_O(i) + O(j) \;\rightarrow\; O(i) + V_O(j)\,,\qquad P_{\mathrm{accept}} \;=\; \min\!\left[ 1,\ e^{-\Delta E / k_B T} \right]$', 'S15'))
R.append(('eq', r'$\mu_V \;=\; -k_B T\, \ln\! \left[ \frac{1}{N_V + 1} \left\langle \sum_{i\ \mathrm{empty}} e^{-\Delta E_i^{\mathrm{ins}} / k_B T} \right\rangle \right]$', 'S16'))
R.append(('note', r'Widom insertion on the sampling slab; verified against exact generating-polynomial ensembles and exact enumerations. Those enumerations build the oxygen sublattice alone, so they certify the interaction bookkeeping of S14, not the compensation of S1.'))

R.append(('head', '8 · CO<sub>2</sub> accessibility (placeholder kinetics)'))
R.append(('eq', r'$k_c \;=\; A_c\ p_{\mathrm{CO_2}}^{\,m}\ e^{-E_{a,c}/k_B T}\,,\qquad P_c(t) \;=\; 1 - e^{-k_c\, t}$', 'S17'))
R.append(('boxeq', r'$V_O^{\mathrm{acc}} \;=\; \sum_c N_c\, \theta_c\, P_c\,,\qquad f_{\mathrm{rec}} \;=\; V_O^{\mathrm{acc}} \,/\, V_O^{\mathrm{tot}}$', 'S18'))
R.append(('note', r'First-order refill per class; $A_c$ are effective prefactors (deeper classes are transport-limited), to be replaced by fitted or DFT values.'))

R.append(('head', '9 · Dielectric correlation (empirical)'))
R.append(('eq', r"$\varepsilon'' \;=\; 10^{\,a}\ \left(V_O^{\mathrm{tot}}\right)^{b}\,,\qquad a = -3.070 \pm 0.695\,,\ \ b = 1.224 \pm 0.230\,,\qquad \tan\delta \;=\; \varepsilon'' / \varepsilon'$", 'S19'))
R.append(('note', r'Linear regression of $\log_{10}\varepsilon''$ on $\log_{10}V_O$ measured in this work ($N = 4$, $R^2 = 0.934$); an extrapolation below the fitted 95-5000 $\mu\mathrm{mol\text{-}O/g}$ range.'))

R.append(('head', '10 · Reoxidation transport — coarse-grained KMC'))
R.append(('eq', r'$R(k \to l) \;=\; \Lambda_{kl}\ \eta_k\ (q_l - \eta_l)\ e^{\,\beta \epsilon_k}\,,\qquad \Lambda_{kl} \;=\; \frac{\nu_0\ a^2 A_{kl}\ e^{-\beta E_{TS}}}{d_{kl}\ v_O\ q_k\ q_l}$', 'S20'))
R.append(('note', r'Radial coarse cells over the particle, occupancies $\eta_k \leq q_k$ (generalised exclusion); framework: Katsoulakis and Vlachos, J. Chem. Phys. 119, 9412 (2003). $E_{TS} = E_m + \epsilon_b$ and $E_{\mathrm{act}}(k \to l) = E_{TS} - \epsilon_k$.'))
R.append(('boxeq', r'$R(k \to l)\,\mu(\eta) \;=\; R(l \to k)\,\mu(\eta + \delta_l - \delta_k)\,,\qquad \mu(\eta) \propto \prod_k \binom{q_k}{\eta_k}\, e^{-\beta \epsilon_k \eta_k}$', 'S21'))
R.append(('note', r'Detailed balance against a product-binomial measure on the oxygen sublattice, with the surface refill as the one irreversible channel. OPEN ITEM: that measure is the NEUTRAL-vacancy equilibrium. S1 is not, so the long-time limit of this kinetic layer and the equilibrium partition above no longer agree, and reconciling them is tracked in docs/roadmap.md.'))
