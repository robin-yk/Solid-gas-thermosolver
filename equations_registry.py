"""The equation registry: the complete system the solver solves, (1)-(30).

Single source for both renderings - build_equations.py (PNG/PDF sheet) and
export_equations.py (per-equation SVG fragment inlined into the page's
Method tab). Row kinds: title, sub, head, eq (math, number), boxeq, note,
gap.
"""

R = []
R.append(('title', 'The complete system the solver solves'))
R.append(('sub', 'Ti–O gas–solid Gibbs equilibrium · solidgas/activeset.py · every output is this one minimisation'))
R.append(('gap', 0.10))

R.append(('head', '1 · Inputs → element inventory'))
R.append(('eq', r'$n_{g,0} \;=\; \frac{P\,V}{R_{\mathrm{atm}}\,T_{\mathrm{charge}}}$', '1'))
R.append(('eq', r'$y_{i,0} \;=\; \frac{r_i}{\sum_k r_k}\,,\qquad n_{i,0} \;=\; y_{i,0}\,n_{g,0}$', '2'))
R.append(('eq', r'$n_{s,0} \;=\; \frac{m_{\mathrm{solid}}}{t_h W_{\mathrm{Ti}} + o_h W_{\mathrm{O}}}$', '3'))
R.append(('eq', r'$b_e \;=\; \sum_i a_{ie}\,n_{i,0} \;+\; a_{he}\,n_{s,0}\,,\qquad e \in \{\mathrm{C,H,N,O,Ti}\}$', '4'))
R.append(('note', r'$i \in \{\mathrm{CO_2,\ H_2,\ CO,\ H_2O,\ N_2,\ O_2}\}\,;\qquad j \in \{\mathrm{TiO_2,\ Ti_{10}O_{19},\ Ti_8O_{15},\ Ti_6O_{11},\ Ti_5O_9,\ Ti_4O_7,\ Ti_3O_5,\ Ti_2O_3}\} = \mathrm{Ti}_{t_j}\mathrm{O}_{o_j}$'))

R.append(('head', '2 · Standard chemical potentials'))
R.append(('note', 'gases: NIST-JANAF Shomate,  t = T/1000 :'))
R.append(('eq', r'$H^{\circ}(T) = \Delta_f H^{\circ}_{298} + A t + \frac{B t^{2}}{2} + \frac{C t^{3}}{3} + \frac{D t^{4}}{4} - \frac{E}{t} + F - H$', '5'))
R.append(('eq', r'$S^{\circ}(T) = A\,\ln t + B t + \frac{C t^{2}}{2} + \frac{D t^{3}}{3} - \frac{E}{2 t^{2}} + G\,,\qquad \mu^{\circ}(T) = H^{\circ}(T) - \frac{T\,S^{\circ}(T)}{1000}$', '6'))
R.append(('note', 'condensed phases: Waldner & Eriksson CALPHAD, integrated piecewise from 298.15 K :'))
R.append(('eq', r'$\mu^{\circ}(T) = \frac{1}{1000}\left[\, \Delta_f H^{\circ}_{298} + \sum_{\mathrm{seg}}\left(\Delta H_{tr} + \int C_p\, dT\right) \;-\; T\left( S^{\circ}_{298} + \sum_{\mathrm{seg}}\left(\frac{\Delta H_{tr}}{T_{tr}} + \int \frac{C_p}{T}\, dT\right)\right) \right]$', '7'))
R.append(('eq', r'$C_p = a + bT + cT^{2} + d\,T^{-2} \ \ \mathrm{(poly)}\,;\qquad C_p = k_0 + k_1 T^{-1/2} + k_2 T^{-2} + k_3 T^{-3} \ \ \mathrm{(magn\acute{e}li)}$', '8'))

R.append(('head', '3 · The master problem'))
R.append(('boxeq', r'$\min_{n,\,m}\ \ G(n,m) \;=\; \sum_i n_i\left[\, \mu^{\circ}_i + RT\,\ln\!\left( \frac{n_i}{N_g}\cdot\frac{P}{P^{\circ}} \right) \right] \;+\; \sum_j m_j\, \mu^{\circ}_j$', '9'))
R.append(('boxeq', r'$\mathrm{s.t.}\ \ \sum_i a_{ie}\, n_i + \sum_j a_{je}\, m_j = b_e \ \ \forall e\,,\qquad n_i \geq 0\,,\ \ m_j \geq 0\,,\qquad N_g = \sum_i n_i$', '10'))
R.append(('note', r'ideal gas mixture; each condensed phase a pure line compound, $a_j = 1$;  $P^{\circ} = 1\ \mathrm{atm}$,  $R = 8.314472\times10^{-3}\ \mathrm{kJ\,mol^{-1}K^{-1}}$'))

R.append(('head', '4 · KKT conditions, solved per candidate active set  A'))
R.append(('eq', r'$n_i > 0: \qquad \mu^{\circ}_i + RT\,\ln\!\left( y_i P / P^{\circ} \right) \;=\; \sum_e a_{ie}\, \lambda_e$', '11'))
R.append(('eq', r'$j \in A: \qquad \mu^{\circ}_j \;=\; t_j\, \lambda_{\mathrm{Ti}} + o_j\, \lambda_{\mathrm{O}}\,,\qquad m_j > 0$', '12'))
R.append(('eq', r'$j \notin A: \qquad m_j = 0\,,\qquad r_j \;\equiv\; \mu^{\circ}_j - t_j \lambda_{\mathrm{Ti}} - o_j \lambda_{\mathrm{O}} \;\geq\; 0$', '13'))
R.append(('eq', r'$\sum_i a_{ie}\, y_i N_g + \sum_j a_{je}\, m_j = b_e\,, \qquad \sum_i y_i = 1$', '14'))

R.append(('head', '5 · Decomposed solve'))
R.append(('note', r'$|A| = 2$  (two-phase buffer; the pair pins both solid potentials linearly):'))
R.append(('eq', r'$t_1 \lambda_{\mathrm{Ti}} + o_1 \lambda_{\mathrm{O}} = \mu^{\circ}_1 \qquad \wedge \qquad t_2 \lambda_{\mathrm{Ti}} + o_2 \lambda_{\mathrm{O}} = \mu^{\circ}_2$', '15'))
R.append(('note', r'gas Newton over  $x = (\theta_e,\ \ln N_g)$,  $\theta_e = \lambda_e / RT$ :'))
R.append(('eq', r'$\ln y_i \;=\; \sum_e a_{ie}\, \theta_e \;-\; \frac{\mu^{\circ}_i}{RT} \;-\; \ln\frac{P}{P^{\circ}}$', '16'))
R.append(('eq', r'$F_e = \sum_i a_{ie}\, y_i N_g - \hat b_e = 0\,, \qquad F_0 = \sum_i y_i - 1 = 0$', '17'))
R.append(('eq', r'$\frac{\partial F_e}{\partial \theta_f} = \sum_i a_{ie} a_{if}\, y_i N_g\,, \qquad \frac{\partial F_e}{\partial \ln N_g} = \sum_i a_{ie}\, y_i N_g\,, \qquad x \leftarrow x - \alpha\, J^{-1} F\,,\ \ \alpha: \Vert F \Vert \downarrow$', '18'))
R.append(('note', r'condensed amounts recovered linearly; the combination  $\mathrm{O}^{\prime} = \mathrm{O} - \rho_h\,\mathrm{Ti}$  removes the host oxygen symbolically:'))
R.append(('eq', r'$\sum_j t_j\, m_j = b_{\mathrm{Ti}}\,, \qquad \sum_i a_{i\mathrm{O}}\, y_i N_g \;=\; b^{\,\mathrm{gas}}_{\mathrm{O}} + \sum_j d_j\, m_j\,, \qquad d_j \equiv \rho_h t_j - o_j\,,\ \ \rho_h = o_h / t_h$', '19'))
R.append(('eq', r'$|A| = 1: \qquad m_1 = b_{\mathrm{Ti}} / t_1\,, \qquad \lambda_{\mathrm{O}} = RT\,\theta_{\mathrm{O}}\,, \qquad \lambda_{\mathrm{Ti}} = \left( \mu^{\circ}_1 - o_1 \lambda_{\mathrm{O}} \right) / t_1$', '20'))
R.append(('eq', r'$|A| = 3\ \mathrm{(probe)}: \qquad \mathrm{res}_3 = \mu^{\circ}_3 - t_3 \lambda_{\mathrm{Ti}} - o_3 \lambda_{\mathrm{O}}\,; \qquad |\mathrm{res}_3| < \varepsilon_{\mathrm{deg}} \ \Leftrightarrow\ \mathrm{invariant\ point}$', '21'))

R.append(('head', '6 · Decision rules'))
R.append(('eq', r'$G\ \mathrm{convex} \ \Rightarrow\ \left[\, m_j > 0\ \ \forall j \in A \ \ \wedge\ \ r_j \geq 0\ \ \forall j \notin A \,\right] \ \Leftrightarrow\ \mathrm{global\ minimum}$', '22'))
R.append(('eq', r'$|r_j| / d_j < \varepsilon_{\mathrm{deg}} \ \Rightarrow\ \mathrm{degenerate\ boundary\ (reported,\ never\ resolved\ arbitrarily)}$', '23'))
R.append(('eq', r'$b_e = 0 \ \Rightarrow\ n_i = 0\ \ \forall\, i \ni e \quad \mathrm{(removed\ from\ the\ variables:\ exactly\ zero)}$', '24'))
R.append(('eq', r'$\mathrm{avail}_{\mathrm{O}} = 0 \ \wedge\ \exists\, j \notin A: d_j > 0 \ \Rightarrow\ \frac{\partial G}{\partial \xi}|_{0^{+}} \sim \frac{RT}{2} \ln \xi \to -\infty \ \Rightarrow\ \mathrm{reject}\ A$', '25'))

R.append(('head', '7 · Reduced cost, two bases'))
R.append(('eq', r'$r_j^{\mathrm{formula}} = \mu^{\circ}_j - t_j \lambda_{\mathrm{Ti}} - o_j \lambda_{\mathrm{O}} \ \ \left[\mathrm{kJ/mol\ compound}\right]\,, \qquad r_j^{\mathrm{O}} = r_j^{\mathrm{formula}} / d_j \ \ \left[\mathrm{kJ/mol\ O}\right]$', '26'))
R.append(('note', r'every $\mathrm{TiO_2} \rightarrow$ Magnéli step has $d_j = 1$, so the two bases coincide for these steps'))

R.append(('head', '8 · Derived outputs'))
R.append(('eq', r'$\mathrm{Ti}^{3+}\% = 100\, \frac{\sum_j c_j\, m_j}{b_{\mathrm{Ti}}}\,,\ \ c_j = 2\ (j \neq \mathrm{TiO_2})\,;\qquad x = 2 - \frac{\sum_j o_j m_j}{\sum_j t_j m_j}$', '27'))
R.append(('eq', r'$\Delta O_{\mathrm{gas}} = \sum_i a_{i\mathrm{O}} ( n_i - n_{i,0} )\,,\qquad \Delta n_i = n_i - n_{i,0}\,,\qquad X_{\mathrm{CO_2}} = \frac{n_{\mathrm{CO_2},0} - n_{\mathrm{CO_2}}}{n_{\mathrm{CO_2},0}}$', '28'))
R.append(('eq', r'$G_{\mathrm{rel}} = G(n,m) - G(n_0, n_{s,0})$', '29'))
R.append(('eq', r'$\lambda_{\mathrm{O}} = \mathrm{oxygen\ potential\ (a\ Lagrange\ multiplier\ of\ the\ minimisation,\ not\ an\ input)}\,,\qquad p_{\mathrm{O_2}} = \exp\!\left[ \frac{2 \lambda_{\mathrm{O}} - \mu^{\circ}_{\mathrm{O_2}}}{RT} \right]$', '30'))
