"""What every manuscript figure claims, and the vocabulary it may use.

The registry is metadata only: page geometry, the colour-role table, and one
entry per plate carrying its manuscript section, the claim it discharges, and
a caption written with {{slots}} instead of numbers. Slots are resolved at
build time against data/figure_data.json, which is extracted from the
certified engines, so a caption cannot drift away from the model. The gate in
tests/test_plates.py enforces both halves of that rule: every slot must
resolve, and no decimal number may appear in caption prose outside a slot.

Two solvers, two independent tracks. They share the page and the drawing kit;
they do not share a figure.
"""

# ------------------------------------------------------------------ page

SPEC = {
    'journal': 'Nature',
    # Nature artwork: 89 mm single column, 183 mm double column,
    # 247 mm maximum height, sans-serif, nothing below 5 pt at print size.
    'width_mm': {'single': 89, 'onehalf': 120, 'double': 183},
    'max_height_mm': 247,
    'type_floor_pt': 5.0,
    'type_body_pt': 7.0,
    'font': 'Arial, Helvetica, sans-serif',
    'mm_per_pt': 25.4 / 72.0,
}

# One hue per physical role, identical in every plate. Okabe-Ito, so the set
# survives colour-blind readers and greyscale printing; structure is never
# coloured, which is what keeps a hue meaningful.
ROLES = [
    ('surface', '#0072B2',
     'The rutile (110) bridging-oxygen class and its (1x2) line phase.'),
    ('subsurface', '#D55E00',
     'The declared subsurface shell: the oxygen layers below the bridging '
     'plane.'),
    ('bulk', '#CC79A7',
     'Bulk oxygen sites and the bulk-class occupancy.'),
    ('extended', '#E69F00',
     'Extended defects: the crystallographic-shear precipitate beyond the '
     'solubility ceiling.'),
    ('gas', '#009E73',
     'The gas charge and everything it carries: feed composition, oxygen '
     'potential, CO2 refill.'),
    ('structure', '#6E6E6E',
     'Axes, boundaries, inactive phases, and anything held fixed.'),
]

# ------------------------------------------------------------- figures

FIGURES = [
    # ---------------------------------------------------- defect track
    {
        'id': 'd1', 'track': 'defect', 'label': 'Figure 1',
        'width': 'double',
        'section': 'Defect thermodynamics - site capacity and the phase '
                   'ladder',
        'claim': 'Particle geometry, not site energy, is what forces the '
                 'inventory below the surface; the ladder that follows is a '
                 'set of certified numbers.',
        'caption':
            '**Site capacity and the phase ladder of reduced rutile.** '
            '**a**, Oxygen-site capacity of the three classes for a '
            '{{defect.geometry.d_um}} um particle '
            '({{defect.geometry.area_m2_g}} m2 g-1): the bridging-oxygen row '
            'holds {{defect.geometry.N_s}} umol-O g-1, '
            '{{defect.geometry.surface_pct_of_O}}% of all oxygen sites, and the '
            'subsurface shell {{defect.geometry.N_ss}}, so both saturate '
            'long before a measured inventory is accommodated. '
            '**b**, Class inventory against total vacancy content at '
            '{{defect.flagship.T_C}} C. The (1x2) reconstruction enters as a '
            'first-order surface phase at '
            '{{defect.boundaries.VO_onset_umol_g}} umol-O g-1 and completes '
            'at {{defect.boundaries.VO_recon_complete_umol_g}}; beyond the '
            'estimated shear-plane ceiling '
            '({{defect.boundaries.VO_cs_umol_g}}) the vacancy chemical '
            'potential pins and the excess precipitates as extended '
            'defects. Dashed verticals are the closed-form boundaries; the '
            'marker is the measured inventory '
            '({{defect.flagship.VO_total_umol_g}} umol-O g-1).',
        'symbols': ['N_s', 'N_ss', 'N_b', 'mu_V', 'theta_c'],
        'methods': ['site-exclusion isotherm', 'lever rule',
                    'closed-form phase boundaries'],
    },
    {
        'id': 'd2', 'track': 'defect', 'label': 'Figure 2',
        'width': 'double',
        'section': 'Defect thermodynamics - the measured inventory',
        'claim': 'At the measured inventory the surface is reconstructed and '
                 'saturated, and the near-surface shell - not the bulk - '
                 'carries the majority of the vacancies.',
        'caption':
            '**Where a measured inventory of '
            '{{defect.flagship.VO_total_umol_g}} umol-O g-1 sits.** '
            '**a**, Cross-section of the outer shell, drawn to scale: one '
            'dot is one vacancy on one oxygen site. The bridging plane '
            'carries one oxygen per (1x1) column and each '
            '{{defect.geometry.layer_nm}} nm layer below it carries four, so '
            'the subsurface class is a slab and not an atomic row. '
            '**b**, Class occupancy against depth. '
            '**c**, The partition: surface {{defect.flagship.umol.surface}}, '
            'subsurface {{defect.flagship.umol.subsurface}}, bulk '
            '{{defect.flagship.umol.bulk}} umol-O g-1 at a common vacancy '
            'chemical potential mu_V = {{defect.flagship.mu_V_eV}} eV. The '
            'shell holds {{defect.flagship.shell_pct}}% of the inventory at '
            '{{defect.flagship.shell_depletion_pct}}% local depletion, just '
            'below the estimated shear-plane ceiling.',
        'symbols': ['mu_V', 'theta_c', 'N_c'],
        'methods': ['site-exclusion isotherm', 'canonical Monte Carlo'],
    },
    {
        'id': 'd3', 'track': 'defect', 'label': 'Figure 3',
        'width': 'onehalf',
        'section': 'Defect thermodynamics - phase construction',
        'claim': 'The reconstruction is a competing phase handled by the '
                 'lever rule, and the high-coverage branch of the ideal '
                 'lattice gas is excluded deliberately, not overlooked.',
        'caption':
            '**Grand-potential construction of the surface phase set.** '
            'Surface grand potential per bridging site against the vacancy '
            'chemical potential. The (1x1) lattice gas (solid) crosses the '
            'added-row Ti2O3 line phase (straight, slope '
            '-{{defect.omega.theta_eff}}) at mu_t = '
            '{{defect.boundaries.mu_t_eV}} eV, which is the first-order '
            'transition; between the branches the surface is a two-phase '
            'mixture at pinned mu_V. Because the ideal branch keeps a '
            'steeper slope it re-crosses the line phase at '
            '{{defect.omega.mu_recross_eV}} eV (dashed): that high-coverage '
            'branch is excluded from the phase set on physical grounds - '
            'repulsive ordering penalises it and the observed strongly '
            'reduced surface is reconstructed, not bridging-stripped.',
        'symbols': ['mu_V', 'omega_s', 'theta_t'],
        'methods': ['grand-potential common tangent', 'lever rule'],
    },
    {
        'id': 'd4', 'track': 'defect', 'label': 'Figure 4',
        'width': 'double',
        'section': 'Defect kinetics - reoxidation transport',
        'claim': 'Coarse-grained kinetic Monte Carlo turns a measured '
                 'recovery fraction into an effective migration barrier, and '
                 'the neutral-vacancy barrier cannot reproduce the '
                 'measurement.',
        'caption':
            '**Coarse-grained reoxidation transport.** '
            '**a**, Recovered fraction of the inventory against exposure '
            'time at {{defect.cgmc.T_C}} C for the neutral-vacancy bulk '
            'barrier ({{defect.cgmc.E_m_eV}} eV) and for the barrier fitted '
            'to a measured recovery of {{defect.cgmc.target_pct}}% at '
            '{{defect.cgmc.target_s}} s. '
            '**b**, Radial vacancy occupancy at the same times. Transport at '
            'the neutral barrier is not rate-limiting, so a partial measured '
            'recovery reads back as an effective barrier of '
            '{{defect.cgmc.E_m_fit_eV}} eV - the signature of a channel the '
            'neutral-vacancy picture does not contain. Refill rate constants '
            'are placeholders; the barrier fit is the quantity this panel '
            'reports.',
        'symbols': ['E_m', 'k_fill', 'theta_k'],
        'methods': ['coarse-grained KMC', 'detailed balance',
                    'matrix exponential'],
    },
    # ----------------------------------------------- equilibrium track
    {
        'id': 'e1', 'track': 'equilibrium', 'label': 'Figure 5',
        'width': 'double',
        'section': 'Gas-solid equilibrium - stability of rutile under the '
                   'working feed',
        'claim': 'An equilibrated CO2/H2 feed does not reduce rutile at any '
                 'temperature in the working range, and the margin is a KKT '
                 'reduced cost rather than a read-off from a phase diagram.',
        'caption':
            '**Rutile stability under three gas charges.** '
            '**a**, Kuhn-Tucker reduced cost of the nearest reduced phase '
            '({{equilibrium.margin.phase}}) for the '
            '{{equilibrium.margin.feed_label}} feed. Positive is unstable, '
            'so rutile is the equilibrium solid throughout: the margin falls '
            'from {{equilibrium.margin.first_kJ}} to '
            '{{equilibrium.margin.last_kJ}} kJ per mol O between '
            '{{equilibrium.margin.first_T_C}} and '
            '{{equilibrium.margin.last_T_C}} C without changing sign. '
            '**b**, Reduced titanium fraction for three charges over the '
            'same range; the CO2/H2 feed returns exact zero at every '
            'temperature, pure hydrogen sits on the '
            '{{equilibrium.feeds.h2_phase}} buffer, and the sealed inert '
            'case is a genuine two-phase-buffer trace rather than an '
            'optimiser floor. Closed finite charge: '
            '{{equilibrium.conditions.mass_mg}} mg TiO2, '
            '{{equilibrium.conditions.volume_mL}} mL, '
            '{{equilibrium.conditions.pressure_atm}} atm.',
        'symbols': ['lambda_e', 'r_j', 'x_red'],
        'methods': ['active-set Gibbs minimisation', 'KKT reduced cost'],
    },
    {
        'id': 'e2', 'track': 'equilibrium', 'label': 'Figure 6',
        'width': 'onehalf',
        'section': 'Gas-solid equilibrium - method',
        'claim': 'The assemblage is selected by enumerating candidate active '
                 'sets and testing optimality, not by assuming which phases '
                 'are present.',
        'caption':
            '**Active-set selection of the condensed assemblage.** Every '
            'candidate subset of the {{equilibrium.method.n_solids}} '
            'implemented Ti-O phases is solved against the elemental '
            'balance; a candidate survives only if it is feasible and every '
            'excluded phase has a non-negative reduced cost. The smallest '
            'surviving set is the equilibrium assemblage, and excluded '
            'phases are returned as exact zero rather than as an optimiser '
            'residual. Shown for the {{equilibrium.margin.feed_label}} feed '
            'at {{equilibrium.method.T_C}} C.',
        'symbols': ['lambda_e', 'r_j', 'm_j'],
        'methods': ['active-set enumeration', 'KKT reduced cost',
                    'element-potential Newton'],
    },
    # ------------------------------------------------------- SI plates
    {
        'id': 's1', 'track': 'defect', 'label': 'Supplementary Figure 1',
        'width': 'onehalf',
        'section': 'Supplementary - declared shell thickness',
        'claim': 'The shell/bulk split is sensitive to the declared '
                 'subsurface thickness, and the conclusion that survives it '
                 'is the shell share, not the bulk remainder.',
        'caption':
            '**Sensitivity to the declared subsurface shell thickness.** '
            'Class inventory at {{defect.flagship.VO_total_umol_g}} umol-O '
            'g-1 and {{defect.flagship.T_C}} C against the number of oxygen '
            'layers assigned to the subsurface class. The surface holds at '
            'the line-phase deficiency throughout, so the ladder trades '
            'shell against bulk only: the subsurface class takes '
            '{{defect.shell_ladder_text.subsurface}} umol-O g-1 and the bulk '
            '{{defect.shell_ladder_text.bulk}} for one to four layers, while '
            'the shell share rises from {{defect.shell_ladder_text.first_pct}} '
            'to {{defect.shell_ladder_text.last_pct}}%. The default of one '
            'layer follows the formation-energy set, which resolves a single '
            'first-subsurface energy; real energies relax smoothly towards '
            'the bulk value, so the physical answer lies between the first '
            'two columns.',
        'symbols': ['n_lay', 'N_ss', 'theta_ss'],
        'methods': ['site-exclusion isotherm'],
    },
    {
        'id': 's2', 'track': 'defect', 'label': 'Supplementary Figure 2',
        'width': 'onehalf',
        'section': 'Supplementary - surface ordering',
        'claim': 'Surface arrangement statistics are a separate calculation '
                 'from the layer partition and are reported as such.',
        'caption':
            '**Surface vacancy arrangement.** Sampled bridging-row '
            'configuration and the distribution of vacancy separations along '
            'the row at the coverage of the surviving (1x1) patches. '
            'Effective pair interactions are constrained to reproduce '
            'published vacancy-correlation statistics; they parameterise '
            'arrangement, never layer energetics, and do not enter the '
            'partition of the preceding figures.',
        'symbols': ['V_1', 'theta_s'],
        'methods': ['canonical swap Metropolis', 'Widom insertion'],
    },
    {
        'id': 's3', 'track': 'defect', 'label': 'Supplementary Figure 3',
        'width': 'single',
        'section': 'Supplementary - particle geometry',
        'claim': 'The measured BET area and the sphere model agree on the '
                 'ladder, so the geometric input is not carrying the result.',
        'caption':
            '**Geometry sensitivity.** Phase-boundary ladder computed from '
            'the measured BET area ({{defect.bet.area_m2_g}} m2 g-1) against '
            'the equivalent-sphere model ({{defect.geometry.d_um}} um). The '
            'two agree to {{defect.bet.max_shift_pct}}% on every boundary.',
        'symbols': ['A_BET', 'd', 'N_c'],
        'methods': ['equivalent-sphere geometry'],
    },
    {
        'id': 's4', 'track': 'defect', 'label': 'Supplementary Figure 4',
        'width': 'onehalf',
        'section': 'Supplementary - thermodynamic bridge',
        'claim': 'The defect calculation is conditional on rutile being the '
                 'stable solid, and the condition is evaluated rather than '
                 'assumed.',
        'caption':
            '**Validity of the conditional-equilibrium partition.** The gas '
            'charge is re-solved at the defect temperature and the '
            'assemblage classified: rutile alone, rutile with a reduced '
            'phase, or an assemblage without rutile, in which case the '
            'partition is outside its own scope. The partition is a '
            'conditional equilibrium at fixed inventory and carries no '
            'statement about how that inventory was created.',
        'symbols': ['r_j', 'mu_O'],
        'methods': ['active-set Gibbs minimisation'],
    },
    {
        'id': 's5', 'track': 'verification', 'label': 'Supplementary Figure 5',
        'width': 'double',
        'section': 'Supplementary - numerical verification',
        'claim': 'Production values are reproduced by independent '
                 'high-precision references and by a second implementation.',
        'caption':
            '**Verification.** **a**, Production engine against the '
            '{{verification.statmech_dps}}-digit defect oracle at the '
            'shipped configuration and across the phase ladder. '
            '**b**, JavaScript engine against the Python engine, same seed, '
            'compared as integer states where the trajectory is discrete. '
            '**c**, Coarse-grained transport against a high-precision matrix '
            'exponential of the same generator. '
            '**d**, Worst deviation measured in each layer, against the '
            'tolerance the test suite enforces. The suite is '
            '{{verification.n_tests}} tests with none skipped.',
        'symbols': ['mu_V', 'theta_c', 'p_k'],
        'methods': ['mpmath reference', 'matrix exponential',
                    'same-seed trajectory comparison'],
    },
    {
        'id': 's6', 'track': 'equilibrium',
        'label': 'Supplementary Figure 6',
        'width': 'double',
        'section': 'Supplementary - equilibrium verification',
        'claim': 'The active-set solution is the constrained minimum, and '
                 'trace amounts are resolved rather than truncated.',
        'caption':
            '**Equilibrium verification.** **a**, Total Gibbs energy of '
            'every candidate assemblage relative to the winner, showing the '
            'reported solution is the minimum and not a local stationary '
            'point. **b**, Trace solid amounts under the sealed inert charge '
            'compared logarithmically against the '
            '{{verification.thermo_dps}}-digit reference, down to '
            '{{verification.trace_floor_log10}} in log10 mol. '
            '**c**, Elemental balance residual per element. '
            '**d**, Deviation ledger against the enforced tolerance.',
        'symbols': ['G', 'lambda_e', 'm_j'],
        'methods': ['mpmath reference', 'active-set enumeration'],
    },
]

FIG_BY_ID = {f['id']: f for f in FIGURES}

# Symbols and numerical methods named in the plates, defined once so a plate
# can be read without the paper and without the code.
SYMBOLS = {
    'N_s': ('N_s', 'Bridging-oxygen site capacity of the particle, '
                   'umol-O g-1.'),
    'N_ss': ('N_ss', 'Subsurface site capacity: four oxygens per (1x1) cell '
                     'per declared layer, umol-O g-1.'),
    'N_b': ('N_b', 'Bulk oxygen site capacity, umol-O g-1.'),
    'N_c': ('N_c', 'Site capacity of class c.'),
    'n_lay': ('n_lay', 'Declared subsurface thickness in oxygen layers.'),
    'mu_V': ('mu_V', 'Vacancy chemical potential, common to all classes at '
                     'the conditional equilibrium, eV.'),
    'omega_s': ('omega_s', 'Surface grand potential per bridging site, eV.'),
    'theta_c': ('theta_c', 'Fractional vacancy occupancy of class c.'),
    'theta_ss': ('theta_ss', 'Occupancy of the subsurface class.'),
    'theta_s': ('theta_s', 'Occupancy of the bridging-oxygen row.'),
    'theta_t': ('theta_t', 'Bridging-row coverage at the reconstruction '
                           'onset.'),
    'E_m': ('E_m', 'Vacancy migration barrier used by the transport layer, '
                   'eV.'),
    'k_fill': ('k_fill', 'First-order CO2 refill rate constant of a site '
                         'class, s-1.'),
    'theta_k': ('theta_k', 'Vacancy occupancy of radial coarse cell k.'),
    'p_k': ('p_k', 'Occupancy vector of the coarse-grained master equation.'),
    'V_1': ('V_1, V_2, V_3', 'Effective pair repulsion along the bridging '
                             'row at one, two and three site separations, '
                             'eV.'),
    'lambda_e': ('lambda_e', 'Element potential of element e, kJ mol-1.'),
    'r_j': ('r_j', 'Kuhn-Tucker reduced cost of excluded condensed phase j; '
                   'positive means unstable.'),
    'm_j': ('m_j', 'Moles of condensed phase j.'),
    'x_red': ('x_red', 'Fraction of titanium in reduced phases.'),
    'mu_O': ('mu_O', 'Oxygen chemical potential set by the gas charge.'),
    'A_BET': ('A_BET', 'Measured specific surface area, m2 g-1.'),
    'd': ('d', 'Equivalent-sphere particle diameter, um.'),
    'G': ('G', 'Total Gibbs energy of the charge, kJ.'),
}

METHODS = {
    'site-exclusion isotherm':
        'Exact occupancy of a class of non-interacting exclusive sites at a '
        'given chemical potential.',
    'lever rule':
        'Partition of a two-phase coexistence at a pinned chemical '
        'potential.',
    'closed-form phase boundaries':
        'Analytic inventory at which each phase transition begins and ends.',
    'grand-potential common tangent':
        'First-order transition located where two branches of the grand '
        'potential cross.',
    'canonical Monte Carlo':
        'Fixed-population Metropolis sampling of the interacting surface '
        'lattice.',
    'canonical swap Metropolis':
        'Metropolis moves that exchange an occupied and an empty site, '
        'conserving the vacancy count.',
    'Widom insertion':
        'Chemical potential from the mean Boltzmann factor of a test '
        'insertion.',
    'coarse-grained KMC':
        'Kinetic Monte Carlo on coarse cells with exchange rates in detailed '
        'balance with the cell Hamiltonian.',
    'detailed balance':
        'Rate pairs constrained so the closed system relaxes to the exact '
        'equilibrium distribution.',
    'matrix exponential':
        'High-precision propagation of the linear master equation, used as '
        'an independent reference.',
    'equivalent-sphere geometry':
        'Site counts from a specific surface area, expressed as a sphere of '
        'equal area per gram.',
    'active-set Gibbs minimisation':
        'Gibbs energy minimised over the full species set at fixed elemental '
        'inventory, with the condensed assemblage selected by enumeration.',
    'active-set enumeration':
        'Every subset of condensed phases solved and tested, the smallest '
        'optimal set retained.',
    'KKT reduced cost':
        'Optimality test for an excluded phase: its formation energy '
        'measured against the element potentials.',
    'element-potential Newton':
        'Newton solution for the element potentials that satisfy the '
        'elemental balance.',
    'mpmath reference':
        'Independent arbitrary-precision implementation sharing only the '
        'input data.',
    'same-seed trajectory comparison':
        'Two implementations driven by one random stream and compared as '
        'integer states.',
}
