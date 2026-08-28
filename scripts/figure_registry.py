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
    # The house drawing style, in points. Every panel is a closed
    # rectangle rather than an open L, ticks point into the panel, and
    # tick labels appear on the bottom and left only. No grid, white
    # ground, one stroke weight for data and one for everything else.
    'axis_pt': 0.9,
    'tick_major_pt': 4.0,
    'tick_minor_pt': 2.0,
    'curve_pt': 2.2,
    'marker_pt': 5.0,
    'axis_box': 'four spines, ticks inward, labels bottom and left',
    'grid': 'none',
    'export': 'print SVG at the declared width, plus 600 dpi PNG',
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
    ('bulk', '#6A51A3',
     'Bulk oxygen sites and the bulk-class occupancy.'),
    ('extended', '#777777',
     'Extended defects: the crystallographic-shear precipitate beyond the '
     'solubility ceiling. The one neutral role, so structure is lighter '
     'than it is.'),
    ('gas', '#009E73',
     'The gas charge and everything it carries: feed composition, oxygen '
     'potential, CO2 refill.'),
    ('structure', '#AAAAAA',
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
        'claim': 'Particle geometry decides where the inventory goes. The '
                 'surface classes are too small to hold it.',
        'caption':
            '**Site capacity and the phase ladder.** '
            '**a**, Oxygen-site capacity of the three classes for a '
            '{{defect.geometry.d_um}} um particle, '
            '{{defect.geometry.area_m2_g}} m2 g-1. The bridging-oxygen row '
            'holds {{defect.geometry.N_s}} umol-O g-1. That is '
            '{{defect.geometry.surface_pct_of_O}}% of all oxygen sites. The '
            'subsurface shell holds {{defect.geometry.N_ss}}. Both fill '
            'before the measured inventory is placed. '
            '**b**, Class inventory against total vacancy content at '
            '{{defect.flagship.T_C}} C. The (1x2) reconstruction starts at '
            '{{defect.boundaries.VO_onset_umol_g}} umol-O g-1 and is '
            'complete at '
            '{{defect.boundaries.VO_recon_complete_umol_g}}. Above the '
            'shear-plane ceiling at {{defect.boundaries.VO_cs_umol_g}} the '
            'vacancy chemical potential is pinned and the excess forms '
            'extended defects. That ceiling is an estimate. Dashed lines '
            'mark the boundaries. The marker is the measured inventory, '
            '{{defect.flagship.VO_total_umol_g}} umol-O g-1.',
        'symbols': ['N_s', 'N_ss', 'N_b', 'mu_V', 'theta_c'],
        'methods': ['site-exclusion isotherm', 'lever rule',
                    'closed-form phase boundaries'],
    },
    {
        'id': 'd2', 'track': 'defect', 'label': 'Figure 2',
        'width': 'double',
        'section': 'Defect thermodynamics - the measured inventory',
        'claim': 'At the measured inventory the surface is reconstructed and '
                 'full. The near-surface shell carries most of the '
                 'vacancies.',
        'caption':
            '**The measured inventory of '
            '{{defect.flagship.VO_total_umol_g}} umol-O g-1.** '
            '**a**, Cross-section of the outer shell, drawn to scale. One '
            'dot is one vacancy on one oxygen site. The bridging plane '
            'carries one oxygen per (1x1) column. Each '
            '{{defect.geometry.layer_nm}} nm layer below it carries four. '
            'The subsurface class is a slab, not an atomic row. '
            '**b**, Class occupancy against depth. '
            '**c**, The partition. Surface '
            '{{defect.flagship.umol.surface}}, subsurface '
            '{{defect.flagship.umol.subsurface}}, bulk '
            '{{defect.flagship.umol.bulk}} umol-O g-1, at a common vacancy '
            'chemical potential of {{defect.flagship.mu_V_eV}} eV. The '
            'shell holds {{defect.flagship.shell_pct}}% of the inventory at '
            '{{defect.flagship.shell_depletion_pct}}% local depletion. This '
            'sits just below the shear-plane ceiling.',
        'symbols': ['mu_V', 'theta_c', 'N_c'],
        'methods': ['site-exclusion isotherm', 'canonical Monte Carlo'],
    },
    {
        'id': 'd3', 'track': 'defect', 'label': 'Figure 3',
        'width': 'onehalf',
        'section': 'Defect thermodynamics - phase construction',
        'claim': 'The reconstruction is a competing phase, handled by the '
                 'lever rule. The high-coverage branch is excluded on '
                 'purpose.',
        'caption':
            '**The surface phase set.** Grand potential per bridging site '
            'against vacancy chemical potential. The (1x1) lattice gas '
            'crosses the added-row Ti2O3 line phase at '
            '{{defect.boundaries.mu_t_eV}} eV. That crossing is the '
            'first-order transition. The line phase has slope '
            '-{{defect.omega.theta_eff}}. Between the branches the surface '
            'is a two-phase mixture at fixed chemical potential. The ideal '
            'branch is steeper, so it crosses back at '
            '{{defect.omega.mu_recross_eV}} eV. That high-coverage branch is '
            'excluded from the phase set. Repulsive ordering penalises it, '
            'and the observed reduced surface is reconstructed rather than '
            'stripped of bridging oxygen.',
        'symbols': ['mu_V', 'omega_s', 'theta_t'],
        'methods': ['grand-potential common tangent', 'lever rule'],
    },
    {
        'id': 'd4', 'track': 'defect', 'label': 'Figure 4',
        'width': 'double',
        'section': 'Defect kinetics - reoxidation transport',
        'claim': 'A measured recovery fraction reads back as an effective '
                 'migration barrier. The neutral-vacancy barrier does not '
                 'reproduce the measurement.',
        'caption':
            '**Reoxidation transport.** '
            '**a**, Recovered fraction of the inventory against exposure '
            'time at {{defect.cgmc.T_C}} C. One curve uses the '
            'neutral-vacancy bulk barrier of {{defect.cgmc.E_m_eV}} eV. The '
            'other uses the barrier fitted to a measured recovery of '
            '{{defect.cgmc.target_pct}}% at {{defect.cgmc.target_s}} s. '
            '**b**, Recovered fraction at {{defect.cgmc.target_s}} s '
            'against the migration barrier. Transport at the neutral '
            'barrier is fast enough that it does not limit the rate, so the '
            'curve saturates there. Reading the measured recovery across '
            'gives an effective barrier of {{defect.cgmc.E_m_fit_eV}} eV. '
            'Refill rate constants are placeholders. The barrier is the '
            'result this panel reports.',
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
                 'temperature in the working range.',
        'caption':
            '**Rutile stability under three gas charges.** '
            '**a**, Kuhn-Tucker reduced cost of the nearest reduced phase, '
            '{{equilibrium.margin.phase}}, for the '
            '{{equilibrium.margin.feed_label}} feed. Positive means the '
            'phase is unstable, so rutile is the equilibrium solid. The '
            'margin falls from {{equilibrium.margin.first_kJ}} to '
            '{{equilibrium.margin.last_kJ}} kJ per mol O between '
            '{{equilibrium.margin.first_T_C}} and '
            '{{equilibrium.margin.last_T_C}} C and does not change sign. '
            '**b**, Reduced titanium fraction for three charges over the '
            'same range. The CO2/H2 feed returns exact zero at every '
            'temperature. Pure hydrogen sits on the '
            '{{equilibrium.feeds.h2_phase}} buffer. The sealed inert case is '
            'a real two-phase-buffer trace, not an optimiser floor. The '
            'charge is closed: {{equilibrium.conditions.mass_mg}} mg TiO2 in '
            '{{equilibrium.conditions.volume_mL}} mL at '
            '{{equilibrium.conditions.pressure_atm}} atm.',
        'symbols': ['lambda_e', 'r_j', 'x_red'],
        'methods': ['active-set Gibbs minimisation', 'KKT reduced cost'],
    },
    {
        'id': 'e2', 'track': 'equilibrium', 'label': 'Figure 6',
        'width': 'onehalf',
        'section': 'Gas-solid equilibrium - method',
        'claim': 'The assemblage is selected by testing every candidate, not '
                 'by assuming which phases are present.',
        'caption':
            '**Active-set selection.** Every subset of the '
            '{{equilibrium.method.n_solids}} implemented Ti-O phases is '
            'solved against the elemental balance. A subset survives if it '
            'is feasible and every excluded phase has a non-negative reduced '
            'cost. The smallest surviving subset is the equilibrium '
            'assemblage. The bars are that optimality test for the selected '
            'assemblage: the reduced cost of each excluded phase, all '
            'positive. Excluded phases are returned as exact zero, not as an '
            'optimiser residual. Shown for the '
            '{{equilibrium.margin.feed_label}} feed at '
            '{{equilibrium.method.T_C}} C.',
        'symbols': ['lambda_e', 'r_j', 'm_j'],
        'methods': ['active-set enumeration', 'KKT reduced cost',
                    'element-potential Newton'],
    },
    # ------------------------------------------------------- SI plates
    {
        'id': 's1', 'track': 'defect', 'label': 'Supplementary Figure 1',
        'width': 'onehalf',
        'section': 'Supplementary - declared shell thickness',
        'claim': 'The shell and bulk split depends on the declared '
                 'subsurface thickness. The shell share does not.',
        'caption':
            '**Declared subsurface shell thickness.** Class inventory at '
            '{{defect.flagship.VO_total_umol_g}} umol-O g-1 and '
            '{{defect.flagship.T_C}} C against the number of oxygen layers '
            'assigned to the subsurface class. The surface stays at the '
            'line-phase deficiency, so the ladder moves inventory between '
            'shell and bulk only. The subsurface class takes '
            '{{defect.shell_ladder_text.subsurface}} umol-O g-1 for one to '
            'four layers. The bulk takes '
            '{{defect.shell_ladder_text.bulk}}. The shell share rises from '
            '{{defect.shell_ladder_text.first_pct}} to '
            '{{defect.shell_ladder_text.last_pct}}%. One layer is the '
            'default because the formation-energy set resolves a single '
            'first-subsurface energy. Real energies relax towards the bulk '
            'value over the first few layers, so the physical answer lies '
            'between the first two columns.',
        'symbols': ['n_lay', 'N_ss', 'theta_ss'],
        'methods': ['site-exclusion isotherm'],
    },
    {
        'id': 's2', 'track': 'defect', 'label': 'Supplementary Figure 2',
        'width': 'onehalf',
        'section': 'Supplementary - surface ordering',
        'claim': 'Surface arrangement is a separate calculation from the '
                 'layer partition.',
        'caption':
            '**Surface vacancy arrangement.** Sampled bridging-row '
            'configuration and the distribution of vacancy separations along '
            'the row, taken at the coverage of the surviving (1x1) patches. '
            'Effective pair interactions are constrained to published '
            'vacancy-correlation statistics. They set arrangement only. They '
            'do not enter the layer partition of the preceding figures.',
        'symbols': ['V_1', 'theta_s'],
        'methods': ['canonical swap Metropolis', 'Widom insertion'],
    },
    {
        'id': 's3', 'track': 'defect', 'label': 'Supplementary Figure 3',
        'width': 'single',
        'section': 'Supplementary - particle geometry',
        'claim': 'The measured BET area and the sphere model give the same '
                 'ladder.',
        'caption':
            '**Geometry sensitivity.** Phase-boundary ladder from the '
            'measured BET area of {{defect.bet.area_m2_g}} m2 g-1 against '
            'the equivalent-sphere model at {{defect.geometry.d_um}} um. The '
            'two agree to {{defect.bet.max_shift_pct}}% on every boundary.',
        'symbols': ['A_BET', 'd', 'N_c'],
        'methods': ['equivalent-sphere geometry'],
    },
    {
        'id': 's4', 'track': 'defect', 'label': 'Supplementary Figure 4',
        'width': 'onehalf',
        'section': 'Supplementary - thermodynamic bridge',
        'claim': 'The defect calculation holds only where rutile is the '
                 'stable solid. That condition is evaluated.',
        'caption':
            '**Validity of the partition.** The gas charge is re-solved at '
            'the defect temperature and the assemblage is classified: rutile '
            'alone, rutile with a reduced phase, or an assemblage without '
            'rutile. The third case is outside the scope of the partition. '
            'The partition is a conditional equilibrium at fixed inventory. '
            'It says nothing about how the inventory was created.',
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
            '**Verification.** **a**, Deviation of the production engine '
            'from the {{verification.statmech_dps}}-digit defect reference '
            'in the vacancy chemical potential, over the whole matching '
            'grid. Each point is one temperature and one inventory. The line '
            'is the tolerance the test suite enforces. '
            '**b**, Worst deviation measured in each verification layer, '
            'against the tolerance for that layer. The layers are the '
            'analytic partition, the declared shell ladder, a finite lattice '
            'against its exact ensemble, and coarse-grained transport '
            'against a high-precision matrix exponential of the same '
            'generator. The suite holds {{verification.n_tests}} tests and '
            'skips none.',
        'symbols': ['mu_V', 'theta_c', 'p_k'],
        'methods': ['mpmath reference', 'matrix exponential',
                    'same-seed trajectory comparison'],
    },
    {
        'id': 's6', 'track': 'equilibrium',
        'label': 'Supplementary Figure 6',
        'width': 'double',
        'section': 'Supplementary - equilibrium verification',
        'claim': 'The active-set solution is the constrained minimum. Trace '
                 'amounts are resolved, not truncated.',
        'caption':
            '**Equilibrium verification.** **a**, Trace amounts under the '
            'sealed inert charge, resolved logarithmically against the '
            '{{verification.thermo_dps}}-digit reference down to '
            '{{verification.trace_floor_log10}} in log10 mol. A minimiser '
            'that truncated instead of resolving would floor these points. '
            '**b**, Worst elemental balance residual per temperature for the '
            'same charge. Closure is at the working precision of the '
            'reference, not at the tolerance of the production solver.',
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
