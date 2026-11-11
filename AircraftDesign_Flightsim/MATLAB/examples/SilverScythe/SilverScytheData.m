function ss = SilverScytheData()
% SILVERSCYTHEDATA  Authoritative input deck for the Silver Scythe AAMF case.
%
% Silver Scythe is the carrier strike-fighter concept of the Virginia Tech
% Gobbler Ghost Program (AIAA final report, 2026).  Every number in this
% file is traceable to one of three sources, tagged in the comments:
%
%   [RPT]  the AIAA final report (table/section noted)
%   [SD]   the senior-design MATLAB deck (aircraft_params('SilverScythe'),
%          Params.m, calculateWeightBalance.m) that produced the report
%   [SIM]  values required only by the 6-DOF simulator and not part of the
%          report's steady-state analysis (frame layout, box inertias,
%          rate-limit and ground-contact settings)
%
% Nothing in this file is tuned to make an answer come out right.  The
% published drag, speed and turn-rate numbers in ss.targets are used for
% comparison only, by SilverScytheValidation.
%
% Unit convention: SI unless the field name ends in _ft, _ft2, _lb, _lbf,
% _deg, _kt, _nm, _fpm or _per_deg.

%% ------------------------------------------------------------------
%  Unit conversions
%% ------------------------------------------------------------------
ss.units.ft_to_m    = 0.3048;
ss.units.ft2_to_m2  = 0.09290304;
ss.units.lbm_to_kg  = 0.45359237;
ss.units.lbf_to_N   = 4.4482216152605;
ss.units.kt_to_mps  = 0.514444;
ss.units.fpm_to_mps = 0.00508;
ss.units.slug_ft2_to_kg_m2 = 1.3558179619;

ss.name    = "Silver Scythe";
ss.program = "Gobbler Ghost Program, Virginia Tech, 2026";

%% ------------------------------------------------------------------
%  Planform and reference geometry                       [RPT Table 1, Fig 18]
%% ------------------------------------------------------------------
ss.length_ft      = 50.0;
ss.height_ft      = 12.0;
ss.span_ft        = 42.0;
ss.folded_span_ft = 35.0;
ss.S_ref_ft2      = 798.0;   % wing + LEX
ss.root_chord_ft  = 32.0;    % [SD] aircraft_params.cr
ss.tip_chord_ft   = 6.0;     % [SD] aircraft_params.ct
ss.LE_sweep_deg   = 38.0;
ss.TE_sweep_deg   = -24.55;  % [SD] forward sweep carried as negative
ss.tc_wing        = 0.06;

ss.taper_ratio = ss.tip_chord_ft / ss.root_chord_ft;
ss.AR          = ss.span_ft^2 / ss.S_ref_ft2;
ss.cbar_ft     = (2/3) * ss.root_chord_ft * ...
                 (1 + ss.taper_ratio + ss.taper_ratio^2) / (1 + ss.taper_ratio);
ss.mean_chord_ft = ss.S_ref_ft2 / ss.span_ft;   % c_a used by LIDRAG spanloads

% Mid-chord sweep, as used by the senior-design stall_solver.m
ss.midchord_sweep_deg = atand(tand(ss.LE_sweep_deg) - ...
    (2/ss.AR)*((1-ss.taper_ratio)/(1+ss.taper_ratio)));

ss.S_ref = ss.S_ref_ft2 * ss.units.ft2_to_m2;
ss.b_ref = ss.span_ft   * ss.units.ft_to_m;
ss.cbar  = ss.cbar_ft   * ss.units.ft_to_m;

%% ------------------------------------------------------------------
%  Airfoil                                          [RPT 11.2] / [SD] cl_alpha
%% ------------------------------------------------------------------
ss.airfoil.name              = "NACA 64-206";
ss.airfoil.polar_file        = "NACA_64206_xfoil.txt";  % XFOIL, Re 12e6, M 0
ss.airfoil.alpha_L0_report_deg = -1.05;   % [RPT 11.2]
ss.airfoil.alpha_max_deg     = 21.0;      % [SD] stall_solver max_aoa, clean
ss.airfoil.alpha_max_land_deg = 8.1;      % [SD] takeoff_landing landing case
ss.airfoil.nominal_cl_alpha_per_deg = 0.104954;  % [SD] aircraft_params
ss.airfoil.nominal_alpha_L0_deg     = -1.0469;   % [SD] aircraft_params

%% ------------------------------------------------------------------
%  FRICTION component breakdown                     [RPT Table 4] / [SD] SUMAC
%
%  Report Table 4 lumps the aft lifting surface into the fuselage
%  (854.610 = 815.462 + 2*19.574).  The senior-design deck keeps them apart
%  because FRICTION needs a reference length and t/c per component.  Both
%  sum to the report's 2380.520 ft^2.
%% ------------------------------------------------------------------
ss.friction.component_names = ["aft_lifting_surface";"fuselage";"LEX"; ...
                               "canopy";"ruddervators";"wings"];
ss.friction.S_wet_ft2 = [2*19.574; 815.462; 2*26.957; 47.938; ...
                         2*153.566; 2*558.463];
ss.friction.l_ref_ft  = [3.50; ss.length_ft; 11.70; 10.43; 9.65; 21.41];
ss.friction.geom_type = [0; 1; 0; 1; 0; 0];   % 0 = plate/lifting, 1 = body
ss.friction.Q         = [1.05; 1.00; 1.05; 1.00; 1.03; 1.05];  % [RPT Table 5]
ss.friction.transition = 0.075 * ones(6,1);   % laminar run fraction
ss.friction.TwTaw      = 1.0;                 % adiabatic wall
ss.friction.leakage_multiplier = 1.05;        % [RPT 11.4] stealth, flat 5%

ss.max_cross_section_ft2   = 28.591;   % [SD] aircraft_params.Amax
ss.wave_reference_length_ft = ss.length_ft;
ss.S_wet_ft2 = sum(ss.friction.S_wet_ft2);

% Fuselage fineness enters FRICTION as the body "t/c" slot, and the canopy
% t/c is the senior-design equivalent-diameter expression.
ss.friction.t_c = [0.10; ...
                   ss.max_cross_section_ft2/ss.length_ft; ...
                   0.06; ...
                   (3.28/2)^2*pi/10.43; ...
                   0.06; ...
                   ss.tc_wing];

ss.friction.stealth_type = "Stealth";

%% ------------------------------------------------------------------
%  Wave drag                                            [RPT 11.5, Eq 15-16]
%% ------------------------------------------------------------------
ss.wave.EWD          = 1.8;   % [SD] aircraft_params.EWD (Raymer 1.8-2.2)
ss.wave.EWD_min      = 1.8;
ss.wave.EWD_max      = 2.2;

%% ------------------------------------------------------------------
%  Weights                                       [RPT Tables 12, 13, 14]
%% ------------------------------------------------------------------
ss.weights.empty_lb            = 33970;
ss.weights.useful_strike_lb    = 27048;
ss.weights.TOGW_strike_lb      = 61018;
ss.weights.max_landing_lb      = 41900;   % Table 14
ss.weights.landing_analysis_lb = 41400;   % [RPT 16.3] value used for Vapp
ss.weights.usable_fuel_lb      = 21000;
ss.weights.strike_weapons_lb   = 4428;
ss.weights.combat_lb           = 44418.6019837056;  % [SD] mission A2A combat
ss.weights.cruise_lb           = 48276.3126782551;  % [SD] mission cruise1
ss.weights.empty_cg_station_ft = 31.31;
ss.weights.takeoff_cg_station_ft = 30.23;

ss.takeoff_mass_kg = ss.weights.TOGW_strike_lb * ss.units.lbm_to_kg;
ss.empty_mass_kg   = ss.weights.empty_lb       * ss.units.lbm_to_kg;

% Weight at which the senior-design aerodynamic deck was generated
% ([SD] PerformanceCall: W_0 = mission.strike.start).  Used only when
% reproducing report Figure 24.
ss.weights.aero_deck_lb = ss.weights.TOGW_strike_lb;

%% ------------------------------------------------------------------
%  Propulsion                                    [RPT Tables 10, 11]
%% ------------------------------------------------------------------
ss.engine.name              = "F414-GE-400";
ss.engine.count             = 2;
ss.engine.weight_lb         = 2445;
ss.engine.T_SL_dry_lbf      = 13000;
ss.engine.T_SL_AB_lbf       = 22000;
ss.engine.TSFC_SL_dry       = 0.84;
ss.engine.TSFC_SL_AB        = 1.85;
ss.engine.bypass_ratio      = 0.25;
ss.engine.mass_flow_lbm_s   = 170;
ss.engine.npss_deck_file    = "F414_NPSS_deck.csv";
% NPSS design point [RPT Table 11]
ss.engine.npss.mass_flow_lbm_s = 168.8;
ss.engine.npss.BPR             = 0.26995;
ss.engine.npss.T_dry_lbf       = 13700;
ss.engine.npss.TSFC_dry        = 0.8552;
ss.engine.npss.T_wet_lbf       = 22000;
ss.engine.npss.TSFC_wet        = 1.8439;
% Mattingly deck constants [RPT 13.5]
ss.engine.mattingly.TR         = 1.0;
ss.engine.mattingly.gamma      = 1.4;

%% ------------------------------------------------------------------
%  Empennage and control surfaces               [RPT Tables 6, 7]
%% ------------------------------------------------------------------
ss.empennage.airfoil       = "NACA 0006";
ss.empennage.span_ft       = 13.1;
ss.empennage.area_ft2      = 158.386;
ss.empennage.LE_sweep_deg  = 38.0;
ss.empennage.TE_sweep_deg  = -24.6;
ss.empennage.dihedral_deg  = 34.5;

ss.controls.ruddervator = struct('area_ft2',125.8,'station_ft',41.9, ...
    'y_ft',12.4,'max_deg',25);
ss.controls.nozzle = struct('area_ft2',25.6,'station_ft',47.9, ...
    'y_ft',2.8,'max_deg',20,'optimum_turn_deg',11);
ss.controls.aileron = struct('area_ft2',22.1,'station_ft',39.8, ...
    'y_ft',14.5,'max_deg',25);
ss.controls.slats = struct('area_ft2',42.5,'station_ft',20.4, ...
    'y_ft',11.01,'max_deg',25);
ss.controls.flaps = struct('area_ft2',22.1,'station_ft',42.1, ...
    'y_ft',9.68,'max_deg',30);

% High-lift multipliers on the clean 3-D CLmax, as used by the
% senior-design analyses.
ss.high_lift.clean_multiplier   = 1.0;
ss.high_lift.landing_multiplier = 2.0;   % [SD] takeoff_landing: CLmax*2
ss.high_lift.takeoff_multiplier = 1.6;   % [SD] takeoff_landing: CLmax*1.6
ss.high_lift.LEX_vortex_multiplier = 1.4;% [SD] turn_performance: CLmax*1.4

%% ------------------------------------------------------------------
%  Static stability                              [RPT Tables 8, 9]
%% ------------------------------------------------------------------
ss.static.cg_station_ft            = 27.9;
ss.static.neutral_point_station_ft = 22.7;
ss.static.static_margin_MAC        = -0.23;  % negative: relaxed stability
ss.static.g_limit                  = 8.0;    % [RPT 16.6] symmetric limit

% Moment reference station of each published derivative source.  Report
% Table 8 quotes its static margin about station 27.9 ft, and the DATCOM
% model was written with XCG = 28.3053 ft ([SD] call_write_datcom.m).  The
% coefficient lookup is referenced to the matching station for whichever
% source is selected, and the framework transfers the moment to the
% instantaneous mass CG.
ss.static.textbook_reference_station_ft = 27.9;
ss.static.datcom_reference_station_ft   = 28.3053;

% Report Table 9 DATCOM derivatives, stored verbatim for reporting only.
% The report does not state the derivative units, and the two possible
% readings (per degree / per radian) differ by 57.3x, so these values are
% NOT used to drive the simulation.  SilverScytheStabDerivatives builds the
% flight-model derivatives from geometry and prints this table alongside
% them for comparison.
ss.datcom_report.CL_alpha  = 0.003;
ss.datcom_report.Cm_alpha  = 0.001;
ss.datcom_report.Cm_deltae = -0.062;
ss.datcom_report.Cn_beta   = 0.007;
ss.datcom_report.Cl_beta   = -0.006;
ss.datcom_report.Cn_p      = 0.0004;
ss.datcom_report.Cl_p      = -0.003;
ss.datcom_report.CY_p      = 0.0001;
ss.datcom_report.Cn_r      = 0.0002;
ss.datcom_report.Cl_r      = 0.0008;
ss.datcom_report.Cm_q      = -0.0009;

% Lateral control derivative set actually exercised by the senior-design
% engine-out trim solution ([SD] takeoff_landing.m), per radian.
ss.sd_control.CY_beta  = -3.909e-3;
ss.sd_control.Cl_beta  = -0.006;
ss.sd_control.Cn_beta  =  0.007;
ss.sd_control.CY_da    =  0.00;
ss.sd_control.Cl_da    = -0.12;
ss.sd_control.Cn_da    = -0.03;
ss.sd_control.CY_dr    =  0.30;
ss.sd_control.Cl_dr    =  0.0008;
ss.sd_control.Cn_dr    =  0.0002;

%% ------------------------------------------------------------------
%  Carrier suitability constants                        [RPT 16.2, Eq 17]
%% ------------------------------------------------------------------
ss.carrier.catapult_stroke_ft   = 302;    % C13-2 power stroke
ss.carrier.catapult_equiv_wt_lb = 6680;
ss.carrier.deadload_weight_lb   = [40000 50000 62400];
ss.carrier.deadload_endspeed_kt = [150 142 134];
ss.carrier.tropical_day_K       = 305.2611;   % 89.8 degF
ss.carrier.launch_alpha_deg     = 7;
ss.carrier.Vapp_over_Vs         = 1.10;
ss.carrier.Varrest_over_Vapp    = 1.05;
ss.carrier.nozzle_launch_deg    = 12;         % [SD] lift credit at launch
ss.carrier.thrust_lift_per_engine_lbf = 9746; % [SD] takeoff_landing
ss.carrier.min_SEROC_launch_fpm   = 200;
ss.carrier.min_SEROC_approach_fpm = 500;

%% ------------------------------------------------------------------
%  Published validation targets                  [RPT Table 1, 16.x]
%% ------------------------------------------------------------------
ss.targets.max_mach_30000_AB      = 1.64;
ss.targets.max_mach_SL_mil        = 0.93;   % [RPT 16.5]; Table 1 lists 0.94
ss.targets.sustained_turn_20000_deg_s     = 12.0;
ss.targets.instantaneous_turn_20000_deg_s = 22.0;
ss.targets.corner_speed_20000_kt  = 400;
ss.targets.catapult_end_speed_kt  = 133.5;
ss.targets.approach_speed_kt      = 133.8;
ss.targets.arresting_speed_kt     = 140.0;
ss.targets.SEROC_approach_fpm     = 3320;
ss.targets.SEROC_launch_fpm       = 2308;
ss.targets.combat_radius_strike_nm = 706;
ss.targets.combat_radius_A2A_nm    = 740;
ss.targets.best_cruise_mach       = 0.9;
ss.targets.best_loiter_mach       = 0.8;
ss.targets.optimum_nozzle_deg     = 11.0;

% Report Figure 24 drag components at 30,000 ft, digitised from the figure.
ss.targets.drag.altitude_ft = 30000;
ss.targets.drag.Mach        = [0.80, 1.64];
ss.targets.drag.CD_profile  = [0.0117084819216830, 0.0092506209444690];
ss.targets.drag.CD_induced  = [0.0185671200638900, 0.0043985854955960];
ss.targets.drag.CD_wave     = [0.000511131519113305, 0.0161335031431980];
ss.targets.drag.CD_total    = ss.targets.drag.CD_profile + ...
                              ss.targets.drag.CD_induced + ...
                              ss.targets.drag.CD_wave;

% LIDRAG span efficiency from the team's archived run
% (reference/lidrag_silverscythe.out).
ss.targets.lidrag_e_span = 0.89127;

%% ------------------------------------------------------------------
%  Low-fidelity (textbook) inputs
%
%  These reproduce the Raymer / Roskam / Nikolai drag-polar method used in
%  the senior-design textbook script drag_calc_sol.m.
%% ------------------------------------------------------------------
ss.low.Cfe                  = 0.0035;  % Raymer Tab 12.3, jet fighter
ss.low.airfoil_a0_per_deg   = 2*pi*pi/180;   % thin-airfoil 2*pi per rad
ss.low.alpha_L0_deg         = ss.airfoil.alpha_L0_report_deg;
ss.low.EWD                  = 2.0;     % midpoint of the Raymer range
ss.low.transonic_blend      = [0.95 1.05];
ss.low.e_alpha_floor        = 0.60;

%% ------------------------------------------------------------------
%  Aero-deck build settings                                        [SIM]
%% ------------------------------------------------------------------
ss.deck.mach_grid = unique([0:0.05:0.90, 0.92:0.02:1.10, 1.15:0.05:1.30, ...
                            1.4:0.02:1.80, 1.9:0.1:2.2, ...
                            ss.targets.drag.Mach]);
ss.deck.altitude_ft_grid = [0 5000 10000 15000 20000 25000 30000 ...
                            35000 40000 45000 50000 55000 60000];
ss.deck.CL_grid = -0.4:0.1:1.8;
ss.deck.spanload_stations = 50;   % [RPT 11.4] 50 spanwise stations
ss.deck.cache_file = "SilverScytheAeroDeck.mat";

%% ------------------------------------------------------------------
%  Simulator-only settings                                         [SIM]
%% ------------------------------------------------------------------
ss.sim.ground_k = 3.0e5;
ss.sim.ground_c = 5.0e4;
ss.sim.aero_mach_max     = 2.20;
ss.sim.aero_beta_limit_deg = 15;
ss.sim.apply_compressibility = true;   % Prandtl-Glauert on the lift slope

% Body inertia about the takeoff CG.  Built by SilverScytheBuild from the
% component box inertias below; these overrides stay empty unless a CAD
% inertia tensor becomes available.
ss.sim.inertia_override_kgm2 = [];

%% ------------------------------------------------------------------
%  Reference stations, aft of the nose, ft                     [RPT Tab 12/13]
%% ------------------------------------------------------------------
ss.station.cg_ft       = ss.weights.takeoff_cg_station_ft;
ss.station.empty_cg_ft = ss.weights.empty_cg_station_ft;
ss.station.engine_ft   = 41.0;
ss.station.nozzle_ft   = 47.9;
ss.station.aero_ref_ft = ss.weights.takeoff_cg_station_ft;
ss.engine_y_ft         = 2.8;   % [RPT Table 7] nozzle spanwise location

%% ------------------------------------------------------------------
%  Takeoff / landing parameters consumed by TakeoffAnalysis and
%  LandingAnalysis through the coefficient lookup.
%% ------------------------------------------------------------------
ss.takeoff_params = struct( ...
    'mu_ground',                0.025, ...
    'mu_rolling',               0.025, ...
    'mu_braking',               0.40, ...
    'safety_factor',            1.15, ...
    'V2_to_Vs_ratio',           1.20, ...
    'VR_to_Vs_ratio',           1.10, ...
    'V1_to_VR_ratio',           0.92, ...
    'screen_height_takeoff_ft', 35, ...
    'rotation_alpha_deg',       10, ...
    'initial_climb_angle_deg',  10, ...
    'reaction_time_s',          1.0, ...
    'continue_time_after_VR_s', 2.0, ...
    'min_accel_mps2',           0.5, ...
    'min_reduced_accel_mps2',   0.3, ...
    'min_brake_decel_mps2',     1.0);

ss.landing_params = struct( ...
    'mu_braking',               0.45, ...
    'mu_spoiler',               0.00, ...
    'mu_reverser',              0.00, ...
    'approach_angle_deg',       3.5, ...
    'safety_factor',            1.0, ...
    'CD_spoiler',               0.04, ...
    'CD_reverser',              0.00, ...
    'CL_landing_touchdown',     0.05, ...
    'screen_height_landing_ft', 50, ...
    'flare_height_m',           5, ...
    'Vapp_to_Vs_ratio',         ss.carrier.Vapp_over_Vs, ...
    'Vtd_to_Vapp_ratio',        ss.carrier.Varrest_over_Vapp, ...
    'idle_throttle',            0.05, ...
    'use_idle_thrust',          true, ...
    'min_brake_decel_mps2',     1.5);

%% ------------------------------------------------------------------
%  Component mass build-up                        [RPT Tables 12 and 13]
%  {name, weight_lb, station_ft, y_ft, z_ft, [L W H] ft, type}
%  Box dimensions are [SIM] packaging estimates used only to give each
%  component a finite inertia tensor.
%% ------------------------------------------------------------------
ss.components = { ...
    'wing',                 6750, 29.0,  0.0, 0.0, [18 42 0.6], 'structure'; ...
    'ruddervator_pair',     3076, 42.0,  0.0, 0.0, [ 8 25 0.5], 'structure'; ...
    'fuselage_strakes',     5827, 29.0,  0.0, 0.0, [50  8 8.0], 'structure'; ...
    'main_landing_gear',    3008, 33.0,  0.0, 4.0, [ 4 12 5.0], 'landing_gear'; ...
    'nose_landing_gear',     532, 12.0,  0.0, 4.0, [ 3  4 5.0], 'landing_gear'; ...
    'arresting_gear',        330, 48.0,  0.0, 1.0, [ 3  6 2.0], 'carrier'; ...
    'catapult_gear',         201,  4.0,  0.0, 1.0, [ 3  4 2.0], 'carrier'; ...
    'nacelle',               200, 41.0,  0.0, 0.0, [ 8  8 4.0], 'structure'; ...
    'engine_left',          2445, 41.0, -2.8, 0.0, [13  3 3.0], 'engine'; ...
    'engine_right',         2445, 41.0,  2.8, 0.0, [13  3 3.0], 'engine'; ...
    'accessory_gearbox',     345, 41.0,  0.0, 0.0, [ 2  3 2.0], 'propulsion'; ...
    'exhaust_system',        782, 48.0,  0.0, 0.0, [ 5  6 3.0], 'propulsion'; ...
    'cooling_provisions',    343, 41.0,  0.0, 0.0, [ 4  4 2.0], 'propulsion'; ...
    'engine_controls',       234, 42.5,  0.0, 0.0, [ 2  2 1.0], 'propulsion'; ...
    'starting_system',       213, 42.5,  0.0, 0.0, [ 2  2 1.0], 'propulsion'; ...
    'fuel_system_pumps',    1006, 37.0,  0.0, 0.0, [ 6  4 2.0], 'propulsion'; ...
    'flight_controls',       312, 37.0,  0.0, 0.0, [ 5  4 2.0], 'equipment'; ...
    'apu',                   250, 41.0,  0.0, 0.0, [ 3  3 2.0], 'equipment'; ...
    'instruments',           150, 10.0,  0.0, 0.0, [ 3  3 2.0], 'equipment'; ...
    'hydraulic_system',      860, 37.0,  0.0, 0.0, [ 6  4 2.0], 'equipment'; ...
    'pneumatic_system',       61, 37.0,  0.0, 0.0, [ 3  2 1.0], 'equipment'; ...
    'electrical_system',     732, 10.0,  0.0, 0.0, [ 8  3 2.0], 'equipment'; ...
    'avionics',             2500, 27.0,  0.0, 0.0, [ 6  5 2.0], 'equipment'; ...
    'armament_provisions',   450, 11.0,  0.0, 0.0, [ 5  3 2.0], 'equipment'; ...
    'furnishings',           350, 37.0,  0.0, 0.0, [ 6  4 2.0], 'equipment'; ...
    'ecs',                   568, 29.0,  0.0, 0.0, [ 5  4 2.0], 'equipment'; ...
    'crew',                  400, 14.0,  0.0, 0.0, [ 3  2 4.0], 'useful'; ...
    'usable_fuel',         21000, 29.0,  0.0, 0.0, [18 24 3.0], 'fuel'; ...
    'trapped_fuel',          200, 29.0,  0.0, 0.0, [ 6  6 2.0], 'fuel'; ...
    'oil',                    70, 42.5,  0.0, 0.0, [ 2  2 1.0], 'useful'; ...
    'gun_system',            200, 15.0,  0.0, 0.0, [ 6  2 2.0], 'useful'; ...
    'ammunition',            200, 18.0,  0.0, 0.0, [ 4  2 2.0], 'useful'; ...
    'pylons_racks',          400, 34.5,  0.0, 0.0, [ 6  4 2.0], 'useful'; ...
    'strike_weapons',       4428, 30.0,  0.0, 0.0, [15  4 3.0], 'payload'; ...
    'flares_chaff',          120, 42.5,  0.0, 0.0, [ 2  2 1.0], 'useful' ...
};

%% ------------------------------------------------------------------
%  Paths
%% ------------------------------------------------------------------
ss.paths.package_dir   = fileparts(mfilename('fullpath'));
ss.paths.data_dir      = fullfile(ss.paths.package_dir,'data');
ss.paths.results_dir   = fullfile(ss.paths.package_dir,'results');
ss.paths.reference_dir = fullfile(ss.paths.package_dir,'reference');

end
