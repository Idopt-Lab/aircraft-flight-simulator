%% F-16 COMPLETE TRIM-CONSTRAINED PERFORMANCE ANALYSIS
% Aircraft definition and analysis configuration are separated below with
% MATLAB sections. All aerodynamic, propulsion, trim and performance
% quantities are returned by Aircraft and PerformanceAnalysis; this driver
% only selects analyses, reports their outputs and plots them.

clear;
clc;
close all;

%% PATH BOOTSTRAP
% Ensures the shared class library and example helper functions are on the
% MATLAB path regardless of the current working directory or session state.
this_script_dir = fileparts(mfilename('fullpath'));
matlab_root_dir = fileparts(fileparts(this_script_dir));
addpath(genpath(fullfile(matlab_root_dir,'classes')));
addpath(genpath(fullfile(matlab_root_dir,'examples')));

%% AIRCRAFT DEFINITION
% Geometry, mass properties, frames, controls, propulsion and load sources.
[aircraft,engine] = defineF16Aircraft();

%% PERFORMANCE ANALYSIS CONFIGURATION
% Flight conditions, grids, limits, solver options and trim bounds.

setup = F16PerformanceConfig(aircraft);
performance = aircraft.get_performance();

condition = setup.condition;
cfg = setup.performance;
solver = setup.solver;
fast_solver = setup.fast_solver;
bounds = setup.bounds;

%% DRY-POWER PERFORMANCE
% Every result below is returned directly by PerformanceAnalysis.

engine.set_rating_mode("dry");

results.dry.stall = performance.optimize_stall_speed(condition,solver,cfg,bounds.stall);
results.dry.best_range = performance.optimize_paper_best_range(condition,solver,cfg,bounds.range);
results.dry.best_endurance = performance.optimize_paper_best_endurance(condition,solver,cfg,bounds.endurance);
results.glide = performance.analyze_glide_performance(condition,solver,cfg,bounds.glide);
results.dry.velocity_sweep = performance.evaluate_velocity_sweep(condition,setup.velocity_grid_mps,setup.level_spec,fast_solver,cfg);
results.dry.required_available = performance.compute_paper_required_available_curve(condition,setup.velocity_grid_mps,fast_solver,cfg,bounds.range);

%% AFTERBURNER PERFORMANCE
% Climb, maximum speed, required/available curves and sustained turns.

engine.set_rating_mode("afterburner");

results.afterburner.best_angle_climb = performance.optimize_best_angle_climb(condition,solver,cfg,bounds.climb);
results.afterburner.best_rate_climb = performance.optimize_best_rate_climb(condition,solver,cfg,bounds.climb);
results.afterburner.maximum_level_speed = performance.optimize_max_level_speed(condition,solver,cfg,bounds.maximum_speed);
results.afterburner.velocity_sweep = performance.evaluate_velocity_sweep(condition,setup.velocity_grid_mps,setup.level_spec,fast_solver,cfg);
results.afterburner.required_available = performance.compute_paper_required_available_curve(condition,setup.velocity_grid_mps,fast_solver,cfg,bounds.range);
results.turn.sustained_curve = performance.compute_sustained_turn_curve(condition,setup.turn_velocity_grid_mps,fast_solver,cfg,bounds.turn);
results.turn.sustained_optimum = performance.optimize_sustained_turn(condition,solver,cfg,bounds.turn);
results.turn.coordinated = performance.optimize_coordinated_turn(setup.turn_condition,setup.coordinated_turn_bank_deg, solver,cfg,bounds.turn);

%% MANEUVER AND GUST ENVELOPES

results.envelope.accelerated_stall = performance.compute_accelerated_stall(condition,setup.accelerated_stall_load_factors,cfg);
results.envelope.vn = performance.compute_vn_diagram(condition,setup.envelope_velocity_grid_mps,cfg);
results.envelope.gust = performance.compute_gust_envelope(condition,setup.envelope_velocity_grid_mps,setup.gust);
results.turn.instantaneous = performance.compute_instantaneous_turn_curve(condition,setup.turn_velocity_grid_mps,cfg);

%% ENERGY-MANEUVERABILITY AND ACCELERATION

results.energy.Ps_map = performance.compute_specific_excess_power_map(condition,setup.Ps_velocity_grid_mps,setup.Ps_load_factors, fast_solver,cfg,bounds.turn);
results.energy.acceleration = performance.compute_acceleration_schedule(condition,setup.velocity_grid_mps,fast_solver,cfg,bounds.maximum_speed);

%% ALTITUDE ENVELOPE, CLIMB SCHEDULE AND CEILINGS

results.altitude.mach_envelope = performance.compute_mach_altitude_envelope(condition,setup.altitude_grid_m,solver,cfg,bounds.envelope);
results.altitude.climb = performance.optimize_climb_schedule(condition,setup.altitude_grid_m,solver,cfg,bounds.climb);

%% COLLECT COMPLETE FRAMEWORK OUTPUTS
% Trim states, control vectors and control summaries remain inside point_opt.

results.aircraft = aircraft;
results.engine = engine;
results.setup = setup;

%% PERFORMANCE SUMMARY
% Values and control vectors are extracted directly from point_opt results.

printF16PerformanceSummary(results);

%% PERFORMANCE PLOTS
% Plot only fields already calculated and returned by PerformanceAnalysis.

plotF16Performance(results);

F16Performance = results;

%% LOCAL AIRCRAFT-DEFINITION FUNCTION
function [aircraft,engine] = defineF16Aircraft()

%% Reference geometry and mass properties

ft_to_m = 0.3048;
ft2_to_m2 = 0.09290304;
lbf_to_N = 4.4482216152605;
g = 9.80665;

aircraft = Aircraft();
aircraft.geometry.set_reference_geometry(300*ft2_to_m2,30*ft_to_m,11.3201786951274*ft_to_m);
aircraft.geometry.set_reference_point([0 0 0]);
aircraft.geometry.aero_mach_max = 2;

analysis_weight_lbf = 0.899666962714079*31377;
empty_weight_lbf = 19980.7005781593;
permanent_payload_lbf = 700;
expendable_payload_lbf = 4400;
engine_weight_lbf = 0.199*23770;
fuel_weight_lbf = max(analysis_weight_lbf-empty_weight_lbf- permanent_payload_lbf-expendable_payload_lbf,0);
airframe_weight_lbf = empty_weight_lbf-engine_weight_lbf;

airframe_mass_kg = airframe_weight_lbf*lbf_to_N/g;
fuel_mass_kg = fuel_weight_lbf*lbf_to_N/g;
engine_mass_kg = engine_weight_lbf*lbf_to_N/g;
permanent_payload_mass_kg = permanent_payload_lbf*lbf_to_N/g;
expendable_payload_mass_kg = expendable_payload_lbf*lbf_to_N/g;

inertia_scale = analysis_weight_lbf/20500;
inertia_correction = 1.3558179483314;
airframe_inertia_kgm2 = inertia_scale*[9496*inertia_correction,0,982*inertia_correction; 0,55814*inertia_correction,0; 982*inertia_correction,0,63100*inertia_correction];

%% Reference frames and component hierarchy

aircraft.add_frame("fuel_tank","body",[0;0;0],@(x) eye(3));
aircraft.add_frame("engine","body",[0;0;0],@(x) eye(3));
aircraft.add_frame("cg","body",[0;0;0],@(x) eye(3));
aircraft.add_frame("gravity_cg","body",[0;0;0], @(x) ReferenceFrame.ned_to_body_dcm(x));

airframe = Component("airframe",airframe_mass_kg,[0;0;0], airframe_inertia_kgm2,aircraft.get_frame("body"),"airframe");
fuel = Component("fuel",fuel_mass_kg,[0;0;0],zeros(3), aircraft.get_frame("fuel_tank"),"fuel");
engine_component = Component("engine",engine_mass_kg,[0;0;0],zeros(3), aircraft.get_frame("engine"),"engine");
permanent_payload = Component("permanent_payload", permanent_payload_mass_kg,[0;0;0],zeros(3), aircraft.get_frame("body"),"payload");
expendable_payload = Component("expendable_payload", expendable_payload_mass_kg,[0;0;0],zeros(3), aircraft.get_frame("body"),"payload");

aircraft.add_component(airframe);
aircraft.add_component(fuel);
aircraft.add_component(engine_component);
aircraft.add_component(permanent_payload);
aircraft.add_component(expendable_payload);

%% Control surfaces

aircraft.add_control_surface(ControlSurface("aileron","aileron","primary",[1 0 0], deg2rad(21.5),deg2rad(-21.5),0,0,0));
aircraft.add_control_surface(ControlSurface("stabilator","elevator","primary",[0 1 0], deg2rad(25),deg2rad(-25),0,0,0));
aircraft.add_control_surface(ControlSurface("rudder","rudder","primary",[0 0 1], deg2rad(30),deg2rad(-30),0,0,0));

%% Propulsion model

engine = BrandtAfterburningEngine("F16A_engine",aircraft.get_frame("engine"),[1;0;0]);
engine.set_rating_mode("dry");
aircraft.add_propulsive_element(engine);
engine_component.add_load_source(PropulsionLoadSolver(engine,aircraft.get_frame("engine")));

%% Aerodynamic and gravity load sources

aerodynamics = CoefficientAerodynamics(@F16ABrandtLookup);
aircraft.add_load_source(AeroLoadSolver(aerodynamics,aircraft.geometry,aircraft,aircraft.get_frame("body")));
aircraft.add_load_source(GravityLoadSolver(aircraft,aircraft.get_frame("gravity_cg")));

%% Final mass-property reference

[~,cg_body_m,~] = aircraft.compute_total_mass_properties();
aircraft.update_frame_position("cg",cg_body_m);
aircraft.update_frame_position("gravity_cg",cg_body_m);
aircraft.set_reference_frame("cg");

end

%% LOCAL PERFORMANCE-CONFIGURATION FUNCTION
function setup = F16PerformanceConfig(aircraft)

%% Conditions and analysis grids

setup.condition = struct('altitude_m',9000);
setup.turn_condition = struct('altitude_m',9000,'velocity_mps',250);
setup.coordinated_turn_bank_deg = 60;

setup.velocity_grid_mps = (140:10:550).';
setup.turn_velocity_grid_mps = (150:15:330).';
setup.Ps_velocity_grid_mps = (150:20:330).';
setup.Ps_load_factors = (1:1:6).';
% The ceiling interpolator needs optimized ROC values on both sides of the
% requested threshold. This is only the search domain; PerformanceAnalysis
% determines the actual ceiling from the resulting ROC schedule.
setup.altitude_grid_m = (0:2000:30000).';
setup.envelope_velocity_grid_mps = linspace(0,550,500).';
setup.accelerated_stall_load_factors = (1:0.5:9).';

%% Aircraft performance limits and propulsion settings

setup.performance = struct( ...
    'available_throttle',1, ...
    'continuous_throttle',1, ...
    'takeoff_throttle',1, ...
    'CL_max',0.9840156989811455, ...
    'CL_min',-0.8, ...
    'max_Mach',2, ...
    'propulsion_type',"thrust", ...
    'range_endurance_metric',"actual_flow", ...
    'fuel_LHV_Jpkg',43e6, ...
    'n_limit_pos',9, ...
    'n_limit_neg',-3, ...
    'max_dynamic_pressure_Pa',2133*47.88025898033584, ...
    'service_ceiling_threshold_fpm',100, ...
    'enforce_propulsion_constraints',true);

%% Nonlinear solver settings

setup.solver = struct( ...
    'residual_tolerance',1e-5, ...
    'inequality_tolerance',1e-6, ...
    'equality_tolerance',1e-5, ...
    'optimality_tolerance',1e-5, ...
    'fmincon_options',optimoptions("fmincon", ...
        "Algorithm","sqp", ...
        "Display","none", ...
        "MaxIterations",1200, ...
        "MaxFunctionEvaluations",20000, ...
        "OptimalityTolerance",1e-8, ...
        "ConstraintTolerance",1e-8, ...
        "StepTolerance",1e-10, ...
        "FiniteDifferenceType","central", ...
        "FiniteDifferenceStepSize",1e-4));

setup.fast_solver = struct( ...
    'residual_tolerance',1e-5, ...
    'inequality_tolerance',1e-6, ...
    'equality_tolerance',1e-5, ...
    'optimality_tolerance',1e-5, ...
    'fmincon_options',optimoptions("fmincon", ...
        "Algorithm","sqp", ...
        "Display","none", ...
        "MaxIterations",350, ...
        "MaxFunctionEvaluations",6000, ...
        "OptimalityTolerance",1e-6, ...
        "ConstraintTolerance",1e-6, ...
        "StepTolerance",1e-8, ...
        "FiniteDifferenceType","forward", ...
        "FiniteDifferenceStepSize",1e-4));

%% Level-flight trim specification

setup.level_spec = struct( ...
    'mode',"level", ...
    'variables',["alpha";"control_pitch";"throttle"], ...
    'initial_guess',[deg2rad(8);deg2rad(1.5);0.4], ...
    'lb',[deg2rad(-5);deg2rad(-25);0], ...
    'ub',[deg2rad(20);deg2rad(25);1], ...
    'fixed',struct('beta',0,'phi',0,'psi',0,'gamma',0, ...
        'control_roll',0,'control_yaw',0), ...
    'reference_frame_name',"cg", ...
    'residual_scale',[aircraft.g;aircraft.g;1]);

%% Optimization bounds

setup.bounds.stall = struct( ...
    'V_lb',40,'V_ub',250,'V0',120, ...
    'alpha0_deg',15,'alpha_min_deg',-5,'alpha_max_deg',20, ...
    'elevator0_deg',1,'elevator_min_deg',-25,'elevator_max_deg',25, ...
    'throttle0',0.4,'throttle_min',0,'throttle_max',1);

setup.bounds.climb = struct( ...
    'V_lb',60,'V_ub',600,'V0',350, ...
    'alpha0_deg',1.5,'alpha_min_deg',-5,'alpha_max_deg',20, ...
    'elevator0_deg',1.7,'elevator_min_deg',-25,'elevator_max_deg',25, ...
    'gamma0_deg',12,'gamma_min_deg',0,'gamma_max_deg',45, ...
    'ceiling_gamma_min_deg',-10,'ceiling_descent_seed_deg',-3, ...
    'throttle0',1);

setup.bounds.range = struct( ...
    'V_lb',60,'V_ub',400,'V0',230, ...
    'alpha0_deg',5,'alpha_min_deg',-5,'alpha_max_deg',20, ...
    'elevator0_deg',1,'elevator_min_deg',-25,'elevator_max_deg',25, ...
    'throttle0',0.4,'throttle_min',0,'throttle_max',1);

setup.bounds.endurance = struct( ...
    'V_lb',60,'V_ub',400,'V0',190, ...
    'alpha0_deg',6,'alpha_min_deg',-5,'alpha_max_deg',20, ...
    'elevator0_deg',1,'elevator_min_deg',-25,'elevator_max_deg',25, ...
    'throttle0',0.35,'throttle_min',0,'throttle_max',1);

setup.bounds.glide = struct( ...
    'V_lb',50,'V_ub',350,'V0',190, ...
    'alpha0_deg',6,'alpha_min_deg',-5,'alpha_max_deg',20, ...
    'elevator0_deg',1,'elevator_min_deg',-25,'elevator_max_deg',25, ...
    'gamma0_deg',-5,'gamma_min_deg',-30,'gamma_max_deg',-0.01);

setup.bounds.maximum_speed = struct( ...
    'V_lb',100,'V_ub',600,'V0',450, ...
    'alpha0_deg',1,'alpha_min_deg',-5,'alpha_max_deg',20, ...
    'elevator0_deg',1,'elevator_min_deg',-25,'elevator_max_deg',25, ...
    'throttle0',1,'throttle_min',0,'throttle_max',1);

setup.bounds.turn = struct( ...
    'V_lb',100,'V_ub',350,'V0',240, ...
    'bank0_deg',45,'bank_min_deg',1,'bank_max_deg',80, ...
    'throttle0',1,'throttle_min',0,'throttle_max',1, ...
    'turn_rate_ub_radps',1.2, ...
    'alpha0_deg',8,'alpha_min_deg',-5,'alpha_max_deg',20, ...
    'elevator0_deg',0,'elevator_min_deg',-25,'elevator_max_deg',25, ...
    'aileron0_deg',0,'aileron_min_deg',-21.5,'aileron_max_deg',21.5, ...
    'rudder0_deg',0,'rudder_min_deg',-30,'rudder_max_deg',30, ...
    'Ps_trim_tolerance_mps',0.05);

setup.bounds.envelope = struct( ...
    'V_lb',40,'V_ub',650,'V0',450, ...
    'max_Mach',2,'stall_speed_factor',1.001, ...
    'alpha0_deg',2,'alpha_min_deg',-5,'alpha_max_deg',20, ...
    'elevator0_deg',1,'elevator_min_deg',-25,'elevator_max_deg',25);

%% Gust-model inputs

x_reference = zeros(12,1);
x_reference(3) = -setup.condition.altitude_m;
x_reference(4) = setup.turn_condition.velocity_mps;
reference_coefficients = F16ABrandtLookup(x_reference,zeros(4,1),aircraft.geometry);
setup.gust = struct('CL_alpha_per_rad',reference_coefficients.derivatives.CLa, 'Ude_mps',15.24, 'Kg',0.88);

end

%% LOCAL PLOTTING FUNCTION
function plotF16Performance(results)

%% Extract valid velocity-sweep points

dry = sweepSeries(results.dry.velocity_sweep);
afterburner = sweepSeries(results.afterburner.velocity_sweep);
dry_curve = results.dry.required_available;
afterburner_curve = results.afterburner.required_available;
instantaneous = results.turn.instantaneous;
sustained = results.turn.sustained_curve;
vn = results.envelope.vn;
gust = results.envelope.gust;
Ps_map = results.energy.Ps_map;
acceleration = results.energy.acceleration;
mach_envelope = results.altitude.mach_envelope;
climb = results.altitude.climb;

%% Aerodynamic plots

figure("Name","Aerodynamics");
tiledlayout(2,2,"TileSpacing","compact");
nexttile;
plot(dry.V_mps,dry.L_N/1000,"LineWidth",1.5);
hold on;
yline(dry.weight_N(1)/1000,"--");
grid on;
xlabel("V [m/s]");
ylabel("Force [kN]");
legend("Lift","Weight","Location","best");
nexttile;
plot(dry.V_mps,dry.D_N/1000,"LineWidth",1.5);
grid on;
xlabel("V [m/s]");
ylabel("Drag [kN]");
nexttile;
plot(dry.V_mps,dry.CL,"LineWidth",1.5);
hold on;
plot(dry.V_mps,dry.CD,"LineWidth",1.5);
grid on;
xlabel("V [m/s]");
ylabel("Coefficient");
legend("C_L","C_D","Location","best");
nexttile;
plot(dry.V_mps,dry.L_over_D,"LineWidth",1.5);
grid on;
xlabel("V [m/s]");
ylabel("L/D");

%% Trim-state and control plots

figure("Name","Trim States and Controls");
tiledlayout(3,2,"TileSpacing","compact");
nexttile;
plot(dry.V_mps,dry.alpha_deg,"LineWidth",1.5);
grid on;
xlabel("V [m/s]");
ylabel("\alpha [deg]");
nexttile;
plot(dry.V_mps,dry.theta_deg,"LineWidth",1.5);
grid on;
xlabel("V [m/s]");
ylabel("\theta [deg]");
nexttile;
plot(dry.V_mps,dry.elevator_deg,"LineWidth",1.5);
grid on;
xlabel("V [m/s]");
ylabel("Stabilator [deg]");
nexttile;
plot(dry.V_mps,dry.throttle,"LineWidth",1.5);
hold on;
plot(afterburner.V_mps,afterburner.throttle,"LineWidth",1.5);
grid on;
xlabel("V [m/s]");
ylabel("Throttle");
legend("Dry","Afterburner","Location","best");
nexttile;
plot(dry.V_mps,dry.fuel_flow_trim,"LineWidth",1.5);
hold on;
plot(afterburner.V_mps,afterburner.fuel_flow_trim,"LineWidth",1.5);
grid on;
xlabel("V [m/s]");
ylabel("Fuel flow [kg/s]");
legend("Dry","Afterburner","Location","best");
nexttile;
semilogy(dry.V_mps,dry.residual_inf,"LineWidth",1.5);
grid on;
xlabel("V [m/s]");
ylabel("Residual infinity norm");

%% Propulsion and power plots

figure("Name","Thrust and Power");
tiledlayout(2,2,"TileSpacing","compact");
nexttile;
plot(dry.V_mps,dry.T_required_N/1000,"LineWidth",1.5);
hold on;
plot(dry.V_mps,dry.T_available_N/1000,"LineWidth",1.5);
plot(afterburner.V_mps,afterburner.T_available_N/1000,"LineWidth",1.5);
grid on;
xlabel("V [m/s]");
ylabel("Thrust [kN]");
legend("Required","Dry available","AB available","Location","best");
nexttile;
plot(dry.V_mps,dry.T_excess_N/1000,"LineWidth",1.5);
hold on;
plot(afterburner.V_mps,afterburner.T_excess_N/1000,"LineWidth",1.5);
yline(0,"--");
grid on;
xlabel("V [m/s]");
ylabel("Excess thrust [kN]");
legend("Dry","Afterburner","Location","best");
nexttile;
plot(dry.V_mps,dry.P_thrust_required_W/1e6,"LineWidth",1.5);
hold on;
plot(dry.V_mps,dry.P_thrust_available_W/1e6,"LineWidth",1.5);
plot(afterburner.V_mps,afterburner.P_thrust_available_W/1e6, "LineWidth",1.5);
grid on;
xlabel("V [m/s]");
ylabel("Power [MW]");
legend("Required","Dry available","AB available","Location","best");
nexttile;
plot(dry.V_mps,dry.P_excess_W/1e6,"LineWidth",1.5);
hold on;
plot(afterburner.V_mps,afterburner.P_excess_W/1e6,"LineWidth",1.5);
yline(0,"--");
grid on;
xlabel("V [m/s]");
ylabel("Excess power [MW]");
legend("Dry","Afterburner","Location","best");

%% Energy, range and endurance plots

figure("Name","Energy, Range and Endurance");
tiledlayout(2,3,"TileSpacing","compact");
nexttile;
plot(dry.V_mps,dry.Ps_mps,"LineWidth",1.5);
hold on;
plot(afterburner.V_mps,afterburner.Ps_mps,"LineWidth",1.5);
yline(0,"--");
grid on;
xlabel("V [m/s]");
ylabel("P_s [m/s]");
legend("Dry","Afterburner","Location","best");
nexttile;
plot(dry.V_mps,dry.ROC_available_mps,"LineWidth",1.5);
hold on;
plot(afterburner.V_mps,afterburner.ROC_available_mps,"LineWidth",1.5);
yline(0,"--");
grid on;
xlabel("V [m/s]");
ylabel("Available ROC [m/s]");
legend("Dry","Afterburner","Location","best");
nexttile;
plot(dry.V_mps,dry.speed_acceleration_mps2,"LineWidth",1.5);
hold on;
plot(afterburner.V_mps,afterburner.speed_acceleration_mps2,"LineWidth",1.5);
yline(0,"--");
grid on;
xlabel("V [m/s]");
ylabel("Acceleration [m/s^2]");
legend("Dry","Afterburner","Location","best");
nexttile;
plot(dry.V_mps,dry.fuel_specific_range_m_per_kg/1000,"LineWidth",1.5);
grid on;
xlabel("V [m/s]");
ylabel("Specific range [km/kg]");
nexttile;
plot(dry.V_mps,dry.fuel_endurance_s_per_kg,"LineWidth",1.5);
grid on;
xlabel("V [m/s]");
ylabel("Endurance [s/kg]");
nexttile;
plot(dry.V_mps,dry.energy_rate_trim_W/1e6,"LineWidth",1.5);
grid on;
xlabel("V [m/s]");
ylabel("Energy rate [MW]");

%% Required and available performance plots

figure("Name","Required and Available Performance");
tiledlayout(2,2,"TileSpacing","compact");
nexttile;
plot(dry_curve.V_mps,dry_curve.required.thrust_N/1000,"LineWidth",1.5);
hold on;
plot(dry_curve.V_mps,dry_curve.available.thrust_N/1000,"LineWidth",1.5);
plot(afterburner_curve.V_mps, afterburner_curve.available.thrust_N/1000,"LineWidth",1.5);
grid on;
xlabel("V [m/s]");
ylabel("Thrust [kN]");
legend("Required","Dry available","AB available","Location","best");
nexttile;
plot(dry_curve.V_mps,dry_curve.required.power_W/1e6,"LineWidth",1.5);
hold on;
plot(dry_curve.V_mps,dry_curve.available.power_W/1e6,"LineWidth",1.5);
plot(afterburner_curve.V_mps, afterburner_curve.available.power_W/1e6,"LineWidth",1.5);
grid on;
xlabel("V [m/s]");
ylabel("Power [MW]");
legend("Required","Dry available","AB available","Location","best");
nexttile;
plot(dry_curve.V_mps,dry_curve.excess_thrust_N/1000,"LineWidth",1.5);
hold on;
plot(afterburner_curve.V_mps, afterburner_curve.excess_thrust_N/1000,"LineWidth",1.5);
yline(0,"--");
grid on;
xlabel("V [m/s]");
ylabel("Excess thrust [kN]");
legend("Dry","Afterburner","Location","best");
nexttile;
plot(dry_curve.V_mps,dry_curve.excess_power_W/1e6,"LineWidth",1.5);
hold on;
plot(afterburner_curve.V_mps, afterburner_curve.excess_power_W/1e6,"LineWidth",1.5);
yline(0,"--");
grid on;
xlabel("V [m/s]");
ylabel("Excess power [MW]");
legend("Dry","Afterburner","Location","best");

%% Flight-envelope plots

figure("Name","V-n and Gust Envelope");
fill([vn.V_mps;flipud(vn.V_mps)], [vn.n_pos;flipud(vn.n_neg)],[0.9 0.9 0.9],"EdgeColor","none");
hold on;
plot(vn.V_mps,vn.n_pos,"LineWidth",2);
plot(vn.V_mps,vn.n_neg,"LineWidth",2);
plot(gust.V_mps,gust.n_gust_pos,"--","LineWidth",1.5);
plot(gust.V_mps,gust.n_gust_neg,"--","LineWidth",1.5);
yline(1,":");
xline(vn.Vs_mps,"--","V_s");
xline(vn.Va_mps,"--","V_A");
grid on;
xlabel("V_{TAS} [m/s]");
ylabel("Load factor");
legend("Envelope","Positive maneuver","Negative maneuver", "Positive gust","Negative gust","Location","best");

%% Turn-performance plots

figure("Name","Turn Performance");
tiledlayout(2,3,"TileSpacing","compact");
iv = instantaneous.valid;
sv = sustained.valid;
nexttile;
plot(instantaneous.V_mps(iv),instantaneous.turn_rate_degps(iv), "LineWidth",1.5);
hold on;
plot(sustained.V_mps(sv),sustained.turn_rate_degps(sv), "o-","LineWidth",1.5);
grid on;
xlabel("V [m/s]");
ylabel("Turn rate [deg/s]");
legend("Instantaneous","Sustained","Location","best");
nexttile;
plot(instantaneous.V_mps(iv),instantaneous.turn_radius_m(iv), "LineWidth",1.5);
hold on;
plot(sustained.V_mps(sv),sustained.turn_radius_m(sv), "o-","LineWidth",1.5);
grid on;
xlabel("V [m/s]");
ylabel("Turn radius [m]");
legend("Instantaneous","Sustained","Location","best");
nexttile;
plot(sustained.V_mps(sv),sustained.bank_deg(sv),"o-","LineWidth",1.5);
grid on;
xlabel("V [m/s]");
ylabel("Bank [deg]");
nexttile;
plot(sustained.V_mps(sv),sustained.load_factor_n(sv), "o-","LineWidth",1.5);
grid on;
xlabel("V [m/s]");
ylabel("Load factor");
nexttile;
plot(sustained.V_mps(sv),sustained.throttle(sv),"o-","LineWidth",1.5);
grid on;
xlabel("V [m/s]");
ylabel("Throttle");
nexttile;
plot(sustained.V_mps(sv),sustained.elevator_deg(sv), "o-","LineWidth",1.5);
hold on;
plot(sustained.V_mps(sv),sustained.aileron_deg(sv), "o-","LineWidth",1.5);
plot(sustained.V_mps(sv),sustained.rudder_deg(sv), "o-","LineWidth",1.5);
grid on;
xlabel("V [m/s]");
ylabel("Control [deg]");
legend("Elevator","Aileron","Rudder","Location","best");

%% Energy-maneuverability plots

figure("Name","Specific Excess Power Map");
contourf(Ps_map.V_mps,Ps_map.n,Ps_map.Ps_mps,20);
colorbar;
grid on;
xlabel("V [m/s]");
ylabel("Load factor");
title("P_s [m/s]");

figure("Name","Acceleration Schedule");
tiledlayout(2,2,"TileSpacing","compact");
nexttile;
plot(acceleration.V_mps,acceleration.acceleration_mps2,"LineWidth",1.5);
yline(0,"--");
grid on;
xlabel("V [m/s]");
ylabel("Acceleration [m/s^2]");
nexttile;
plot(acceleration.V_mps,acceleration.Ps_mps,"LineWidth",1.5);
yline(0,"--");
grid on;
xlabel("V [m/s]");
ylabel("P_s [m/s]");
nexttile;
plot(acceleration.V_mps,acceleration.time_s,"LineWidth",1.5);
grid on;
xlabel("V [m/s]");
ylabel("Cumulative time [s]");
nexttile;
plot(acceleration.V_mps,acceleration.fuel_used_kg,"LineWidth",1.5);
grid on;
xlabel("V [m/s]");
ylabel("Cumulative fuel [kg]");

%% Altitude, climb and ceiling plots

figure("Name","Mach-Altitude Envelope and Climb");
tiledlayout(2,2,"TileSpacing","compact");
nexttile;
plot(mach_envelope.stall_Mach,mach_envelope.altitude_m/1000, "LineWidth",1.5);
hold on;
plot(mach_envelope.max_level_Mach,mach_envelope.altitude_m/1000, "LineWidth",1.5);
grid on;
xlabel("Mach");
ylabel("Altitude [km]");
legend("Stall","Maximum level","Location","best");
nexttile;
plot(climb.best_ROC_mps,climb.altitude_m/1000,"o-","LineWidth",1.5);
xline(0,"--");
grid on;
xlabel("ROC [m/s]");
ylabel("Altitude [km]");
nexttile;
plot(climb.best_V_mps,climb.altitude_m/1000,"o-","LineWidth",1.5);
grid on;
xlabel("Best climb V [m/s]");
ylabel("Altitude [km]");
nexttile;
plot(climb.best_gamma_deg,climb.altitude_m/1000,"o-","LineWidth",1.5);
xline(0,"--");
grid on;
xlabel("Climb angle [deg]");
ylabel("Altitude [km]");

end

function data = sweepSeries(sweep)

indices = find(sweep.valid);
if isempty(indices)
    error('plotF16Performance:NoValidSweepPoints', 'The velocity sweep contains no valid performance points.');
end
points = [sweep.point{indices}];

fields = [ ...
    "V_mps";"weight_N";"L_N";"D_N";"CL";"CD";"L_over_D"; ...
    "alpha_deg";"theta_deg";"elevator_deg";"throttle"; ...
    "fuel_flow_trim";"T_required_N"; ...
    "T_available_N";"T_excess_N";"P_thrust_required_W"; ...
    "P_thrust_available_W";"P_excess_W";"Ps_mps"; ...
    "ROC_available_mps";"speed_acceleration_mps2"; ...
    "fuel_specific_range_m_per_kg";"fuel_endurance_s_per_kg"; ...
    "energy_rate_trim_W"];

for k = 1:numel(fields)
    name = fields(k);
    data.(name) = reshape([points.(name)],[],1);
end
data.residual_inf = sweep.residual_norm(indices);

end

%% LOCAL RESULTS-SUMMARY FUNCTION
function printF16PerformanceSummary(results)

rows = [ ...
    pointTableRow("Stall",results.dry.stall); ...
    pointTableRow("Best range",results.dry.best_range); ...
    pointTableRow("Best endurance",results.dry.best_endurance); ...
    pointTableRow("Best glide",results.glide.best_glide); ...
    pointTableRow("Minimum sink",results.glide.min_sink); ...
    pointTableRow("Best angle climb", ...
        results.afterburner.best_angle_climb); ...
    pointTableRow("Best rate climb", ...
        results.afterburner.best_rate_climb); ...
    pointTableRow("Maximum level speed", ...
        results.afterburner.maximum_level_speed); ...
    pointTableRow("Maximum sustained turn", ...
        results.turn.sustained_optimum); ...
    pointTableRow("Coordinated turn",results.turn.coordinated)];

fprintf('\n=== F-16 PERFORMANCE POINTS AND TRIM CONTROLS ===\n');
disp(rows);

climb = results.altitude.climb;
fprintf('\n=== CEILING RESULTS FROM OPTIMIZED ROC SCHEDULE ===\n');
fprintf('Service ceiling  : %.3f m (%.3f ft)\n', climb.service_ceiling_m,climb.service_ceiling_m*3.280839895);
fprintf('Absolute ceiling : %.3f m (%.3f ft)\n', climb.absolute_ceiling_m,climb.absolute_ceiling_m*3.280839895);

climb_table = table( ...
    climb.altitude_m,climb.best_V_mps,climb.best_Mach, ...
    climb.best_gamma_deg,climb.best_ROC_mps,climb.best_ROC_fpm, ...
    climb.best_fuel_flow,climb.point_valid,climb.solution_status, ...
    'VariableNames',{'Altitude_m','V_y_mps','Mach','Gamma_deg', ...
    'ROC_mps','ROC_fpm','FuelFlow_kgps','Valid','Status'});
disp(climb_table);

end

function row = pointTableRow(metric,opt)

p = opt.point_opt;
row = table( ...
    string(metric),logical(opt.converged),p.V_mps,p.Mach,p.alpha_deg, ...
    p.beta_deg,p.gamma_trim_deg,p.phi_deg,p.theta_deg, ...
    p.aileron_deg,p.elevator_deg,p.rudder_deg,p.throttle, ...
    p.ROC_trim_mps,p.fuel_flow_trim,{p.u_trim(:).'}, ...
    'VariableNames',{'Metric','Converged','V_mps','Mach','Alpha_deg', 'Beta_deg','Gamma_deg','Phi_deg','Theta_deg','Aileron_deg', 'Stabilator_deg','Rudder_deg','Throttle','ROC_mps', 'FuelFlow_kgps','ControlVector'});

end
