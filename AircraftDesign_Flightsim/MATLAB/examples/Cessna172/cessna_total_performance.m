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

%% USER CONFIGURATION
% DATCOM files, reference flight condition and analysis grids.

data_directory = string(fileparts(mfilename('fullpath')));
datcom_files = fullfile(data_directory,["cessna010.out"; "cessna012.out"; "cessna014.out"; "cessna016.out"; "cessna018.out"; "cessna020.out"; "cessna022.out"; "cessna024.out"]);

reference_altitude_m = 1500;
reference_velocity_mps = 50.173;

velocity_grid_mps = linspace(25,80,70).';
turn_velocity_grid_mps = linspace(25,80,18).';
Ps_velocity_grid_mps = linspace(25,80,14).';
Ps_load_factors = [1;1.5;2;2.5;3;3.5];
acceleration_velocity_grid_mps = linspace(25,80,24).';
envelope_velocity_grid_mps = linspace(0,80,220).';
altitude_grid_m = (0:500:9000).';

%% AIRCRAFT DEFINITION
% Geometry, mass properties, frames, controls, aerodynamics and propulsion.

aircraft = Aircraft();

%% Geometry and reference frames

cg_from_nose_m = 1.90;
wing_ac_from_nose_m = 1.57;
propeller_from_nose_m = 1.68;
fuel_from_nose_m = 1.70;

S_ref_m2 = 16.17;
b_ref_m = 10.98;
c_ref_m = 1.50;

wing_ac_body_m = [wing_ac_from_nose_m-cg_from_nose_m;0;0];
propeller_body_m = [propeller_from_nose_m-cg_from_nose_m;0;0];
fuel_body_m = [fuel_from_nose_m-cg_from_nose_m;0;0];

aircraft.geometry.set_reference_geometry(S_ref_m2,b_ref_m,c_ref_m);
aircraft.geometry.set_reference_point(wing_ac_body_m);
aircraft.set_body_frame("body");
aircraft.set_reference_frame("body");

aircraft.add_frame("aero_ref","body",wing_ac_body_m,@(x) eye(3));
aircraft.add_frame("fuel_tank","body",fuel_body_m,@(x) eye(3));
aircraft.add_frame("propeller","body",propeller_body_m,@(x) eye(3));
aircraft.add_frame("cg","body",[0;0;0],@(x) eye(3));
aircraft.add_frame("gravity_cg","body",[0;0;0], @(x) ReferenceFrame.ned_to_body_dcm(x));

%% Mass properties and component hierarchy

target_gross_mass_kg = 2550*0.45359237;
empty_mass_kg = 767;
fuel_mass_kg = 144;
payload_mass_kg = target_gross_mass_kg-empty_mass_kg-fuel_mass_kg;

airframe_inertia_kgm2 = diag([1285 1824 2666]);

airframe = Component("airframe",empty_mass_kg,[0;0;0], airframe_inertia_kgm2,aircraft.get_frame("body"),"airframe");
fuel = Component("fuel",fuel_mass_kg,[0;0;0],zeros(3), aircraft.get_frame("fuel_tank"),"fuel");
payload = Component("payload",payload_mass_kg,[0.20;0;0],zeros(3), aircraft.get_frame("body"),"payload");

aircraft.add_component(airframe);
aircraft.add_component(fuel);
aircraft.add_component(payload);

%% Control surfaces

aircraft.add_control_surface(ControlSurface("aileron","aileron","primary",[1 0 0], deg2rad(20),deg2rad(-20),0,0,0));
aircraft.add_control_surface(ControlSurface("elevator","elevator","primary",[0 1 0], deg2rad(28),deg2rad(-26),0,0,0));
aircraft.add_control_surface(ControlSurface("rudder","rudder","primary",[0 0 1], deg2rad(30),deg2rad(-30),0,0,0));

%% DATCOM aerodynamic model

datcom_options = struct('CD0_add',0.01681, 'CD_alpha_extra',0, 'CL_max_abs',1.9);

[datcom_lookup,datcom_data] = c172datcom(datcom_files,true,datcom_options);

aerodynamic_model = CoefficientAerodynamics(datcom_lookup);
aircraft.add_load_source(AeroLoadSolver(aerodynamic_model,aircraft.geometry,aircraft, aircraft.get_frame("aero_ref")));

aero_mach_limit = max(datcom_data.mach,[],'all');
alpha_data_min_deg = rad2deg(min(datcom_data.alpha,[],'all'));
alpha_data_max_deg = rad2deg(max(datcom_data.alpha,[],'all'));
alpha_trim_min_deg = max(-8,alpha_data_min_deg);
alpha_trim_max_deg = min(18,alpha_data_max_deg);

if alpha_trim_max_deg <= alpha_trim_min_deg
    error('C172PerformanceAnalysis:InvalidDATCOMAlphaRange', 'The DATCOM database does not contain a usable trim-alpha range.');
end

CL_max_analysis = min(1.60,max(abs(datcom_data.CL),[],'all'));
CL_min_analysis = max(-0.80,min(datcom_data.CL,[],'all'));

%% Piston-engine and fixed-pitch propeller models

avgas_density_kg_per_usgal = 2.72;
maximum_engine_power_W = 180*745.7;

engine_model = PistonEngineModel(maximum_engine_power_W,2700,"Lycoming_O_360_A4M");
engine_model.minimum_rpm = 600;
engine_model.maximum_rpm = 2700;
engine_model.density_lapse_exponent = 1;
engine_model.installation_power_factor = 1;
engine_model.set_rpm_power_map([0.00,0.25,0.40,0.55,0.70,0.85,1.00], [0.00,0.18,0.38,0.60,0.78,0.91,1.00]);
engine_model.set_throttle_map([0.00,0.25,0.50,0.75,1.00], [0.00,0.25,0.50,0.75,1.00]);
engine_model.set_bsfc_map([0.10,0.40,0.47,0.53,0.59,0.65,0.72,0.76,1.00], [0.42,0.32,0.289,0.279,0.275,0.27435,0.270,0.269,0.290]);
engine_model.idle_fuel_flow_kgps = 1.2*avgas_density_kg_per_usgal/3600;

propeller_scale = 0.84091;
propeller_model = PropellerPerformanceModel( ...
    1.9304, ...
    [0.00,0.20,0.40,0.60,0.80,1.00,1.20,1.40], ...
    propeller_scale*[0.115,0.105,0.088,0.065,0.050,0.033,0.015,0], ...
    propeller_scale*[0.060,0.059,0.057,0.054,0.051,0.048,0.045,0.043], ...
    2,"Sensenich_76EM8S14_0_60");
propeller_model.minimum_rpm = 600;
propeller_model.maximum_rpm = 2700;
propeller_model.set_map_metadata("calibrated provisional coefficient shape", "calibrated_provisional",false, "Ct/Cp shape calibrated to the 6000-ft, 112-KTAS cruise condition.");

propulsion = PropellerPropulsion("O360",aircraft.get_frame("propeller"),[1;0;0], maximum_engine_power_W,0,[]);
propulsion.set_models(engine_model,propeller_model);
propulsion.set_air_velocity_frame(aircraft.get_body_frame());
propulsion.set_rotation_direction(1);
propulsion.rpm_grid_size = 160;

aircraft.add_propulsive_element(propulsion);

engine_component = Component("engine",0,[0;0;0],zeros(3), aircraft.get_frame("propeller"),"engine");
engine_component.add_load_source(PropulsionLoadSolver(propulsion,aircraft.get_frame("propeller")));
aircraft.add_component(engine_component);

%% Gravity and final CG reference

aircraft.ensure_gravity_source();
[~,cg_body_m,~] = aircraft.compute_total_mass_properties([]);
aircraft.update_frame_position("cg",cg_body_m);
aircraft.update_frame_position("gravity_cg",cg_body_m);
aircraft.set_reference_frame("cg");

%% PERFORMANCE ANALYSIS CONFIGURATION
% These values configure PerformanceAnalysis; they do not duplicate its
% aerodynamic, propulsion or performance calculations.

performance = aircraft.get_performance();
condition = struct('altitude_m',reference_altitude_m);

performance_cfg = struct( ...
    'available_throttle',1, ...
    'takeoff_throttle',1, ...
    'continuous_throttle',1, ...
    'propulsion_type',"power", ...
    'range_endurance_metric',"classical", ...
    'CL_max',CL_max_analysis, ...
    'CL_min',CL_min_analysis, ...
    'n_limit_pos',3.8, ...
    'n_limit_neg',-1.52, ...
    'max_load_factor',3.8, ...
    'min_load_factor',-1.52, ...
    'max_Mach',aero_mach_limit, ...
    'service_ceiling_threshold_fpm',100, ...
    'enforce_propulsion_constraints',true);

solver_cfg = struct( ...
    'residual_tolerance',1e-5, ...
    'inequality_tolerance',1e-6, ...
    'equality_tolerance',1e-5, ...
    'optimality_tolerance',1e-5, ...
    'fmincon_options',optimoptions("fmincon", ...
        "Algorithm","sqp", ...
        "Display","none", ...
        "MaxIterations",1500, ...
        "MaxFunctionEvaluations",30000, ...
        "OptimalityTolerance",1e-8, ...
        "ConstraintTolerance",1e-8, ...
        "StepTolerance",1e-10, ...
        "FiniteDifferenceType","central", ...
        "FiniteDifferenceStepSize",1e-4));

fast_solver_cfg = struct( ...
    'residual_tolerance',1e-5, ...
    'inequality_tolerance',1e-6, ...
    'equality_tolerance',1e-5, ...
    'optimality_tolerance',1e-5, ...
    'fmincon_options',optimoptions("fmincon", ...
        "Algorithm","sqp", ...
        "Display","none", ...
        "MaxIterations",500, ...
        "MaxFunctionEvaluations",9000, ...
        "OptimalityTolerance",1e-6, ...
        "ConstraintTolerance",1e-6, ...
        "StepTolerance",1e-8, ...
        "FiniteDifferenceType","forward", ...
        "FiniteDifferenceStepSize",1e-4));

level_spec = struct( ...
    'mode',"level", ...
    'variables',["alpha";"control_pitch";"throttle"], ...
    'initial_guess',[deg2rad(3);deg2rad(-4);0.45], ...
    'lb',[deg2rad(alpha_trim_min_deg);deg2rad(-26);0], ...
    'ub',[deg2rad(alpha_trim_max_deg);deg2rad(28);1], ...
    'fixed',struct('beta',0,'phi',0,'psi',0,'gamma',0, ...
        'control_roll',0,'control_yaw',0,'p',0,'q',0,'r',0), ...
    'reference_frame_name',"cg", ...
    'residual_scale',[aircraft.g;aircraft.g;1]);

[~,speed_of_sound_mps,~,~] = aircraft.get_atmosphere(reference_altitude_m);
speed_lower_mps = 25;
speed_upper_mps = min(80,aero_mach_limit*speed_of_sound_mps);

paper_bounds = struct();
paper_bounds.stall = struct( ...
    'V_lb',20,'V_ub',50,'V0',30, ...
    'alpha0_deg',12,'alpha_min_deg',alpha_trim_min_deg, ...
    'alpha_max_deg',alpha_trim_max_deg, ...
    'elevator0_deg',-5,'elevator_min_deg',-26,'elevator_max_deg',28, ...
    'throttle0',0.35,'throttle_min',0,'throttle_max',1);
paper_bounds.climb = struct( ...
    'V_lb',speed_lower_mps,'V_ub',speed_upper_mps,'V0',42, ...
    'alpha0_deg',5,'alpha_min_deg',alpha_trim_min_deg, ...
    'alpha_max_deg',alpha_trim_max_deg, ...
    'elevator0_deg',-4,'elevator_min_deg',-26,'elevator_max_deg',28, ...
    'gamma0_deg',5,'gamma_min_deg',0,'gamma_max_deg',25, ...
    'ceiling_gamma_min_deg',-10,'ceiling_descent_seed_deg',-3, ...
    'throttle0',1,'throttle_min',0,'throttle_max',1);
paper_bounds.range = struct( ...
    'V_lb',speed_lower_mps,'V_ub',speed_upper_mps,'V0',50, ...
    'alpha0_deg',4,'alpha_min_deg',alpha_trim_min_deg, ...
    'alpha_max_deg',alpha_trim_max_deg, ...
    'elevator0_deg',-4,'elevator_min_deg',-26,'elevator_max_deg',28, ...
    'throttle0',0.45,'throttle_min',0,'throttle_max',1);
paper_bounds.endurance = struct( ...
    'V_lb',speed_lower_mps,'V_ub',speed_upper_mps,'V0',40, ...
    'alpha0_deg',6,'alpha_min_deg',alpha_trim_min_deg, ...
    'alpha_max_deg',alpha_trim_max_deg, ...
    'elevator0_deg',-5,'elevator_min_deg',-26,'elevator_max_deg',28, ...
    'throttle0',0.40,'throttle_min',0,'throttle_max',1);
paper_bounds.glide = struct( ...
    'V_lb',speed_lower_mps,'V_ub',min(65,speed_upper_mps),'V0',50, ...
    'alpha0_deg',4,'alpha_min_deg',alpha_trim_min_deg, ...
    'alpha_max_deg',alpha_trim_max_deg, ...
    'elevator0_deg',-4,'elevator_min_deg',-26,'elevator_max_deg',28, ...
    'gamma0_deg',-4,'gamma_min_deg',-20,'gamma_max_deg',-0.01);
paper_bounds.max_speed = struct( ...
    'V_lb',45,'V_ub',speed_upper_mps,'V0',0.9*speed_upper_mps, ...
    'alpha0_deg',1,'alpha_min_deg',alpha_trim_min_deg, ...
    'alpha_max_deg',alpha_trim_max_deg, ...
    'elevator0_deg',-3,'elevator_min_deg',-26,'elevator_max_deg',28, ...
    'throttle0',1,'throttle_min',0,'throttle_max',1);
paper_bounds.V_sweep_mps = velocity_grid_mps;
paper_bounds.alpha0_deg = 5;
paper_bounds.alpha_min_deg = alpha_trim_min_deg;
paper_bounds.alpha_max_deg = alpha_trim_max_deg;
paper_bounds.elevator0_deg = -4;
paper_bounds.elevator_min_deg = -26;
paper_bounds.elevator_max_deg = 28;
paper_bounds.throttle0 = 0.5;
paper_bounds.throttle_min = 0;
paper_bounds.throttle_max = 1;
paper_bounds.gamma0_deg = 5;
paper_bounds.gamma_min_deg = -20;
paper_bounds.gamma_max_deg = 25;

turn_bounds = struct( ...
    'V_lb',speed_lower_mps,'V_ub',speed_upper_mps,'V0',42, ...
    'bank0_deg',30,'bank_min_deg',0.1,'bank_max_deg',60, ...
    'alpha0_deg',5,'alpha_min_deg',alpha_trim_min_deg, ...
    'alpha_max_deg',alpha_trim_max_deg, ...
    'elevator0_deg',-4,'elevator_min_deg',-26,'elevator_max_deg',28, ...
    'aileron0_deg',0,'aileron_min_deg',-20,'aileron_max_deg',20, ...
    'rudder0_deg',0,'rudder_min_deg',-30,'rudder_max_deg',30, ...
    'throttle0',0.65,'throttle_min',0,'throttle_max',1, ...
    'turn_rate_ub_radps',0.8,'Ps_trim_tolerance_mps',0.05);

envelope_bounds = struct( ...
    'V_lb',20,'V_ub',speed_upper_mps,'V0',0.9*speed_upper_mps, ...
    'max_Mach',aero_mach_limit,'stall_speed_factor',1, ...
    'alpha0_deg',3,'alpha_min_deg',alpha_trim_min_deg, ...
    'alpha_max_deg',alpha_trim_max_deg, ...
    'elevator0_deg',-4,'elevator_min_deg',-26,'elevator_max_deg',28);

%% BASELINE CRUISE TRIM
% The trim point and all reported quantities come from PerformanceAnalysis.

baseline_condition = struct('altitude_m',reference_altitude_m, 'velocity_mps',reference_velocity_mps);

[baseline_x,baseline_u,baseline_converged,baseline_info] = performance.solve_trim(baseline_condition,level_spec,solver_cfg);
baseline = performance.evaluate_trim_point(baseline_x,baseline_u,baseline_condition, performance_cfg,baseline_info.extras);
baseline.trim_converged = baseline_converged;
baseline.trim_info = baseline_info;

%% PAPER-ALIGNED AIRBORNE PERFORMANCE
% VS, Vx, Vy, best range, best endurance, best glide, Vh and curves.

paper_suite = performance.analyze_paper_trim_suite(condition,solver_cfg,performance_cfg,paper_bounds);

%% GLIDE PERFORMANCE
% Adds minimum sink and glide range/time to the paper-suite best glide.

glide = performance.analyze_glide_performance(condition,solver_cfg,performance_cfg,paper_bounds.glide);

%% MANEUVER ENVELOPE AND TURN PERFORMANCE

vn = performance.compute_vn_diagram(condition,envelope_velocity_grid_mps,performance_cfg);
instantaneous_turn = performance.compute_instantaneous_turn_curve(condition,turn_velocity_grid_mps,performance_cfg);
sustained_turn = performance.compute_sustained_turn_curve(condition,turn_velocity_grid_mps,fast_solver_cfg, performance_cfg,turn_bounds);
sustained_turn_optimum = performance.optimize_sustained_turn(condition,solver_cfg,performance_cfg,turn_bounds);

%% ENERGY-MANEUVERABILITY AND ACCELERATION

Ps_map = performance.compute_specific_excess_power_map(condition,Ps_velocity_grid_mps,Ps_load_factors, fast_solver_cfg,performance_cfg,turn_bounds);
acceleration = performance.compute_acceleration_schedule(condition,acceleration_velocity_grid_mps,fast_solver_cfg, performance_cfg,envelope_bounds);

%% ALTITUDE ENVELOPE, CLIMB SCHEDULE AND CEILINGS

climb = performance.optimize_climb_schedule(condition,altitude_grid_m,solver_cfg, performance_cfg,paper_bounds.climb);
mach_altitude_envelope = performance.compute_mach_altitude_envelope(condition,altitude_grid_m,solver_cfg, performance_cfg,envelope_bounds);

%% COLLECT COMPLETE FRAMEWORK OUTPUTS
% point_opt contains x_trim, u_trim, controls, forces, moments and metrics.

C172Performance = struct();
C172Performance.aircraft = aircraft;
C172Performance.propulsion = propulsion;
C172Performance.configuration = struct('condition',condition, 'performance',performance_cfg, 'paper_bounds',paper_bounds, 'turn_bounds',turn_bounds, 'envelope_bounds',envelope_bounds);
C172Performance.baseline = baseline;
C172Performance.paper_suite = paper_suite;
C172Performance.glide = glide;
C172Performance.vn = vn;
C172Performance.instantaneous_turn = instantaneous_turn;
C172Performance.sustained_turn = sustained_turn;
C172Performance.sustained_turn_optimum = sustained_turn_optimum;
C172Performance.Ps_map = Ps_map;
C172Performance.acceleration = acceleration;
C172Performance.climb = climb;
C172Performance.mach_altitude_envelope = mach_altitude_envelope;

%% PERFORMANCE RESULTS

fprintf('\n=== C172 BASELINE TRIM ===\n');
performance.print_point(baseline);

fprintf('\n=== C172 CHARACTERISTIC PERFORMANCE POINTS ===\n');
performance.print_point(paper_suite.stall_speed.point_opt);
performance.print_point(paper_suite.best_angle_climb.point_opt);
performance.print_point(paper_suite.best_rate_climb.point_opt);
performance.print_point(paper_suite.best_range.point_opt);
performance.print_point(paper_suite.best_endurance.point_opt);
performance.print_point(paper_suite.best_glide.point_opt);
performance.print_point(paper_suite.maximum_level_speed.point_opt);

%% PERFORMANCE PLOTS
% Plots use fields already returned by PerformanceAnalysis.

plotC172Performance(C172Performance);

%% LOCAL PLOTTING FUNCTION

function plotC172Performance(results)

curve = results.paper_suite.required_available_curve;
required = pointSeries(curve.required.point,curve.required.valid);
instantaneous = results.instantaneous_turn;
sustained = results.sustained_turn;
vn = results.vn;
Ps_map = results.Ps_map;
acceleration = results.acceleration;
climb = results.climb;
mach_envelope = results.mach_altitude_envelope;

%% Aerodynamics and trim controls

figure("Name","C172 Aerodynamics and Trim");
tiledlayout(2,3,"TileSpacing","compact");
nexttile;
plot(required.V_kts,required.L_N/1000,"LineWidth",1.5);
hold on;
yline(required.weight_N(1)/1000,"--");
grid on;
xlabel("True airspeed [kt]");
ylabel("Force [kN]");
legend("Lift","Weight","Location","best");
nexttile;
plot(required.V_kts,required.D_N,"LineWidth",1.5);
grid on;
xlabel("True airspeed [kt]");
ylabel("Drag [N]");
nexttile;
plot(required.V_kts,required.CL,"LineWidth",1.5);
hold on;
plot(required.V_kts,required.CD,"LineWidth",1.5);
grid on;
xlabel("True airspeed [kt]");
ylabel("Coefficient");
legend("C_L","C_D","Location","best");
nexttile;
plot(required.V_kts,required.L_over_D,"LineWidth",1.5);
grid on;
xlabel("True airspeed [kt]");
ylabel("L/D");
nexttile;
plot(required.V_kts,required.alpha_deg,"LineWidth",1.5);
hold on;
plot(required.V_kts,required.theta_deg,"LineWidth",1.5);
grid on;
xlabel("True airspeed [kt]");
ylabel("Angle [deg]");
legend("Alpha","Theta","Location","best");
nexttile;
plot(required.V_kts,required.elevator_deg,"LineWidth",1.5);
hold on;
plot(required.V_kts,required.throttle,"LineWidth",1.5);
grid on;
xlabel("True airspeed [kt]");
ylabel("Control");
legend("Elevator [deg]","Throttle","Location","best");

%% Required and available propulsion performance

figure("Name","C172 Required and Available Performance");
tiledlayout(2,2,"TileSpacing","compact");
nexttile;
plot(curve.V_kts,curve.required.thrust_N/1000,"LineWidth",1.5);
hold on;
plot(curve.V_kts,curve.available.thrust_N/1000,"LineWidth",1.5);
grid on;
xlabel("True airspeed [kt]");
ylabel("Thrust [kN]");
legend("Required","Available","Location","best");
nexttile;
plot(curve.V_kts,curve.required.power_W/1000,"LineWidth",1.5);
hold on;
plot(curve.V_kts,curve.available.power_W/1000,"LineWidth",1.5);
grid on;
xlabel("True airspeed [kt]");
ylabel("Thrust power [kW]");
legend("Required","Available","Location","best");
nexttile;
plot(curve.V_kts,curve.excess_thrust_N,"LineWidth",1.5);
yline(0,"--");
grid on;
xlabel("True airspeed [kt]");
ylabel("Excess thrust [N]");
nexttile;
plot(curve.V_kts,curve.excess_power_W/1000,"LineWidth",1.5);
yline(0,"--");
grid on;
xlabel("True airspeed [kt]");
ylabel("Excess thrust power [kW]");

%% Range, endurance and fuel flow

figure("Name","C172 Range and Endurance");
tiledlayout(1,3,"TileSpacing","compact");
nexttile;
plot(required.V_kts,required.fuel_specific_range_m_per_kg/1852, "LineWidth",1.5);
grid on;
xlabel("True airspeed [kt]");
ylabel("Specific range [nm/kg]");
nexttile;
plot(required.V_kts,required.fuel_endurance_s_per_kg/3600, "LineWidth",1.5);
grid on;
xlabel("True airspeed [kt]");
ylabel("Endurance [hr/kg]");
nexttile;
plot(required.V_kts,required.fuel_flow_trim*3600,"LineWidth",1.5);
grid on;
xlabel("True airspeed [kt]");
ylabel("Fuel flow [kg/hr]");

%% Maneuver envelope and turns

figure("Name","C172 Maneuver Envelope");
fill([vn.V_EAS_kts;flipud(vn.V_EAS_kts)], [vn.n_pos;flipud(vn.n_neg)],[0.9 0.9 0.9],"EdgeColor","none");
hold on;
plot(vn.V_EAS_kts,vn.n_pos,"LineWidth",2);
plot(vn.V_EAS_kts,vn.n_neg,"LineWidth",2);
yline(1,"--");
grid on;
xlabel("Equivalent airspeed [kt]");
ylabel("Load factor");
legend("Envelope","Positive limit","Negative limit","Location","best");

figure("Name","C172 Turn Performance");
tiledlayout(1,2,"TileSpacing","compact");
iv = instantaneous.valid;
sv = sustained.valid;
nexttile;
plot(instantaneous.V_kts(iv),instantaneous.turn_rate_degps(iv), "LineWidth",1.5);
hold on;
plot(sustained.V_kts(sv),sustained.turn_rate_degps(sv), "o-","LineWidth",1.5);
grid on;
xlabel("True airspeed [kt]");
ylabel("Turn rate [deg/s]");
legend("Instantaneous","Sustained","Location","best");
nexttile;
plot(instantaneous.V_kts(iv),instantaneous.turn_radius_m(iv), "LineWidth",1.5);
hold on;
plot(sustained.V_kts(sv),sustained.turn_radius_m(sv), "o-","LineWidth",1.5);
grid on;
xlabel("True airspeed [kt]");
ylabel("Turn radius [m]");
legend("Instantaneous","Sustained","Location","best");

%% Specific excess power and acceleration

figure("Name","C172 Specific Excess Power");
contourf(Ps_map.V_kts,Ps_map.n,Ps_map.Ps_mps,18);
hold on;
contour(Ps_map.V_kts,Ps_map.n,Ps_map.Ps_mps,[0 0], "LineColor","k","LineWidth",2);
colorbar;
grid on;
xlabel("True airspeed [kt]");
ylabel("Load factor");
title("P_s [m/s]");

figure("Name","C172 Acceleration");
tiledlayout(1,2,"TileSpacing","compact");
nexttile;
plot(acceleration.V_kts,acceleration.acceleration_mps2,"LineWidth",1.5);
yline(0,"--");
grid on;
xlabel("True airspeed [kt]");
ylabel("Acceleration [m/s^2]");
nexttile;
plot(acceleration.V_kts,acceleration.Ps_mps,"LineWidth",1.5);
yline(0,"--");
grid on;
xlabel("True airspeed [kt]");
ylabel("P_s [m/s]");

%% Altitude envelope, climb and ceilings

figure("Name","C172 Altitude Performance");
tiledlayout(2,2,"TileSpacing","compact");
nexttile;
plot(mach_envelope.stall_Mach,mach_envelope.altitude_ft,"LineWidth",1.5);
hold on;
plot(mach_envelope.max_level_Mach,mach_envelope.altitude_ft, "LineWidth",1.5);
grid on;
xlabel("Mach");
ylabel("Altitude [ft]");
legend("Stall","Maximum level","Location","best");
nexttile;
plot(climb.best_ROC_fpm,climb.altitude_ft,"o-","LineWidth",1.5);
xline(0,"--");
grid on;
xlabel("Best ROC [ft/min]");
ylabel("Altitude [ft]");
nexttile;
plot(climb.best_V_mps*1.943844492,climb.altitude_ft, "o-","LineWidth",1.5);
grid on;
xlabel("Best climb speed [kt]");
ylabel("Altitude [ft]");
nexttile;
plot(climb.best_gamma_deg,climb.altitude_ft,"o-","LineWidth",1.5);
xline(0,"--");
grid on;
xlabel("Climb angle [deg]");
ylabel("Altitude [ft]");

end

function data = pointSeries(point_cells,valid)

indices = find(valid);
if isempty(indices)
    error('C172PerformanceAnalysis:NoValidCurvePoints', 'The performance curve contains no valid points.');
end
points = [point_cells{indices}];

fields = ["V_kts";"weight_N";"L_N";"D_N";"CL";"CD";"L_over_D"; "alpha_deg";"theta_deg";"elevator_deg";"throttle"; "fuel_specific_range_m_per_kg";"fuel_endurance_s_per_kg"; "fuel_flow_trim"];

for k = 1:numel(fields)
    name = fields(k);
    data.(name) = reshape([points.(name)],[],1);
end

end
