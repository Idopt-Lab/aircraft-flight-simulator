if bdIsLoaded('flightsim7')
    set_param('flightsim7','SimulationCommand','stop');
end

clear; clc; close all;

%% PATH BOOTSTRAP
% Ensures the shared class library and example helper functions are on the
% MATLAB path regardless of the current working directory or session state.
this_script_dir = fileparts(mfilename('fullpath'));
matlab_root_dir = fileparts(fileparts(this_script_dir));
addpath(genpath(fullfile(matlab_root_dir,'classes')));
addpath(genpath(fullfile(matlab_root_dir,'examples')));

%% User settings

sim_time = 100.0;
dt_fc = 0.01;
disable_rate_limiting = true;

data_directory = string(fileparts(mfilename('fullpath')));
datcom_files = fullfile(data_directory,["cessna010.out"; "cessna012.out"; "cessna014.out"; "cessna016.out"; "cessna018.out"; "cessna020.out"; "cessna022.out"; "cessna024.out"]);

datcom_options = struct('CD0_add',0.01681, 'CD_alpha_extra',0, 'CL_max_abs',1.9);

avgas_density_kg_per_usgal = 2.72;

%% Aircraft object

ac = Aircraft();

%% Geometry / reference positions

CG_FROM_NOSE = 1.90;

cg_to_nose    = -CG_FROM_NOSE;
cg_to_wing_le = 1.20 - CG_FROM_NOSE;
cg_to_wing_ac = 1.57 - CG_FROM_NOSE;
cg_to_prop    = 1.68 - CG_FROM_NOSE;
cg_to_fuel    = 1.70 - CG_FROM_NOSE;

S_ref = 16.17;
b_ref = 10.98;
c_ref = 1.50;

ac.geometry.set_reference_geometry(S_ref,b_ref,c_ref);
ac.geometry.set_reference_point([cg_to_wing_ac 0 0]);

% Keep aliases populated for older code paths / lookup functions.
ac.geometry.wing_area = S_ref;
ac.geometry.wing_span = b_ref;
ac.geometry.mean_aerodynamic_chord = c_ref;
ac.geometry.ref_area = S_ref;
ac.geometry.ref_span = b_ref;
ac.geometry.ref_chord = c_ref;
ac.geometry.ref_point = [cg_to_wing_ac 0 0];

%% Frames

ac.set_body_frame("body");
ac.set_reference_frame("body");

add_or_update_frame(ac,"aero_ref",   "body",[cg_to_wing_ac;0;0], @(x) eye(3));
add_or_update_frame(ac,"fuel_tank",  "body",[cg_to_fuel;0;0],    @(x) eye(3));
add_or_update_frame(ac,"propeller",  "body",[cg_to_prop;0;1.26], @(x) eye(3));
add_or_update_frame(ac,"cg",         "body",[0;0;0],             @(x) eye(3));
add_or_update_frame(ac,"gravity_cg", "body",[0;0;0],             @(x) ned_to_body_dcm(x));

%% Mass properties

empty_mass = 767;
fuel_mass  = 90;

Ixx = 1285;
Iyy = 1824;
Izz = 2666;

I_body = [Ixx 0 0;
          0 Iyy 0;
          0 0 Izz];

airframe = Component("airframe",empty_mass,[0;0;0],I_body,ac.get_frame("body"),"airframe");

% cg_local is relative to fuel_tank frame.
% fuel_tank frame is already located at tank CG, so cg_local = zero.
fuel = Component("fuel",fuel_mass,[0;0;0],zeros(3,3),ac.get_frame("fuel_tank"),"fuel");

ac.add_component(airframe);
ac.add_component(fuel);

%% Controls

aileron = ControlSurface("aileron","aileron","primary", [1 0 0], deg2rad(20),deg2rad(-20), 0,0,0);

elevator = ControlSurface("elevator","elevator","primary", [0 1 0], deg2rad(28),deg2rad(-26), 0,0,0);

rudder = ControlSurface("rudder","rudder","primary", [0 0 1], deg2rad(30),deg2rad(-30), 0,0,0);

ac.add_control_surface(aileron);
ac.add_control_surface(elevator);
ac.add_control_surface(rudder);

%% Aerodynamics

[datcom_lookup,datcom_data] = c172datcom(datcom_files,true,datcom_options);
assignin("base","datcom_data",datcom_data);

aero_model = CoefficientAerodynamics(datcom_lookup);

aero_solver = AeroLoadSolver(aero_model, ac.geometry, ac, ac.get_frame("aero_ref"));

ac.add_load_source(aero_solver);

%% Propulsion

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

prop_pe = PropellerPropulsion("O360", ac.get_frame("propeller"), [1;0;0], maximum_engine_power_W, 0, []);
prop_pe.set_models(engine_model,propeller_model);
prop_pe.set_air_velocity_frame(ac.get_body_frame());
prop_pe.set_rotation_direction(1);
prop_pe.rpm_grid_size = 160;

ac.add_propulsive_element(prop_pe);

engine_c = Component("engine", 0, [0;0;0], zeros(3,3), ac.get_frame("propeller"), "engine");

prop_solver = PropulsionLoadSolver(prop_pe,ac.get_frame("propeller"));
engine_c.add_load_source(prop_solver);

ac.add_component(engine_c);

%% Gravity load

gravity_solver = GravityLoadSolver(ac,ac.get_frame("gravity_cg"));
ac.add_load_source(gravity_solver);

%% Sync mass / CG

n_cs = numel(ac.control_surfaces);
n_pe = numel(ac.propulsive_elements);
n_total = n_cs + n_pe;

[m,cg,I_total] = ac.compute_total_mass_properties([]);

ac.update_frame_position("cg",cg);
ac.update_frame_position("gravity_cg",cg);

W = m * ac.g;
cbar = c_ref;

fprintf('\n=== CESSNA SETUP ===\n');
fprintf('Reference frame : %s\n', ac.reference_frame_name);
fprintf('Empty mass      : %.1f kg\n', empty_mass);
fprintf('Fuel mass       : %.1f kg\n', fuel.mass);
fprintf('Total mass      : %.1f kg\n', m);
fprintf('Weight          : %.1f N\n', W);
fprintf('CG position     : [% .6f % .6f % .6f] m\n', cg);
fprintf('Nose            : [% .2f 0 0] m from CG\n', cg_to_nose);
fprintf('Wing LE         : [% .2f 0 0] m from CG\n', cg_to_wing_le);
fprintf('Wing AC         : [% .2f 0 0] m from CG\n', cg_to_wing_ac);
fprintf('Propeller       : [% .2f 0 1.26] m from CG\n', cg_to_prop);
fprintf('Fuel tank       : [% .2f 0 0] m from CG\n', cg_to_fuel);

%% Test point

x_test = zeros(12,1);
x_test(3) = -1000;
x_test(4) = 50*cosd(5);
x_test(5) = 0;
x_test(6) = 50*sind(5);
x_test(7) = 0;
x_test(8) = deg2rad(5);
x_test(9) = 0;
x_test(10:12) = 0;

u_test = zeros(n_total,1);
u_test(1) = 0;
u_test(2) = 0;
u_test(3) = 0;
u_test(n_cs+1) = 0.65;

ac.state.set_full_state(x_test);
ac.set_controls_from_vector(u_test);

[F_cg,M_cg] = ac.compute_total_loads_about_cg(x_test,u_test);

fprintf('\n=== TEST POINT ===\n');
fprintf('At alpha = %.1f deg, elev = 0 deg, thr = 0.65\n', rad2deg(atan2(x_test(6),x_test(4))));
fprintf('F(CG)    : [% .1f % .1f % .1f] N\n', F_cg);
fprintf('M(CG)    : [% .1f % .1f % .1f] N-m\n', M_cg);
fprintf('Fx/W     : %.4f\n', F_cg(1)/W);
fprintf('Fz/W     : %.4f\n', F_cg(3)/W);
fprintf('My/(Wc)  : %.4f\n', M_cg(2)/(W*cbar));

%% Trim setup

alt_trim_m = 1500;
mach_trim = 0.15;

[~,a_trim,~,~] = ac.get_atmosphere(alt_trim_m);
V_trim = mach_trim * a_trim;

condition = struct();
condition.altitude_m = alt_trim_m;
condition.mach = mach_trim;

trim_cfg = struct();
trim_cfg.variables = ["alpha","elevator","throttle"];
trim_cfg.residuals = ["Fx","Fz","My"];

trim_cfg.initial_guess = [deg2rad(3); deg2rad(0); 0.45];

trim_cfg.lb = [deg2rad(-8); deg2rad(-26); 0.0];
trim_cfg.ub = [deg2rad(15); deg2rad(28); 1.0];

trim_cfg.weights = [1; 1; 1];

trim_cfg.fixed = struct();
trim_cfg.fixed.beta = 0;
trim_cfg.fixed.phi = 0;
trim_cfg.fixed.psi = 0;
trim_cfg.fixed.gamma = 0;
trim_cfg.fixed.aileron = 0;
trim_cfg.fixed.rudder = 0;
trim_cfg.fixed.p = 0;
trim_cfg.fixed.q = 0;
trim_cfg.fixed.r = 0;

trim_cfg.reference_frame_name = "cg";

solver_cfg = struct();
solver_cfg.residual_tolerance = 1e-5;
solver_cfg.fmincon_options = optimoptions( ...
    "fmincon", ...
    "Algorithm","sqp", ...
    "Display","iter", ...
    "MaxIterations",5000, ...
    "MaxFunctionEvaluations",50000, ...
    "OptimalityTolerance",1e-10, ...
    "StepTolerance",1e-10, ...
    "ConstraintTolerance",1e-10, ...
    "FunctionTolerance",1e-10, ...
    "FiniteDifferenceStepSize",1e-8);

%% Run trim

fprintf('\n=== RUNNING CESSNA CRUISE TRIM ===\n');
fprintf('Altitude : %.1f m\n', alt_trim_m);
fprintf('Mach     : %.3f\n', mach_trim);
fprintf('V        : %.3f m/s\n', V_trim);

ac.set_reference_frame("cg");

solver = ac.get_trim_solver();
[x_trim,u_trim,converged,info] = solver.solve_trim(condition,trim_cfg,solver_cfg); %#ok<NASGU>

if converged
    fprintf('\n=== TRIM CONVERGED ===\n');
else
    fprintf('\n=== TRIM NOT CONVERGED: using returned best solution ===\n');
end

solver.print_summary();

%% Final trim verification about CG

ac.state.set_full_state(x_trim);
ac.set_controls_from_vector(u_trim);

[m_trim,cg_trim,I_trim] = ac.compute_total_mass_properties(x_trim);
ac.update_frame_position("cg",cg_trim);
ac.update_frame_position("gravity_cg",cg_trim);
ac.set_reference_frame("cg");

[F_trim,M_trim,fuel_flow_trim] = ac.compute_total_loads(x_trim,u_trim);

W_trim = m_trim * ac.g;
cbar_trim = c_ref;

fprintf('\n=== TRIM CHECK ABOUT CG ===\n');
fprintf('Mass     : %.2f kg\n', m_trim);
fprintf('CG       : [% .4f % .4f % .4f] m\n', cg_trim);
fprintf('Alpha    : %.4f deg\n', rad2deg(atan2(x_trim(6),x_trim(4))));
fprintf('Theta    : %.4f deg\n', rad2deg(x_trim(8)));
fprintf('Aileron  : %.4f deg\n', rad2deg(u_trim(1)));
fprintf('Elevator : %.4f deg\n', rad2deg(u_trim(2)));
fprintf('Rudder   : %.4f deg\n', rad2deg(u_trim(3)));
fprintf('Throttle : %.6f\n', u_trim(n_cs+1));
fprintf('F        : [% .3e % .3e % .3e] N\n', F_trim);
fprintf('M        : [% .3e % .3e % .3e] N-m\n', M_trim);
fprintf('Fx/W     : %.6e\n', F_trim(1)/W_trim);
fprintf('Fz/W     : %.6e\n', F_trim(3)/W_trim);
fprintf('My/(Wc)  : %.6e\n', M_trim(2)/(W_trim*cbar_trim));
fprintf('Fuel flow: %.6e kg/s\n', fuel_flow_trim);

%% Load contributions using new solver path

fprintf('\n=== LOAD CONTRIBUTIONS ABOUT CURRENT REF: %s ===\n', ac.reference_frame_name);

body_frame = ac.get_body_frame();
ref_frame = ac.get_reference_frame();

[F_aero,M_aero,~,src_aero] = ac.compute_loads_by_solver(x_trim,u_trim,"AeroLoadSolver",body_frame,ref_frame);

[F_prop,M_prop,ff_prop,src_prop] = ac.compute_loads_by_solver(x_trim,u_trim,"PropulsionLoadSolver",body_frame,ref_frame);

[F_grav,M_grav,~,src_grav] = ac.compute_loads_by_solver(x_trim,u_trim,"GravityLoadSolver",body_frame,ref_frame);

fprintf('Aero sources      : %d\n', numel(src_aero));
fprintf('Aero              F=[% .3e % .3e % .3e] M=[% .3e % .3e % .3e]\n', F_aero, M_aero);

fprintf('Propulsion sources: %d\n', numel(src_prop));
fprintf('Propulsion        F=[% .3e % .3e % .3e] M=[% .3e % .3e % .3e] fuel=%.6e\n', F_prop, M_prop, ff_prop);

fprintf('Gravity sources   : %d\n', numel(src_grav));
fprintf('Gravity           F=[% .3e % .3e % .3e] M=[% .3e % .3e % .3e]\n', F_grav, M_grav);

fprintf('Total             F=[% .3e % .3e % .3e] M=[% .3e % .3e % .3e]\n', F_trim, M_trim);

%% Direct aero / prop diagnostic checks

fprintf('\n=== AERO DIRECT CHECK ===\n');
[F_aero_local,M_aero_local,coeff] = aero_model.get_FM(x_trim,u_trim,ac.geometry,ac);
fprintf('F aero local/wind [N]   : [% .3e % .3e % .3e]\n', F_aero_local);
fprintf('M aero local/wind [N-m] : [% .3e % .3e % .3e]\n', M_aero_local);
disp(coeff);

fprintf('\n=== PROPULSION DIRECT CHECK ===\n');
[F_prop_local,M_prop_local,fuel_flow_direct] = ac.propulsive_elements{1}.get_FM(x_trim,u_trim);
fprintf('F prop local direct [N]  : [% .3e % .3e % .3e]\n', F_prop_local);
fprintf('M prop local direct [N-m]: [% .3e % .3e % .3e]\n', M_prop_local);
fprintf('Fuel flow direct         : %.6e\n', fuel_flow_direct);

%% Performance analysis

fprintf('\n=== RUNNING PERFORMANCE ANALYSIS ===\n');

perf = ac.get_performance();

perf_cfg = struct('available_throttle',1, 'propulsion_type',"power", 'reference_frame_name',"cg");

perf_point = perf.evaluate_trim_point(x_trim,u_trim,condition,perf_cfg);
perf_point.trim_converged = converged;

perf.print_point(perf_point);

%% Stability analysis

fprintf('\n=== RUNNING STABILITY ANALYSIS ===\n');

stab = ac.get_stability();
stab.set_trim(x_trim,u_trim);

dx_lin = 1e-5;
du_lin = 1e-5;

[A,B] = stab.linearize(dx_lin,du_lin);

stab.analyze_modes();
stab.print_modes();

alpha_range = deg2rad(-5:0.5:15);
beta_range  = deg2rad(-10:0.5:10);

[alpha_vec,Cm_vec,Cm_alpha] = stab.compute_pitch_static_stability(alpha_range);
[beta_vec,Cn_vec,Cn_beta]   = stab.compute_yaw_static_stability(beta_range);

fprintf('\n=== STATIC STABILITY ===\n');
fprintf('Cm_alpha : %.6f 1/rad\n', Cm_alpha);
fprintf('Cn_beta  : %.6f 1/rad\n', Cn_beta);

if Cm_alpha < 0
    fprintf('Pitch static stability       : stable tendency\n');
else
    fprintf('Pitch static stability       : unstable/neutral tendency\n');
end

if Cn_beta > 0
    fprintf('Directional static stability : stable tendency\n');
else
    fprintf('Directional static stability : unstable/neutral tendency\n');
end

%% Takeoff analysis

fprintf('\n=== RUNNING TAKEOFF ANALYSIS ===\n');

to = ac.get_takeoff();

takeoff_altitude_m = 0;
takeoff_runway_slope_deg = 0;
takeoff_runway_length_m = 1000;

[TO_m,to_res] = to.calculate_takeoff(takeoff_altitude_m, takeoff_runway_slope_deg, takeoff_runway_length_m);

to.print_takeoff_summary(to_res);

%% Landing analysis

fprintf('\n=== RUNNING LANDING ANALYSIS ===\n');

ld = ac.get_landing();

landing_altitude_m = 0;
landing_temp_offset_K = 0;
landing_runway_length_m = 1000;

[LD_m,ld_res] = ld.calculate_landing(landing_altitude_m, landing_temp_offset_K, landing_runway_length_m);

ld.print_landing_summary(ld_res);

%% CG migration test

fprintf('\n=== CG MIGRATION TEST ===\n');

fuel_levels = [180,142,110,85,50,20,0];

for f = fuel_levels
    fuel.set_mass_properties(f,[0;0;0],zeros(3,3));

    [~,cg_current,~] = ac.compute_total_mass_properties([]);
    cg_current = cg_current(:);

    ac.update_frame_position("cg",cg_current);
    ac.update_frame_position("gravity_cg",cg_current);

    fprintf('Fuel %3.0f kg -> CG = [%+.4f %+.4f %+.4f] m\n', f,cg_current(1),cg_current(2),cg_current(3));
end

fuel.set_mass_properties(fuel_mass,[0;0;0],zeros(3,3));

[m_final,cg_final,~] = ac.compute_total_mass_properties([]);
ac.update_frame_position("cg",cg_final);
ac.update_frame_position("gravity_cg",cg_final);

%% Simulink variables

t_vec = (0:dt_fc:sim_time).';

x0 = x_trim(:);
u0 = u_trim(:);

Initialpos = x0(1:3);
InitialVel = x0(4:6);
InitialOri = x0(7:9);
InitialRot = x0(10:12);

u_trim_mat = repmat(u0.',numel(t_vec),1);

control_input_data = struct();
control_input_data.time = t_vec;
control_input_data.signals.values = u_trim_mat;
control_input_data.signals.dimensions = n_total;

ground_k = ac.ground_k;
ground_c = ac.ground_c;

assignin("base","ac",ac);
assignin("base","n_cs",n_cs);
assignin("base","n_pe",n_pe);
assignin("base","n_total",n_total);

assignin("base","initial_state",x0);
assignin("base","Initialpos",Initialpos);
assignin("base","InitialVel",InitialVel);
assignin("base","InitialOri",InitialOri);
assignin("base","InitialRot",InitialRot);

assignin("base","u_trim",u0);
assignin("base","control_input_data",control_input_data);

assignin("base","sim_stop_time",sim_time);
assignin("base","dt_fc",dt_fc);
assignin("base","disable_rate_limiting",disable_rate_limiting);
assignin("base","ground_k",ground_k);
assignin("base","ground_c",ground_c);

assignin("base","perf",perf);
assignin("base","perf_point",perf_point);

assignin("base","stab",stab);
assignin("base","A_matrix",A);
assignin("base","B_matrix",B);
assignin("base","eigenvalues",stab.eigenvalues);
assignin("base","alpha_static_vec",alpha_vec);
assignin("base","Cm_static_vec",Cm_vec);
assignin("base","Cm_alpha",Cm_alpha);
assignin("base","beta_static_vec",beta_vec);
assignin("base","Cn_static_vec",Cn_vec);
assignin("base","Cn_beta",Cn_beta);

assignin("base","to",to);
assignin("base","TO_m",TO_m);
assignin("base","takeoff_results",to_res);

assignin("base","ld",ld);
assignin("base","LD_m",LD_m);
assignin("base","landing_results",ld_res);

fprintf('\n=== SIMULINK VARIABLES READY ===\n');
fprintf('n_cs             = %d\n', n_cs);
fprintf('n_pe             = %d\n', n_pe);
fprintf('n_total          = %d\n', n_total);
fprintf('sim_stop_time    = %.2f s\n', sim_time);
fprintf('dt_fc            = %.4f s\n', dt_fc);
fprintf('rate limiting off= %d\n', disable_rate_limiting);
fprintf('x0 alpha         = %.4f deg\n', rad2deg(atan2(x0(6),x0(4))));
fprintf('x0 theta         = %.4f deg\n', rad2deg(x0(8)));
fprintf('u_trim elevator  = %.4f deg\n', rad2deg(u0(2)));
fprintf('u_trim throttle  = %.6f\n', u0(n_cs+1));

%% Local functions

function add_or_update_frame(ac,name,parent_name,r_parent,dcm_fn)
    if ac.has_frame(name)
        ac.update_frame_position(name,r_parent);
        ac.update_frame_orientation(name,dcm_fn);
    else
        ac.add_frame(name,parent_name,r_parent,dcm_fn);
    end
end

function C = ned_to_body_dcm(x)
    phi = x(7);
    theta = x(8);
    psi = x(9);

    cp = cos(phi);
    sp = sin(phi);
    ct = cos(theta);
    st = sin(theta);
    cs = cos(psi);
    ss = sin(psi);

    C = [ct*cs,              ct*ss,              -st;
         sp*st*cs-cp*ss,     sp*st*ss+cp*cs,     sp*ct;
         cp*st*cs+sp*ss,     cp*st*ss-sp*cs,     cp*ct];
end