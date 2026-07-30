clear; clc

%% Aircraft object

ac = Aircraft();

%% Geometry
% Direct property assignment (not set_reference_geometry) since that
% setter requires strictly positive values -- fill in real numbers below.

ac.geometry.wing_area = 0;
ac.geometry.wing_span = 0;
ac.geometry.mean_aerodynamic_chord = 0;
ac.geometry.ref_area = 0;
ac.geometry.ref_span = 0;
ac.geometry.ref_chord = 0;
ac.geometry.set_reference_point([0 0 0]);

cfg = ac.get_configurator();

%% Mass properties

cfg.add_component('name',"airframe",'type',"airframe", ...
    'mass',0,'position',[0 0 0],'inertia',zeros(3));

%% Aerodynamics
% Replace @YourLookup with your own coefficient-lookup function handle,
% e.g. one built by ExampleLookup.m or DATCOMLookup.m.

cfg.add_aero_solver(CoefficientAerodynamics(@YourLookup),'name',"aero");

%% Controls

cfg.add_control_surface('name',"aileron",'surface_type',"aileron",'classification',"primary",'axis',[1 0 0], ...
    'max_deflection',0,'min_deflection',0,'dCl',0,'dCm',0,'dCn',0);
cfg.add_control_surface('name',"elevator",'surface_type',"elevator",'classification',"primary",'axis',[0 1 0], ...
    'max_deflection',0,'min_deflection',0,'dCl',0,'dCm',0,'dCn',0);
cfg.add_control_surface('name',"rudder",'surface_type',"rudder",'classification',"primary",'axis',[0 0 1], ...
    'max_deflection',0,'min_deflection',0,'dCl',0,'dCm',0,'dCn',0);

%% Propulsion

cfg.add_propulsive_element('name',"engine",'element_type',"generic",'max_output',0, ...
    'position',[0 0 0],'direction',[1 0 0],'fuel_rate',0,'thrust_model',[]);

n_cs = numel(ac.control_surfaces);
n_pe = numel(ac.propulsive_elements);
n_total = n_cs + n_pe;

%% Trim setup

alt_cruise  = 0;
mach_cruise = 0;
dur_cruise  = 0;

[~, a, ~, ~] = ac.get_atmosphere(alt_cruise);
V_cruise = mach_cruise * a;

condition = struct();
condition.altitude_m = alt_cruise;
condition.mach = mach_cruise;

% Fill in real variable bounds/initial guess before running -- lb = ub = 0
% below is a degenerate (zero-width) trim search and will not converge.
trim_cfg = struct();
trim_cfg.variables = ["alpha","elevator","throttle"];
trim_cfg.residuals = ["Fx","Fz","My"];
trim_cfg.initial_guess = [0; 0; 0];
trim_cfg.lb = [0; 0; 0];
trim_cfg.ub = [0; 0; 0];
trim_cfg.weights = [1; 1; 1];
trim_cfg.fixed = struct('beta',0,'phi',0,'psi',0,'gamma',0, ...
    'aileron',0,'rudder',0,'p',0,'q',0,'r',0);
trim_cfg.reference_frame_name = "body";

solver_cfg = struct();
solver_cfg.residual_tolerance = 1e-5;
solver_cfg.fmincon_options = optimoptions("fmincon","Display","none");

%% Run trim

solver = ac.get_trim_solver();
[x_trim, u_trim, converged] = solver.solve_trim(condition,trim_cfg,solver_cfg);

if ~converged
    x_trim = zeros(12,1);
    x_trim(3) = -alt_cruise;
    x_trim(4) = V_cruise;
    u_trim = zeros(n_total,1);
end

%% Mission export

N_points = 0;
time_vector = linspace(0, dur_cruise, max(N_points,2));

state_ref   = repmat(x_trim, 1, numel(time_vector));
control_ref = repmat(u_trim, 1, numel(time_vector));

mission = struct();
mission.waypoints = struct('name',"Cruise",'type',"cruise",'altitude',alt_cruise,'mach',mach_cruise,'duration',dur_cruise, ...
    'state',x_trim,'controls',u_trim,'converged',converged,'velocity',V_cruise);
mission.timeline = struct('phase',"Cruise",'t_start',0,'duration',dur_cruise,'wp_start',1,'wp_end',1);
mission.time_vector = time_vector;
mission.altitude_profile = alt_cruise * ones(size(time_vector));
mission.velocity_profile = V_cruise * ones(size(time_vector));
mission.state_profile = state_ref;
mission.control_profile = control_ref;
mission.phase_profile = ones(size(time_vector));
mission.total_duration = dur_cruise;

control_input_data = [time_vector(:), control_ref.'];
sim_stop_time = dur_cruise;

assignin('base','ac',ac);
assignin('base','mission',mission);
assignin('base','waypoints',mission.waypoints);
assignin('base','timeline',mission.timeline);
assignin('base','initial_state',x_trim);
assignin('base','initial_controls',u_trim);
assignin('base','n_cs',n_cs);
assignin('base','n_pe',n_pe);
assignin('base','n_total',n_total);
assignin('base','control_input_data',control_input_data);
assignin('base','sim_stop_time',sim_stop_time);
assignin('base','autopilot',[]);
assignin('base','autopilot_enabled',0);
assignin('base','autopilot_mode',"off");
