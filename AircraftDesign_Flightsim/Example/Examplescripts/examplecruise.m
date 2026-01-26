clear; clc

ac = Aircraft();

ac.mass.set_mass_properties(0,0,0,0,0,0,0,[0 0 0]);
ac.mass.set_fuel_capacity(0,0);
ac.mass.set_fuel_mass(0);

ac.geometry.wing_area = 0;
ac.geometry.wing_span = 0;
ac.geometry.mean_aerodynamic_chord = 0;
ac.geometry.ref_area = 0;
ac.geometry.ref_span = 0;
ac.geometry.ref_chord = 0;

ac.aero.set_lookup(@YourLookup);

cfg = ac.get_configurator();
cfg.add_control_surface('name',"aileron",'surface_type',"aileron",'classification',"primary",'axis',[1 0 0], ...
    'max_deflection',0,'min_deflection',0,'dCl',0,'dCm',0,'dCn',0);
cfg.add_control_surface('name',"elevator",'surface_type',"elevator",'classification',"primary",'axis',[0 1 0], ...
    'max_deflection',0,'min_deflection',0,'dCl',0,'dCm',0,'dCn',0);
cfg.add_control_surface('name',"rudder",'surface_type',"rudder",'classification',"primary",'axis',[0 0 1], ...
    'max_deflection',0,'min_deflection',0,'dCl',0,'dCm',0,'dCn',0);

cfg.add_propulsive_element('name',"engine",'element_type',"generic",'max_output',0, ...
    'position',[0 0 0],'direction',[1 0 0],'fuel_rate',0,'thrust_model',[]);

n_cs = numel(ac.control_surfaces);
n_pe = numel(ac.propulsive_elements);
n_total = n_cs + n_pe;

alt_cruise  = 0;
mach_cruise = 0;
dur_cruise  = 0;

[~, a, ~, ~] = atmosisa(alt_cruise);
V_cruise = mach_cruise * a;

solver = ac.get_trim_solver();
solver.trim_tolerance = 0;
[x_trim, u_trim, converged] = solver.solve_cruise_trim(alt_cruise, mach_cruise);

if ~converged
    x_trim = zeros(12,1);
    x_trim(3) = -alt_cruise;
    x_trim(4) = V_cruise;
    u_trim = zeros(n_total,1);
end

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
