clear; clc

ac = Aircraft();

empty_mass = 767;
Ixx = 1100; Iyy = 1500; Izz = 2400;
ac.mass.set_mass_properties(empty_mass, Ixx, Iyy, Izz, 0, 0, 0, [0 0 0]);
ac.mass.set_fuel_capacity(200, 1111);
ac.mass.set_fuel_mass(120);

ac.geometry.wing_area = 16.2;
ac.geometry.wing_span = 11.0;
ac.geometry.mean_aerodynamic_chord = 1.49;
ac.geometry.ref_area = ac.geometry.wing_area;
ac.geometry.ref_span = ac.geometry.wing_span;
ac.geometry.ref_chord = ac.geometry.mean_aerodynamic_chord;

ac.aero.set_lookup(@C172Lookup);

cfg = ac.get_configurator();
cfg.add_control_surface('name','aileron','surface_type','aileron','classification','primary','axis',[1 0 0], ...
    'max_deflection',deg2rad(20),'min_deflection',deg2rad(-20),'dCl',0.08,'dCm',0,'dCn',0.02);
cfg.add_control_surface('name','elevator','surface_type','elevator','classification','primary','axis',[0 1 0], ...
    'max_deflection',deg2rad(25),'min_deflection',deg2rad(-25),'dCl',0,'dCm',-1.05,'dCn',0);
cfg.add_control_surface('name','rudder','surface_type','rudder','classification','primary','axis',[0 0 1], ...
    'max_deflection',deg2rad(30),'min_deflection',deg2rad(-30),'dCl',0.02,'dCm',0,'dCn',-0.08);

max_thrust_sl = 2500;
prop_pos = [1.6 0 0];
prop_dir = [1 0 0];
cfg.add_propulsive_element('name','O320','element_type','prop','max_output',max_thrust_sl, ...
    'position',prop_pos,'direction',prop_dir,'fuel_rate',0,'thrust_model',[]);
ac.propulsive_elements{1}.thrust_model = @(thr,M,alt,V,rho) C172PropModel(thr,M,alt,V,rho,max_thrust_sl);


n_cs = numel(ac.control_surfaces);
n_pe = numel(ac.propulsive_elements);
n_total = n_cs + n_pe;

alt_cruise  = 2000;
mach_cruise = 0.15;
dur_cruise  = 600;

[~, a, ~, ~] = atmosisa(alt_cruise);
V_cruise = mach_cruise * a;

solver = ac.get_trim_solver();
solver.trim_tolerance = 1e-4;
[x_trim, u_trim, converged] = solver.solve_cruise_trim(alt_cruise, mach_cruise);

if ~converged
    x_trim = zeros(12,1);
    x_trim(3) = -alt_cruise;
    x_trim(4) = V_cruise;
    u_trim = zeros(n_total,1);
    if n_pe > 0
        u_trim(n_cs+1:end) = 0.55;
    end
end

N_points = 1000;
time_vector = linspace(0, dur_cruise, N_points);

mission = struct();
mission.waypoints = struct('name',"Cruise",'type',"cruise",'altitude',alt_cruise,'mach',mach_cruise,'duration',dur_cruise, ...
    'state',x_trim,'controls',u_trim,'converged',converged,'velocity',V_cruise);
mission.timeline = struct('phase',"Cruise",'t_start',0,'duration',dur_cruise,'wp_start',1,'wp_end',1);
mission.time_vector = time_vector;
mission.altitude_profile = alt_cruise * ones(size(time_vector));
mission.velocity_profile = V_cruise * ones(size(time_vector));
mission.state_profile = repmat(x_trim, 1, numel(time_vector));
mission.control_profile = repmat(u_trim, 1, numel(time_vector));
mission.phase_profile = ones(size(time_vector));
mission.total_duration = dur_cruise;

control_input_data = [time_vector(:), mission.control_profile.'];
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
