clear; clc;

ac = Aircraft();

empty_mass = 37648;
Ixx = 7.62e5; Iyy = 2.00e6; Izz = 2.57e6;
Ixz = 1.08e4; Ixy = 0; Iyz = 0;
ac.mass.set_mass_properties(empty_mass, Ixx, Iyy, Izz, Ixy, Ixz, Iyz, [0 0 0]);

ac.mass.add_fuel_tank([0 -5 -1], 5000, 'distributed');
ac.mass.add_fuel_tank([0  5 -1], 5000, 'distributed');
ac.mass.add_fuel_tank([0  0 -1], 5000, 'distributed');
ac.mass.set_fuel_mass(10000);
ac.mass.set_fuel_capacity(15000, 70000);

wing_area = 108.77;
wing_span = 28.86;
wing_chord = 3.75;

ac.geometry.wing_area = wing_area;
ac.geometry.wing_span = wing_span;
ac.geometry.mean_aerodynamic_chord = wing_chord;
ac.geometry.ref_area = wing_area;
ac.geometry.ref_span = wing_span;
ac.geometry.ref_chord = wing_chord;

ac.aero.set_lookup(@Boeing737Lookup);

cfg = ac.get_configurator();

cfg.add_control_surface('name','aileron','surface_type','aileron','classification','primary','axis',[1 0 0], ...
    'max_deflection',deg2rad(20),'min_deflection',deg2rad(-20),'dCl',0.10,'dCm',0,'dCn',0);

cfg.add_control_surface('name','elevator','surface_type','elevator','classification','primary','axis',[0 1 0], ...
    'max_deflection',deg2rad(17),'min_deflection',deg2rad(-17),'dCl',0,'dCm',-1.20,'dCn',0);

cfg.add_control_surface('name','rudder','surface_type','rudder','classification','primary','axis',[0 0 1], ...
    'max_deflection',deg2rad(20),'min_deflection',deg2rad(-20),'dCl',0,'dCm',0,'dCn',-0.20);

max_thrust_per_engine = 89000 * 4.44822;

cfg.add_propulsive_element('name','CFM56_L','element_type','turbofan_blockset','max_output',max_thrust_per_engine, ...
    'position',[-3.5 -4.9 -1.0],'direction',[1;0;0],'fuel_rate',0.55,'thrust_model',[]);

cfg.add_propulsive_element('name','CFM56_R','element_type','turbofan_blockset','max_output',max_thrust_per_engine, ...
    'position',[-3.5  4.9 -1.0],'direction',[1;0;0],'fuel_rate',0.55,'thrust_model',[]);

bypass_ratio = 5.5;
fan_pressure_ratio = 1.65;
compressor_pressure_ratio = 27.5;
turbine_inlet_temp = 1600;
mass_flow_corrected = 450;
fuel_heating_value = 43e6;
fan_efficiency = 0.89;
compressor_efficiency = 0.87;
turbine_efficiency = 0.89;
nozzle_efficiency = 0.97;

ac.propulsive_elements{1}.set_turbofan_params(bypass_ratio, fan_pressure_ratio, compressor_pressure_ratio, ...
    turbine_inlet_temp, mass_flow_corrected, fuel_heating_value, ...
    fan_efficiency, compressor_efficiency, turbine_efficiency, nozzle_efficiency);

ac.propulsive_elements{2}.set_turbofan_params(bypass_ratio, fan_pressure_ratio, compressor_pressure_ratio, ...
    turbine_inlet_temp, mass_flow_corrected, fuel_heating_value, ...
    fan_efficiency, compressor_efficiency, turbine_efficiency, nozzle_efficiency);

n_cs = numel(ac.control_surfaces);
n_pe = numel(ac.propulsive_elements);
n_total = n_cs + n_pe;

alt_cruise = 9144;
mach_cruise = 0.62;
dur_cruise = 600;

[~, a, ~, ~] = atmosisa(alt_cruise);
V_cruise = mach_cruise * a;

fprintf('\n=== BOEING 737-800 CRUISE SETUP ===\n');
fprintf('Altitude: %.0f m (%.0f ft)\n', alt_cruise, alt_cruise/0.3048);
fprintf('Mach: %.2f\n', mach_cruise);
fprintf('Velocity: %.2f m/s\n', V_cruise);
fprintf('Mass: %.1f kg\n', ac.mass.get_total_mass());

solver = ac.get_trim_solver();
solver.trim_tolerance = 1e-3;
solver.max_iterations = 15000;
solver.use_fmincon = true;
solver.initial_guess = [deg2rad(2.5); deg2rad(0.0); 0.65];

[x_trim, u_trim, converged] = solver.solve_cruise_trim(alt_cruise, mach_cruise);

if converged
    fprintf('\nTrim converged!\n');
    solver.print_summary();
else
    fprintf('\nTrim failed, using fallback\n');
    x_trim = zeros(12,1);
    x_trim(3) = -alt_cruise;
    x_trim(4) = V_cruise;
    u_trim = zeros(n_total,1);
    if n_pe > 0
        u_trim(n_cs+1:end) = 0.65;
    end
end

N_points = 1000;
time_vector = linspace(0, dur_cruise, N_points);

state_ref = repmat(x_trim, 1, N_points);
control_ref = repmat(u_trim, 1, N_points);

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

initial_state = x_trim;
initial_controls = u_trim;

Initialpos = initial_state(1:3).';
InitialVel = initial_state(4:6).';
InitialOri = initial_state(7:9).';
InitialRot = initial_state(10:12).';

autopilot = [];
autopilot_enabled = 0;
autopilot_mode = "off";
autopilot_enable_time = inf;
autopilot_min_alt = 0;

ground_k = 1e6;
ground_c = 1e5;

control_input_data = [time_vector(:), control_ref.'];
sim_stop_time = dur_cruise;

assignin('base','ac',ac);
assignin('base','mission',mission);
assignin('base','waypoints',mission.waypoints);
assignin('base','timeline',mission.timeline);
assignin('base','initial_state',initial_state);
assignin('base','initial_controls',initial_controls);
assignin('base','Initialpos',Initialpos);
assignin('base','InitialVel',InitialVel);
assignin('base','InitialOri',InitialOri);
assignin('base','InitialRot',InitialRot);
assignin('base','autopilot',autopilot);
assignin('base','autopilot_enabled',autopilot_enabled);
assignin('base','autopilot_mode',autopilot_mode);
assignin('base','autopilot_enable_time',autopilot_enable_time);
assignin('base','autopilot_min_alt',autopilot_min_alt);
assignin('base','n_cs',n_cs);
assignin('base','n_pe',n_pe);
assignin('base','n_total',n_total);
assignin('base','control_input_data',control_input_data);
assignin('base','sim_stop_time',sim_stop_time);
assignin('base','ground_k',ground_k);
assignin('base','ground_c',ground_c);


ac.state.set_full_state(x_trim);
ac.set_controls_from_vector(u_trim);
[F_check, M_check, ~] = ac.calculate_total_forces_moments_with_gravity();
fprintf('Total forces (should be ~0): [%.2f, %.2f, %.2f] N\n', F_check(1), F_check(2), F_check(3));
fprintf('Total moments (should be ~0): [%.2f, %.2f, %.2f] N-m\n', M_check(1), M_check(2), M_check(3));
fprintf('Weight: %.2f N\n', ac.mass.get_total_mass() * 9.81);

