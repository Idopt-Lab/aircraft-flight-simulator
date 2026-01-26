clear; clc; close all

modelName = 'flightsim7';
load_system(modelName);

ac = Aircraft();

empty_mass = 767;
Ixx = 1285;
Iyy = 1824;
Izz = 2666;
Ixz = 0;

cg = [2.11, 0.00, 1.26];
ac.mass.set_mass_properties(empty_mass, Ixx, Iyy, Izz, 0, Ixz, 0, cg);
ac.mass.set_fuel_capacity(200, 1111);
ac.mass.set_fuel_mass(120);

wing_span = 2 * 5.49;
wing_root_chord = 1.626;
wing_tip_chord = 1.13;
wing_mac = (2/3) * wing_root_chord * (1 + wing_tip_chord/wing_root_chord + (wing_tip_chord/wing_root_chord)^2) / ...
           (1 + wing_tip_chord/wing_root_chord);
wing_area = (wing_root_chord + wing_tip_chord) * 5.49;

ac.geometry.wing_area = wing_area;
ac.geometry.wing_span = wing_span;
ac.geometry.mean_aerodynamic_chord = wing_mac;
ac.geometry.ref_area = wing_area;
ac.geometry.ref_span = wing_span;
ac.geometry.ref_chord = wing_mac;

datcom_fp = "C:\Users\naman\Downloads\datcom_pack\datcom_pack\cessna.out";
C172Lookup = c172datcom(datcom_fp);
ac.aero.set_lookup(C172Lookup);

cfg = ac.get_configurator();

cfg.add_control_surface('name','aileron','surface_type','aileron','classification','primary','axis',[1 0 0], ...
    'max_deflection',deg2rad(20),'min_deflection',deg2rad(-20),'dCl',0,'dCm',0,'dCn',0);

cfg.add_control_surface('name','elevator','surface_type','elevator','classification','primary','axis',[0 1 0], ...
    'max_deflection',deg2rad(28),'min_deflection',deg2rad(-26),'dCl',0.30,'dCm',-1.20,'dCn',0);

cfg.add_control_surface('name','rudder','surface_type','rudder','classification','primary','axis',[0 0 1], ...
    'max_deflection',deg2rad(30),'min_deflection',deg2rad(-30),'dCl',0,'dCm',0,'dCn',-0.075);

cfg.add_propulsive_element('name','O320','element_type','propeller','max_output',119310, ...
    'position',[1.68 0 1.26],'direction',[1 0 0],'fuel_rate',0,'thrust_model',@C172PropModel);

n_cs = numel(ac.control_surfaces);
n_pe = numel(ac.propulsive_elements);
n_total = n_cs + n_pe;

alt_cruise = 1000;
mach_cruise = 0.15;
dur_cruise = 120;

[~, a, ~, rho] = atmosisa(alt_cruise);
V_cruise = mach_cruise * max(a,1e-6);

fprintf('\n=== CESSNA 172 CRUISE SETUP ===\n');
fprintf('Altitude: %.0f m\n', alt_cruise);
fprintf('Velocity: %.2f m/s\n', V_cruise);
fprintf('Mass: %.1f kg\n', ac.mass.get_total_mass());

solver = ac.get_trim_solver();
solver.trim_tolerance = 1e-3;
solver.max_iterations = 15000;
solver.use_fmincon = true;
solver.initial_guess = [deg2rad(2.0); deg2rad(-1.0); 0.65];

[x_trim, u_trim, converged] = solver.solve_cruise_trim(alt_cruise, mach_cruise);

if converged
    fprintf('\nTrim converged!\n');
    solver.print_summary();
else
    fprintf('\nTrim failed!\n');
    x_trim = zeros(12,1);
    x_trim(3) = -alt_cruise;
    x_trim(4) = V_cruise;
    u_trim = zeros(n_total,1);
    u_trim(2) = deg2rad(-1.0);
    u_trim(n_cs+1) = 0.65;
end

dt_fc = 0.01;
sim_time = dur_cruise;

N_points = 1000;
time_vector = linspace(0, dur_cruise, N_points);

altitude_ref = alt_cruise * ones(1, N_points);
velocity_ref = V_cruise * ones(1, N_points);
state_ref = repmat(x_trim, 1, N_points);
control_ref = repmat(u_trim, 1, N_points);
phase_ref = ones(1, N_points);

initial_state = x_trim;
initial_controls = u_trim;

Initialpos = initial_state(1:3).';
InitialVel = initial_state(4:6).';
InitialOri = initial_state(7:9).';
InitialRot = initial_state(10:12).';

fprintf('\n=== INITIAL CONDITIONS ===\n');
fprintf('Position (NED): [%.2f, %.2f, %.2f] m\n', Initialpos(1), Initialpos(2), Initialpos(3));
fprintf('Velocity (BODY): [%.4f, %.4f, %.4f] m/s\n', InitialVel(1), InitialVel(2), InitialVel(3));
fprintf('  u (forward) = %.4f m/s\n', InitialVel(1));
fprintf('  v (right) = %.4f m/s\n', InitialVel(2));
fprintf('  w (down) = %.4f m/s\n', InitialVel(3));
fprintf('  |V| = %.4f m/s\n', norm(InitialVel));
fprintf('Orientation: [%.4f, %.4f, %.4f] deg\n', rad2deg(InitialOri(1)), rad2deg(InitialOri(2)), rad2deg(InitialOri(3)));
fprintf('Angular rates: [%.6f, %.6f, %.6f] rad/s\n', InitialRot(1), InitialRot(2), InitialRot(3));

mission = struct();
mission.waypoints = struct('name',"Cruise",'type',"cruise",'altitude',alt_cruise, ...
    'mach',mach_cruise,'duration',dur_cruise,'state',x_trim,'controls',u_trim, ...
    'converged',converged,'velocity',V_cruise);
mission.timeline = struct('phase',"Cruise",'t_start',0,'duration',dur_cruise,'wp_start',1,'wp_end',1);
mission.time_vector = time_vector;
mission.altitude_profile = altitude_ref;
mission.velocity_profile = velocity_ref;
mission.state_profile = state_ref;
mission.control_profile = control_ref;
mission.phase_profile = phase_ref;
mission.total_duration = dur_cruise;
mission.takeoff_distance = 0;
mission.landing_distance = 0;

autopilot = [];
autopilot_enabled = 0;
autopilot_mode = "off";
autopilot_enable_time = inf;
autopilot_min_alt = 0;

ground_k = 5e4;
ground_c = 5e3;

time_vec = time_vector(:);
control_data = repmat(u_trim(:)', length(time_vec), 1);
control_input_data = [time_vec, control_data];

fprintf('\n=== SIMULATION CONFIGURATION ===\n');
fprintf('Fixed step: %.4f s\n', dt_fc);
fprintf('Duration: %.0f s\n', sim_time);

assignin('base','ac',ac);
assignin('base','mission',mission);
assignin('base','autopilot',autopilot);
assignin('base','waypoints',mission.waypoints);
assignin('base','timeline',mission.timeline);
assignin('base','initial_state',initial_state);
assignin('base','initial_controls',initial_controls);
assignin('base','Initialpos',Initialpos);
assignin('base','InitialVel',InitialVel);
assignin('base','InitialOri',InitialOri);
assignin('base','InitialRot',InitialRot);
assignin('base','autopilot_enabled',autopilot_enabled);
assignin('base','autopilot_mode',autopilot_mode);
assignin('base','autopilot_enable_time',autopilot_enable_time);
assignin('base','autopilot_min_alt',autopilot_min_alt);
assignin('base','n_cs',n_cs);
assignin('base','n_pe',n_pe);
assignin('base','n_total',n_total);
assignin('base','control_input_data',control_input_data);
assignin('base','sim_stop_time',sim_time);
assignin('base','ground_k',ground_k);
assignin('base','ground_c',ground_c);

fprintf('\n=== TRIM VERIFICATION ===\n');
ac.state.set_full_state(x_trim);
ac.set_controls_from_vector(u_trim);
[F_check, M_check, ~] = ac.calculate_total_forces_moments_with_gravity();
fprintf('Total forces (should be ~0): [%.2f, %.2f, %.2f] N\n', F_check(1), F_check(2), F_check(3));
fprintf('Total moments (should be ~0): [%.2f, %.2f, %.2f] N-m\n', M_check(1), M_check(2), M_check(3));
fprintf('Weight: %.2f N\n', ac.mass.get_total_mass() * 9.81);

fprintf('\n=== UPDATING 6DOF BLOCK ===\n');

blockPath = [modelName '/6DOF (Euler Angles)'];

try
    set_param(blockPath, 'Xe_0', mat2str(Initialpos(1)));
    set_param(blockPath, 'Ye_0', mat2str(Initialpos(2)));
    set_param(blockPath, 'Ze_0', mat2str(Initialpos(3)));
    
    set_param(blockPath, 'U_0', mat2str(InitialVel(1)));
    set_param(blockPath, 'V_0', mat2str(InitialVel(2)));
    set_param(blockPath, 'W_0', mat2str(InitialVel(3)));
    
    set_param(blockPath, 'phi_0', mat2str(InitialOri(1)));
    set_param(blockPath, 'theta_0', mat2str(InitialOri(2)));
    set_param(blockPath, 'psi_0', mat2str(InitialOri(3)));
    
    set_param(blockPath, 'p_0', mat2str(InitialRot(1)));
    set_param(blockPath, 'q_0', mat2str(InitialRot(2)));
    set_param(blockPath, 'r_0', mat2str(InitialRot(3)));
    
    fprintf('6DOF block updated:\n');
    fprintf('  Position: [%.2f, %.2f, %.2f]\n', Initialpos(1), Initialpos(2), Initialpos(3));
    fprintf('  Velocity: [%.4f, %.4f, %.4f]\n', InitialVel(1), InitialVel(2), InitialVel(3));
    fprintf('  Euler: [%.6f, %.6f, %.6f]\n', InitialOri(1), InitialOri(2), InitialOri(3));
catch ME
    warning('Could not update 6DOF block automatically: %s', ME.message);
    fprintf('Please manually set 6DOF initial conditions to use workspace variables\n');
end

fprintf('\n=== READY TO SIMULATE ===\n');