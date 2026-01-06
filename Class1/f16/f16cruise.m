clear; clc

ac = Aircraft();

empty_mass = 9300;
Ixx = 55800; Iyy = 63100; Izz = 118000;
ac.mass.set_mass_properties(empty_mass, Ixx, Iyy, Izz, 0, 0, 0, [0 0 0]);
ac.mass.add_fuel_tank([0 0 0], 3000, 'distributed');
ac.mass.add_fuel_tank([-2 0 0], 1500, 'distributed');
ac.mass.set_fuel_mass(2400);
ac.mass.set_fuel_capacity(12500, 12500);

ac.geometry.wing_area = 27.87;
ac.geometry.wing_span = 9.14;
ac.geometry.mean_aerodynamic_chord = 3.45;
ac.geometry.ref_area = 27.87;
ac.geometry.ref_span = 9.14;
ac.geometry.ref_chord = 3.45;

ac.aero.set_lookup(@F16Lookup);

cfg = ac.get_configurator();
cfg.add_control_surface('name','aileron','surface_type','aileron','classification','primary','axis',[1 0 0], ...
    'max_deflection',deg2rad(20),'min_deflection',deg2rad(-20),'dCl',0.05,'dCm',0,'dCn',0);
cfg.add_control_surface('name','elevator','surface_type','elevator','classification','primary','axis',[0 1 0], ...
    'max_deflection',deg2rad(25),'min_deflection',deg2rad(-25),'dCl',0,'dCm',-0.10,'dCn',0);
cfg.add_control_surface('name','rudder','surface_type','rudder','classification','primary','axis',[0 0 1], ...
    'max_deflection',deg2rad(30),'min_deflection',deg2rad(-30),'dCl',0,'dCm',0,'dCn',-0.08);

cfg.add_propulsive_element('name','F100_PW_229','element_type','turbofan_afterburning','max_output',128992, ...
    'position',[-1.5 0 0],'direction',[1 0 0],'fuel_rate',2.5,'thrust_model',[]);
ac.propulsive_elements{1}.thrust_model = @(thr,M,alt,V,rho) F100ThrustModel(thr, M, alt, V, rho, 128992);

n_cs = numel(ac.control_surfaces);
n_pe = numel(ac.propulsive_elements);
n_total = n_cs + n_pe;

alt_cruise = 9144;
mach_cruise = 0.60;
dur_cruise = 600;

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
        u_trim(n_cs+1:end) = 0.7;
    end
end

N_points = 1000;
time_vector = linspace(0, dur_cruise, N_points);

altitude_ref = alt_cruise * ones(1, N_points);
velocity_ref = V_cruise * ones(1, N_points);
state_ref = repmat(x_trim, 1, N_points);
control_ref = repmat(u_trim, 1, N_points);
phase_ref = ones(1, N_points);

autopilot = Autopilot_v5(ac);
autopilot.enabled = true;
autopilot.mode = "altitude+speed";
autopilot.Kp_h = 0.003;
autopilot.Ki_h = 0.0005;
autopilot.Kd_h = 0.02;
autopilot.Kp_theta = 1.2;
autopilot.Kd_q = 0.25;
autopilot.Kp_V = 0.01;
autopilot.Ki_V = 0.001;
autopilot.theta_up_lim = deg2rad(18);
autopilot.theta_dn_lim = deg2rad(12);
autopilot.max_de_rate = deg2rad(60);
autopilot.max_thr_rate = 0.5;
autopilot.min_alt_protect = 2;

autopilot_enabled = 1;
autopilot_mode = "altitude+speed";
autopilot_enable_time = 0;
autopilot_min_alt = 20;

ground_k = 3e5;
ground_c = 5e4;

initial_state = x_trim;
initial_controls = u_trim;

Initialpos = initial_state(1:3).';
InitialVel = initial_state(4:6).';
InitialOri = initial_state(7:9).';
InitialRot = initial_state(10:12).';

mission = struct();
mission.waypoints = struct('name',"Cruise",'type',"cruise",'altitude',alt_cruise,'mach',mach_cruise,'duration',dur_cruise,'state',x_trim,'controls',u_trim,'converged',converged,'velocity',V_cruise);
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

time_vec = time_vector(:);
control_data = control_ref.';
control_input_data = [time_vec, control_data];

sim_stop_time = dur_cruise;

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
assignin('base','sim_stop_time',sim_stop_time);
assignin('base','ground_k',ground_k);
assignin('base','ground_c',ground_c);
