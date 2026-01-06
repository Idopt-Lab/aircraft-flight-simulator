clear; clc;

ac = Aircraft();

empty_mass = 9300;
Ixx = 55800; Iyy = 63100; Izz = 118000;
ac.mass.set_mass_properties(empty_mass, Ixx, Iyy, Izz, 0, 0, 0, [0 0 0]);
ac.mass.add_fuel_tank([0, 0, 0], 3000, 'distributed');
ac.mass.add_fuel_tank([-2, 0, 0], 1500, 'distributed');
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

waypoints = struct();

waypoints(1).name = 'Takeoff';
waypoints(1).type = 'takeoff';
waypoints(1).altitude = 0;
waypoints(1).mach = 0.30;
waypoints(1).gamma = deg2rad(8);

waypoints(2).name = 'Climb Start';
waypoints(2).type = 'climb';
waypoints(2).altitude = 500;
waypoints(2).mach = 0.50;
waypoints(2).gamma = deg2rad(8);

waypoints(3).name = 'Climb End';
waypoints(3).type = 'climb';
waypoints(3).altitude = 9144;
waypoints(3).mach = 0.60;
waypoints(3).gamma = deg2rad(5);

waypoints(4).name = 'Cruise';
waypoints(4).type = 'cruise';
waypoints(4).altitude = 9144;
waypoints(4).mach = 0.60;
waypoints(4).duration = 300;

waypoints(5).name = 'Dash Altitude';
waypoints(5).type = 'climb';
waypoints(5).altitude = 12000;
waypoints(5).mach = 0.60;
waypoints(5).gamma = deg2rad(3);

waypoints(6).name = 'Dash';
waypoints(6).type = 'dash';
waypoints(6).altitude = 12000;
waypoints(6).mach = 1.2;
waypoints(6).duration = 300;

waypoints(7).name = 'Approach';
waypoints(7).type = 'descent';
waypoints(7).altitude = 300;
waypoints(7).mach = 0.25;
waypoints(7).gamma = deg2rad(-3);

waypoints(8).name = 'Landing';
waypoints(8).type = 'landing';
waypoints(8).altitude = 0;
waypoints(8).mach = 0.22;
waypoints(8).gamma = deg2rad(-3);

solver = ac.get_trim_solver();
solver.trim_tolerance = 1e-4;

for i = 1:length(waypoints)
    [~, a, ~, ~] = atmosisa(waypoints(i).altitude);
    V = waypoints(i).mach * a;

    switch waypoints(i).type
        case 'takeoff'
            [x_trim, u_trim, converged, info] = solver.solve_takeoff_rotation_trim(waypoints(i).altitude, V);
        case 'climb'
            gamma = waypoints(i).gamma;
            [x_trim, u_trim, converged, info] = solver.solve_climb_trim(waypoints(i).altitude, waypoints(i).mach, gamma);
        case 'cruise'
            [x_trim, u_trim, converged, info] = solver.solve_cruise_trim(waypoints(i).altitude, waypoints(i).mach);
        case 'dash'
            [x_trim, u_trim, converged, info] = solver.solve_dash_trim(waypoints(i).altitude, waypoints(i).mach);
        case 'descent'
            gamma = waypoints(i).gamma;
            [x_trim, u_trim, converged, info] = solver.solve_descent_trim(waypoints(i).altitude, waypoints(i).mach, gamma);
        case 'landing'
            [x_trim, u_trim, converged, info] = solver.solve_landing_approach_trim(waypoints(i).altitude, waypoints(i).mach);
        otherwise
            [x_trim, u_trim, converged, info] = solver.solve_cruise_trim(waypoints(i).altitude, waypoints(i).mach);
    end

    waypoints(i).state = x_trim;
    waypoints(i).controls = u_trim;
    waypoints(i).converged = converged;
    waypoints(i).velocity = V;

    if ~converged
        waypoints(i).state = zeros(12,1);
        waypoints(i).state(3) = -waypoints(i).altitude;
        waypoints(i).state(4) = V;
        waypoints(i).controls = zeros(n_total,1);
        if n_pe > 0
            waypoints(i).controls(n_cs+1:end) = 0.7;
        end
    end
end

to = TakeoffAnalysis(ac);
[~, to_res] = to.calculate_takeoff(0, 0, inf);
takeoff_distance = to_res.bfl_results.BFL_m;

perf = ac.get_performance();
alt_grid = (0:1000:15000).';
V_grid = (80:5:220).';
climb = perf.calculate_full_throttle_climb_schedule(alt_grid, V_grid);
climb_rate_avg = mean(climb.ROC_opt(1:10));

ld = LandingAnalysis(ac);
[landing_distance, ~] = ld.calculate_landing(0, 0, inf);

timeline = struct();
t = 0;

timeline(1).phase = 'Takeoff';
timeline(1).t_start = t;
timeline(1).duration = 30;
timeline(1).wp_start = 1;
timeline(1).wp_end = 1;
t = t + timeline(1).duration;

timeline(2).phase = 'Initial Climb';
timeline(2).t_start = t;
timeline(2).duration = max((waypoints(2).altitude - waypoints(1).altitude) / max(climb_rate_avg,1), 10);
timeline(2).wp_start = 1;
timeline(2).wp_end = 2;
t = t + timeline(2).duration;

timeline(3).phase = 'Climb to Cruise';
timeline(3).t_start = t;
timeline(3).duration = max((waypoints(3).altitude - waypoints(2).altitude) / max(climb_rate_avg,1), 10);
timeline(3).wp_start = 2;
timeline(3).wp_end = 3;
t = t + timeline(3).duration;

timeline(4).phase = 'Cruise';
timeline(4).t_start = t;
timeline(4).duration = waypoints(4).duration;
timeline(4).wp_start = 3;
timeline(4).wp_end = 4;
t = t + timeline(4).duration;

timeline(5).phase = 'Climb to Dash';
timeline(5).t_start = t;
timeline(5).duration = abs(waypoints(5).altitude - waypoints(4).altitude) / 50;
timeline(5).wp_start = 4;
timeline(5).wp_end = 5;
t = t + timeline(5).duration;

timeline(6).phase = 'Dash';
timeline(6).t_start = t;
timeline(6).duration = waypoints(6).duration;
timeline(6).wp_start = 5;
timeline(6).wp_end = 6;
t = t + timeline(6).duration;

timeline(7).phase = 'Descent to Approach';
timeline(7).t_start = t;
timeline(7).duration = abs(waypoints(7).altitude - waypoints(6).altitude) / 30;
timeline(7).wp_start = 6;
timeline(7).wp_end = 7;
t = t + timeline(7).duration;

timeline(8).phase = 'Landing';
timeline(8).t_start = t;
timeline(8).duration = 60;
timeline(8).wp_start = 7;
timeline(8).wp_end = 8;
t = t + timeline(8).duration;

total_duration = t;

N_points = 2000;
time_vector = linspace(0, total_duration, N_points);

altitude_ref = zeros(size(time_vector));
velocity_ref = zeros(size(time_vector));
state_ref = zeros(12, length(time_vector));
control_ref = zeros(n_total, length(time_vector));
phase_ref = zeros(size(time_vector));

for i = 1:length(time_vector)
    t_curr = time_vector(i);
    phase_idx = find([timeline.t_start] <= t_curr, 1, 'last');
    if isempty(phase_idx), phase_idx = 1; end
    if phase_idx > length(timeline), phase_idx = length(timeline); end

    phase_ref(i) = phase_idx;

    t_phase = t_curr - timeline(phase_idx).t_start;
    tau = t_phase / max(timeline(phase_idx).duration, 1e-6);
    tau = min(max(tau, 0), 1);

    wp1 = timeline(phase_idx).wp_start;
    wp2 = timeline(phase_idx).wp_end;

    altitude_ref(i) = waypoints(wp1).altitude * (1-tau) + waypoints(wp2).altitude * tau;
    velocity_ref(i) = waypoints(wp1).velocity * (1-tau) + waypoints(wp2).velocity * tau;
    state_ref(:,i) = waypoints(wp1).state * (1-tau) + waypoints(wp2).state * tau;
    control_ref(:,i) = waypoints(wp1).controls * (1-tau) + waypoints(wp2).controls * tau;
end

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
autopilot_enable_time = 35;
autopilot_min_alt = 20;

ground_k = 3e5;
ground_c = 5e4;

initial_state = waypoints(1).state;
initial_controls = waypoints(1).controls;

Initialpos = initial_state(1:3).';
InitialVel = initial_state(4:6).';
InitialOri = initial_state(7:9).';
InitialRot = initial_state(10:12).';

mission = struct();
mission.waypoints = waypoints;
mission.timeline = timeline;
mission.time_vector = time_vector;
mission.altitude_profile = altitude_ref;
mission.velocity_profile = velocity_ref;
mission.state_profile = state_ref;
mission.control_profile = control_ref;
mission.phase_profile = phase_ref;
mission.total_duration = total_duration;
mission.takeoff_distance = takeoff_distance;
mission.landing_distance = landing_distance;

time_vec = time_vector(:);
control_data = control_ref.';
control_input_data = [time_vec, control_data];

sim_stop_time = total_duration;

assignin('base','ac',ac);
assignin('base','mission',mission);
assignin('base','autopilot',autopilot);
assignin('base','waypoints',waypoints);
assignin('base','timeline',timeline);
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
