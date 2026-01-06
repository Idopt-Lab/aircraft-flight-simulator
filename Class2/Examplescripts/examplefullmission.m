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
cfg.add_control_surface('name',"cs1",'surface_type',"generic",'classification',"primary",'axis',[0 0 0], ...
    'max_deflection',0,'min_deflection',0,'dCl',0,'dCm',0,'dCn',0);
cfg.add_propulsive_element('name',"pe1",'element_type',"generic",'max_output',0, ...
    'position',[0 0 0],'direction',[0 0 0],'fuel_rate',0,'thrust_model',[]);

n_cs = numel(ac.control_surfaces);
n_pe = numel(ac.propulsive_elements);
n_total = n_cs + n_pe;

waypoints = struct();

waypoints(1).name = "Takeoff";
waypoints(1).type = "takeoff";
waypoints(1).altitude = 0;
waypoints(1).mach = 0;
waypoints(1).gamma = 0;

waypoints(2).name = "Climb";
waypoints(2).type = "climb";
waypoints(2).altitude = 0;
waypoints(2).mach = 0;
waypoints(2).gamma = 0;

waypoints(3).name = "Cruise";
waypoints(3).type = "cruise";
waypoints(3).altitude = 0;
waypoints(3).mach = 0;
waypoints(3).duration = 0;

waypoints(4).name = "Dash";
waypoints(4).type = "dash";
waypoints(4).altitude = 0;
waypoints(4).mach = 0;
waypoints(4).duration = 0;

waypoints(5).name = "Descent";
waypoints(5).type = "descent";
waypoints(5).altitude = 0;
waypoints(5).mach = 0;
waypoints(5).gamma = 0;

waypoints(6).name = "Landing";
waypoints(6).type = "landing";
waypoints(6).altitude = 0;
waypoints(6).mach = 0;

solver = ac.get_trim_solver();
solver.trim_tolerance = 0;

for i = 1:numel(waypoints)
    [~, a, ~, ~] = atmosisa(waypoints(i).altitude);
    V = waypoints(i).mach * a;

    switch string(waypoints(i).type)
        case "takeoff"
            [x_trim,u_trim,converged] = solver.solve_takeoff_rotation_trim(waypoints(i).altitude, V);
        case "climb"
            [x_trim,u_trim,converged] = solver.solve_climb_trim(waypoints(i).altitude, waypoints(i).mach, waypoints(i).gamma);
        case "cruise"
            [x_trim,u_trim,converged] = solver.solve_cruise_trim(waypoints(i).altitude, waypoints(i).mach);
        case "dash"
            [x_trim,u_trim,converged] = solver.solve_dash_trim(waypoints(i).altitude, waypoints(i).mach);
        case "descent"
            [x_trim,u_trim,converged] = solver.solve_descent_trim(waypoints(i).altitude, waypoints(i).mach, waypoints(i).gamma);
        case "landing"
            [x_trim,u_trim,converged] = solver.solve_landing_approach_trim(waypoints(i).altitude, waypoints(i).mach);
        otherwise
            [x_trim,u_trim,converged] = solver.solve_cruise_trim(waypoints(i).altitude, waypoints(i).mach);
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
    end
end

timeline = struct();
t = 0;

timeline(1).phase = "Takeoff";
timeline(1).t_start = t;
timeline(1).duration = 0;
timeline(1).wp_start = 1;
timeline(1).wp_end = 1;
t = t + timeline(1).duration;

timeline(2).phase = "Climb";
timeline(2).t_start = t;
timeline(2).duration = 0;
timeline(2).wp_start = 1;
timeline(2).wp_end = 2;
t = t + timeline(2).duration;

timeline(3).phase = "Cruise";
timeline(3).t_start = t;
timeline(3).duration = waypoints(3).duration;
timeline(3).wp_start = 2;
timeline(3).wp_end = 3;
t = t + timeline(3).duration;

timeline(4).phase = "Dash";
timeline(4).t_start = t;
timeline(4).duration = waypoints(4).duration;
timeline(4).wp_start = 3;
timeline(4).wp_end = 4;
t = t + timeline(4).duration;

timeline(5).phase = "Descent";
timeline(5).t_start = t;
timeline(5).duration = 0;
timeline(5).wp_start = 4;
timeline(5).wp_end = 5;
t = t + timeline(5).duration;

timeline(6).phase = "Landing";
timeline(6).t_start = t;
timeline(6).duration = 0;
timeline(6).wp_start = 5;
timeline(6).wp_end = 6;
t = t + timeline(6).duration;

total_duration = t;

N_points = 0;
time_vector = linspace(0, total_duration, max(N_points,2));

altitude_ref = zeros(size(time_vector));
velocity_ref = zeros(size(time_vector));
state_ref = zeros(12, numel(time_vector));
control_ref = zeros(n_total, numel(time_vector));
phase_ref = zeros(size(time_vector));

for k = 1:numel(time_vector)
    t_curr = time_vector(k);
    phase_idx = find([timeline.t_start] <= t_curr, 1, 'last');
    if isempty(phase_idx), phase_idx = 1; end
    if phase_idx > numel(timeline), phase_idx = numel(timeline); end

    phase_ref(k) = phase_idx;

    t_phase = t_curr - timeline(phase_idx).t_start;
    tau = t_phase / max(timeline(phase_idx).duration, 1e-6);
    tau = min(max(tau,0),1);

    wp1 = timeline(phase_idx).wp_start;
    wp2 = timeline(phase_idx).wp_end;

    altitude_ref(k) = waypoints(wp1).altitude*(1-tau) + waypoints(wp2).altitude*tau;
    velocity_ref(k) = waypoints(wp1).velocity*(1-tau) + waypoints(wp2).velocity*tau;
    state_ref(:,k) = waypoints(wp1).state*(1-tau) + waypoints(wp2).state*tau;
    control_ref(:,k) = waypoints(wp1).controls*(1-tau) + waypoints(wp2).controls*tau;
end

autopilot = [];
autopilot_enabled = 0;
autopilot_mode = "off";

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

time_vec = time_vector(:);
control_input_data = [time_vec, control_ref.'];

sim_stop_time = total_duration;

assignin('base','ac',ac);
assignin('base','mission',mission);
assignin('base','autopilot',autopilot);
assignin('base','autopilot_enabled',autopilot_enabled);
assignin('base','autopilot_mode',autopilot_mode);
assignin('base','initial_state',initial_state);
assignin('base','initial_controls',initial_controls);
assignin('base','Initialpos',Initialpos);
assignin('base','InitialVel',InitialVel);
assignin('base','InitialOri',InitialOri);
assignin('base','InitialRot',InitialRot);
assignin('base','n_cs',n_cs);
assignin('base','n_pe',n_pe);
assignin('base','n_total',n_total);
assignin('base','control_input_data',control_input_data);
assignin('base','sim_stop_time',sim_stop_time);
