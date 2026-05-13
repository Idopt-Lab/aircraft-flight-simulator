%% F-16 Mission Using MissionPlanner
clear; clc; close all

model_name = 'flightsim_aircraft';
if bdIsLoaded(model_name)
    set_param(model_name, 'SimulationCommand', 'stop');
end

%% Aircraft Configuration
ac = Aircraft();

empty_mass = 9300;
Ixx = 55800;
Iyy = 63100;
Izz = 118000;

ac.mass.set_mass_properties(empty_mass, Ixx, Iyy, Izz, 0, 0, 0, [0 0 0]);
ac.mass.add_fuel_tank([0, 0, 0], 3000, 'distributed');
ac.mass.add_fuel_tank([-2, 0, 0], 1500, 'distributed');
ac.mass.set_fuel_mass(2400);

ac.geometry.wing_area = 27.87;
ac.geometry.wing_span = 9.14;
ac.geometry.mean_aerodynamic_chord = 3.45;

ac.set_aerodynamics(CoefficientAerodynamics(@F16Lookup));

%% Controls
cfg = ac.get_configurator();

cfg.add_control_surface( ...
    'name','aileron', ...
    'surface_type','aileron', ...
    'classification','primary', ...
    'axis',[1 0 0], ...
    'max_deflection',deg2rad(20), ...
    'min_deflection',deg2rad(-20), ...
    'dCl',0.05,'dCm',0,'dCn',0);

cfg.add_control_surface( ...
    'name','elevator', ...
    'surface_type','elevator', ...
    'classification','primary', ...
    'axis',[0 1 0], ...
    'max_deflection',deg2rad(25), ...
    'min_deflection',deg2rad(-25), ...
    'dCl',0,'dCm',-0.10,'dCn',0);

cfg.add_control_surface( ...
    'name','rudder', ...
    'surface_type','rudder', ...
    'classification','primary', ...
    'axis',[0 0 1], ...
    'max_deflection',deg2rad(30), ...
    'min_deflection',deg2rad(-30), ...
    'dCl',0,'dCm',0,'dCn',-0.08);

%% Propulsion
cfg.add_propulsive_element( ...
    'name','F100_PW_229', ...
    'element_type','turbofan_afterburning', ...
    'max_output',128992, ...
    'position',[-1.5 0 0], ...
    'direction',[1 0 0], ...
    'fuel_rate',2.5, ...
    'thrust_model',@(thr,M,alt,V,rho) F100ThrustModel(thr, M, alt, V, rho, 128992));

n_cs = numel(ac.control_surfaces);
n_pe = numel(ac.propulsive_elements);
n_total = n_cs + n_pe;

%% Mission Segments
segments = struct([]);

segments(1).type     = 'hold_time';
segments(1).altitude = 0;
segments(1).mach     = 0.30;
segments(1).duration = 30;
segments(1).gamma    = deg2rad(8);

segments(2).type     = 'climb';
segments(2).h_start  = 0;
segments(2).h_end    = 500;
segments(2).mach     = 0.30;
segments(2).mach_end = 0.50;
segments(2).gamma    = deg2rad(8);
segments(2).gamma_end= deg2rad(8);

segments(3).type     = 'climb';
segments(3).h_start  = 500;
segments(3).h_end    = 9144;
segments(3).mach     = 0.50;
segments(3).mach_end = 0.60;
segments(3).gamma    = deg2rad(8);
segments(3).gamma_end= deg2rad(5);

segments(4).type     = 'cruise';
segments(4).altitude = 9144;
segments(4).mach     = 0.60;
segments(4).duration = 300;
segments(4).gamma    = 0;

segments(5).type     = 'descent';
segments(5).h_start  = 9144;
segments(5).h_end    = 300;
segments(5).mach     = 0.60;
segments(5).mach_end = 0.25;
segments(5).gamma    = deg2rad(-3);
segments(5).gamma_end= deg2rad(-3);

segments(6).type     = 'approach';
segments(6).h_start  = 300;
segments(6).h_end    = 0;
segments(6).mach     = 0.25;
segments(6).mach_end = 0.22;
segments(6).gamma    = deg2rad(-3);
segments(6).gamma_end= deg2rad(-3);

%% Build + Trim Mission
mp = ac.get_mission_planner(0.25);

mp.build_mission(segments);
mp.compute_trim_schedule('max_iters',15000,'tol',1e-4,'dh_m',300);
mp.export_to_workspace('blend_fraction',0.25);

mission = mp.mission;

%% Useful local variables
initial_state = mission.state_profile(:,1);
initial_controls = mission.control_profile(:,1);

control_input_data = [mission.time_vector(:), mission.control_profile'];

assignin('base','mission_segments',segments);
assignin('base','mission',mission);
assignin('base','initial_state',initial_state);
assignin('base','initial_controls',initial_controls);
assignin('base','control_input_data',control_input_data);
assignin('base','n_cs',n_cs);
assignin('base','n_pe',n_pe);
assignin('base','n_total',n_total);
assignin('base','sim_stop_time',mission.total_duration);

fprintf('Mission duration: %.1f s (%.2f min)\n', mission.total_duration, mission.total_duration/60);
fprintf('Mission points  : %d\n', numel(mission.time_vector));
fprintf('Trim success    : %d / %d points\n', sum(mission.converged_profile), numel(mission.converged_profile));
