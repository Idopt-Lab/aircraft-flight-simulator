%% F-16 Mission Using MissionPlanner
clear; clc; close all

%% PATH BOOTSTRAP
% Ensures the shared class library and example helper functions are on the
% MATLAB path regardless of the current working directory or session state.
this_script_dir = fileparts(mfilename('fullpath'));
matlab_root_dir = fileparts(fileparts(this_script_dir));
addpath(genpath(fullfile(matlab_root_dir,'classes')));
addpath(genpath(fullfile(matlab_root_dir,'examples')));

model_name = 'flightsim_aircraft';
if bdIsLoaded(model_name)
    set_param(model_name, 'SimulationCommand', 'stop');
end

%% Aircraft object

ac = Aircraft();

%% Geometry / reference positions

S_ref = 27.87;
b_ref = 9.14;
c_ref = 3.45;
ac.geometry.set_reference_geometry(S_ref,b_ref,c_ref);
ac.geometry.set_reference_point([0;0;0]);

%% Frames

ac.set_body_frame("body");
ac.set_reference_frame("body");

cfg = ac.get_configurator();

%% Mass properties

empty_mass = 9300;
Ixx = 55800;
Iyy = 63100;
Izz = 118000;
I_body = diag([Ixx Iyy Izz]);

cfg.add_component('name','airframe', 'type','airframe', 'mass',empty_mass, 'position',[0 0 0], 'inertia',I_body, 'parent_frame','body');

% Two fuel tanks (nose and -2 m station), current load split by capacity.
fuel_tank_1_capacity = 3000;
fuel_tank_2_capacity = 1500;
total_fuel_mass = 2400;
fuel_tank_1_mass = total_fuel_mass* fuel_tank_1_capacity/(fuel_tank_1_capacity+fuel_tank_2_capacity);
fuel_tank_2_mass = total_fuel_mass-fuel_tank_1_mass;

cfg.add_component('name','fuel_tank_1', 'type','fuel', 'mass',fuel_tank_1_mass, 'position',[0 0 0], 'parent_frame','body');

cfg.add_component('name','fuel_tank_2', 'type','fuel', 'mass',fuel_tank_2_mass, 'position',[-2 0 0], 'parent_frame','body');

%% Aerodynamics

cfg.add_aero_solver(CoefficientAerodynamics(@F16Lookup), 'name','f16_aero', 'position',[0;0;0], 'geom',ac.geometry, 'parent_frame','body');

%% Controls

cfg.add_control_surface('name','aileron', 'surface_type','aileron', 'classification','primary', 'axis',[1 0 0], 'max_deflection',deg2rad(20), 'min_deflection',deg2rad(-20), 'dCl',0.05,'dCm',0,'dCn',0);

cfg.add_control_surface('name','elevator', 'surface_type','elevator', 'classification','primary', 'axis',[0 1 0], 'max_deflection',deg2rad(25), 'min_deflection',deg2rad(-25), 'dCl',0,'dCm',-0.10,'dCn',0);

cfg.add_control_surface('name','rudder', 'surface_type','rudder', 'classification','primary', 'axis',[0 0 1], 'max_deflection',deg2rad(30), 'min_deflection',deg2rad(-30), 'dCl',0,'dCm',0,'dCn',-0.08);

%% Propulsion
cfg.add_propulsive_element( ...
    'name','F100_PW_229', ...
    'element_type','turbofan_afterburning', ...
    'max_output',128992, ...
    'position',[-1.5 0 0], ...
    'direction',[1 0 0], ...
    'fuel_rate',2.5, ...
    'thrust_model',@(thr,M,alt,V) F100ThrustModel(thr, M, alt, V));

n_cs = numel(ac.control_surfaces);
n_pe = numel(ac.propulsive_elements);
n_total = n_cs + n_pe;

%% Mission Segments
segments = struct([]);

segments(1).type     = 'hold_time';
segments(1).altitude = 0;
segments(1).mach     = 0.30;
segments(1).duration = 30;
segments(1).gamma    = 0;

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
