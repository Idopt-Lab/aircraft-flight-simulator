%% Quadcopter Hover Simulation - Multirotor Physics
clear; clc; close all;

if bdIsLoaded('flightsimquad')
    set_param('flightsimquad', 'SimulationCommand', 'stop');
end

%% Aircraft Configuration
quad = Aircraft();

m_total = 1.5;
Ixx = 0.02;
Iyy = 0.02;
Izz = 0.04;
cg = [0 0 0];

quad.mass.set_mass_properties(m_total, Ixx, Iyy, Izz, 0, 0, 0, cg);
quad.mass.set_fuel_capacity(0, m_total);
quad.mass.set_fuel_mass(0);

quad.geometry.wing_area = 0.1;
quad.geometry.wing_span = 0.5;
quad.geometry.mean_aerodynamic_chord = 0.1;

quad.aero.set_lookup(@QuadAeroLookup);

%% Multirotor Configuration
cfg = quad.get_configurator();

L_arm = 0.25;
rotor_diameter = 0.254;
rotor_speed_hover = 5000;
CT = 0.0048;
CP = 0.0002;

rotor_positions = [
    L_arm,  L_arm,  0;
   -L_arm,  L_arm,  0;
   -L_arm, -L_arm,  0;
    L_arm, -L_arm,  0
];

rotor_directions = [1; -1; 1; -1];

cfg.add_propulsive_element('name','multirotor','element_type','multirotor','max_output',24.0, ...
    'position',[0 0 0],'direction',[0 0 -1],'fuel_rate',0);

quad.propulsive_elements{1}.set_multirotor_params(4, rotor_diameter, rotor_speed_hover, ...
    CT, CP, rotor_positions, rotor_directions);

%% Flight Parameters
h_target = -10.0;
sim_time = 45;
dt_fc = 0.01;

%% Hover Controller Gains
Kp_h = 0.8;
Ki_h = 0.03;
Kd_h = 0.5;

%% Calculate Hover Throttle
mass = quad.mass.get_total_mass();
I_mat = quad.mass.get_inertia_matrix();
g = 9.80665;

T_total_hover = mass * g;
T_max_total = quad.propulsive_elements{1}.max_output;
u_hover_single = T_total_hover / T_max_total;

u_hover = u_hover_single;

fprintf('\n=== QUADCOPTER MULTIROTOR SETUP ===\n');
fprintf('Mass: %.2f kg\n', mass);
fprintf('Weight: %.2f N\n', mass * g);
fprintf('Total max thrust: %.2f N\n', T_max_total);
fprintf('Hover throttle: %.1f%%\n', u_hover*100);
fprintf('Target altitude: %.1f m\n', h_target);
fprintf('Rotor diameter: %.3f m\n', rotor_diameter);
fprintf('Hover RPM: %.0f\n', rotor_speed_hover);

%% Initial Conditions
x0 = zeros(12,1);
x0(3) = -0.1;
x0(4) = 1e-3;

Initialpos = x0(1:3)';
InitialVel = x0(4:6)';
InitialOri = x0(7:9)';
InitialRot = x0(10:12)';

%% Control Schedule
time_vector = (0:dt_fc:sim_time)';
control_input_data = [time_vector, repmat(u_hover, length(time_vector), 1)];

%% Workspace Assignment
assignin('base','quad', quad);
assignin('base','ac', quad);
assignin('base','AIRCRAFT_OBJ', quad);
assignin('base','initial_state', x0);
assignin('base','Initialpos', Initialpos);
assignin('base','InitialVel', InitialVel);
assignin('base','InitialOri', InitialOri);
assignin('base','InitialRot', InitialRot);
assignin('base','mass', double(mass));
assignin('base','Inertia', double(I_mat));
assignin('base','inertia', double(I_mat));
assignin('base','Kp_h', Kp_h);
assignin('base','Ki_h', Ki_h);
assignin('base','Kd_h', Kd_h);
assignin('base','h_cmd', h_target);
assignin('base','u_hover', u_hover);
assignin('base','n_total', 1);
assignin('base','n_cs', 0);
assignin('base','n_pe', 1);
assignin('base','control_input_data', control_input_data);
assignin('base','sim_stop_time', sim_time);
assignin('base','ground_k', 10000);
assignin('base','ground_c', 1000);
assignin('base','dt_fc', dt_fc);

%% Update S-function Input Port
modelName = 'flightsimquad';
if ~bdIsLoaded(modelName)
    load_system(modelName);
end

sfunc_block = [modelName '/Level-2 MATLAB S-Function'];
try
    set_param(sfunc_block, 'InputPort(2).Dimensions', '1');
catch
    warning('Could not update S-function input dimensions automatically');
end

set_param(modelName, 'Solver', 'ode4');
set_param(modelName, 'FixedStep', num2str(dt_fc));
set_param(modelName, 'StopTime', num2str(sim_time));
set_param(modelName, 'SaveOutput', 'on');
set_param(modelName, 'OutputSaveName', 'out');

fprintf('\n=== READY TO SIMULATE ===\n');
fprintf('Simulation time: %.1f s\n', sim_time);
fprintf('NOTE: S-function now expects 1 control input (collective throttle)\n');