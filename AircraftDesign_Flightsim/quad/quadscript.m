clear; clc; close all;

quad = Aircraft();

h0 = 5;
x0 = zeros(12,1);
x0(3) = -h0;
x0(4) = 1e-3;

m_total = 1.5;
Ixx = 0.02; Iyy = 0.02; Izz = 0.04;
Ixy = 0; Ixz = 0; Iyz = 0;
cg = [0 0 0];

quad.mass.set_mass_properties(m_total, Ixx, Iyy, Izz, Ixy, Ixz, Iyz, cg);
quad.mass.set_fuel_capacity(0, m_total);
quad.mass.set_fuel_mass(0);
quad.mass.set_fuel_burn_rate(0);

S_ref = 0.1;
b_ref = 0.5;
c_ref = 0.1;

quad.geometry.wing_area = S_ref;
quad.geometry.wing_span = b_ref;
quad.geometry.mean_aerodynamic_chord = c_ref;
quad.geometry.ref_area = S_ref;
quad.geometry.ref_span = b_ref;
quad.geometry.ref_chord = c_ref;

quad.aero.set_lookup(@QuadAeroLookup);

cfg = quad.get_configurator();
L_arm = 0.25;
T_max_per_motor = 6.0;

cfg.add_propulsive_element('name','motor_1','element_type','electric', ...
    'max_output',T_max_per_motor,'position',[ L_arm,  L_arm, 0],'direction',[0 0 -1],'fuel_rate',0);
cfg.add_propulsive_element('name','motor_2','element_type','electric', ...
    'max_output',T_max_per_motor,'position',[-L_arm,  L_arm, 0],'direction',[0 0 -1],'fuel_rate',0);
cfg.add_propulsive_element('name','motor_3','element_type','electric', ...
    'max_output',T_max_per_motor,'position',[-L_arm, -L_arm, 0],'direction',[0 0 -1],'fuel_rate',0);
cfg.add_propulsive_element('name','motor_4','element_type','electric', ...
    'max_output',T_max_per_motor,'position',[ L_arm, -L_arm, 0],'direction',[0 0 -1],'fuel_rate',0);

dt_fc = 0.01;
sim_time = 45;

Kp_h = 1.0;
Ki_h = 0.05;
Kd_h = 0.8;
h_cmd = 5.0;

mass = quad.mass.get_total_mass();
I_mat = quad.mass.get_inertia_matrix();

mass = double(mass);
Inertia = double(I_mat);

g = 9.80665;
n_mot = numel(quad.propulsive_elements);
T_total_hover = mass * g;
T_each = T_total_hover / n_mot;

u_hover = zeros(1, n_mot);
for k = 1:n_mot
    T_max_k = quad.propulsive_elements{k}.max_output;
    u_hover(k) = min(1.0, max(0.0, T_each / T_max_k));
end

t_vec = (0:dt_fc:sim_time)';
control_input_data = [t_vec, repmat(u_hover, numel(t_vec), 1)];

Initialpos = x0(1:3)';
InitialVel = x0(4:6)';
InitialOri = x0(7:9)';
InitialRot = x0(10:12)';

assignin('base', 'AIRCRAFT_OBJ', quad);
assignin('base', 'Initialpos', Initialpos);
assignin('base', 'InitialVel', InitialVel);
assignin('base', 'InitialOri', InitialOri);
assignin('base', 'InitialRot', InitialRot);
assignin('base', 'InitialState', x0);
assignin('base', 'dt_fc', dt_fc);
assignin('base', 'mass', mass);
assignin('base', 'Inertia', Inertia);
assignin('base', 'inertia', Inertia);
assignin('base', 'Kp_h', Kp_h);
assignin('base', 'Ki_h', Ki_h);
assignin('base', 'Kd_h', Kd_h);
assignin('base', 'h_cmd', h_cmd);
assignin('base', 'u_hover', u_hover);
assignin('base', 'control_input_data', control_input_data);
assignin('base', 'n_total', n_mot);
assignin('base', 'ground_k', 10000);
assignin('base', 'ground_c', 1000);

modelName = 'flightsimquad';

if ~bdIsLoaded(modelName)
    load_system(modelName);
end

set_param(modelName, 'Solver', 'ode4');
set_param(modelName, 'FixedStep', num2str(dt_fc));
set_param(modelName, 'StopTime', num2str(sim_time));
set_param(modelName, 'SaveOutput', 'on');
set_param(modelName, 'OutputSaveName', 'out');

out = sim(modelName);