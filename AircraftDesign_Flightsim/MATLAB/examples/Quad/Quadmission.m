clear; clc; close all

%% PATH BOOTSTRAP
% Ensures the shared class library and example helper functions are on the
% MATLAB path regardless of the current working directory or session state.
this_script_dir = fileparts(mfilename('fullpath'));
matlab_root_dir = fileparts(fileparts(this_script_dir));
addpath(genpath(fullfile(matlab_root_dir,'classes')));
addpath(genpath(fullfile(matlab_root_dir,'examples')));

modelName = 'FlightSim_Aircraft';

if bdIsLoaded(modelName)
    set_param(modelName,'SimulationCommand','stop');
end

%% Aircraft object

quad = Aircraft();
quad.geometry.set_reference_point([0 0 0]);

cfg = quad.get_configurator();

%% Geometry%%Dummy

quad.geometry.wing_area = 0.10;
quad.geometry.wing_span = 0.10;
quad.geometry.mean_aerodynamic_chord = 0.10;
quad.geometry.ref_area = quad.geometry.wing_area;
quad.geometry.ref_span = quad.geometry.wing_span;
quad.geometry.ref_chord = quad.geometry.mean_aerodynamic_chord;
quad.geometry.ref_point = quad.geometry.get_reference_point().';

%% Mass properties

I_body = diag([0.02 0.02 0.04]);
cfg.add_component('name','airframe', 'type','airframe', 'mass',1.5, 'position',[0 0 0], 'inertia',I_body);

%% Aerodynamics

cfg.add_aero_solver(CoefficientAerodynamics(@QuadAeroLookup), 'name','quad_aero');

%% Ground contact

quad.ground_k = 1.0e4;
quad.ground_c = 1.0e3;
ground_k = quad.ground_k;
ground_c = quad.ground_c;

%% Rotor properties

d_rotor = 0.254;
CT = 0.02186;
CQ = CT*0.016;
omega_max = 5000*(2*pi/60);
L_arm = 0.25;

rotor_config = {'Omega_FR',[ L_arm, L_arm,0], 1; 'Omega_RR',[-L_arm, L_arm,0],-1; 'Omega_RL',[-L_arm,-L_arm,0], 1; 'Omega_FL',[ L_arm,-L_arm,0],-1};

[~,~,~,rho0] = atmosisa(0);
A_rotor = pi*(d_rotor/2)^2;
R_rotor = d_rotor/2;
T_max = CT*rho0*A_rotor*(omega_max*R_rotor)^2;

for k = 1:size(rotor_config,1)
    cfg.add_propulsive_element('name',rotor_config{k,1}, 'element_type','rotor_variable', 'max_output',T_max, 'position',rotor_config{k,2}, 'mount_euler',[0 0 0], 'direction',[0 0 -1], 'fuel_rate',0);
    quad.propulsive_elements{k}.set_variable_rotor_params(d_rotor,CT,CQ,omega_max,rotor_config{k,3});
end

%% Derived parameters

[mass,~,I_mat] = quad.compute_total_mass_properties([]);
W = mass*9.80665;

n_cs = numel(quad.control_surfaces);
n_pe = numel(quad.propulsive_elements);
n_total = n_cs + n_pe;

%% Flight profile

dt_fc = 0.01;
t_takeoff = 8;
t_hover = 5;
t_land = 5;
sim_time = t_takeoff + t_hover + t_land;

[~,~,~,rho10] = atmosisa(10);
K_thr10 = CT*rho10*A_rotor*R_rotor^2;
u_hover10 = max(0,min(1,sqrt((W/4)/K_thr10)/omega_max));

time_vector = (0:dt_fc:sim_time)';
z_ref = zeros(size(time_vector));

i1 = time_vector < t_takeoff;
i2 = time_vector >= t_takeoff & time_vector < (t_takeoff + t_hover);
i3 = time_vector >= (t_takeoff + t_hover);

z_ref(i1) = -10*(time_vector(i1)/t_takeoff);
z_ref(i2) = -10;
tau = (time_vector(i3) - t_takeoff - t_hover)/t_land;
z_ref(i3) = -10*(1 - tau);

z_ref_data = [time_vector z_ref];

%% Initial state

x0 = zeros(12,1);
Initialpos = [0 0 0];
InitialVel = [0 0 0];
InitialOri = [0 0 0];
InitialRot = [0 0 0];

%% Workspace assignment

assignin('base','Initialpos',Initialpos);
assignin('base','InitialVel',InitialVel);
assignin('base','InitialOri',InitialOri);
assignin('base','InitialRot',InitialRot);
assignin('base','ac',quad);
assignin('base','n_total',n_total);
assignin('base','n_cs',n_cs);
assignin('base','n_pe',n_pe);
assignin('base','sim_stop_time',sim_time);
assignin('base','dt_fc',dt_fc);
assignin('base','disable_rate_limiting',true);
assignin('base','initial_state',x0);
assignin('base','u_hover10',u_hover10);
assignin('base','z_ref_data',z_ref_data);
assignin('base','ground_k',ground_k);
assignin('base','ground_c',ground_c);

%% Simulink setup

if ~bdIsLoaded(modelName)
    load_system(modelName);
end

try
    set_param([modelName '/Level-2 MATLAB S-Function'],'MFunctionName','aircraft_sfunc_generalized');
catch
end

try
    set_param([modelName '/Level-2 MATLAB S-Function'],'InputPort(2).Dimensions',num2str(n_total));
catch
end

try
    set_param([modelName '/6DOF (Euler Angles)'],'mass',num2str(mass));
catch
end

try
    set_param([modelName '/6DOF (Euler Angles)'],'inertia',mat2str(I_mat));
catch
end

try
    set_param([modelName '/6DOF (Euler Angles)'],'Xe_0','[0 0 0]');
catch
end

set_param(modelName,'Solver','ode4','FixedStep',num2str(dt_fc),'StopTime',num2str(sim_time),'SaveOutput','on','OutputSaveName','out');
