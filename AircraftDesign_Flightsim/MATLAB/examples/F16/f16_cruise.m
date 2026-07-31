%% F16_Steady_Level_Trim.m
% Finds one steady, wings-level, constant-speed trim condition.
%
% Trim unknowns:
%   alpha
%   stabilator deflection
%   throttle
%
% Trim equations:
%   Fx = 0
%   Fz = 0
%   My = 0
%
% Steady level assumptions:
%   beta  = 0
%   gamma = 0
%   phi   = 0
%   p     = 0
%   q     = 0
%   r     = 0
%   aileron = 0
%   rudder  = 0

clear;
clc;
close all;
clear classes;
clear functions;
rehash toolboxcache;

%% PATH BOOTSTRAP
% Ensures the shared class library and example helper functions are on the
% MATLAB path regardless of the current working directory or session state.
this_script_dir = fileparts(mfilename('fullpath'));
matlab_root_dir = fileparts(fileparts(this_script_dir));
addpath(genpath(fullfile(matlab_root_dir,'classes')));
addpath(genpath(fullfile(matlab_root_dir,'examples')));

%% Aircraft geometry

S=300*0.09290304;
b=30*0.3048;
c=11.3201786951274*0.3048;

ac=Aircraft();

ac.geometry.set_reference_geometry(S,b,c);
ac.geometry.set_reference_point([0 0 0]);

%% Aircraft mass properties

lbfN=4.4482216152605;
g=9.80665;

Wlbf=0.899666962714079*31377;
Wfuel=max(Wlbf-20680.7005781593,0);
Weng=0.199*23770;
Wair=Wlbf-Wfuel-Weng;

mair=Wair*lbfN/g;
mfuel=Wfuel*lbfN/g;
meng=Weng*lbfN/g;

scale=Wlbf/20500;
k=1.3558179483314;

I=scale*[9496*k,0,982*k; 0,55814*k,0; 982*k,0,63100*k];

%% Aircraft reference frames

ac.add_frame("fuel_tank", "body", [0;0;0], @(x) eye(3));

ac.add_frame("engine", "body", [0;0;0], @(x) eye(3));

ac.add_frame("cg", "body", [0;0;0], @(x) eye(3));

ac.add_frame("gravity_cg", "body", [0;0;0], @(x) ReferenceFrame.ned_to_body_dcm(x));

%% Aircraft components

airframe=Component("airframe", mair, [0;0;0], I, ac.get_frame("body"), "airframe");

fuel=Component("fuel", mfuel, [0;0;0], zeros(3), ac.get_frame("fuel_tank"), "fuel");

engine_component=Component("engine", meng, [0;0;0], zeros(3), ac.get_frame("engine"), "engine");

ac.add_component(airframe);
ac.add_component(fuel);
ac.add_component(engine_component);

%% Control surfaces

aileron=ControlSurface("aileron", "aileron", "primary", [1 0 0], deg2rad(21.5), deg2rad(-21.5), 0, 0, 0);

stabilator=ControlSurface("stabilator", "elevator", "primary", [0 1 0], deg2rad(25), deg2rad(-25), 0, 0, 0);

rudder=ControlSurface("rudder", "rudder", "primary", [0 0 1], deg2rad(30), deg2rad(-30), 0, 0, 0);

ac.add_control_surface(aileron);
ac.add_control_surface(stabilator);
ac.add_control_surface(rudder);

%% Propulsion model

engine=BrandtAfterburningEngine("F16A_engine", ac.get_frame("engine"), [1;0;0]);

engine.set_rating_mode("dry");

ac.add_propulsive_element(engine);

engine_load_solver=PropulsionLoadSolver(engine, ac.get_frame("engine"));

engine_component.add_load_source(engine_load_solver);

%% Aerodynamic model

aero=CoefficientAerodynamics(@F16ABrandtLookup);

aero_load_solver=AeroLoadSolver(aero, ac.geometry, ac, ac.get_frame("body"));

ac.add_load_source(aero_load_solver);

%% Gravity model

gravity_load_solver=GravityLoadSolver(ac, ac.get_frame("gravity_cg"));

ac.add_load_source(gravity_load_solver);

%% Update CG reference

[mass,cg,Icg]=ac.compute_total_mass_properties();

ac.update_frame_position("cg",cg);
ac.update_frame_position("gravity_cg",cg);
ac.set_reference_frame("cg");

weight_N=mass*ac.g;

fprintf("\n=== F-16 MODEL ===\n");
fprintf("Mass   : %.2f kg\n",mass);
fprintf("Weight : %.2f N\n",weight_N);
fprintf("CG     : [% .4f % .4f % .4f] m\n",cg);

%% Flight condition

altitude_m=9000;

% Select one speed-input mode:
speed_input_mode="mach";

mach=0.80;
velocity_mps=250;

condition=struct();
condition.altitude_m=altitude_m;

if speed_input_mode=="mach"

    condition.mach=mach;

    [~,speed_of_sound_mps,~,~]= atmosisa(altitude_m);

    commanded_velocity_mps= mach*speed_of_sound_mps;

elseif speed_input_mode=="velocity"

    condition.velocity_mps=velocity_mps;
    commanded_velocity_mps=velocity_mps;

    [~,speed_of_sound_mps,~,~]= atmosisa(altitude_m);

    mach=velocity_mps/speed_of_sound_mps;

else

    error("speed_input_mode must be ""mach"" or ""velocity"".");

end

fprintf("\n=== REQUESTED STEADY LEVEL CONDITION ===\n");
fprintf("Altitude : %.2f m\n",altitude_m);
fprintf("Velocity : %.3f m/s\n",commanded_velocity_mps);
fprintf("Mach     : %.4f\n",mach);
fprintf("Engine   : dry\n");

%% Trim variables and equations
% Same "level" formulation as brand4pointtrim.m's F-16 baseline trim:
% alpha/stabilator/throttle solved against the exact EOM equality
% constraints via PerformanceAnalysis (GenericTrimSolver's weighted-
% penalty formulation is no longer used by any example script).

level_spec=struct();

level_spec.mode="level";

level_spec.variables=["alpha"; "stabilator"; "throttle"];

level_spec.initial_guess=[deg2rad(5); deg2rad(1.5); 0.45];

level_spec.lb=[deg2rad(-5); deg2rad(-25); 0];

level_spec.ub=[deg2rad(20); deg2rad(25); 1];

%% Fixed steady-level values

level_spec.fixed=struct();

level_spec.fixed.beta=0;
level_spec.fixed.gamma=0;

level_spec.fixed.phi=0;
level_spec.fixed.psi=0;

level_spec.fixed.control_roll=0;
level_spec.fixed.control_yaw=0;

level_spec.reference_frame_name="cg";

level_spec.residual_scale=[ac.g;ac.g;1];

%% Solver settings

solver_cfg=struct();

solver_cfg.residual_tolerance=1e-5;
solver_cfg.inequality_tolerance=1e-6;
solver_cfg.equality_tolerance=1e-5;
solver_cfg.optimality_tolerance=1e-5;

solver_cfg.fmincon_options=optimoptions( ...
    "fmincon", ...
    "Algorithm","sqp", ...
    "Display","iter", ...
    "MaxIterations",500, ...
    "MaxFunctionEvaluations",10000, ...
    "OptimalityTolerance",1e-10, ...
    "ConstraintTolerance",1e-10, ...
    "StepTolerance",1e-12, ...
    "FunctionTolerance",1e-12, ...
    "FiniteDifferenceType","central", ...
    "FiniteDifferenceStepSize",1e-5);

%% Solve steady-level trim

performance=ac.get_performance();

fprintf("\n=== SOLVING STEADY LEVEL TRIM ===\n");

[x_trim,u_trim,trim_converged,trim_info]= performance.solve_trim(condition, level_spec, solver_cfg);

%% Extract trim variables

trim_alpha_rad=x_trim(6);
trim_alpha_rad=atan2(x_trim(6),x_trim(4));

trim_alpha_deg=rad2deg(trim_alpha_rad);

trim_theta_rad=x_trim(8);
trim_theta_deg=rad2deg(trim_theta_rad);

trim_beta_rad=asin(max(-1,min(1,x_trim(5)/norm(x_trim(4:6)))));

trim_beta_deg=rad2deg(trim_beta_rad);

stabilator_index= ac.control.get_index_by_name("stabilator");

trim_stabilator_rad=u_trim(stabilator_index);
trim_stabilator_deg=rad2deg(trim_stabilator_rad);

trim_throttle=engine.throttle;

trim_velocity_mps=norm(x_trim(4:6));

[~,trim_speed_of_sound_mps,~,~]= atmosisa(altitude_m);

trim_mach= trim_velocity_mps/trim_speed_of_sound_mps;

%% Extract force and moment residuals

trim_force_N=trim_info.loads.F_total_body_N;
trim_moment_Nm=trim_info.loads.M_total_cg_body_Nm;

Fx_N=trim_force_N(1);
Fy_N=trim_force_N(2);
Fz_N=trim_force_N(3);

Mx_Nm=trim_moment_Nm(1);
My_Nm=trim_moment_Nm(2);
Mz_Nm=trim_moment_Nm(3);

Fx_normalized=Fx_N/weight_N;
Fz_normalized=Fz_N/weight_N;
My_normalized=My_Nm/(weight_N*c);

%% Display trim solution

fprintf("\n=== STEADY LEVEL TRIM RESULT ===\n");

fprintf("Converged       : %d\n",trim_converged);
fprintf("Exit flag       : %d\n",trim_info.exitflag);
fprintf("Residual norm   : %.6e\n",trim_info.residual_norm);

fprintf("\nFlight condition\n");
fprintf("Altitude        : %.3f m\n",altitude_m);
fprintf("Velocity        : %.6f m/s\n",trim_velocity_mps);
fprintf("Mach            : %.6f\n",trim_mach);

fprintf("\nTrim state\n");
fprintf("Alpha           : %.6f deg\n",trim_alpha_deg);
fprintf("Beta            : %.6f deg\n",trim_beta_deg);
fprintf("Theta           : %.6f deg\n",trim_theta_deg);
fprintf("Phi             : %.6f deg\n",rad2deg(x_trim(7)));
fprintf("Psi             : %.6f deg\n",rad2deg(x_trim(9)));
fprintf("p               : %.6e rad/s\n",x_trim(10));
fprintf("q               : %.6e rad/s\n",x_trim(11));
fprintf("r               : %.6e rad/s\n",x_trim(12));

fprintf("\nTrim controls\n");
fprintf("Aileron         : %.6f deg\n",rad2deg(aileron.deflection));
fprintf("Stabilator      : %.6f deg\n",trim_stabilator_deg);
fprintf("Rudder          : %.6f deg\n",rad2deg(rudder.deflection));
fprintf("Throttle        : %.8f\n",trim_throttle);

fprintf("\nTotal loads about CG\n");
fprintf("Fx              : % .6e N\n",Fx_N);
fprintf("Fy              : % .6e N\n",Fy_N);
fprintf("Fz              : % .6e N\n",Fz_N);
fprintf("Mx              : % .6e N-m\n",Mx_Nm);
fprintf("My              : % .6e N-m\n",My_Nm);
fprintf("Mz              : % .6e N-m\n",Mz_Nm);

fprintf("\nNormalized trim residuals\n");
fprintf("Fx/W            : % .6e\n",Fx_normalized);
fprintf("Fz/W            : % .6e\n",Fz_normalized);
fprintf("My/(W*c)        : % .6e\n",My_normalized);

%% Result tables

TrimCondition=table( ...
    altitude_m, ...
    trim_velocity_mps, ...
    trim_mach, ...
    trim_alpha_deg, ...
    trim_beta_deg, ...
    trim_theta_deg, ...
    trim_stabilator_deg, ...
    trim_throttle, ...
    trim_converged, ...
    trim_info.residual_norm, ...
    'VariableNames',{ ...
    'Altitude_m', ...
    'Velocity_mps', ...
    'Mach', ...
    'Alpha_deg', ...
    'Beta_deg', ...
    'Theta_deg', ...
    'Stabilator_deg', ...
    'Throttle', ...
    'Converged', ...
    'ResidualNorm'});

TrimLoads=table(Fx_N, Fy_N, Fz_N, Mx_Nm, My_Nm, Mz_Nm, Fx_normalized, Fz_normalized, My_normalized, 'VariableNames',{'Fx_N', 'Fy_N', 'Fz_N', 'Mx_Nm', 'My_Nm', 'Mz_Nm', 'Fx_over_W', 'Fz_over_W', 'My_over_Wc'});

fprintf("\n=== TRIM CONDITION TABLE ===\n");
disp(TrimCondition);

fprintf("\n=== TRIM LOAD TABLE ===\n");
disp(TrimLoads);

%% Save useful workspace outputs

SteadyLevelTrim=struct();

SteadyLevelTrim.condition=condition;
SteadyLevelTrim.x_trim=x_trim;
SteadyLevelTrim.u_trim=u_trim;
SteadyLevelTrim.converged=trim_converged;
SteadyLevelTrim.info=trim_info;
SteadyLevelTrim.condition_table=TrimCondition;
SteadyLevelTrim.load_table=TrimLoads;
SteadyLevelTrim.mass_kg=mass;
SteadyLevelTrim.weight_N=weight_N;
SteadyLevelTrim.cg_m=cg;
SteadyLevelTrim.inertia_kgm2=Icg;

fprintf("\n=== STEADY LEVEL TRIM COMPLETE ===\n");
fprintf("Workspace outputs:\n");
fprintf("  SteadyLevelTrim\n");
fprintf("  TrimCondition\n");
fprintf("  TrimLoads\n");
fprintf("  x_trim\n");
fprintf("  u_trim\n");
