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

%% ========================================================================
% F-16A FORMULA MODEL
%
% Frame/load architecture:
%   Component -> LoadSolver -> ForceMoment -> ReferenceFrame
%
% Conventions:
%   x = [X;Y;Z;u;v;w;phi;theta;psi;p;q;r]
%   NED position, so altitude = -x(3)
%   Body axes: x forward, y right, z down
%
% Important:
%   Propulsion models return local force and intrinsic moment only.
%   ReferenceFrame.transform_FM_to() adds r x F exactly once.
% ========================================================================

%% ========================================================================
% 1. AIRCRAFT REFERENCE GEOMETRY
% ========================================================================

aircraft_name = "F-16A Formula Model";

S_ref = 300.0 * 0.09290304;
b_ref = 30.0  * 0.3048;
c_ref = 11.3201786951274 * 0.3048;

ac = Aircraft();

ac.geometry.set_reference_geometry(S_ref, b_ref, c_ref);
ac.geometry.set_reference_point([0 0 0]);

ac.set_body_frame("body");
ac.set_reference_frame("body");

%% ========================================================================
% 2. MASS AND INERTIA
% ========================================================================

lbf_to_N = 4.4482216152605;
g0 = 9.80665;

takeoff_weight_lbf = 31377.0;
landing_weight_lbf = 20680.7005781593;
weight_fraction = 0.899666962714079;

aircraft_weight_lbf = weight_fraction * takeoff_weight_lbf;

fuel_weight_lbf = max(aircraft_weight_lbf - landing_weight_lbf, 0);

engine_weight_lbf = 0.199 * 23770.0;

airframe_weight_lbf = aircraft_weight_lbf - fuel_weight_lbf - engine_weight_lbf;

if airframe_weight_lbf <= 0
    error("F16ABrandt:InvalidMassBreakdown", "The selected mass breakdown is not physical.");
end

airframe_mass_kg = airframe_weight_lbf * lbf_to_N / g0;
fuel_mass_kg     = fuel_weight_lbf     * lbf_to_N / g0;
engine_mass_kg   = engine_weight_lbf   * lbf_to_N / g0;
target_mass_kg   = aircraft_weight_lbf * lbf_to_N / g0;

% Reference F-16 inertia tensor scaled to selected aircraft mass.
reference_inertia_weight_lbf = 20500.0;
inertia_scale = aircraft_weight_lbf / reference_inertia_weight_lbf;

slug_ft2_to_kg_m2 = 1.3558179483314;

Ixx_ref = 9496  * slug_ft2_to_kg_m2;
Iyy_ref = 55814 * slug_ft2_to_kg_m2;
Izz_ref = 63100 * slug_ft2_to_kg_m2;
Ixz_ref = -982  * slug_ft2_to_kg_m2;

I_body = inertia_scale * [Ixx_ref, 0,      -Ixz_ref; 0,       Iyy_ref, 0; -Ixz_ref, 0,      Izz_ref];

%% ========================================================================
% 3. PHYSICAL LOCATIONS AND REFERENCE FRAMES
% ========================================================================

% All positions are vectors from the body-frame origin to the relevant
% frame origin, expressed in body coordinates.
%
% Update these values when your aircraft datum and component locations are
% known. Keeping them at zero is valid but produces zero position moment.
fuel_tank_position_body_m = [0;0;0];
engine_position_body_m    = [0;0;0];

add_or_update_frame(ac, "fuel_tank", "body", fuel_tank_position_body_m, @(x) eye(3));

add_or_update_frame(ac, "engine", "body", engine_position_body_m, @(x) eye(3));

add_or_update_frame(ac, "cg", "body", [0;0;0], @(x) eye(3));

% This frame is located at the CG but oriented like NED.
% Its child-to-parent DCM is therefore NED -> body.
add_or_update_frame(ac, "gravity_cg", "body", [0;0;0], @(x) ned_to_body_dcm(x));

%% ========================================================================
% 4. AIRCRAFT COMPONENTS
% ========================================================================

airframe = Component("airframe", airframe_mass_kg, [0;0;0], I_body, ac.get_frame("body"), "airframe");

fuel = Component("fuel", fuel_mass_kg, [0;0;0], zeros(3,3), ac.get_frame("fuel_tank"), "fuel");

engine_component = Component("engine", engine_mass_kg, [0;0;0], zeros(3,3), ac.get_frame("engine"), "engine");

airframe.set_metadata("description", "Formula-based F-16A airframe");

fuel.set_metadata("fuel_mass_kg", fuel_mass_kg);

engine_component.set_metadata("engine_type", "Dry/afterburning formula model");

ac.add_component(airframe);
ac.add_component(fuel);
ac.add_component(engine_component);

%% ========================================================================
% 5. CONTROL SURFACES
% ========================================================================

% F16ABrandtLookup evaluates aerodynamic control effects. Keep dCl/dCm/dCn
% zero here to avoid duplicate control moment increments.

aileron = ControlSurface("aileron", "aileron", "primary", [1 0 0], deg2rad(21.5), deg2rad(-21.5), 0, 0, 0);

stabilator = ControlSurface("stabilator", "elevator", "primary", [0 1 0], deg2rad(25), deg2rad(-25), 0, 0, 0);

rudder = ControlSurface("rudder", "rudder", "primary", [0 0 1], deg2rad(30), deg2rad(-30), 0, 0, 0);

ac.add_control_surface(aileron);
ac.add_control_surface(stabilator);
ac.add_control_surface(rudder);

%% ========================================================================
% 6. PROPULSION
% ========================================================================

engine = BrandtAfterburningEngine("F16A_engine", ac.get_frame("engine"), [1;0;0]);

engine.set_rating_mode("dry");

% Registers throttle in the global control vector.
ac.add_propulsive_element(engine);

% Attaches the actual propulsion load source to the engine component.
% The solver frame is exactly the propulsion application frame.
engine_solver = PropulsionLoadSolver(engine, ac.get_frame("engine"));

engine_component.add_load_source(engine_solver);

%% ========================================================================
% 7. AERODYNAMICS AND GRAVITY
% ========================================================================

% This self-contained adapter converts coefficient data into BODY-AXIS
% force and BODY-AXIS moment. It avoids mixing wind-axis force with
% body-axis moments in one ForceMoment object.
aero_model = BodyAerodynamics(@f16_body_load_lookup);

aero_solver = AeroLoadSolver(aero_model, ac.geometry, ac, ac.get_frame("body"));

gravity_solver = GravityLoadSolver(ac, ac.get_frame("gravity_cg"));

ac.add_load_source(aero_solver);
ac.add_load_source(gravity_solver);

%% ========================================================================
% 8. TOTAL MASS PROPERTIES AND CG FRAMES
% ========================================================================

[m_total, cg_total, I_total] = ac.compute_total_mass_properties();

ac.update_frame_position("cg", cg_total);
ac.update_frame_position("gravity_cg", cg_total);
ac.set_reference_frame("cg");

fprintf("\n============================================================\n");
fprintf("F-16A AIRCRAFT MODEL\n");
fprintf("============================================================\n");
fprintf("Mass            : %.6f kg\n", m_total);
fprintf("Target mass     : %.6f kg\n", target_mass_kg);
fprintf("Mass error      : %.6e kg\n", m_total-target_mass_kg);
fprintf("CG              : [% .6f % .6f % .6f] m\n", cg_total);
fprintf("Ixx / Iyy / Izz : %.3f / %.3f / %.3f kg-m^2\n", I_total(1,1), I_total(2,2), I_total(3,3));
fprintf("S / b / c       : %.6f / %.6f / %.6f m\n", S_ref, b_ref, c_ref);
fprintf("Controls        : %d\n", numel(ac.control_surfaces));
fprintf("Engines         : %d\n", numel(ac.propulsive_elements));

%% ========================================================================
% 9. FRAME AND MOMENT-ARM VERIFICATION
% ========================================================================

verify_gravity_frame(ac, m_total);
verify_propulsion_transport(ac, engine);

fprintf("Frame and moment-arm verification passed.\n");

%% ========================================================================
% 10. LONGITUDINAL LEVEL-FLIGHT TRIM
% ========================================================================

perf = PerformanceAnalysis(ac);

trim_condition = struct();
trim_condition.altitude_m = 9144;
trim_condition.mach = 0.60;

trim_spec = struct();
trim_spec.mode = "level";

trim_spec.variables = ["alpha"; "control_pitch"; "throttle"];

trim_spec.initial_guess = [deg2rad(5); deg2rad(1.7); 0.45];

trim_spec.lb = [deg2rad(-5); deg2rad(-25); 0];

trim_spec.ub = [deg2rad(20); deg2rad(25); 1];

trim_spec.fixed = struct('beta', 0, 'phi', 0, 'psi', 0, 'gamma', 0, 'p', 0, 'q', 0, 'r', 0, 'control_roll', 0, 'control_yaw', 0);

trim_spec.reference_frame_name = "cg";

% Longitudinal rigid-body residual scaling:
%   udot/g, wdot/g, qdot/(1 rad/s^2)
trim_spec.residual_scale = [ac.g; ac.g; 1];

solver_cfg = struct();
solver_cfg.residual_tolerance = 1e-5;
solver_cfg.inequality_tolerance = 1e-6;
solver_cfg.equality_tolerance = 1e-5;
solver_cfg.optimality_tolerance = 1e-5;

solver_cfg.fmincon_options = optimoptions( ...
    "fmincon", ...
    "Algorithm", "sqp", ...
    "Display", "iter", ...
    "MaxIterations", 1500, ...
    "MaxFunctionEvaluations", 30000, ...
    "OptimalityTolerance", 1e-8, ...
    "ConstraintTolerance", 1e-8, ...
    "StepTolerance", 1e-10, ...
    "FiniteDifferenceType", "central", ...
    "FiniteDifferenceStepSize", 1e-4);

[x_trim, u_trim, trim_converged, trim_info] = perf.solve_trim(trim_condition, trim_spec, solver_cfg);

perf_cfg = struct();
perf_cfg.available_throttle = 1;
perf_cfg.CL_max = 0.9840156989811455;
perf_cfg.CL_min = -0.80;
perf_cfg.max_Mach = 2.0;
perf_cfg.propulsion_type = "thrust";

trim_point = perf.evaluate_trim_point(x_trim, u_trim, trim_condition, perf_cfg, trim_info.extras);

C_trim = F16ABrandtLookup(x_trim, u_trim, ac.geometry);

% Restore solved operating point after performance evaluation.
ac.state.set_full_state(x_trim);
ac.set_controls_from_vector(u_trim);

fprintf("\n============================================================\n");
fprintf("LONGITUDINAL TRIM RESULT\n");
fprintf("============================================================\n");
fprintf("Converged       : %d\n", trim_converged);
fprintf("Residual norm   : %.6e\n", trim_info.residual_norm);
fprintf("Residual inf    : %.6e\n", trim_info.residual_inf);
fprintf("Altitude        : %.3f m\n", trim_point.altitude_m);
fprintf("Mach            : %.8f\n", trim_point.Mach);
fprintf("Velocity        : %.6f m/s\n", trim_point.V_mps);
fprintf("Alpha           : %.6f deg\n", trim_point.alpha_deg);
fprintf("Pitch angle     : %.6f deg\n", trim_point.theta_deg);
fprintf("Stabilator      : %.6f deg\n", trim_point.elevator_deg);
fprintf("Throttle        : %.8f\n", u_trim(end));
fprintf("CL              : %.8f\n", trim_point.CL);
fprintf("CD              : %.8f\n", trim_point.CD);
fprintf("L/D             : %.8f\n", trim_point.L_over_D);
fprintf("Fuel flow       : %.8f kg/s\n", trim_point.fuel_flow_trim);

aero_valid = get_struct_field_or(C_trim, "valid", true);
fprintf("Aero valid      : %d\n", logical(aero_valid));

if isstruct(engine.last_debug) && isfield(engine.last_debug, "thrust_N")
    fprintf("Engine thrust   : %.6f kN\n", engine.last_debug.thrust_N/1000);
end

if isfield(trim_info, "residual")
    fprintf("\nLongitudinal residual vector:\n");
    disp(trim_info.residual(:));
end

if ~trim_converged
    warning("F16ABrandt:TrimNotConverged", "The longitudinal trim solution did not meet the tolerance.");
end

if ~logical(aero_valid)
    warning("F16ABrandt:AeroModelInvalid", "The trim point lies outside the aerodynamic validity envelope.");
end

%% ========================================================================
% 11. SIMULATION INITIALIZATION
% ========================================================================

F16A_Simulation = struct();

F16A_Simulation.aircraft = ac;
F16A_Simulation.engine = engine;

F16A_Simulation.x0 = x_trim;
F16A_Simulation.u0 = u_trim;

F16A_Simulation.trim_condition = trim_condition;
F16A_Simulation.trim_point = trim_point;
F16A_Simulation.trim_converged = trim_converged;

F16A_Simulation.engine_mode = engine.rating_mode;
F16A_Simulation.t_start_s = 0;
F16A_Simulation.t_end_s = 30;

assignin("base", "ac", ac);
assignin("base", "perf", perf);
assignin("base", "engine", engine);

assignin("base", "x_trim", x_trim);
assignin("base", "u_trim", u_trim);
assignin("base", "trim_point", trim_point);
assignin("base", "trim_info", trim_info);

assignin("base", "F16A_Simulation", F16A_Simulation);

%% ========================================================================
% 12. SAVE
% ========================================================================

F16A_Trim_Results = struct();

F16A_Trim_Results.aircraft_name = aircraft_name;
F16A_Trim_Results.mass_kg = m_total;
F16A_Trim_Results.cg_m = cg_total;
F16A_Trim_Results.inertia_kgm2 = I_total;

F16A_Trim_Results.condition = trim_condition;
F16A_Trim_Results.converged = trim_converged;
F16A_Trim_Results.x_trim = x_trim;
F16A_Trim_Results.u_trim = u_trim;
F16A_Trim_Results.trim_info = trim_info;
F16A_Trim_Results.trim_point = trim_point;
F16A_Trim_Results.aero_coefficients = C_trim;

script_path = mfilename("fullpath");
script_dir = fileparts(script_path);

if strlength(string(script_dir)) == 0 || ~isfolder(script_dir)
    script_dir = pwd;
end

results_dir = fullfile(script_dir, "results");

if ~isfolder(results_dir)
    mkdir(results_dir);
end

results_file = fullfile(results_dir, "F16A_Longitudinal_Trim.mat");

save(results_file, "F16A_Trim_Results", "-v7.3");

assignin("base", "F16A_Trim_Results", F16A_Trim_Results);
assignin("base", "F16A_results_file", string(results_file));

fprintf("\n=== F-16A MODEL READY ===\n");
fprintf("Aircraft, trim state and controls exported to workspace.\n");
fprintf("Saved: %s\n", results_file);

%% ========================================================================
% LOCAL FUNCTIONS
% ========================================================================

function add_or_update_frame(ac, name, parent_name, translation, dcm_handle)
% ADD_OR_UPDATE_FRAME Create or consistently update a frame.

translation = translation(:);

if numel(translation) ~= 3 || ~isreal(translation) || any(~isfinite(translation))

    error("F16ABrandt:InvalidFrameTranslation", "Frame translation must be a finite real three-vector.");
end

if ~isa(dcm_handle, "function_handle")
    error("F16ABrandt:InvalidFrameDCM", "Frame DCM definition must be a function handle.");
end

if isempty(parent_name)
    requested_parent = [];
else
    requested_parent = ac.get_frame(parent_name);
end

if ac.has_frame(name)

    frame = ac.get_frame(name);

    if ~isequal(frame.parent, requested_parent)
        error("F16ABrandt:FrameParentMismatch", ['Existing frame "%s" has a different parent. ' 'Delete and recreate the frame to change its parent.'], char(name));
    end

    ac.update_frame_position(name, translation);
    ac.update_frame_orientation(name, dcm_handle);

else

    ac.add_frame(name, parent_name, translation, dcm_handle);
end
end

function C_ned_to_body = ned_to_body_dcm(x)
% NED_TO_BODY_DCM 3-2-1 Euler DCM.
%
%   v_body = C_ned_to_body * v_ned

x = x(:);

if numel(x) < 9
    error("F16ABrandt:InvalidStateForGravityDCM", "State vector must contain Euler angles x(7:9).");
end

phi   = x(7);
theta = x(8);
psi   = x(9);

cphi = cos(phi);
sphi = sin(phi);

cth = cos(theta);
sth = sin(theta);

cpsi = cos(psi);
spsi = sin(psi);

C_ned_to_body = [cth*cpsi, cth*spsi, -sth; sphi*sth*cpsi-cphi*spsi, sphi*sth*spsi+cphi*cpsi, sphi*cth; cphi*sth*cpsi+sphi*spsi, cphi*sth*spsi-sphi*cpsi, cphi*cth];
end

function d = f16_body_load_lookup(x, u, geom)
% F16_BODY_LOAD_LOOKUP Convert F16 coefficients to body-axis loads.
%
% The returned structure retains all configuration metadata from
% F16ABrandtLookup, including takeoff_params and landing_params.

x = x(:);
u = u(:);

if numel(x) < 12
    error("F16ABrandt:InvalidAeroState", "Aerodynamic state vector must contain 12 elements.");
end

C = F16ABrandtLookup(x, u, geom);
d = C;

u_b = x(4);
v_b = x(5);
w_b = x(6);

p = x(10);
q = x(11);
r = x(12);

V = norm([u_b;v_b;w_b]);

if V < 1e-9
    d.Fx = 0;
    d.Fy = 0;
    d.Fz = 0;
    d.Mx = 0;
    d.My = 0;
    d.Mz = 0;
    return;
end

alpha = atan2(w_b, u_b);
beta = atan2(v_b, sqrt(u_b^2+w_b^2));

altitude_m = max(-x(3), 0);
[~, ~, ~, rho] = isa1976(altitude_m);

qbar = 0.5*rho*V^2;

[S, b, cbar] = geom.get_reference_geometry();

S = max(S, 1e-9);
b = max(b, 1e-9);
cbar = max(cbar, 1e-9);

CL = get_struct_field_or(C, "CL", 0);
CD = get_struct_field_or(C, "CD", 0);
CY = get_struct_field_or(C, "CY", 0);
Cl = get_struct_field_or(C, "Cl", 0);
Cm = get_struct_field_or(C, "Cm", 0);
Cn = get_struct_field_or(C, "Cn", 0);

% Optional static derivative fields.
CY = CY + get_struct_field_or(C, "CYb", 0)*beta;
Cl = Cl + get_struct_field_or(C, "Clb", 0)*beta;
Cm = Cm + get_struct_field_or(C, "Cma", 0)*alpha;
Cn = Cn + get_struct_field_or(C, "Cnb", 0)*beta;

% Optional nondimensional rate derivatives.
p_hat = p*b/(2*V);
q_hat = q*cbar/(2*V);
r_hat = r*b/(2*V);

CL = CL + get_struct_field_or(C, "CLq", 0)*q_hat;
CD = CD + get_struct_field_or(C, "CDq", 0)*q_hat;

CY = CY + get_struct_field_or(C, "CYp", 0)*p_hat + get_struct_field_or(C, "CYr", 0)*r_hat;

Cl = Cl + get_struct_field_or(C, "Clp", 0)*p_hat + get_struct_field_or(C, "Clr", 0)*r_hat;

Cm = Cm + get_struct_field_or(C, "Cmq", 0)*q_hat;

Cn = Cn + get_struct_field_or(C, "Cnp", 0)*p_hat + get_struct_field_or(C, "Cnr", 0)*r_hat;

L = CL*qbar*S;
D = CD*qbar*S;
Y = CY*qbar*S;

Mx = Cl*qbar*S*b;
My = Cm*qbar*S*cbar;
Mz = Cn*qbar*S*b;

% Wind-axis aerodynamic force.
F_wind = [-D;Y;-L];

ca = cos(alpha);
sa = sin(alpha);
cb = cos(beta);
sb = sin(beta);

% Wind -> body DCM:
%   F_body = C_wind_to_body * F_wind
C_wind_to_body = [ca*cb, -ca*sb, -sa; sb,     cb,      0; sa*cb, -sa*sb,  ca];

F_body = C_wind_to_body*F_wind;

d.CL = CL;
d.CD = CD;
d.CY = CY;
d.Cl = Cl;
d.Cm = Cm;
d.Cn = Cn;

d.Fx = F_body(1);
d.Fy = F_body(2);
d.Fz = F_body(3);

% Aerodynamic moment coefficients are interpreted in body axes.
d.Mx = Mx;
d.My = My;
d.Mz = Mz;
end

function verify_gravity_frame(ac, mass_kg)
% VERIFY_GRAVITY_FRAME Confirm NED gravity transforms correctly to body.

x_test = zeros(12,1);
x_test(7) = deg2rad(4);
x_test(8) = deg2rad(10);
x_test(9) = deg2rad(-7);

gravity_frame = ac.get_frame("gravity_cg");
body_frame = ac.get_body_frame();

F_ned = [0;0;mass_kg*ac.g];

F_body_actual = gravity_frame.transform_vector_to(body_frame, F_ned, x_test);

F_body_expected = ned_to_body_dcm(x_test)*F_ned;

tol = 1e-9*max(1, norm(F_body_expected));

assert(norm(F_body_actual-F_body_expected) <= tol, "Gravity frame transformation is incorrect.");
end

function verify_propulsion_transport(ac, engine)
% VERIFY_PROPULSION_TRANSPORT Confirm M = C*M0 + r x F exactly once.

x_test = zeros(12,1);
x_test(4) = 100;

body_frame = ac.get_body_frame();
cg_frame = ac.get_frame("cg");
engine_frame = engine.frame;

old_throttle = engine.throttle;
cleanup = onCleanup(@() engine.set_throttle(old_throttle)); %#ok<NASGU>

engine.set_throttle(0.5);

[F_local, M_intrinsic_local, ~] = engine.get_FM(x_test, []);

[F_body_actual, M_cg_actual] = engine_frame.transform_FM_to(body_frame, F_local, M_intrinsic_local, x_test, cg_frame);

C_engine_to_body = engine_frame.get_dcm_to(body_frame, x_test);

F_body_expected = C_engine_to_body*F_local;
M_intrinsic_body = C_engine_to_body*M_intrinsic_local;

r_body_to_engine = engine_frame.get_position_to(body_frame, x_test);
r_body_to_cg = cg_frame.get_position_to(body_frame, x_test);
r_cg_to_engine = r_body_to_engine-r_body_to_cg;

M_cg_expected = M_intrinsic_body + cross(r_cg_to_engine, F_body_expected);

tol_F = 1e-9*max(1, norm(F_body_expected));
tol_M = 1e-9*max(1, norm(M_cg_expected));

assert(norm(F_body_actual-F_body_expected) <= tol_F, "Propulsion force transformation is incorrect.");

assert(norm(M_cg_actual-M_cg_expected) <= tol_M, "Propulsion moment-arm transport is incorrect.");
end

function value = get_struct_field_or(s, field_name, default_value)
% GET_STRUCT_FIELD_OR Read a nonempty structure field or use a default.

field_name = char(field_name);

if isstruct(s) && isfield(s, field_name) && ~isempty(s.(field_name))

    value = s.(field_name);
else
    value = default_value;
end
end
