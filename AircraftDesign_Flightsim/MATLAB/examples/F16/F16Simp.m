clearvars; clc;
clear aircraft_sfunc_generalized

alt_trim_m       = 9144;      % 30,000 ft
mach_trim        = 0.60;
sim_time         = 100.0;
dt_fc            = 0.01;

disable_rate_limiting = true;

%% Aircraft object

ac = Aircraft();

%% Reference geometry

S_ref = 27.87;
b_ref = 9.14;
c_ref = 3.45;

ac.geometry.set_reference_geometry(S_ref, b_ref, c_ref);

nozzle_cg       = 193.0;
nozzle_aerorp   = 189.5;
nozzle_fueltank = 174.4;
nozzle_engine   = 0.0;

cg_to_aerorp = (nozzle_aerorp  - nozzle_cg) * 0.0254 * [1;0;0];
cg_to_tank   = (nozzle_fueltank - nozzle_cg) * 0.0254 * [1;0;0];
cg_to_engine = (nozzle_engine   - nozzle_cg) * 0.0254 * [1;0;0];

ac.geometry.set_reference_point(cg_to_aerorp.');

%% Frames

ac.set_body_frame("body");
ac.set_reference_frame("body");

add_or_update_frame(ac, "aero_ref", "body", cg_to_aerorp, @(x) wind_to_body_dcm(x));
add_or_update_frame(ac, "fuel_tank", "body", cg_to_tank, @(x) eye(3));
add_or_update_frame(ac, "engine", "body", cg_to_engine, @(x) eye(3));
add_or_update_frame(ac, "cg", "body", [0;0;0], @(x) eye(3));
add_or_update_frame(ac, "gravity_cg", "body", [0;0;0], @(x) ned_to_body_dcm(x));

%% Mass components

slug_ft2_to_kg_m2 = 1.35582;

Ixx = 9496  * slug_ft2_to_kg_m2;
Iyy = 55814 * slug_ft2_to_kg_m2;
Izz = 63100 * slug_ft2_to_kg_m2;
Ixz = -982  * slug_ft2_to_kg_m2;

I_body = [ Ixx  0   -Ixz;
           0    Iyy  0;
          -Ixz  0    Izz ];

empty_mass = 9300;
fuel_mass  = 2400;
tank_mass  = fuel_mass / 2;

airframe = Component( ...
    "airframe", ...
    empty_mass, ...
    [0;0;0], ...
    I_body, ...
    ac.get_frame("body"), ...
    "airframe");

tank1 = Component( ...
    "tank1", ...
    tank_mass, ...
    [0;0;0], ...
    zeros(3,3), ...
    ac.get_frame("fuel_tank"), ...
    "fuel");

tank2 = Component( ...
    "tank2", ...
    tank_mass, ...
    [0;0;0], ...
    zeros(3,3), ...
    ac.get_frame("fuel_tank"), ...
    "fuel");

engine_comp = Component( ...
    "engine", ...
    0, ...
    [0;0;0], ...
    zeros(3,3), ...
    ac.get_frame("engine"), ...
    "engine");

airframe.set_metadata("description", "F16 airframe");
tank1.set_metadata("fuel_mass", tank_mass);
tank2.set_metadata("fuel_mass", tank_mass);
engine_comp.set_metadata("max_thrust", 128992);
engine_comp.set_metadata("engine_type", "F100-PW-229");

ac.add_component(airframe);
ac.add_component(tank1);
ac.add_component(tank2);
ac.add_component(engine_comp);

%% Controls

aileron = ControlSurface( ...
    "aileron", ...
    "aileron", ...
    "primary", ...
    [1 0 0], ...
    deg2rad(21.5), ...
    deg2rad(-21.5), ...
    0.05, 0, 0);

elevator = ControlSurface( ...
    "elevator", ...
    "elevator", ...
    "primary", ...
    [0 1 0], ...
    deg2rad(25), ...
    deg2rad(-25), ...
    0, -0.10, 0);

rudder = ControlSurface( ...
    "rudder", ...
    "rudder", ...
    "primary", ...
    [0 0 1], ...
    deg2rad(30), ...
    deg2rad(-30), ...
    0, 0, -0.08);

ac.add_control_surface(aileron);
ac.add_control_surface(elevator);
ac.add_control_surface(rudder);

%% Propulsion

engine_pe = TurbofanPropulsion( ...
    "F100_PW_229", ...
    ac.get_frame("engine"), ...
    [1;0;0], ...
    128992, ...
    2.5, ...
    @(thr,M,alt,V) F100ThrustModel(thr, M, alt, V));

ac.add_propulsive_element(engine_pe);

prop_solver = PropulsionLoadSolver(engine_pe, ac.get_frame("engine"));
engine_comp.add_load_source(prop_solver);

%% Aerodynamics and gravity loads

aero_model = CoefficientAerodynamics(@F16Lookup);

aero_solver = AeroLoadSolver( ...
    aero_model, ...
    ac.geometry, ...
    ac, ...
    ac.get_frame("aero_ref"));

gravity_solver = GravityLoadSolver( ...
    ac, ...
    ac.get_frame("gravity_cg"));

ac.add_load_source(aero_solver);
ac.add_load_source(gravity_solver);

%% Sync mass and CG frames

[m_total, cg_total, I_total] = ac.compute_total_mass_properties(); %#ok<ASGLU>
ac.update_frame_position("cg", cg_total);
ac.update_frame_position("gravity_cg", cg_total);

W = m_total * ac.g;

fprintf("\n=== AIRCRAFT READY ===\n");
fprintf("Mass    : %.2f kg\n", m_total);
fprintf("Weight  : %.2f N\n", W);
fprintf("CG body : [% .4f % .4f % .4f] m\n", cg_total);

%% Initial guess

[~, a_trim, ~, ~] = atmosisa(alt_trim_m);
V_trim = mach_trim * a_trim;

alpha0 = deg2rad(7.0);
elev0  = deg2rad(-0.5);
thr0   = 0.05;

x_guess = zeros(12,1);
x_guess(3) = -alt_trim_m;
x_guess(4) = V_trim * cos(alpha0);
x_guess(5) = 0;
x_guess(6) = V_trim * sin(alpha0);
x_guess(7) = 0;
x_guess(8) = alpha0;
x_guess(9) = 0;

ac.state.set_full_state(x_guess);
ac.control_surfaces(1).set_deflection(0);
ac.control_surfaces(2).set_deflection(elev0);
ac.control_surfaces(3).set_deflection(0);
ac.propulsive_elements{1}.set_throttle(thr0);
ac.sync_control_vector_from_components();

%% Trim

fprintf("\n=== RUNNING CRUISE TRIM ===\n");
fprintf("Altitude : %.1f m\n", alt_trim_m);
fprintf("Mach     : %.3f\n", mach_trim);

ac.set_reference_frame("cg");

solver = ac.get_trim_solver();
solver.trim_tolerance = 1e-5;
solver.max_iterations = 15000;
solver.initial_guess  = [alpha0; elev0; thr0];
solver.debug_failures = true;
solver.use_fmincon = exist("fmincon", "file") == 2;

[x_trim, u_trim, converged, info] = solver.solve_cruise_trim(alt_trim_m, mach_trim); %#ok<ASGLU>

if converged
    fprintf("\n=== TRIM CONVERGED ===\n");
else
    fprintf("\n=== TRIM NOT CONVERGED: using returned best solution ===\n");
end

solver.print_summary();

%% Final verification about CG

ac.state.set_full_state(x_trim);
ac.set_controls_from_vector(u_trim);

[m_trim, cg_trim, I_trim] = ac.compute_total_mass_properties(x_trim); %#ok<ASGLU>
ac.update_frame_position("cg", cg_trim);
ac.update_frame_position("gravity_cg", cg_trim);
ac.set_reference_frame("cg");

[F_trim, M_trim, fuel_flow_trim] = ac.compute_total_loads(x_trim, u_trim); %#ok<ASGLU>

W_trim = m_trim * ac.g;
cbar_trim = ac.geometry.mean_aerodynamic_chord;

fprintf("\n=== TRIM VERIFICATION ABOUT CG ===\n");
fprintf("Mass       : %.2f kg\n", m_trim);
fprintf("CG         : [% .4f % .4f % .4f] m\n", cg_trim);
fprintf("Alpha      : %.4f deg\n", rad2deg(atan2(x_trim(6), x_trim(4))));
fprintf("Theta      : %.4f deg\n", rad2deg(x_trim(8)));
fprintf("Aileron    : %.4f deg\n", rad2deg(u_trim(1)));
fprintf("Elevator   : %.4f deg\n", rad2deg(u_trim(2)));
fprintf("Rudder     : %.4f deg\n", rad2deg(u_trim(3)));
fprintf("Throttle   : %.6f\n", u_trim(4));
fprintf("F [N]      : [% .4e % .4e % .4e]\n", F_trim);
fprintf("M [Nm]     : [% .4e % .4e % .4e]\n", M_trim);
fprintf("Fx/W       : %.6e\n", F_trim(1) / W_trim);
fprintf("Fz/W       : %.6e\n", F_trim(3) / W_trim);
fprintf("My/Wc      : %.6e\n", M_trim(2) / (W_trim * cbar_trim));

%% Simulink 

n_cs    = numel(ac.control_surfaces);
n_pe    = numel(ac.propulsive_elements);
n_total = n_cs + n_pe;

t_vec = (0:dt_fc:sim_time).';

x0 = x_trim(:);
u0 = u_trim(:);

Initialpos = x0(1:3);
InitialVel = x0(4:6);
InitialOri = x0(7:9);
InitialRot = x0(10:12);

u_trim_mat = repmat(u0.', numel(t_vec), 1);

control_input_data = struct();
control_input_data.time = t_vec;
control_input_data.signals.values = u_trim_mat;
control_input_data.signals.dimensions = n_total;

ground_k = ac.ground_k;
ground_c = ac.ground_c;

assignin("base", "ac", ac);
assignin("base", "n_cs", n_cs);
assignin("base", "n_pe", n_pe);
assignin("base", "n_total", n_total);

assignin("base", "initial_state", x0);
assignin("base", "Initialpos", Initialpos);
assignin("base", "InitialVel", InitialVel);
assignin("base", "InitialOri", InitialOri);
assignin("base", "InitialRot", InitialRot);

assignin("base", "u_trim", u0);
assignin("base", "control_input_data", control_input_data);

assignin("base", "sim_stop_time", sim_time);
assignin("base", "dt_fc", dt_fc);
assignin("base", "disable_rate_limiting", disable_rate_limiting);
assignin("base", "ground_k", ground_k);

fprintf("n_cs             = %d\n", n_cs);
fprintf("n_pe             = %d\n", n_pe);
fprintf("n_total          = %d\n", n_total);
fprintf("sim_stop_time    = %.2f s\n", sim_time);
fprintf("dt_fc            = %.4f s\n", dt_fc);
fprintf("rate limiting off= %d\n", disable_rate_limiting);
fprintf("x0 alpha         = %.4f deg\n", rad2deg(atan2(x0(6), x0(4))));

%% Local helper functions

function add_or_update_frame(ac, name, parent_name, r_parent, dcm_fn)
    if ~ac.has_frame(name)
        ac.add_frame(name, parent_name, r_parent(:), dcm_fn);
    else
        ac.update_frame_position(name, r_parent(:));
        ac.update_frame_orientation(name, dcm_fn);
    end
end

function C = wind_to_body_dcm(x)
    u = x(4);
    v = x(5);
    w = x(6);

    V = sqrt(u^2 + v^2 + w^2);

    if V < 1e-9
        C = eye(3);
        return;
    end

    alpha = atan2(w, u);
    beta  = asin(max(-1, min(1, v / V)));

    ca = cos(alpha);
    sa = sin(alpha);
    cb = cos(beta);
    sb = sin(beta);

    C = [ca*cb, -ca*sb, -sa;
         sb,     cb,      0;
         sa*cb, -sa*sb,  ca];
end

function C = ned_to_body_dcm(x)
    phi   = x(7);
    theta = x(8);
    psi   = x(9);

    cp = cos(phi);
    sp = sin(phi);

    ct = cos(theta);
    st = sin(theta);

    cs = cos(psi);
    ss = sin(psi);

    C = [ ct*cs,              ct*ss,             -st;
          sp*st*cs-cp*ss,     sp*st*ss+cp*cs,    sp*ct;
          cp*st*cs+sp*ss,     cp*st*ss-sp*cs,    cp*ct ];
end