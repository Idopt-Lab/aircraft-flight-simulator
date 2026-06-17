if bdIsLoaded('flightsim7')
    set_param('flightsim7','SimulationCommand','stop');
end

clear; clc

ac = Aircraft();

%% Geometry / reference positions

CG_FROM_NOSE = 1.90;

cg_to_nose    = -CG_FROM_NOSE;
cg_to_wing_le = 1.20 - CG_FROM_NOSE;
cg_to_wing_ac = 1.57 - CG_FROM_NOSE;
cg_to_prop    = 1.68 - CG_FROM_NOSE;
cg_to_fuel    = 1.70 - CG_FROM_NOSE;

ac.geometry.wing_area = 16.17;
ac.geometry.wing_span = 10.98;
ac.geometry.mean_aerodynamic_chord = 1.50;
ac.geometry.ref_area = 16.17;
ac.geometry.ref_span = 10.98;
ac.geometry.ref_chord = 1.50;
ac.geometry.ref_point = [cg_to_wing_ac,0,0];

%% Frames

ac.add_frame("earth", [], [0;0;0], @(x) eye(3));
ac.add_frame("body", "earth", [0;0;0], @(x) eye(3));

ac.set_body_frame("body");
ac.set_reference_frame("body");

ac.add_frame("aero_ref",  "body", [cg_to_wing_ac;0;0], @(x) wind_to_body_dcm(x));
ac.add_frame("fuel_tank", "body", [cg_to_fuel;0;0],    @(x) eye(3));
ac.add_frame("propeller", "body", [cg_to_prop;0;1.26], @(x) eye(3));
ac.add_frame("cg",        "body", [0;0;0],             @(x) eye(3));

%% Mass properties

empty_mass = 767;
fuel_mass  = 90;

Ixx = 1285;
Iyy = 1824;
Izz = 2666;

I_body = [Ixx 0 0;
          0 Iyy 0;
          0 0 Izz];

airframe = Component("airframe", empty_mass, [0;0;0], I_body, ac.get_frame("body"));

% Important: cg_local is relative to fuel_tank frame.
% Since fuel_tank frame is already located at tank CG, cg_local should be zero.
fuel = Component("fuel", fuel_mass, [0;0;0], zeros(3,3), ac.get_frame("fuel_tank"));

ac.add_component(airframe);
ac.add_component(fuel);

%% Aerodynamics

datcom_fp = "C:\Users\naman\Downloads\Aircrfatdummy\dummy3\aircraft-flight-simulator\AircraftDesign_Flightsim\Matlab\Example\cessna\cessna.out";

if isfile(datcom_fp)
    aero_model = CoefficientAerodynamics(c172datcom(datcom_fp));
else
    error('DATCOM file not found: %s', datcom_fp);
end

aero_solver = AeroLoadSolver(aero_model, ac.geometry, ac, ac.get_frame("aero_ref"));
ac.add_load_source(aero_solver);

%% Controls

cfg = ac.get_configurator();

cfg.add_control_surface( ...
    'name','aileron', ...
    'surface_type','aileron', ...
    'classification','primary', ...
    'axis',[1 0 0], ...
    'max_deflection',deg2rad(20), ...
    'min_deflection',deg2rad(-20), ...
    'dCl',0,'dCm',0,'dCn',0);

cfg.add_control_surface( ...
    'name','elevator', ...
    'surface_type','elevator', ...
    'classification','primary', ...
    'axis',[0 1 0], ...
    'max_deflection',deg2rad(28), ...
    'min_deflection',deg2rad(-26), ...
    'dCl',0,'dCm',0,'dCn',0);

cfg.add_control_surface( ...
    'name','rudder', ...
    'surface_type','rudder', ...
    'classification','primary', ...
    'axis',[0 0 1], ...
    'max_deflection',deg2rad(30), ...
    'min_deflection',deg2rad(-30), ...
    'dCl',0,'dCm',0,'dCn',0);

%% Propulsion

prop_pe = PropellerPropulsion( ...
    "O360", ...
    ac.get_frame("propeller"), ...
    [1;0;0], ...
    134228, ...
    0.0089, ...
    []);

prop_pe.set_propeller_params( ...
    1.905, ...
    1.219, ...
    2, ...
    0.80, ...
    [0.10,-0.05], ...
    [0.05,0.02], ...
    0.05);

ac.add_propulsive_element(prop_pe);

engine_c = Component("engine", 0, [0;0;0], zeros(3,3), ac.get_frame("propeller"));

prop_solver = PropulsionLoadSolver(prop_pe, ac.get_frame("propeller"));
engine_c.add_load_source(prop_solver);

ac.add_component(engine_c);

%% Setup mass / CG

n_cs = numel(ac.control_surfaces);
n_pe = numel(ac.propulsive_elements);
n_total = n_cs + n_pe; %#ok<NASGU>

[m, cg, I_total] = ac.compute_total_mass_properties([]);
W = m * 9.80665;
cbar = ac.geometry.mean_aerodynamic_chord;

ac.update_frame_position("cg", cg);

fprintf('\n=== CESSNA SETUP ===\n');
fprintf('Reference frame : %s\n', ac.reference_frame_name);
fprintf('Empty mass      : %.1f kg\n', empty_mass);
fprintf('Fuel mass       : %.1f kg\n', fuel.mass);
fprintf('Total mass      : %.1f kg\n', m);
fprintf('Weight          : %.1f N\n', W);
fprintf('CG position     : [%.6f %.6f %.6f] m\n', cg);
fprintf('Nose            : [%.2f 0 0] m from CG\n', cg_to_nose);
fprintf('Wing LE         : [%.2f 0 0] m from CG\n', cg_to_wing_le);
fprintf('Wing AC         : [%.2f 0 0] m from CG\n', cg_to_wing_ac);
fprintf('Propeller       : [%.2f 0 1.26] m from CG\n', cg_to_prop);
fprintf('Fuel tank       : [%.2f 0 0] m from CG\n', cg_to_fuel);

%% Test point

x_test = zeros(12,1);
x_test(3) = -1000;
x_test(4) = 50*cosd(5);
x_test(5) = 0;
x_test(6) = 50*sind(5);
x_test(7) = 0;
x_test(8) = deg2rad(5);
x_test(9) = 0;

ac.state.set_full_state(x_test);
ac.control_surfaces(1).set_deflection(0);
ac.control_surfaces(2).set_deflection(0);
ac.control_surfaces(3).set_deflection(0);
ac.propulsive_elements{1}.set_throttle(0.65);
ac.sync_control_vector_from_components();

[F_cg, M_cg] = ac.compute_total_loads_about_cg(x_test, ac.get_control_vector());

fprintf('\n=== TEST POINT ===\n');
fprintf('At alpha = %.1f deg, elev = 0 deg, thr = 0.65\n', rad2deg(atan2(x_test(6),x_test(4))));
fprintf('F(CG)    : [%.1f %.1f %.1f] N\n', F_cg);
fprintf('M(CG)    : [%.1f %.1f %.1f] N-m\n', M_cg);
fprintf('Fx/W     : %.4f\n', F_cg(1)/W);
fprintf('Fz/W     : %.4f\n', F_cg(3)/W);
fprintf('My/(Wc)  : %.4f\n', M_cg(2)/(W*cbar));

%% Trim

% If you want trim about CG, set reference to CG before trim.
% If you want trim about body origin, keep reference_frame = body.
ac.set_reference_frame("cg");

solver = ac.get_trim_solver();
solver.trim_tolerance = 1e-6;
solver.max_iterations = 10000;
solver.use_fmincon = true;
solver.initial_guess = [deg2rad(3); deg2rad(0); 0.45];

alt_test = 1500;
mach_test = 0.15;

fprintf('\nRunning Cessna trim...\n');
[x_trim, u_trim, converged] = solver.solve_cruise_trim(alt_test, mach_test);

if converged
    fprintf('\n=== TRIM CONVERGED ===\n');
    solver.print_summary();

    ac.state.set_full_state(x_trim);
    ac.set_controls_from_vector(u_trim);

    [m_trim, cg_trim, ~] = ac.compute_total_mass_properties(x_trim);
    ac.update_frame_position("cg", cg_trim);

    [F_trim, M_trim] = ac.compute_total_loads(x_trim, u_trim);
    W_trim = m_trim * 9.80665;

    fprintf('\n=== TRIM CHECK ABOUT CURRENT REF: %s ===\n', ac.reference_frame_name);
    fprintf('Alpha    : %.2f deg\n', rad2deg(atan2(x_trim(6),x_trim(4))));
    fprintf('Theta    : %.2f deg\n', rad2deg(x_trim(8)));
    fprintf('Elevator : %.2f deg\n', rad2deg(u_trim(2)));
    fprintf('Throttle : %.3f\n', u_trim(n_cs+1));
    fprintf('CG       : [%.4f %.4f %.4f] m\n', cg_trim);
    fprintf('F        : [%.3f %.3f %.3f] N\n', F_trim);
    fprintf('M        : [%.3f %.3f %.3f] N-m\n', M_trim);
    fprintf('Fx/W     : %.6f\n', F_trim(1)/W_trim);
    fprintf('Fz/W     : %.6f\n', F_trim(3)/W_trim);
    fprintf('My/(Wc)  : %.6f\n', M_trim(2)/(W_trim*cbar));
else
    error('TRIM FAILED TO CONVERGE');
end

%% Load contributions

fprintf('\n=== LOAD CONTRIBUTIONS ABOUT CURRENT REF: %s ===\n', ac.reference_frame_name);

body_ref = ac.get_body_frame();
ref_frame = ac.get_reference_frame();

for k = 1:numel(ac.load_sources)
    fm = ac.load_sources{k}.get_force_moment(x_trim,u_trim);
    fm_b = fm.transform_to(body_ref,ref_frame,x_trim);

    fprintf('Aircraft source %d %-25s F=[% .3e % .3e % .3e] M=[% .3e % .3e % .3e]\n', ...
        k, class(ac.load_sources{k}), fm_b.F, fm_b.M);
end

for k = 1:numel(ac.components)
    fm = ac.components{k}.compute_force_moment(x_trim,u_trim,body_ref,ref_frame);

    fprintf('Component %d %-20s F=[% .3e % .3e % .3e] M=[% .3e % .3e % .3e]\n', ...
        k, char(ac.components{k}.name), fm.F, fm.M);
end

fm_g = ac.compute_gravity_force_moment(x_trim, body_ref, ref_frame);
fprintf('Gravity                    F=[% .3e % .3e % .3e] M=[% .3e % .3e % .3e]\n', ...
    fm_g.F, fm_g.M);

%% Direct checks

fprintf('\n=== AERO DIRECT CHECK ===\n');
[F_aero, M_aero, coeff] = aero_model.get_FM(x_trim, u_trim, ac.geometry, ac);
fprintf('F aero local/wind [N]   : [% .3e % .3e % .3e]\n', F_aero);
fprintf('M aero local/wind [N-m] : [% .3e % .3e % .3e]\n', M_aero);
disp(coeff);

fprintf('\n=== PROPULSION DIRECT CHECK ===\n');
[F_prop, M_prop, fuel_flow] = ac.propulsive_elements{1}.get_FM(x_trim, u_trim);
fprintf('F prop local direct [N]  : [% .3e % .3e % .3e]\n', F_prop);
fprintf('M prop local direct [N-m]: [% .3e % .3e % .3e]\n', M_prop);
fprintf('fuel flow                : %.6e\n', fuel_flow);

%% CG migration

fprintf('\n=== CG MIGRATION TEST ===\n');

fuel_levels = [180,142,110,85,50,20,0];

for f = fuel_levels
    fuel.set_mass_properties(f, [0;0;0], zeros(3,3));
    [~, cg_current, ~] = ac.compute_total_mass_properties([]);
    cg_current = cg_current(:);

    fprintf('Fuel %3.0f kg -> CG = [%+.4f %+.4f %+.4f] m\n', ...
        f, cg_current(1), cg_current(2), cg_current(3));
end

fuel.set_mass_properties(fuel_mass, [0;0;0], zeros(3,3));

%% Local functions

function C = wind_to_body_dcm(x)
    u = x(4);
    v = x(5);
    w = x(6);

    V = sqrt(u^2 + v^2 + w^2);

    if V < 1e-9
        C = eye(3);
        return;
    end

    alpha = atan2(w,u);
    beta  = asin(max(-1,min(1,v/V)));

    ca = cos(alpha); sa = sin(alpha);
    cb = cos(beta);  sb = sin(beta);

    C = [ca*cb, -ca*sb, -sa;
         sb,     cb,      0;
         sa*cb, -sa*sb,  ca];
end