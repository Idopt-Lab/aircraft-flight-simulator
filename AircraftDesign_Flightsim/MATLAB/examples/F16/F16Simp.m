
clear; clc

if bdIsLoaded('flightsim7')
    set_param('flightsim7', 'SimulationCommand', 'stop');
end

%% Aircraft Initialization
ac = Aircraft();

% All positions are defined relative to the aircraft reference point.
% Here, the reference point is the initial CG.
ac.set_reference_point([0 0 0]);

%% Reference Frame Definition
% Nose-referenced coordinates:
%   position_cg_ref = position_nose_ref - NOSE_REF_CG
%
% All positions are full 3x1 body-axis vectors:
%   x forward, y right, z down

NOSE_REF_CG = [4.88; 0.00; 0.00];   % CG location in nose-referenced frame [m]

nose_nose_ref      = [0.00; 0.00; 0.00];
wing_apex_nose_ref = [3.80; 0.00; 0.00];
wing_ac_nose_ref   = [4.20; 0.00; 0.00];
tank_1_nose_ref    = [4.00; 0.00; 0.00];
tank_2_nose_ref    = [5.50; 0.00; 0.00];
engine_nose_ref    = [6.00; 0.00; 0.00];

cg_to_nose      = nose_nose_ref      - NOSE_REF_CG;
cg_to_wing_apex = wing_apex_nose_ref - NOSE_REF_CG;
cg_to_wing_ac   = wing_ac_nose_ref   - NOSE_REF_CG;
cg_to_tank_1    = tank_1_nose_ref    - NOSE_REF_CG;
cg_to_tank_2    = tank_2_nose_ref    - NOSE_REF_CG;
cg_to_engine    = engine_nose_ref    - NOSE_REF_CG;

%% Mass Properties
empty_mass = 9300;

Ixx = 12875;
Iyy = 75674;
Izz = 85552;
Ixy = 0;
Ixz = 1331;
Iyz = 0;

ac.mass.set_mass_properties(empty_mass, Ixx, Iyy, Izz, Ixy, Ixz, Iyz, [0 0 0]);

ac.mass.add_fuel_tank(cg_to_tank_1.', 1800, 'distributed');
ac.mass.add_fuel_tank(cg_to_tank_2.', 1200, 'distributed');
ac.mass.set_fuel_mass(2400);
ac.set_reference_point([0 0 0]);

fprintf('Aircraft reference point: [%.3f %.3f %.3f] m\n', ac.get_reference_point());
fprintf('Current CG              : [%.3f %.3f %.3f] m\n', ac.mass.get_cg());
fprintf('Empty mass              : %.1f kg\n', ac.mass.get_empty_mass());
fprintf('Fuel mass               : %.1f kg\n', ac.mass.get_fuel_mass());
fprintf('Total mass              : %.1f kg\n', ac.mass.get_total_mass());

fprintf('Nose       : [%.2f %.2f %.2f] m from ref\n', cg_to_nose);
fprintf('Wing apex  : [%.2f %.2f %.2f] m from ref\n', cg_to_wing_apex);
fprintf('Wing AC    : [%.2f %.2f %.2f] m from ref\n', cg_to_wing_ac);
fprintf('Fuel tank 1: [%.2f %.2f %.2f] m from ref\n', cg_to_tank_1);
fprintf('Fuel tank 2: [%.2f %.2f %.2f] m from ref\n', cg_to_tank_2);
fprintf('Engine     : [%.2f %.2f %.2f] m from ref\n', cg_to_engine);

%% Geometry
ac.geometry.wing_area = 27.87;
ac.geometry.wing_span = 9.14;
ac.geometry.mean_aerodynamic_chord = 3.45;
ac.geometry.ref_point = cg_to_wing_ac.';

%% Aerodynamics
ac.set_aerodynamics(CoefficientAerodynamics(@F16Lookup));

%% Control Surfaces
cfg = ac.get_configurator();

cfg.add_control_surface('name','aileron','surface_type','aileron','classification','primary','axis',[1 0 0], ...
    'max_deflection',deg2rad(20),'min_deflection',deg2rad(-20),'dCl',0.05,'dCm',0,'dCn',0);

cfg.add_control_surface('name','elevator','surface_type','elevator','classification','primary','axis',[0 1 0], ...
    'max_deflection',deg2rad(25),'min_deflection',deg2rad(-25),'dCl',0,'dCm',-0.10,'dCn',0);

cfg.add_control_surface('name','rudder','surface_type','rudder','classification','primary','axis',[0 0 1], ...
    'max_deflection',deg2rad(30),'min_deflection',deg2rad(-30),'dCl',0,'dCm',0,'dCn',-0.08);

%% Propulsion
cfg.add_propulsive_element('name','F100_PW_229','element_type','turbofan_afterburning','max_output',128992, ...
    'position',cg_to_engine.','direction',[1 0 0],'fuel_rate',2.5, ...
    'thrust_model',@(thr,M,alt,V,rho) F100ThrustModel(thr, M, alt, V, rho, 128992));

%% Control Vector Setup
n_cs = numel(ac.control_surfaces);
n_pe = numel(ac.propulsive_elements);
n_total = n_cs + n_pe;

%% Cruise Condition
alt_cruise  = 9144;
mach_cruise = 0.6;
dur_cruise  = 60;

[~, a, ~, ~] = atmosisa(alt_cruise);
V_cruise = mach_cruise * a;


alpha_test = deg2rad(6);

x = zeros(12,1);
x(3) = -alt_cruise;
x(4) = V_cruise * cos(alpha_test);
x(6) = V_cruise * sin(alpha_test);
x(8) = alpha_test;

ac.state.set_full_state(x);
ac.control_surfaces(2).set_deflection(0);
ac.propulsive_elements{1}.set_throttle(0.35);
ac.sync_control_vector_from_components();

[F_ref, M_ref, ~] = ac.calculate_total_forces_moments_with_gravity();

m = ac.mass.get_total_mass();
g = 9.80665;
W = m * g;
cbar = ac.geometry.mean_aerodynamic_chord;

res_check = [F_ref(1)/W; F_ref(3)/W; M_ref(2)/(W*cbar)];

fprintf('Test at alpha=%.1f deg, elevator=0 deg:\n', rad2deg(alpha_test));
fprintf('  F(ref) = [%.1f %.1f %.1f] N\n', F_ref);
fprintf('  M(ref) = [%.1f %.1f %.1f] N-m\n', M_ref);
fprintf('  My/(Wc) = %.4f\n', res_check(3));



%% Trim Solution
solver = ac.get_trim_solver();
solver.trim_tolerance = 1e-6;
solver.max_iterations = 15000;
solver.use_fmincon = true;

alpha0 = deg2rad(6.0);
elev0  = deg2rad(-1.0);
thr0   = 0.35;

solver.initial_guess = [alpha0; elev0; thr0];


[x_trim, u_trim, converged] = solver.solve_cruise_trim(alt_cruise, mach_cruise);

if converged
    fprintf('\n=== TRIM CONVERGED ===\n');
    solver.print_summary();
else
    fprintf('\n=== TRIM FAILED - USING FALLBACK ===\n');

    alpha_fb = solver.initial_guess(1);

    x_trim = zeros(12,1);
    x_trim(3) = -alt_cruise;
    x_trim(4) = V_cruise * cos(alpha_fb);
    x_trim(6) = V_cruise * sin(alpha_fb);
    x_trim(8) = alpha_fb;

    u_trim = zeros(n_total,1);
    u_trim(2) = elev0;
    u_trim(n_cs+1) = 0.5;
end
%% Simulation Setup
time_vector = linspace(0, dur_cruise, 1000);
control_input_data = [time_vector(:), repmat(u_trim.', length(time_vector), 1)];

Initialpos = x_trim(1:3).';
InitialVel = x_trim(4:6).';
InitialOri = x_trim(7:9).';
InitialRot = x_trim(10:12).';

assignin('base','ac', ac);
assignin('base','initial_state', x_trim);
assignin('base','Initialpos', Initialpos);
assignin('base','InitialVel', InitialVel);
assignin('base','InitialOri', InitialOri);
assignin('base','InitialRot', InitialRot);
assignin('base','n_cs', n_cs);
assignin('base','n_pe', n_pe);
assignin('base','n_total', n_total);
assignin('base','control_input_data', control_input_data);
assignin('base','sim_stop_time', dur_cruise);

%% Ground Contact
ac.ground_k = 3e5;
ac.ground_c = 5e4;

%% Verification
ac.state.set_full_state(x_trim);
ac.set_controls_from_vector(u_trim);

[F_check, M_check, ~] = ac.calculate_total_forces_moments_with_gravity();
[F_cg, M_cg, ~] = ac.get_forces_moments_about_cg();


fprintf('Reference pt  : [%.3f %.3f %.3f] m\n', ac.get_reference_point());
fprintf('Total mass    : %.1f kg\n', ac.mass.get_total_mass());
fprintf('Empty mass    : %.1f kg\n', ac.mass.get_empty_mass());
fprintf('Fuel mass     : %.1f kg\n', ac.mass.get_fuel_mass());
fprintf('Current CG    : [%.3f %.3f %.3f] m\n', ac.mass.get_cg());
fprintf('Altitude      : %.1f m\n', ac.state.get_altitude());
fprintf('Airspeed      : %.1f m/s, Mach %.2f\n', ac.state.get_airspeed(), mach_cruise);
fprintf('Alpha / Beta  : %.2f / %.2f deg\n', rad2deg(ac.state.get_alpha()), rad2deg(ac.state.get_beta()));
fprintf('Pitch / Roll  : %.2f / %.2f deg\n', rad2deg(x_trim(8)), rad2deg(x_trim(7)));

fprintf('Aileron       : %.3f deg\n', rad2deg(u_trim(1)));
fprintf('Elevator      : %.3f deg\n', rad2deg(u_trim(2)));
fprintf('Rudder        : %.3f deg\n', rad2deg(u_trim(3)));
fprintf('Throttle      : %.3f\n', u_trim(4));

fprintf('Inertia tensor [kg-m^2]:\n');

I = ac.mass.get_inertia_matrix();

fprintf('  [%10.1f %10.1f %10.1f]\n', I(1,1), I(1,2), I(1,3));
fprintf('  [%10.1f %10.1f %10.1f]\n', I(2,1), I(2,2), I(2,3));
fprintf('  [%10.1f %10.1f %10.1f]\n', I(3,1), I(3,2), I(3,3));

fuel_levels = [2400, 1800, 1200, 600, 0];

for f = fuel_levels

    ac.mass.set_fuel_mass(f);

    cg_current = ac.mass.get_cg();
    cg_shift   = cg_current(:) - ac.get_reference_point();

    fprintf(['Fuel %4.0f kg -> ' ...
             'CG = [%.3f %.3f %.3f] m, ' ...
             'CG-ref shift = [%.2f %.2f %.2f] cm, ' ...
             '|shift| = %.2f cm\n'], ...
        f, ...
        cg_current(1), cg_current(2), cg_current(3), ...
        cg_shift(1)*100, cg_shift(2)*100, cg_shift(3)*100, ...
        norm(cg_shift)*100);

end

ac.mass.set_fuel_mass(2400);
