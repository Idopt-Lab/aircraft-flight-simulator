%% F16_CRUISE_COMPONENTMASS_REFPOINT.m
% F-16 cruise trim using ComponentMass with explicit aircraft reference point

clear; clc

if bdIsLoaded('flightsim7')
    set_param('flightsim7', 'SimulationCommand', 'stop');
end

%% Aircraft Initialization
ac = Aircraft();
ac.set_mass_model(ComponentMass());

% Nose-referenced coordinate system:
%   x forward/aft along fuselage, y right, z down
% Here the aircraft reference point is the nose.
ac.set_reference_point([0 0 0]);

%% Reference Frame Definition
NOSE_REF_POINT = [0.00; 0.00; 0.00];

fuselage_fwd_nose_ref  = [2.00;  0.00;  0.00];
fuselage_mid_nose_ref  = [4.50;  0.00;  0.00];
fuselage_aft_nose_ref  = [7.00;  0.00;  0.10];

wing_center_nose_ref   = [4.20;  0.00;  0.15];
wing_left_nose_ref     = [4.20; -2.50;  0.20];
wing_right_nose_ref    = [4.20;  2.50;  0.20];

h_tail_nose_ref        = [8.50;  0.00;  0.30];
v_tail_nose_ref        = [8.50;  0.00; -0.50];

nose_gear_nose_ref     = [1.50;  0.00;  0.50];
main_gear_left_nose_ref  = [4.80; -1.20; 0.80];
main_gear_right_nose_ref = [4.80;  1.20; 0.80];

engine_core_nose_ref   = [6.00;  0.00;  0.00];
inlet_nose_ref         = [3.50;  0.00;  0.00];
nozzle_nose_ref        = [8.00;  0.00;  0.00];

avionics_nose_ref      = [1.80;  0.00; -0.30];
hydraulics_nose_ref    = [5.50;  0.00;  0.20];
fuel_system_nose_ref   = [4.50;  0.00;  0.30];
flight_controls_nose_ref = [4.00; 0.00; -0.20];

cockpit_nose_ref       = [1.20;  0.00; -0.50];
ejection_seat_nose_ref = [1.50;  0.00; -0.30];
pilot_nose_ref         = [1.50;  0.00; -0.40];

radar_nose_ref         = [0.50;  0.00; -0.20];
targeting_pod_nose_ref = [3.50;  0.00; -0.80];

fuel_tank_1_nose_ref   = [4.00;  0.00;  0.00];
fuel_tank_2_nose_ref   = [5.50;  0.00;  0.00];

wing_ac_nose_ref       = [4.20;  0.00;  0.00];
engine_nose_ref        = [6.00;  0.00;  0.00];

%% Mass Properties
fprintf('=== COMPONENTMASS CONFIGURATION ===\n');

ac.mass.add_component('fuselage_fwd', 1200, fuselage_fwd_nose_ref.', diag([80 300 300]), 'structure');
ac.mass.add_component('fuselage_mid', 1500, fuselage_mid_nose_ref.', diag([100 400 400]), 'structure');
ac.mass.add_component('fuselage_aft', 800, fuselage_aft_nose_ref.', diag([60 250 250]), 'structure');

ac.mass.add_component('wing_center', 1200, wing_center_nose_ref.', diag([150 30 150]), 'structure');
ac.mass.add_component('wing_left', 600, wing_left_nose_ref.', diag([80 15 80]), 'structure');
ac.mass.add_component('wing_right', 600, wing_right_nose_ref.', diag([80 15 80]), 'structure');

ac.mass.add_component('h_tail', 250, h_tail_nose_ref.', diag([20 5 20]), 'structure');
ac.mass.add_component('v_tail', 200, v_tail_nose_ref.', diag([15 15 5]), 'structure');

ac.mass.add_component('nose_gear', 150, nose_gear_nose_ref.', [], 'landing_gear');
ac.mass.add_component('main_gear_left', 200, main_gear_left_nose_ref.', [], 'landing_gear');
ac.mass.add_component('main_gear_right', 200, main_gear_right_nose_ref.', [], 'landing_gear');

ac.mass.add_component('engine_core', 1800, engine_core_nose_ref.', diag([200 100 200]), 'propulsion');
ac.mass.add_component('inlet', 150, inlet_nose_ref.', [], 'propulsion');
ac.mass.add_component('nozzle', 100, nozzle_nose_ref.', [], 'propulsion');

ac.mass.add_component('avionics', 300, avionics_nose_ref.', [], 'avionics');
ac.mass.add_component('hydraulics', 150, hydraulics_nose_ref.', [], 'systems');
ac.mass.add_component('fuel_system', 120, fuel_system_nose_ref.', [], 'systems');
ac.mass.add_component('flight_controls', 180, flight_controls_nose_ref.', [], 'systems');

ac.mass.add_component('cockpit', 250, cockpit_nose_ref.', [], 'cockpit');
ac.mass.add_component('ejection_seat', 100, ejection_seat_nose_ref.', [], 'cockpit');
ac.mass.add_component('pilot', 100, pilot_nose_ref.', [], 'crew');

ac.mass.add_component('radar', 150, radar_nose_ref.', [], 'avionics');
ac.mass.add_component('targeting_pod', 80, targeting_pod_nose_ref.', [], 'avionics');

ac.mass.add_fuel_tank(fuel_tank_1_nose_ref.', 1800, 'distributed');
ac.mass.add_fuel_tank(fuel_tank_2_nose_ref.', 1200, 'distributed');
ac.mass.set_fuel_mass(2400);

ac.mass.list_components();

fprintf('\nMass summary:\n');
fprintf('Aircraft reference point: [%.3f %.3f %.3f] m\n', ac.get_reference_point());
fprintf('Empty mass              : %.1f kg\n', ac.mass.get_empty_mass());
fprintf('Fuel mass               : %.1f kg\n', ac.mass.get_fuel_mass());
fprintf('Total mass              : %.1f kg\n', ac.mass.get_total_mass());
fprintf('Current CG              : [%.3f %.3f %.3f] m\n', ac.mass.get_cg());

%% Geometry
ac.geometry.wing_area = 27.87;
ac.geometry.wing_span = 9.14;
ac.geometry.mean_aerodynamic_chord = 3.45;
ac.geometry.ref_point = wing_ac_nose_ref.';

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
    'position',engine_nose_ref.','direction',[1 0 0],'fuel_rate',2.5, ...
    'thrust_model',@(thr,M,alt,V,rho) F100ThrustModel(thr, M, alt, V, rho, 128992));

%% Control Vector Setup
n_cs = numel(ac.control_surfaces);
n_pe = numel(ac.propulsive_elements);
n_total = n_cs + n_pe;

%% Cruise Condition
alt_cruise = 9144;
mach_cruise = 0.6;
dur_cruise = 60;

[~, a, ~, ~] = atmosisa(alt_cruise);
V_cruise = mach_cruise * a;

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