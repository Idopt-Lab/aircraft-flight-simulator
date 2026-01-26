clear; clc
if bdIsLoaded('flightsim7')
    set_param('flightsim7', 'SimulationCommand', 'stop');
end
ac = Aircraft();
empty_mass = 9300;
Ixx = 55800; 
Iyy = 63100; 
Izz = 118000;
ac.mass.set_mass_properties(empty_mass, Ixx, Iyy, Izz, 0, 0, 0, [0 0 0]);
ac.mass.add_fuel_tank([0 0 0], 3000, 'distributed');
ac.mass.add_fuel_tank([-2 0 0], 1500, 'distributed');
ac.mass.set_fuel_mass(2400);
ac.geometry.wing_area = 27.87;
ac.geometry.wing_span = 9.14;
ac.geometry.mean_aerodynamic_chord = 3.45;
ac.aero.set_lookup(@F16Lookup);
cfg = ac.get_configurator();
cfg.add_control_surface('name','aileron','surface_type','aileron','classification','primary','axis',[1 0 0], ...
    'max_deflection',deg2rad(20),'min_deflection',deg2rad(-20),'dCl',0.05,'dCm',0,'dCn',0);
cfg.add_control_surface('name','elevator','surface_type','elevator','classification','primary','axis',[0 1 0], ...
    'max_deflection',deg2rad(25),'min_deflection',deg2rad(-25),'dCl',0,'dCm',-0.10,'dCn',0);
cfg.add_control_surface('name','rudder','surface_type','rudder','classification','primary','axis',[0 0 1], ...
    'max_deflection',deg2rad(30),'min_deflection',deg2rad(-30),'dCl',0,'dCm',0,'dCn',-0.08);
cfg.add_propulsive_element('name','F100_PW_229','element_type','turbofan_afterburning','max_output',128992, ...
    'position',[-1.5 0 0],'direction',[1 0 0],'fuel_rate',2.5,'thrust_model',@(thr,M,alt,V,rho) F100ThrustModel(thr, M, alt, V, rho, 128992));
n_cs = numel(ac.control_surfaces);
n_pe = numel(ac.propulsive_elements);
n_total = n_cs + n_pe;
alt_cruise = 9144;
mach_cruise = 0.60;
dur_cruise = 60;
[~, a, ~, ~] = atmosisa(alt_cruise);
V_cruise = mach_cruise * a;
solver = ac.get_trim_solver();
solver.trim_tolerance = 1e-3;
solver.max_iterations = 15000;
solver.use_fmincon = true;
solver.initial_guess = [deg2rad(3.0); deg2rad(0.0); 0.50];
[x_trim, u_trim, converged] = solver.solve_cruise_trim(alt_cruise, mach_cruise);
if converged
    fprintf('\nTrim converged!\n');
    solver.print_summary();
else
    fprintf('\nTrim failed, using fallback\n');
    x_trim = zeros(12,1);
    x_trim(3) = -alt_cruise;
    x_trim(4) = V_cruise;
    u_trim = zeros(n_total,1);
    u_trim(n_cs+1) = 0.5;
end
time_vector = linspace(0, dur_cruise, 1000);
control_input_data = [time_vector(:), repmat(u_trim', length(time_vector), 1)];

Initialpos = x_trim(1:3)';
InitialVel = x_trim(4:6)';
InitialOri = x_trim(7:9)';
InitialRot = x_trim(10:12)';

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
assignin('base','ground_k', 3e5);
assignin('base','ground_c', 5e4);
ac.state.set_full_state(x_trim);
ac.set_controls_from_vector(u_trim);
[F_check, M_check, ~] = ac.calculate_total_forces_moments_with_gravity();
