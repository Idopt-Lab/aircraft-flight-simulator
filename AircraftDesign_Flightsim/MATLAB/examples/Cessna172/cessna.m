if bdIsLoaded('flightsim7')
set_param('flightsim7','SimulationCommand','stop');
end
clear; clc

ac = Aircraft();
ac.set_reference_point([0 0 0]);

%% Mass Properties
empty_mass = 767;
Ixx = 1285;
Iyy = 1824;
Izz = 2666;
ac.mass.set_mass_properties(empty_mass,Ixx,Iyy,Izz,0,0,0,[0 0 0]);
ac.mass.set_fuel_mass(90);

%% Reference Definition
CG_FROM_NOSE = 2.00;
cg_to_nose = -CG_FROM_NOSE;
cg_to_wing_le = 1.20 - CG_FROM_NOSE;
cg_to_wing_ac = 1.57 - CG_FROM_NOSE;
cg_to_prop = 1.68 - CG_FROM_NOSE;

%% Geometry
ac.geometry.wing_area = 16.17;
ac.geometry.wing_span = 10.98;
ac.geometry.mean_aerodynamic_chord = 1.50;
ac.geometry.ref_area = 16.17;
ac.geometry.ref_span = 10.98;
ac.geometry.ref_chord = 1.50;
ac.geometry.ref_point = [cg_to_wing_ac,0,0];

%% Aerodynamics
datcom_fp = "C:\Users\naman\Downloads\Aircrfatdummy\dummy3\aircraft-flight-simulator\AircraftDesign_Flightsim\Matlab\Example\cessna\cessna.out";
if isfile(datcom_fp)
ac.set_aerodynamics(CoefficientAerodynamics(c172datcom(datcom_fp)));
else
warning('DATCOM file not found');
end

%% Control Surfaces
cfg = ac.get_configurator();
cfg.add_control_surface('name','aileron','surface_type','aileron','classification','primary','axis',[1 0 0],'max_deflection',deg2rad(20),'min_deflection',deg2rad(-20),'dCl',0,'dCm',0,'dCn',0);
cfg.add_control_surface('name','elevator','surface_type','elevator','classification','primary','axis',[0 1 0],'max_deflection',deg2rad(28),'min_deflection',deg2rad(-26),'dCl',0,'dCm',0,'dCn',0);
cfg.add_control_surface('name','rudder','surface_type','rudder','classification','primary','axis',[0 0 1],'max_deflection',deg2rad(30),'min_deflection',deg2rad(-30),'dCl',0,'dCm',0,'dCn',0);

%% Propulsion
cfg.add_propulsive_element('name','O360','element_type','propeller','max_output',134228,'position',[cg_to_prop,0,1.26],'mount_euler',[0 0 0],'direction',[1 0 0],'fuel_rate',0.0089);
ac.propulsive_elements{1}.set_propeller_params(1.905,1.219,2,0.80,[0.10,-0.05],[0.05,0.02],0.05,0.1);

%% Control Vector
n_cs = numel(ac.control_surfaces);
n_pe = numel(ac.propulsive_elements);
n_total = n_cs + n_pe;

fprintf('Reference frame : CG\n');
fprintf('Empty mass      : %.1f kg\n',empty_mass);
fprintf('Fuel mass       : %.1f kg\n',ac.mass.get_fuel_mass());
fprintf('Total mass      : %.1f kg\n',ac.mass.get_total_mass());
fprintf('CG position     : [%.6f %.6f %.6f] m\n',ac.mass.get_cg());
fprintf('Nose            : [%.2f 0 0] m from CG\n',cg_to_nose);
fprintf('Wing LE         : [%.2f 0 0] m from CG\n',cg_to_wing_le);
fprintf('Wing AC         : [%.2f 0 0] m from CG\n',cg_to_wing_ac);
fprintf('Propeller       : [%.2f 0 1.26] m from CG\n',cg_to_prop);

%% Trim Feasibility Check
fprintf('\n=== TRIM FEASIBILITY CHECK ===\n');

x_test = zeros(12,1);
x_test(3) = -1000;
x_test(4) = 50;
x_test(8) = deg2rad(5);

ac.state.set_full_state(x_test);
ac.control_surfaces(2).set_deflection(0);
ac.propulsive_elements{1}.set_throttle(0.65);
ac.sync_control_vector_from_components();

[F_cg,M_cg,~] = ac.get_forces_moments_about_cg();

m = ac.mass.get_total_mass();
W = m*9.80665;
cbar = ac.geometry.mean_aerodynamic_chord;

fprintf('At alpha = %.1f deg, elev = 0 deg, thr = 0.65\n',rad2deg(ac.state.get_alpha()));
fprintf('M(CG)    : [%.1f %.1f %.1f] N-m\n',M_cg);
fprintf('My/(Wc)  : %.4f\n',M_cg(2)/(W*cbar));

%% Trim Solver Test

solver = ac.get_trim_solver();
solver.trim_tolerance = 1e-6;
solver.max_iterations = 5000;
solver.use_fmincon = true;

alt_test = 1000;
mach_test = 0.15;

[x_trim,u_trim,converged] = solver.solve_cruise_trim(alt_test,mach_test);

if converged
solver.print_summary();
ac.state.set_full_state(x_trim);
ac.set_controls_from_vector(u_trim);
[F_trim,M_trim,~] = ac.get_forces_moments_about_cg();
fprintf('Alpha    : %.2f deg\n',rad2deg(ac.state.get_alpha()));
fprintf('Elevator : %.2f deg\n',rad2deg(u_trim(2)));
fprintf('Throttle : %.3f\n',u_trim(n_cs+1));
fprintf('F(CG)    : [%.3f %.3f %.3f] N\n',F_trim);
fprintf('M(CG)    : [%.3f %.3f %.3f] N-m\n',M_trim);
fprintf('My/(Wc)  : %.6f\n',M_trim(2)/(W*cbar));
else
fprintf('\nTRIM FAILED TO CONVERGE\n');
end

%% CG Migration Test
fuel_levels = [180,142,110,85,50,20,0];

for f = fuel_levels
ac.mass.set_fuel_mass(f);
cg_current = ac.mass.get_cg();
cg_current = cg_current(:);
fprintf('Fuel %3.0f kg -> CG = [%+.3f %+.3f %+.3f] m\n',f,cg_current(1),cg_current(2),cg_current(3));
end

ac.mass.set_fuel_mass(90);

