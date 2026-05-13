clear; clc

ac = Aircraft();
ac.set_reference_point([0 0 0]);

%% Reference Definition
NOSE_REF_CG = [1.95;0.00;1.26];
airframe_nose_ref = [1.95;0.00;1.26];
payload_nose_ref = [2.10;0.00;1.20];
wing_ac_nose_ref = [1.95;0.00;1.26];
engine_nose_ref = [1.68;0.00;1.26];
cg_to_airframe = airframe_nose_ref - NOSE_REF_CG;
cg_to_payload = payload_nose_ref - NOSE_REF_CG;
cg_to_wing_ac = wing_ac_nose_ref - NOSE_REF_CG;
cg_to_engine = engine_nose_ref - NOSE_REF_CG;

%% Mass Properties
target_gross_kg = 2550*0.45359237;
empty_mass_kg = 767;
fuel_mass_kg = 144;
payload_kg = target_gross_kg - empty_mass_kg - fuel_mass_kg;

ac.set_mass_model(ComponentMass());
ac.mass.add_component('empty_airframe',empty_mass_kg,cg_to_airframe.',[1285 0 0;0 1824 0;0 0 2666],'airframe');
ac.mass.add_component('payload',payload_kg,cg_to_payload.',zeros(3,3),'payload');
ac.mass.add_component('fuel_fixed',fuel_mass_kg,[0 0 0],zeros(3,3),'fuel');
ac.mass.set_fuel_mass(0);
ac.mass.set_fuel_capacity(0,target_gross_kg);

%% Geometry
ac.geometry.wing_area = 16.17;
ac.geometry.wing_span = 10.98;
ac.geometry.mean_aerodynamic_chord = 1.50;
ac.geometry.ref_area = 16.17;
ac.geometry.ref_span = 10.98;
ac.geometry.ref_chord = 1.50;
ac.geometry.ref_point = cg_to_wing_ac.';

%% Aerodynamics
datcom_path = "C:\Users\naman\Downloads\Aircrfatdummy\dummy7\AircraftDesign_Flightsim\Matlab\Example\cessna\cessna.out";
c172_lookup = c172datcom(datcom_path);
ac.set_aerodynamics(CoefficientAerodynamics(c172_lookup));

%% Control Surfaces
cfg = ac.get_configurator();
cfg.add_control_surface('name','aileron','surface_type','aileron','classification','primary','axis',[1 0 0],'max_deflection',deg2rad(20),'min_deflection',deg2rad(-20),'dCl',0,'dCm',0,'dCn',0);
cfg.add_control_surface('name','elevator','surface_type','elevator','classification','primary','axis',[0 1 0],'max_deflection',deg2rad(28),'min_deflection',deg2rad(-26),'dCl',0,'dCm',0,'dCn',0);
cfg.add_control_surface('name','rudder','surface_type','rudder','classification','primary','axis',[0 0 1],'max_deflection',deg2rad(30),'min_deflection',deg2rad(-30),'dCl',0,'dCm',0,'dCn',0);

%% Propulsion
cfg.add_propulsive_element('name','O360','element_type','propeller','max_output',180*745.7,'position',cg_to_engine.','mount_euler',[0 0 0],'direction',[1 0 0],'fuel_rate',0.0102);
ac.propulsive_elements{1}.set_propeller_params(1.905,1.219,2,0.80,[0.10,-0.05],[0.05,0.02],0.05,0.1);
ac.propulsive_elements{1}.set_efficiency_params(0.80,0.58);
ac.propulsive_elements{1}.set_sfc_curve([0.00,0.40,0.55,0.65,0.75,1.00],[0.62,0.58,0.53,0.50,0.50,0.54],true);

%% Validation Cases
eta_prop_eff = 0.52;
test_cases = struct('name',{'2000ft cruise','6000ft cruise','10000ft cruise','6000ft economy'},'alt_ft',{2000,6000,10000,6000},'ktas',{113,112,110,103},'gph',{9.2,8.8,8.1,7.9},'bhp_pct',{68,65,59,50},'rpm',{2500,2500,2500,2300});
n_cases = numel(test_cases);
results = struct([]);

fprintf('=== AIR PLAINS C172 (180HP) CRUISE VALIDATION ===\n\n');
fprintf('%-20s %-8s %-8s %-9s %-8s %-8s %-10s %-8s %-8s %-10s\n','Condition','RefKTAS','SimKTAS','SpdErr%','RefGPH','SimGPH','FuelErr%','RefBHP','SimBHP','BHP Err%');
fprintf('%s\n',repmat('-',107,1));

for i = 1:n_cases
tc = test_cases(i);
alt_m = tc.alt_ft*0.3048;
V_ref = tc.ktas*0.514444;
[~,a,~,~] = atmosisa(alt_m);
mach = V_ref/a;
solver = ac.get_trim_solver();
solver.trim_tolerance = 1e-6;
solver.use_fmincon = true;
solver.initial_guess = [deg2rad(3);0;0.55];
[x_trim,u_trim,converged,info] = solver.solve_cruise_trim(alt_m,mach);
if ~converged
fprintf('%-20s %-8.1f %-8s %-9s %-8.1f %-8s %-10s %-8.0f %-8s %-10s\n',tc.name,tc.ktas,'FAIL','N/A',tc.gph,'N/A','N/A',tc.bhp_pct,'N/A','N/A');
continue;
end
ac.state.set_full_state(x_trim);
ac.set_controls_from_vector(u_trim);
[F_aero,~,~] = ac.aero.calculate_forces_moments(x_trim,u_trim,ac.geometry,ac,ac.time_step);
V_sim = norm(x_trim(4:6));
ktas_sim = V_sim/0.514444;
speed_err = 100*(ktas_sim - tc.ktas)/tc.ktas;
alpha = atan2(x_trim(6),x_trim(4));
Vhat = [cos(alpha);0;sin(alpha)];
D = -dot(F_aero,Vhat);
P_thrust_W = D*V_sim;
P_shaft_W = P_thrust_W/eta_prop_eff;
sim_bhp = P_shaft_W/745.7;
sim_bhp_pct = 100*sim_bhp/180;
bhp_err = sim_bhp_pct - tc.bhp_pct;

results(i).name = tc.name;
results(i).alt_ft = tc.alt_ft;
results(i).ref_ktas = tc.ktas;
results(i).sim_ktas = ktas_sim;
results(i).speed_err = speed_err;
results(i).ref_bhp_pct = tc.bhp_pct;
results(i).sim_bhp_pct = sim_bhp_pct;
results(i).bhp_err = bhp_err;
results(i).alpha_deg = rad2deg(alpha);
results(i).theta_deg = rad2deg(x_trim(8));
results(i).elev_deg = rad2deg(u_trim(2));
results(i).thr = u_trim(4);
results(i).drag_N = D;
results(i).P_thrust_hp = P_thrust_W/745.7;
results(i).P_shaft_hp = sim_bhp;
results(i).gross_lb = ac.mass.get_total_mass()/0.45359237;
results(i).converged = converged;

fprintf('%-20s %-8.1f %-8.1f %-9.2f %-8.0f %-8.1f %-10.2f\n',tc.name,tc.ktas,ktas_sim,speed_err,tc.bhp_pct,sim_bhp_pct,bhp_err);
end

idx = find(strcmp({results.name},'6000ft cruise'),1);

if ~isempty(idx)
r = results(idx);

fprintf('Condition:\n');
fprintf('  Altitude     : %.0f ft\n',r.alt_ft);
fprintf('  Target speed : %.0f KTAS\n',r.ref_ktas);
fprintf('  Weight       : %.0f lbs\n',r.gross_lb);
fprintf('\nTrim Solution:\n');
fprintf('  Speed        : %.1f KTAS (%.1f%% error)\n',r.sim_ktas,r.speed_err);
fprintf('  Alpha        : %.2f deg\n',r.alpha_deg);
fprintf('  Pitch        : %.2f deg\n',r.theta_deg);
fprintf('  Elevator     : %.2f deg\n',r.elev_deg);
fprintf('  Throttle     : %.3f\n',r.thr);
fprintf('\nPerformance:\n');

fprintf('  Ref %%BHP     : %.0f\n',r.ref_bhp_pct);
fprintf('  Sim %%BHP     : %.1f\n',r.sim_bhp_pct);
fprintf('  %%BHP Error   : %.1f\n',r.bhp_err);
fprintf('  Thrust Power : %.1f HP\n',r.P_thrust_hp);
fprintf('  Shaft Power  : %.1f HP\n',r.P_shaft_hp);
end

