clc; clear; clear functions

ac = Aircraft();
ac.set_reference_point([0 0 0]);

%% Reference Definition
NOSE_REF_CG = [12.6;0;0];
wing_ac_nose_ref = [12.6;0;0];
engine_L_nose_ref = [12.9;5.5;-1.5];
engine_R_nose_ref = [12.9;-5.5;-1.5];
cg_to_wing_ac = wing_ac_nose_ref - NOSE_REF_CG;
cg_to_engine_L = engine_L_nose_ref - NOSE_REF_CG;
cg_to_engine_R = engine_R_nose_ref - NOSE_REF_CG;

%% Mass Properties
ac.mass.set_mass_properties(65000,3.5e6,4.5e6,7.0e6,0,0,0,[0 0 0]);

%% Geometry
ac.geometry.wing_area = 1329.9*0.092903;
ac.geometry.wing_span = 93.0*0.3048;
ac.geometry.mean_aerodynamic_chord = 14.3*0.3048;
ac.geometry.ref_area = ac.geometry.wing_area;
ac.geometry.ref_span = ac.geometry.wing_span;
ac.geometry.ref_chord = ac.geometry.mean_aerodynamic_chord;
ac.geometry.ref_point = cg_to_wing_ac.';

%% Aerodynamics
lookup = b737datcom('datcom.out');
ac.load_aerodynamics('coefficient',lookup);

%% Control Surfaces
cfg = ac.get_configurator();
cfg.add_control_surface('name','aileron','surface_type','aileron','classification','primary','axis',[1 0 0],'max_deflection',deg2rad(20),'min_deflection',deg2rad(-20),'dCl',0.08,'dCm',0.0,'dCn',0.01);
cfg.add_control_surface('name','elevator','surface_type','elevator','classification','primary','axis',[0 1 0],'max_deflection',deg2rad(25),'min_deflection',deg2rad(-25),'dCl',0.0,'dCm',-1.10,'dCn',0.0);
cfg.add_control_surface('name','rudder','surface_type','rudder','classification','primary','axis',[0 0 1],'max_deflection',deg2rad(30),'min_deflection',deg2rad(-30),'dCl',0.0,'dCm',0.0,'dCn',-0.07);

%% Propulsion
cfg.add_propulsive_element('name','CFM56_L','element_type','turbofan','max_output',104500,'position',cg_to_engine_L.','mount_euler',[0 0 0],'direction',[1 0 0],'fuel_rate',9.58);
cfg.add_propulsive_element('name','CFM56_R','element_type','turbofan','max_output',104500,'position',cg_to_engine_R.','mount_euler',[0 0 0],'direction',[1 0 0],'fuel_rate',9.58);

for k = 1:2
ac.propulsive_elements{k}.set_jet_params(8500,[0.6 1.0],[1.00 -0.10;0.94 -0.08;0.91 -0.04],0.30,1.0,0.30);
end

%% Trim
alt = 7000;
mach = 0.60;
ts = ac.get_trim_solver();
ts.use_fmincon = true;
ts.trim_tolerance = 1e-6;
ts.max_iterations = 5000;

[x_trim,u_trim,converged,info] = ts.solve_cruise_trim(alt,mach);

[~,a_spd,~,rho] = atmosisa(alt);
V = mach*a_spd;
q = 0.5*rho*V^2;
W = ac.mass.get_total_mass()*9.80665;
C = lookup(x_trim,u_trim(1:3),ac.geometry);

fprintf('\n=== TRIM  M%.2f / %dm ===\n',mach,alt);
fprintf('Converged  : %d\n',converged);
fprintf('Alpha      : %.3f deg\n',rad2deg(atan2(x_trim(6),max(x_trim(4),1e-9))));
fprintf('Elevator   : %.3f deg\n',rad2deg(info.delta_pitch));
fprintf('Throttle   : %.4f\n',info.throttle_trim);
fprintf('CL / CD    : %.4f / %.4f\n',C.CL,C.CD);
fprintf('L/D        : %.1f\n',C.CL/max(C.CD,1e-6));
fprintf('q_bar      : %.0f Pa\n',q);

%% Takeoff
ac.mass.set_mass_properties(79015,3.5e6,4.5e6,7.0e6,0,0,0,[0 0 0]);
to = ac.get_takeoff();
[~,res_to] = to.calculate_takeoff(0,0,3200);
to.print_takeoff_summary(res_to);

%% Landing
ac.mass.set_mass_properties(65000,3.5e6,4.5e6,7.0e6,0,0,0,[0 0 0]);
ld = ac.get_landing();
[landing_distance,res_ld] = ld.calculate_landing(0,0,3000);

fprintf('\n=== LANDING SUMMARY ===\n');
fprintf('Distance         : %.1f m\n',landing_distance);
fprintf('Vstall landing   : %.2f m/s\n',res_ld.V_stall);
fprintf('V approach       : %.2f m/s\n',res_ld.V_approach);
fprintf('V touchdown      : %.2f m/s\n',res_ld.V_touchdown);
fprintf('Runway adequate  : %d\n',res_ld.runway_adequate);

%% Stability
stab = ac.get_stability();
stab.set_trim(x_trim,u_trim);

[~,~,Cm_alpha] = stab.compute_pitch_static_stability(deg2rad(-8:0.5:12));
[~,~,Cn_beta] = stab.compute_yaw_static_stability(deg2rad(-6:0.5:6));

stab.linearize(1e-6,1e-6);
stab.analyze_modes();
s = stab.get_modes_summary();

fprintf('\n=== STABILITY ===\n');
fprintf('Cm_alpha : %.4f /rad  (%s)\n',Cm_alpha,tf(Cm_alpha<0,'stable','UNSTABLE'));
fprintf('Cn_beta  : %.4f /rad  (%s)\n',Cn_beta,tf(Cn_beta>0,'stable','UNSTABLE'));

modes = {'short_period','phugoid','dutch_roll'};

for i = 1:3
m = s.(modes{i});
if ~isempty(m)
fprintf('%-14s : wn=%.3f  zeta=%.4f  T=%.2f s  %s\n',modes{i},abs(m.lambda),m.damping,m.period,tf(real(m.lambda)<0,'stable','UNSTABLE'));
end
end

if ~isempty(s.roll_subsidence)
fprintf('roll_subsid    : tau=%.3f s  %s\n',s.roll_subsidence.time_constant,tf(real(s.roll_subsidence.lambda)<0,'stable','UNSTABLE'));
end

if ~isempty(s.spiral)
fprintf('spiral         : tau=%.1f s  %s\n',s.spiral.time_constant,tf(s.spiral.stable,'stable','unstable'));
end

function s = tf(c,a,b)
if c
s = a;
else
s = b;
end
end