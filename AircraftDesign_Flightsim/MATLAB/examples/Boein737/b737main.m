clc;
clear;
clear functions;
close all;

%% PATH BOOTSTRAP
% Ensures the shared class library and example helper functions are on the
% MATLAB path regardless of the current working directory or session state.
this_script_dir = fileparts(mfilename('fullpath'));
matlab_root_dir = fileparts(fileparts(this_script_dir));
addpath(genpath(fullfile(matlab_root_dir,'classes')));
addpath(genpath(fullfile(matlab_root_dir,'examples')));

%% DATCOM aerodynamic database

datcom_dir = fullfile(fileparts(mfilename('fullpath')),"data");

if ~isfolder(datcom_dir)
    error("B737:DATCOMDirectoryNotFound", "Boeing DATCOM directory not found: %s",datcom_dir);
end

[lookup,b737_datcom_data] = b737datcom(datcom_dir,true);
aero_mach_limit = max(b737_datcom_data.grid.mach_vec);

%% Aircraft object

ac = Aircraft();

%% Reference positions

NOSE_REF_CG = b737_datcom_data.reference.moment_reference_m(:);

aero_ref_station_m = b737_datcom_data.reference.moment_reference_m(1);
engine_station_m = 12.9;
engine_lateral_m = 5.5;
engine_vertical_m = 1.5;

cg_to_aero_ref = [NOSE_REF_CG(1)-aero_ref_station_m;0;0];
cg_to_engine_L = [NOSE_REF_CG(1)-engine_station_m; -engine_lateral_m; engine_vertical_m];
cg_to_engine_R = [NOSE_REF_CG(1)-engine_station_m; engine_lateral_m; engine_vertical_m];

%% Geometry

S_ref = b737_datcom_data.reference.sref_m2;
b_ref = b737_datcom_data.reference.bref_m;
c_ref = b737_datcom_data.reference.cbar_m;

ac.geometry.set_reference_geometry(S_ref,b_ref,c_ref);
ac.geometry.set_reference_point(cg_to_aero_ref.');
ac.geometry.aero_mach_max = aero_mach_limit;

%% Frames

ac.set_body_frame("body");
ac.set_reference_frame("body");

add_or_update_frame(ac,"aero_ref",   "body",cg_to_aero_ref,  @(x) eye(3));
add_or_update_frame(ac,"engine_L",   "body",cg_to_engine_L, @(x) eye(3));
add_or_update_frame(ac,"engine_R",   "body",cg_to_engine_R, @(x) eye(3));
add_or_update_frame(ac,"cg",         "body",[0;0;0],        @(x) eye(3));
add_or_update_frame(ac,"gravity_cg", "body",[0;0;0],        @(x) ned_to_body_dcm(x));

%% Mass properties

cruise_mass_kg = 65000;
takeoff_mass_kg = 79015;
landing_mass_kg = 65000;

I_body = [3.5e6 0     0;
          0     4.5e6 0;
          0     0     7.0e6];

airframe = Component("airframe", cruise_mass_kg, [0;0;0], I_body, ac.get_frame("body"), "airframe");

ac.add_component(airframe);

%% Controls
% Important:
% b737datcom already applies control derivatives inside the lookup.
% So these ControlSurface derivative entries are kept zero to avoid double-counting.

aileron = ControlSurface("aileron","aileron","primary", [1 0 0], deg2rad(20),deg2rad(-20), 0,0,0);

elevator = ControlSurface("elevator","elevator","primary", [0 1 0], deg2rad(25),deg2rad(-25), 0,0,0);

rudder = ControlSurface("rudder","rudder","primary", [0 0 1], deg2rad(30),deg2rad(-30), 0,0,0);

ac.add_control_surface(aileron);
ac.add_control_surface(elevator);
ac.add_control_surface(rudder);

%% Aerodynamics

aero_model = CoefficientAerodynamics(lookup);

aero_solver = AeroLoadSolver(aero_model, ac.geometry, ac, ac.get_frame("aero_ref"));

ac.add_load_source(aero_solver);

%% Propulsion: two CFM56-7B26/3 engines

eng_L = CFM56Turbofan("CFM56_L", ac.get_frame("engine_L"), [1;0;0], "CFM56-7B26/3");

eng_R = CFM56Turbofan("CFM56_R", ac.get_frame("engine_R"), [1;0;0], "CFM56-7B26/3");

% Cruise trim, performance and stability all use the continuous rating.
eng_L.set_rating_mode("continuous");
eng_R.set_rating_mode("continuous");

ac.add_propulsive_element(eng_L);
ac.add_propulsive_element(eng_R);

engine_L_comp = Component("engine_L", 0, [0;0;0], zeros(3,3), ac.get_frame("engine_L"), "engine");

engine_R_comp = Component("engine_R", 0, [0;0;0], zeros(3,3), ac.get_frame("engine_R"), "engine");

engine_L_comp.add_load_source(PropulsionLoadSolver(eng_L,ac.get_frame("engine_L")));
engine_R_comp.add_load_source(PropulsionLoadSolver(eng_R,ac.get_frame("engine_R")));

ac.add_component(engine_L_comp);
ac.add_component(engine_R_comp);

%% Gravity

gravity_solver = GravityLoadSolver(ac,ac.get_frame("gravity_cg"));
ac.add_load_source(gravity_solver);

%% Sync mass / CG

[m_total,cg_total,I_total] = ac.compute_total_mass_properties([]);

ac.update_frame_position("cg",cg_total);
ac.update_frame_position("gravity_cg",cg_total);
ac.set_reference_frame("cg");

W = m_total * ac.g;

n_cs = numel(ac.control_surfaces);
n_pe = numel(ac.propulsive_elements);
n_total = n_cs + n_pe;

fprintf("\n=== B737 READY ===\n");
fprintf("Mass    : %.2f kg\n",m_total);
fprintf("Weight  : %.2f N\n",W);
fprintf("CG body : [% .4f % .4f % .4f] m\n",cg_total);
fprintf("Aero Mach validity: %.2f to %.2f\n", min(b737_datcom_data.grid.mach_vec),aero_mach_limit);
fprintf("Engine variant     : %s\n",eng_L.variant);
fprintf("n_cs    : %d\n",n_cs);
fprintf("n_pe    : %d\n",n_pe);

%% Trim setup

alt_trim_m = 7000;
mach_trim = 0.58;
% Kept strictly inside the usable DATCOM Mach grid (max = 0.60): trimming
% exactly on the grid boundary lets floating-point round-trip in the
% Mach<->velocity conversion push the interpolant query a hair past 0.60,
% which falls outside scatteredInterpolant's convex hull and forces a
% non-smooth nearest-neighbor fallback right where the trim solver needs
% a clean gradient.

[~,a_trim,~,rho_trim] = ac.get_atmosphere(alt_trim_m);
V_trim = mach_trim * a_trim;
qbar_trim = 0.5 * rho_trim * V_trim^2;

condition = struct();
condition.altitude_m = alt_trim_m;
condition.mach = mach_trim;

trim_cfg = struct();
trim_cfg.variables = ["alpha","elevator","throttle"];
trim_cfg.residuals = ["Fx","Fz","My"];

trim_cfg.initial_guess = [deg2rad(1); deg2rad(-3); 0.35];

trim_cfg.lb = [deg2rad(-5); deg2rad(-25); 0.0];
trim_cfg.ub = [deg2rad(15); deg2rad(25); 1.0];

trim_cfg.weights = [1;1;1];

trim_cfg.fixed = struct();
trim_cfg.fixed.beta = 0;
trim_cfg.fixed.phi = 0;
trim_cfg.fixed.psi = 0;
trim_cfg.fixed.gamma = 0;
trim_cfg.fixed.aileron = 0;
trim_cfg.fixed.rudder = 0;
trim_cfg.fixed.p = 0;
trim_cfg.fixed.q = 0;
trim_cfg.fixed.r = 0;

trim_cfg.reference_frame_name = "cg";

solver_cfg = struct();
solver_cfg.residual_tolerance = 1e-5;
solver_cfg.fmincon_options = optimoptions("fmincon", ...
    "Algorithm","sqp", ...
    "Display","iter", ...
    "MaxIterations",5000, ...
    "MaxFunctionEvaluations",50000, ...
    "OptimalityTolerance",1e-10, ...
    "StepTolerance",1e-10, ...
    "ConstraintTolerance",1e-10, ...
    "FunctionTolerance",1e-10, ...
    "FiniteDifferenceStepSize",1e-8);

%% Run trim

fprintf("\n=== RUNNING B737 CRUISE TRIM ===\n");
fprintf("Altitude : %.1f m\n",alt_trim_m);
fprintf("Mach     : %.3f\n",mach_trim);
fprintf("V        : %.3f m/s / %.3f kt\n",V_trim,V_trim*1.94384449);

solver = ac.get_trim_solver();

[x_trim,u_trim,converged,info] = solver.solve_trim(condition,trim_cfg,solver_cfg); %#ok<NASGU>

if converged
    fprintf("\n=== TRIM CONVERGED ===\n");
else
    fprintf("\n=== TRIM NOT CONVERGED: using returned best solution ===\n");
end

try
    solver.print_summary();
catch
end

%% Trim verification

ac.state.set_full_state(x_trim);
ac.set_controls_from_vector(u_trim);

[m_trim,cg_trim,I_trim] = ac.compute_total_mass_properties(x_trim); %#ok<NASGU>
ac.update_frame_position("cg",cg_trim);
ac.update_frame_position("gravity_cg",cg_trim);
ac.set_reference_frame("cg");

[F_trim,M_trim,fuel_flow_trim] = ac.compute_total_loads(x_trim,u_trim);

coeff_trim = lookup(x_trim,u_trim,ac.geometry);

W_trim = m_trim * ac.g;
cbar_trim = c_ref;

alpha_trim_deg = rad2deg(atan2(x_trim(6),max(x_trim(4),1e-9)));

if isfield(coeff_trim,'datcom_valid') && ~coeff_trim.datcom_valid
    warning("B737:BaselineOutsideDATCOMCoverage", "The baseline trim point lies outside valid DATCOM coverage.");
end

fprintf("\n=== B737 TRIM CHECK ABOUT CG ===\n");
fprintf("Converged : %d\n",converged);
fprintf("Mass      : %.2f kg\n",m_trim);
fprintf("Alpha     : %.4f deg\n",alpha_trim_deg);
fprintf("Theta     : %.4f deg\n",rad2deg(x_trim(8)));
fprintf("Aileron   : %.4f deg\n",rad2deg(u_trim(1)));
fprintf("Elevator  : %.4f deg\n",rad2deg(u_trim(2)));
fprintf("Rudder    : %.4f deg\n",rad2deg(u_trim(3)));

if numel(u_trim) >= n_cs+n_pe
    fprintf("Throttle L: %.6f\n",u_trim(n_cs+1));
    fprintf("Throttle R: %.6f\n",u_trim(n_cs+2));
elseif numel(u_trim) >= n_cs+1
    fprintf("Throttle  : %.6f\n",u_trim(n_cs+1));
end

fprintf("CL/CD     : %.6f / %.6f\n",coeff_trim.CL,coeff_trim.CD);
fprintf("L/D       : %.3f\n",coeff_trim.CL/max(coeff_trim.CD,1e-9));
fprintf("q_bar     : %.2f Pa\n",qbar_trim);
fprintf("F         : [% .3e % .3e % .3e] N\n",F_trim);
fprintf("M         : [% .3e % .3e % .3e] N-m\n",M_trim);
fprintf("Fx/W      : %.6e\n",F_trim(1)/W_trim);
fprintf("Fz/W      : %.6e\n",F_trim(3)/W_trim);
fprintf("My/Wc     : %.6e\n",M_trim(2)/(W_trim*cbar_trim));
fprintf("Fuel flow : %.6f kg/s\n",fuel_flow_trim);

%% Load contribution check

fprintf("\n=== LOAD CONTRIBUTIONS ABOUT CURRENT REF: %s ===\n",ac.reference_frame_name);

body_frame = ac.get_body_frame();
ref_frame = ac.get_reference_frame();

[F_aero,M_aero,~,src_aero] = ac.compute_loads_by_solver(x_trim,u_trim,"AeroLoadSolver",body_frame,ref_frame);

[F_prop,M_prop,ff_prop,src_prop] = ac.compute_loads_by_solver(x_trim,u_trim,"PropulsionLoadSolver",body_frame,ref_frame);

[F_grav,M_grav,~,src_grav] = ac.compute_loads_by_solver(x_trim,u_trim,"GravityLoadSolver",body_frame,ref_frame);

fprintf("Aero sources      : %d\n",numel(src_aero));
fprintf("Aero              F=[% .3e % .3e % .3e] M=[% .3e % .3e % .3e]\n",F_aero,M_aero);

fprintf("Propulsion sources: %d\n",numel(src_prop));
fprintf("Propulsion        F=[% .3e % .3e % .3e] M=[% .3e % .3e % .3e] fuel=%.6e\n",F_prop,M_prop,ff_prop);

fprintf("Gravity sources   : %d\n",numel(src_grav));
fprintf("Gravity           F=[% .3e % .3e % .3e] M=[% .3e % .3e % .3e]\n",F_grav,M_grav);

fprintf("Total             F=[% .3e % .3e % .3e] M=[% .3e % .3e % .3e]\n",F_trim,M_trim);

%% Performance analysis

fprintf("\n=== RUNNING B737 PERFORMANCE ANALYSIS ===\n");

perf = ac.get_performance();

perf_cfg = struct('available_throttle',1, 'propulsion_type',"thrust", 'reference_frame_name',"cg");

perf_point = perf.evaluate_trim_point(x_trim,u_trim,condition,perf_cfg);
perf_point.trim_converged = converged;

perf.print_point(perf_point);

%% Stability

fprintf("\n=== RUNNING B737 STABILITY ANALYSIS ===\n");

stab = ac.get_stability();
stab.set_trim(x_trim,u_trim);

[~,~,Cm_alpha] = stab.compute_pitch_static_stability(deg2rad(-8:0.5:12));
[~,~,Cn_beta]  = stab.compute_yaw_static_stability(deg2rad(-6:0.5:6));

stab.linearize(1e-5,1e-5);
stab.analyze_modes();

try
    stab.print_modes();
catch
end

s = stab.get_modes_summary();

fprintf("\n=== STATIC STABILITY ===\n");
fprintf("Cm_alpha : %.4f /rad  (%s)\n",Cm_alpha,tf(Cm_alpha<0,"stable","UNSTABLE"));
fprintf("Cn_beta  : %.4f /rad  (%s)\n",Cn_beta,tf(Cn_beta>0,"stable","UNSTABLE"));

modes = ["short_period","phugoid","dutch_roll"];

for i = 1:numel(modes)
    mode_name = modes(i);

    if isfield(s,mode_name)
        m = s.(mode_name);

        if ~isempty(m)
            fprintf("%-14s : wn=%.3f  zeta=%.4f  T=%.2f s  %s\n", mode_name, abs(m.lambda), m.damping, m.period, tf(real(m.lambda)<0,"stable","UNSTABLE"));
        end
    end
end

if isfield(s,"roll_subsidence") && ~isempty(s.roll_subsidence)
    fprintf("roll_subsid    : tau=%.3f s  %s\n", s.roll_subsidence.time_constant, tf(real(s.roll_subsidence.lambda)<0,"stable","UNSTABLE"));
end

if isfield(s,"spiral") && ~isempty(s.spiral)
    fprintf("spiral         : tau=%.1f s  %s\n", s.spiral.time_constant, tf(s.spiral.stable,"stable","unstable"));
end

%% Takeoff

fprintf("\n=== RUNNING B737 TAKEOFF ANALYSIS ===\n");

eng_L.set_rating_mode("takeoff");
eng_R.set_rating_mode("takeoff");

airframe.set_mass_properties(takeoff_mass_kg,[0;0;0],I_body);

[m_to,cg_to,~] = ac.compute_total_mass_properties([]);
ac.update_frame_position("cg",cg_to);
ac.update_frame_position("gravity_cg",cg_to);

to = ac.get_takeoff();

[TO_m,res_to] = to.calculate_takeoff(0,0,3200); %#ok<NASGU>
to.print_takeoff_summary(res_to);

%% Landing

fprintf("\n=== RUNNING B737 LANDING ANALYSIS ===\n");

eng_L.set_rating_mode("idle");
eng_R.set_rating_mode("idle");

airframe.set_mass_properties(landing_mass_kg,[0;0;0],I_body);

[m_ld,cg_ld,~] = ac.compute_total_mass_properties([]);
ac.update_frame_position("cg",cg_ld);
ac.update_frame_position("gravity_cg",cg_ld);

ld = ac.get_landing();

[landing_distance,res_ld] = ld.calculate_landing(0,0,3000);

fprintf("\n=== LANDING SUMMARY ===\n");
fprintf("Distance         : %.1f m\n",landing_distance);
fprintf("Vstall landing   : %.2f m/s\n",res_ld.V_stall);
fprintf("V approach       : %.2f m/s\n",res_ld.V_approach);
fprintf("V touchdown      : %.2f m/s\n",res_ld.V_touchdown);
fprintf("Runway adequate  : %d\n",res_ld.runway_adequate);

%% Restore cruise configuration

eng_L.set_rating_mode("continuous");
eng_R.set_rating_mode("continuous");

airframe.set_mass_properties(cruise_mass_kg,[0;0;0],I_body);

[m_cr,cg_cr,~] = ac.compute_total_mass_properties([]);
ac.update_frame_position("cg",cg_cr);
ac.update_frame_position("gravity_cg",cg_cr);
ac.set_reference_frame("cg");

ac.state.set_full_state(x_trim);
ac.set_controls_from_vector(u_trim);

fprintf("\n=== B737 SCRIPT COMPLETE ===\n");

%% ============================================================
% Local functions
% ============================================================

function add_or_update_frame(ac,name,parent_name,r_parent,dcm_fn)
    if ac.has_frame(name)
        ac.update_frame_position(name,r_parent(:));
        ac.update_frame_orientation(name,dcm_fn);
    else
        ac.add_frame(name,parent_name,r_parent(:),dcm_fn);
    end
end

function C = ned_to_body_dcm(x)
    phi = x(7);
    theta = x(8);
    psi = x(9);

    cp = cos(phi);
    sp = sin(phi);
    ct = cos(theta);
    st = sin(theta);
    cs = cos(psi);
    ss = sin(psi);

    C = [ct*cs,              ct*ss,              -st;
         sp*st*cs-cp*ss,     sp*st*ss+cp*cs,     sp*ct;
         cp*st*cs+sp*ss,     cp*st*ss-sp*cs,     cp*ct];
end

function out = tf(cond,a,b)
    if cond
        out = a;
    else
        out = b;
    end
end
