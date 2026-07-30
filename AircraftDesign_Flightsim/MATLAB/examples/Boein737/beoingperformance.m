%% B737-800 PERFORMANCE ANALYSIS
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

datcom_dir = "D:\aircraft-flight-simulator - Copy\AircraftDesign_Flightsim\MATLAB\examples\Boein737\data";

if ~isfolder(datcom_dir)
    error("B737:DATCOMDirectoryNotFound", "Boeing DATCOM directory not found: %s",datcom_dir);
end

datcom_listing = dir(fullfile(datcom_dir,"*.out"));

if isempty(datcom_listing)
    error("B737:NoDATCOMFiles", "No DATCOM .out files were found in: %s",datcom_dir);
end

datcom_files = string(fullfile({datcom_listing.folder},{datcom_listing.name})).';

[lookup,b737_datcom_data] = b737datcom(datcom_files,true);

if numel(b737_datcom_data.source_files) ~= 14
    error("B737:UnexpectedDATCOMFileCount", "Expected 14 DATCOM files but the parser received %d.", numel(b737_datcom_data.source_files));
end

requested_mach = [0.10 0.15 0.20 0.25 0.30 0.35 0.40 0.50 0.60 0.65 0.70 0.72 0.76 0.82];
loaded_mach = b737_datcom_data.grid.mach_vec(:).';
missing_mach = requested_mach(~arrayfun(@(m) any(abs(loaded_mach-m) < 1e-8),requested_mach));

fprintf("\n=== B737 DATCOM COVERAGE ===\n");
fprintf("Files read          : %d\n",numel(b737_datcom_data.source_files));
fprintf("Usable Mach values  :");
fprintf(" %.2f",loaded_mach);
fprintf("\n");
fprintf("Unavailable stations:");
fprintf(" %.2f",missing_mach);
fprintf("\n");

if isempty(loaded_mach)
    error("B737:EmptyDATCOMMachGrid", "The DATCOM database contains no usable Mach stations.");
end

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
ac.geometry.aero_mach_max = max(loaded_mach);

%% Frames

ac.set_body_frame("body");
ac.set_reference_frame("body");

add_or_update_frame(ac,"aero_ref",   "body",cg_to_aero_ref, @(x) eye(3));
add_or_update_frame(ac,"engine_L",   "body",cg_to_engine_L, @(x) eye(3));
add_or_update_frame(ac,"engine_R",   "body",cg_to_engine_R, @(x) eye(3));
add_or_update_frame(ac,"cg",         "body",[0;0;0],        @(x) eye(3));

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

% Cruise trim and airborne performance use the continuous rating.
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

ac.ensure_gravity_source([]);

if ~ac.has_gravity_source() || ac.count_gravity_sources() ~= 1
    error("B737:GravityInstallation", "Aircraft must contain exactly one automatic gravity source.");
end

%% Sync mass / CG

[m_total,cg_total,I_total] = ac.compute_total_mass_properties([]);

ac.update_frame_position("cg",cg_total);
ac.ensure_gravity_source([]);
ac.set_reference_frame("cg");

W = m_total * ac.g;

n_cs = numel(ac.control_surfaces);
n_pe = numel(ac.propulsive_elements);
n_total = n_cs + n_pe;

fprintf("\n=== B737 READY ===\n");
fprintf("Mass    : %.2f kg\n",m_total);
fprintf("Weight  : %.2f N\n",W);
fprintf("CG body : [% .4f % .4f % .4f] m\n",cg_total);
fprintf("Aero ref: [% .4f % .4f % .4f] m body\n",cg_to_aero_ref);
fprintf("Engine L: [% .4f % .4f % .4f] m body\n",cg_to_engine_L);
fprintf("Engine R: [% .4f % .4f % .4f] m body\n",cg_to_engine_R);
fprintf("Aero Mach validity: %.2f to %.2f\n", min(loaded_mach),max(loaded_mach));
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

[~,a_trim,~,~] = ac.get_atmosphere(alt_trim_m);
V_trim = mach_trim * a_trim;

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
    "Display","off", ...
    "MaxIterations",5000, ...
    "MaxFunctionEvaluations",50000, ...
    "OptimalityTolerance",1e-10, ...
    "StepTolerance",1e-10, ...
    "ConstraintTolerance",1e-10, ...
    "FunctionTolerance",1e-10, ...
    "FiniteDifferenceType","central", ...
    "FiniteDifferenceStepSize",1e-4);

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
ac.ensure_gravity_source(x_trim);
ac.set_reference_frame("cg");

[F_trim,M_trim,fuel_flow_trim] = ac.compute_total_loads_about_cg(x_trim,u_trim);

[~,~,coeff_trim] = aero_model.get_FM(x_trim,u_trim,ac.geometry,ac);

W_trim = m_trim * ac.g;
cbar_trim = c_ref;

air_trim = ac.get_air_data(x_trim);
alpha_trim_deg = rad2deg(air_trim.alpha_rad);
qbar_trim = air_trim.qbar_Pa;

if isfield(coeff_trim,'datcom_valid') && ~coeff_trim.datcom_valid
    error("B737:InvalidBaselineDATCOMPoint", "The baseline trim point lies outside valid DATCOM coverage.");
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

%% Takeoff

fprintf("\n=== RUNNING B737 TAKEOFF ANALYSIS ===\n");

eng_L.set_rating_mode("takeoff");
eng_R.set_rating_mode("takeoff");

airframe.set_mass_properties(takeoff_mass_kg,[0;0;0],I_body);

[m_to,cg_to,~] = ac.compute_total_mass_properties([]);
ac.update_frame_position("cg",cg_to);
ac.ensure_gravity_source([]);

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
ac.ensure_gravity_source([]);

ld = ac.get_landing();

[landing_distance,res_ld] = ld.calculate_landing(0,0,3000);

fprintf("\n=== LANDING SUMMARY ===\n");
fprintf("Distance         : %.1f m\n",landing_distance);
fprintf("Vstall landing   : %.2f m/s\n",res_ld.V_stall);
fprintf("V approach       : %.2f m/s\n",res_ld.V_approach);
fprintf("V touchdown      : %.2f m/s\n",res_ld.V_touchdown);
fprintf("Runway adequate  : %d\n",res_ld.runway_adequate);

%% Restore cruise mass before stability

eng_L.set_rating_mode("continuous");
eng_R.set_rating_mode("continuous");

airframe.set_mass_properties(cruise_mass_kg,[0;0;0],I_body);

[m_cr,cg_cr,~] = ac.compute_total_mass_properties([]);
ac.update_frame_position("cg",cg_cr);
ac.ensure_gravity_source([]);
ac.set_reference_frame("cg");

ac.state.set_full_state(x_trim);
ac.set_controls_from_vector(u_trim);

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



%% ============================================================
% B737 FIXED-ALTITUDE AND EXTENDED PERFORMANCE
% ============================================================

fprintf("\n=== B737 FIXED-ALTITUDE PERFORMANCE ===\n");

% Restore cruise mass and continuous engine rating before performance work.
airframe.set_mass_properties(cruise_mass_kg,[0;0;0],I_body);

eng_L.set_rating_mode("continuous");
eng_R.set_rating_mode("continuous");

[m_perf,cg_perf,~] = ac.compute_total_mass_properties([]);
ac.update_frame_position("cg",cg_perf);
ac.ensure_gravity_source([]);
ac.set_reference_frame("cg");

ac.state.set_full_state(x_trim);
ac.set_controls_from_vector(u_trim);

perf_alt_m = alt_trim_m;
[~,a_perf,~,~] = ac.get_atmosphere(perf_alt_m);
aero_mach_limit = max(loaded_mach);

M_sweep = linspace(0.20,aero_mach_limit,60).';

baseline_z = [air_trim.alpha_rad; u_trim(2); mean(u_trim(n_cs+1:n_cs+n_pe))];

sweep = run_b737_level_sweep(ac,lookup,solver,trim_cfg,solver_cfg, perf_alt_m,M_sweep,mach_trim,baseline_z, n_cs,n_pe,S_ref,eng_L,eng_R);

b737_debug_table = sweep.debug_table;
b737_perf_table = sweep.perf_table;

disp(b737_debug_table);

if isempty(b737_perf_table) || height(b737_perf_table) == 0
    error("B737:NoValidPerformancePoints", "No valid B737 performance points were found.");
end

%% Fixed-altitude key metrics

[drag_min,idx_drag] = min(b737_perf_table.Drag_N);
[ld_max,idx_ld] = max(b737_perf_table.L_D);
[pmin,idx_pmin] = min(b737_perf_table.P_required_kW);
[roc_max,idx_roc] = max(b737_perf_table.ROC_fpm);
[range_max,idx_range] = max(b737_perf_table.FuelRange_nmpkg);
[endurance_max,idx_endurance] = max(b737_perf_table.FuelEndurance_hrpkg);

idx_positive = find(b737_perf_table.T_excess_N >= 0,1,"last");

if isempty(idx_positive)
    idx_positive = 1;
    vmax_status = "no_nonnegative_excess_thrust_point";
elseif b737_perf_table.SourceIndex(idx_positive) == numel(M_sweep)
    vmax_status = "sweep_upper_bound";
else
    vmax_status = "thrust_or_trim_boundary";
end

b737_key_table = table( ...
    ["Minimum drag"; ...
     "Maximum L/D"; ...
     "Minimum power required"; ...
     "Maximum excess-power ROC"; ...
     "Best fuel range"; ...
     "Best fuel endurance"; ...
     "Maximum level speed in sweep"], ...
    [b737_perf_table.Mach(idx_drag); ...
     b737_perf_table.Mach(idx_ld); ...
     b737_perf_table.Mach(idx_pmin); ...
     b737_perf_table.Mach(idx_roc); ...
     b737_perf_table.Mach(idx_range); ...
     b737_perf_table.Mach(idx_endurance); ...
     b737_perf_table.Mach(idx_positive)], ...
    [b737_perf_table.V_kts(idx_drag); ...
     b737_perf_table.V_kts(idx_ld); ...
     b737_perf_table.V_kts(idx_pmin); ...
     b737_perf_table.V_kts(idx_roc); ...
     b737_perf_table.V_kts(idx_range); ...
     b737_perf_table.V_kts(idx_endurance); ...
     b737_perf_table.V_kts(idx_positive)], ...
    [drag_min; ...
     ld_max; ...
     pmin; ...
     roc_max; ...
     range_max; ...
     endurance_max; ...
     b737_perf_table.V_kts(idx_positive)], ...
    ["N";"-";"kW";"ft/min";"nm/kg";"hr/kg";"kt"], ...
    'VariableNames',{'Metric','Mach','Speed_kts','Value','Unit'});

fprintf("\n=== B737 PERFORMANCE TABLE AT %.0f m ===\n",perf_alt_m);
disp(b737_perf_table);

fprintf("\n=== B737 KEY PERFORMANCE METRICS AT %.0f m ===\n",perf_alt_m);
disp(b737_key_table);
fprintf("Sweep maximum-level-speed status: %s\n",vmax_status);

%% Fixed-altitude plots

figure;
plot(b737_perf_table.Mach,b737_perf_table.Drag_N,"LineWidth",1.5);
hold on;
plot(b737_perf_table.Mach,b737_perf_table.T_available_N,"LineWidth",1.5);
grid on;
xlabel("Mach");
ylabel("Force [N]");
legend("Thrust required","Thrust available","Location","best");
title("B737 thrust required and available");

figure;
plot(b737_perf_table.Mach,b737_perf_table.P_required_kW,"LineWidth",1.5);
hold on;
plot(b737_perf_table.Mach,b737_perf_table.P_available_kW,"LineWidth",1.5);
grid on;
xlabel("Mach");
ylabel("Propulsive power [kW]");
legend("Required","Available","Location","best");
title("B737 propulsive power required and available");

figure;
plot(b737_perf_table.Mach,b737_perf_table.ROC_fpm,"LineWidth",1.5);
hold on;
yline(0,"--");
grid on;
xlabel("Mach");
ylabel("Excess-power ROC [ft/min]");
title("B737 excess-power performance");

figure;
plot(b737_perf_table.Mach,b737_perf_table.L_D,"LineWidth",1.5);
grid on;
xlabel("Mach");
ylabel("L/D");
title("B737 lift-to-drag ratio");

figure;
plot(b737_perf_table.Mach,b737_perf_table.FuelFlow_kgph,"LineWidth",1.5);
grid on;
xlabel("Mach");
ylabel("Total fuel flow [kg/hr]");
title("B737 total fuel flow");

figure;
plot(b737_perf_table.Mach,b737_perf_table.FuelRange_nmpkg,"LineWidth",1.5);
grid on;
xlabel("Mach");
ylabel("Specific range [nm/kg]");
title("B737 fuel specific range");

figure;
yyaxis left
plot(b737_perf_table.Mach,b737_perf_table.Elevator_deg,"LineWidth",1.5);
ylabel("Elevator [deg]");

yyaxis right
plot(b737_perf_table.Mach,b737_perf_table.ThrottleMean,"LineWidth",1.5);
ylabel("Mean throttle");

grid on;
xlabel("Mach");
title("B737 trim controls");

%% ============================================================
% OPTIMIZED AND EXTENDED PERFORMANCE
% ============================================================

fprintf("\n=== B737 OPTIMIZED AND EXTENDED PERFORMANCE ===\n");

required_methods = [ ...
    "analyze_paper_trim_suite"; ...
    "compute_vn_diagram"; ...
    "compute_instantaneous_turn_curve"; ...
    "compute_sustained_turn_curve"; ...
    "compute_specific_excess_power_map"; ...
    "compute_acceleration_schedule"; ...
    "analyze_glide_performance"; ...
    "optimize_climb_schedule"; ...
    "compute_mach_altitude_envelope"];

available_methods = string(methods("PerformanceAnalysis"));
missing_methods = required_methods(~ismember(required_methods,available_methods));

if ~isempty(missing_methods)
    error("B737:PerformanceAnalysisVersion", "PerformanceAnalysis is missing: %s", strjoin(missing_methods,", "));
end

perf = PerformanceAnalysis(ac);

perf_condition = struct();
perf_condition.altitude_m = perf_alt_m;

solver_cfg_extended = solver_cfg;
solver_cfg_extended.fmincon_options = optimoptions(solver_cfg.fmincon_options, "Display","none", "MaxIterations",1800, "MaxFunctionEvaluations",40000, "FiniteDifferenceType","central", "FiniteDifferenceStepSize",1e-4);

perf_cfg = struct();
perf_cfg.available_throttle = 1;
perf_cfg.takeoff_throttle = 1;
perf_cfg.continuous_throttle = 1;
perf_cfg.propulsion_type = "thrust";
perf_cfg.range_endurance_metric = "actual_flow";
perf_cfg.CL_max = 1.8;
perf_cfg.CL_min = -1.0;
perf_cfg.n_limit_pos = 2.5;
perf_cfg.n_limit_neg = -1.0;
perf_cfg.max_load_factor = 2.5;
perf_cfg.min_load_factor = -1.0;
perf_cfg.max_Mach = aero_mach_limit;
perf_cfg.service_ceiling_threshold_fpm = 100;

V_exact_lb = min(b737_perf_table.V_mps);
V_exact_ub = min(aero_mach_limit*a_perf, max(b737_perf_table.V_mps));

if V_exact_ub <= V_exact_lb
    error("B737:InvalidPerformanceRange", "No usable optimized-performance speed range remains.");
end

maneuver_bounds = struct();
maneuver_bounds.V_lb = V_exact_lb;
maneuver_bounds.V_ub = V_exact_ub;
maneuver_bounds.V0 = min(max(V_trim,V_exact_lb),V_exact_ub);
maneuver_bounds.bank0_deg = 30;
maneuver_bounds.bank_min_deg = 0.1;
maneuver_bounds.bank_max_deg = 60;
maneuver_bounds.alpha0_deg = alpha_trim_deg;
maneuver_bounds.alpha_min_deg = -5;
maneuver_bounds.alpha_max_deg = 15;
maneuver_bounds.elevator0_deg = rad2deg(u_trim(2));
maneuver_bounds.elevator_min_deg = -25;
maneuver_bounds.elevator_max_deg = 25;
maneuver_bounds.aileron0_deg = 0;
maneuver_bounds.aileron_min_deg = -20;
maneuver_bounds.aileron_max_deg = 20;
maneuver_bounds.rudder0_deg = 0;
maneuver_bounds.rudder_min_deg = -30;
maneuver_bounds.rudder_max_deg = 30;
maneuver_bounds.throttle0 = mean(u_trim(n_cs+1:n_cs+n_pe));
maneuver_bounds.throttle_min = 0;
maneuver_bounds.throttle_max = 1;
maneuver_bounds.turn_rate_ub_radps = 0.5;

V_turn = linspace(V_exact_lb,V_exact_ub,18).';
V_ps = linspace(V_exact_lb,V_exact_ub,14).';
n_ps = [1;1.25;1.5;1.75;2;2.25;2.5];
V_accel = linspace(V_exact_lb,V_exact_ub,26).';
V_vn = linspace(0,V_exact_ub,260).';
altitude_schedule_m = (0:1000:15000).';

b737_extended_results = struct();

%% Optimized airborne metrics

bounds = struct();

bounds.stall = struct( ...
    'V_lb',max(45,0.75*V_exact_lb), ...
    'V_ub',min(140,V_exact_ub), ...
    'V0',max(V_exact_lb,75), ...
    'alpha0_deg',10, ...
    'alpha_min_deg',-5, ...
    'alpha_max_deg',15, ...
    'elevator0_deg',rad2deg(u_trim(2)), ...
    'elevator_min_deg',-25, ...
    'elevator_max_deg',25, ...
    'throttle0',0.35, ...
    'throttle_min',0, ...
    'throttle_max',1);

bounds.climb = struct( ...
    'V_lb',V_exact_lb, ...
    'V_ub',V_exact_ub, ...
    'V0',min(max(V_trim,V_exact_lb),V_exact_ub), ...
    'alpha0_deg',alpha_trim_deg, ...
    'alpha_min_deg',-5, ...
    'alpha_max_deg',15, ...
    'elevator0_deg',rad2deg(u_trim(2)), ...
    'elevator_min_deg',-25, ...
    'elevator_max_deg',25, ...
    'gamma0_deg',4, ...
    'gamma_min_deg',0, ...
    'gamma_max_deg',20, ...
    'throttle0',1, ...
    'throttle_min',0, ...
    'throttle_max',1);

bounds.range = struct( ...
    'V_lb',V_exact_lb, ...
    'V_ub',V_exact_ub, ...
    'V0',b737_perf_table.V_mps(idx_range), ...
    'alpha0_deg',b737_perf_table.Alpha_deg(idx_range), ...
    'alpha_min_deg',-5, ...
    'alpha_max_deg',15, ...
    'elevator0_deg',b737_perf_table.Elevator_deg(idx_range), ...
    'elevator_min_deg',-25, ...
    'elevator_max_deg',25, ...
    'throttle0',b737_perf_table.ThrottleMean(idx_range), ...
    'throttle_min',0, ...
    'throttle_max',1);

bounds.endurance = bounds.range;
bounds.endurance.V0 = b737_perf_table.V_mps(idx_endurance);
bounds.endurance.alpha0_deg = b737_perf_table.Alpha_deg(idx_endurance);
bounds.endurance.elevator0_deg = b737_perf_table.Elevator_deg(idx_endurance);
bounds.endurance.throttle0 = b737_perf_table.ThrottleMean(idx_endurance);

bounds.glide = struct( ...
    'V_lb',V_exact_lb, ...
    'V_ub',V_exact_ub, ...
    'V0',b737_perf_table.V_mps(idx_ld), ...
    'alpha0_deg',b737_perf_table.Alpha_deg(idx_ld), ...
    'alpha_min_deg',-5, ...
    'alpha_max_deg',15, ...
    'elevator0_deg',b737_perf_table.Elevator_deg(idx_ld), ...
    'elevator_min_deg',-25, ...
    'elevator_max_deg',25, ...
    'gamma0_deg',-3, ...
    'gamma_min_deg',-15, ...
    'gamma_max_deg',-0.01);

bounds.max_speed = struct( ...
    'V_lb',max(V_trim,V_exact_lb), ...
    'V_ub',aero_mach_limit*a_perf, ...
    'V0',min(0.95*aero_mach_limit*a_perf, ...
             aero_mach_limit*a_perf-1), ...
    'alpha0_deg',0.5, ...
    'alpha_min_deg',-5, ...
    'alpha_max_deg',15, ...
    'elevator0_deg',rad2deg(u_trim(2)), ...
    'elevator_min_deg',-25, ...
    'elevator_max_deg',25, ...
    'throttle0',0.9, ...
    'throttle_min',0, ...
    'throttle_max',1);

bounds.V_sweep_mps = V_accel;
bounds.alpha0_deg = alpha_trim_deg;
bounds.alpha_min_deg = -5;
bounds.alpha_max_deg = 15;
bounds.elevator0_deg = rad2deg(u_trim(2));
bounds.elevator_min_deg = -25;
bounds.elevator_max_deg = 25;
bounds.throttle0 = mean(u_trim(n_cs+1:n_cs+n_pe));
bounds.throttle_min = 0;
bounds.throttle_max = 1;
bounds.gamma0_deg = 3;
bounds.gamma_min_deg = -15;
bounds.gamma_max_deg = 20;

try
    optimized_suite = perf.analyze_paper_trim_suite(perf_condition,solver_cfg_extended,perf_cfg,bounds);

    optimized_suite.best_range = robust_level_metric_solution(perf,perf_condition,solver_cfg_extended,perf_cfg, bounds.range,"range",b737_perf_table,idx_range, optimized_suite.best_range);

    optimized_suite.best_endurance = robust_level_metric_solution(perf,perf_condition,solver_cfg_extended,perf_cfg, bounds.endurance,"endurance",b737_perf_table,idx_endurance, optimized_suite.best_endurance);

    optimized_suite.best_glide = robust_glide_solution(perf,perf_condition,solver_cfg_extended,perf_cfg, bounds.glide,b737_perf_table,idx_ld, optimized_suite.best_glide);

    b737_extended_results.optimized_suite = optimized_suite;

    fprintf("\n--- OPTIMIZED AIRBORNE PERFORMANCE ---\n");

    metric_names = ["Minimum powered trim speed"; "Best angle climb"; "Best rate climb"; "Best fuel range"; "Best fuel endurance"; "Best glide"; "Maximum level speed"];

    metric_results = { ...
        optimized_suite.stall_speed; ...
        optimized_suite.best_angle_climb; ...
        optimized_suite.best_rate_climb; ...
        optimized_suite.best_range; ...
        optimized_suite.best_endurance; ...
        optimized_suite.best_glide; ...
        optimized_suite.maximum_level_speed};

    for k = 1:numel(metric_results)
        result_k = metric_results{k};

        if isstruct(result_k) && isfield(result_k,'converged') && result_k.converged && isfield(result_k,'point_opt')

            p = result_k.point_opt;
            source = "direct_optimization";

            if isfield(result_k,'solution_source')
                source = string(result_k.solution_source);
            end

            fprintf("%-28s : M %.4f | %7.2f kt | alpha %7.3f deg | gamma %7.3f deg | throttle %.4f | %s\n", metric_names(k),p.Mach,p.V_kts,p.alpha_deg, p.gamma_deg,p.throttle,source);
        else
            fprintf("%-28s : not converged\n",metric_names(k));
        end
    end

catch ME
    b737_extended_results.optimized_suite = struct();
    b737_extended_results.optimized_suite_error = string(ME.identifier)+": "+string(ME.message);

    warning("B737:OptimizedSuite", "Optimized airborne suite failed: %s",ME.message);
end

%% V-n envelope

try
    vn = perf.compute_vn_diagram(perf_condition,V_vn,perf_cfg);
    b737_extended_results.vn = vn;

    fprintf("\n--- V-N ENVELOPE ---\n");
    fprintf("1-g stall speed : %.3f kt TAS\n",vn.Vs_kts);
    fprintf("Positive Va     : %.3f kt TAS\n",vn.Va_kts);

    figure;
    plot(vn.V_EAS_kts,vn.n_pos,"LineWidth",1.5);
    hold on;
    plot(vn.V_EAS_kts,vn.n_neg,"LineWidth",1.5);
    yline(1,"--");
    grid on;
    xlabel("Equivalent airspeed [kt]");
    ylabel("Load factor n");
    legend("Positive envelope","Negative envelope","1 g", "Location","best");
    title("B737 maneuver V-n envelope");

catch ME
    b737_extended_results.vn = struct();
    b737_extended_results.vn_error = string(ME.identifier)+": "+string(ME.message);
    warning("B737:VN","V-n analysis failed: %s",ME.message);
end

%% Instantaneous and sustained turns

try
    turn_inst = perf.compute_instantaneous_turn_curve(perf_condition,V_turn,perf_cfg);

    b737_extended_results.instantaneous_turn = turn_inst;

    fprintf("\n--- INSTANTANEOUS TURN PERFORMANCE ---\n");
    fprintf("Maximum rate   : %.4f deg/s at %.3f kt\n", turn_inst.max_turn_rate_degps, turn_inst.speed_at_max_turn_rate_mps*1.943844492);
    fprintf("Minimum radius : %.3f m at %.3f kt\n", turn_inst.min_turn_radius_m, turn_inst.speed_at_min_radius_mps*1.943844492);

catch ME
    turn_inst = struct();
    b737_extended_results.instantaneous_turn = struct();
    b737_extended_results.instantaneous_turn_error = string(ME.identifier)+": "+string(ME.message);
    warning("B737:InstantaneousTurn", "Instantaneous-turn analysis failed: %s",ME.message);
end

try
    turn_sust = perf.compute_sustained_turn_curve(perf_condition,V_turn,solver_cfg_extended, perf_cfg,maneuver_bounds);

    b737_extended_results.sustained_turn = turn_sust;

    fprintf("\n--- SUSTAINED TURN PERFORMANCE ---\n");
    fprintf("Valid points   : %d / %d\n", nnz(turn_sust.valid),numel(turn_sust.valid));
    fprintf("Maximum rate   : %.4f deg/s at %.3f kt\n", turn_sust.max_turn_rate_degps, turn_sust.speed_at_max_turn_rate_mps*1.943844492);
    fprintf("Minimum radius : %.3f m at %.3f kt\n", turn_sust.min_turn_radius_m, turn_sust.speed_at_min_radius_mps*1.943844492);

    valid_sust = turn_sust.valid & isfinite(turn_sust.turn_rate_degps);

    figure;
    if isfield(turn_inst,'valid')
        valid_inst = turn_inst.valid & isfinite(turn_inst.turn_rate_degps);
        plot(turn_inst.V_kts(valid_inst), turn_inst.turn_rate_degps(valid_inst), "o-","LineWidth",1.5);
        hold on;
    end

    plot(turn_sust.V_kts(valid_sust), turn_sust.turn_rate_degps(valid_sust), "s-","LineWidth",1.5);

    grid on;
    xlabel("True airspeed [kt]");
    ylabel("Turn rate [deg/s]");
    legend("Instantaneous","Sustained","Location","best");
    title("B737 turn performance");

catch ME
    b737_extended_results.sustained_turn = struct();
    b737_extended_results.sustained_turn_error = string(ME.identifier)+": "+string(ME.message);
    warning("B737:SustainedTurn", "Sustained-turn analysis failed: %s",ME.message);
end

%% Specific excess power map

try
    Ps_map = perf.compute_specific_excess_power_map(perf_condition,V_ps,n_ps,solver_cfg_extended, perf_cfg,maneuver_bounds);

    b737_extended_results.Ps_map = Ps_map;

    fprintf("\n--- SPECIFIC EXCESS POWER MAP ---\n");
    fprintf("Valid points : %d / %d\n", nnz(Ps_map.valid),numel(Ps_map.valid));

    Ps_plot = Ps_map.Ps_mps;
    Ps_plot(~Ps_map.valid) = NaN;

    [V_grid_ps,n_grid_ps] = meshgrid(Ps_map.V_kts,Ps_map.n);

    figure;
    contourf(V_grid_ps,n_grid_ps,Ps_plot,18, "LineStyle","none");
    hold on;
    contour(V_grid_ps,n_grid_ps,Ps_plot,[0 0],"LineWidth",2);
    colorbar;
    grid on;
    xlabel("True airspeed [kt]");
    ylabel("Load factor n");
    title("B737 specific excess power P_s [m/s]");

catch ME
    b737_extended_results.Ps_map = struct();
    b737_extended_results.Ps_map_error = string(ME.identifier)+": "+string(ME.message);
    warning("B737:PsMap","P_s map failed: %s",ME.message);
end

%% Accelerated level flight

try
    acceleration = perf.compute_acceleration_schedule(perf_condition,V_accel,solver_cfg_extended, perf_cfg,maneuver_bounds);

    b737_extended_results.acceleration = acceleration;

    fprintf("\n--- LEVEL-FLIGHT ACCELERATION ---\n");
    fprintf("Maximum acceleration : %.5f m/s^2 at %.3f kt\n", acceleration.max_acceleration_mps2, acceleration.speed_at_max_acceleration_mps*1.943844492);

    valid_accel = acceleration.valid & isfinite(acceleration.acceleration_mps2);

    figure;
    plot(acceleration.V_kts(valid_accel), acceleration.acceleration_mps2(valid_accel), "o-","LineWidth",1.5);
    hold on;
    yline(0,"--");
    grid on;
    xlabel("True airspeed [kt]");
    ylabel("Tangential acceleration [m/s^2]");
    title("B737 full-throttle level acceleration");

catch ME
    b737_extended_results.acceleration = struct();
    b737_extended_results.acceleration_error = string(ME.identifier)+": "+string(ME.message);
    warning("B737:Acceleration", "Acceleration analysis failed: %s",ME.message);
end

%% Engine-off glide and minimum sink

try
    glide = perf.analyze_glide_performance(perf_condition,solver_cfg_extended,perf_cfg,bounds.glide);

    if (~glide.best_glide.converged) && isfield(b737_extended_results,'optimized_suite') && isfield(b737_extended_results.optimized_suite,'best_glide') && b737_extended_results.optimized_suite.best_glide.converged
        glide.best_glide = b737_extended_results.optimized_suite.best_glide;
    end

    b737_extended_results.glide = glide;

    fprintf("\n--- ENGINE-OFF GLIDE PERFORMANCE ---\n");

    if glide.best_glide.converged
        p = glide.best_glide.point_opt;
        fprintf("Best glide : %.3f kt | gamma %.3f deg | L/D %.4f\n", p.V_kts,p.gamma_deg,p.L_over_D);
    else
        fprintf("Best glide did not converge.\n");
    end

    if glide.min_sink.converged
        p = glide.min_sink.point_opt;
        fprintf("Min sink   : %.3f kt | %.3f ft/min\n", p.V_kts,p.sink_rate_fpm);
    else
        fprintf("Minimum sink did not converge.\n");
    end

catch ME
    b737_extended_results.glide = struct();
    b737_extended_results.glide_error = string(ME.identifier)+": "+string(ME.message);
    warning("B737:Glide","Glide analysis failed: %s",ME.message);
end

%% Climb schedule and ceilings

try
    climb = perf.optimize_climb_schedule(perf_condition,altitude_schedule_m, solver_cfg_extended,perf_cfg,bounds.climb);

    b737_extended_results.climb_schedule = climb;

    fprintf("\n--- CLIMB SCHEDULE AND CEILINGS ---\n");

    if isfinite(climb.service_ceiling_m)
        fprintf("Service ceiling : %.1f m  %.0f ft\n", climb.service_ceiling_m, climb.service_ceiling_m*3.280839895);
    else
        fprintf("Service ceiling : not crossed through %.0f m\n", max(altitude_schedule_m));
    end

    if isfinite(climb.absolute_ceiling_m)
        fprintf("Absolute ceiling: %.1f m  %.0f ft\n", climb.absolute_ceiling_m, climb.absolute_ceiling_m*3.280839895);
    else
        fprintf("Absolute ceiling: not crossed through %.0f m\n", max(altitude_schedule_m));
    end

    valid_climb = climb.point_valid & isfinite(climb.best_ROC_fpm);

    figure;
    plot(climb.best_ROC_fpm(valid_climb), climb.altitude_ft(valid_climb), "o-","LineWidth",1.5);
    hold on;
    xline(perf_cfg.service_ceiling_threshold_fpm,"--");
    xline(0,"--");
    grid on;
    xlabel("Best ROC [ft/min]");
    ylabel("Altitude [ft]");
    title("B737 climb schedule");

catch ME
    b737_extended_results.climb_schedule = struct();
    b737_extended_results.climb_schedule_error = string(ME.identifier)+": "+string(ME.message);
    warning("B737:ClimbSchedule", "Climb schedule failed: %s",ME.message);
end

%% Mach-altitude envelope

envelope_bounds = struct();
envelope_bounds.V_lb = max(45,0.8*V_exact_lb);
envelope_bounds.V_ub = aero_mach_limit*a_perf;
envelope_bounds.max_Mach = aero_mach_limit;
envelope_bounds.stall_speed_factor = 1;
envelope_bounds.alpha0_deg = alpha_trim_deg;
envelope_bounds.alpha_min_deg = -5;
envelope_bounds.alpha_max_deg = 15;
envelope_bounds.elevator0_deg = rad2deg(u_trim(2));
envelope_bounds.elevator_min_deg = -25;
envelope_bounds.elevator_max_deg = 25;

try
    envelope = perf.compute_mach_altitude_envelope(perf_condition,altitude_schedule_m, solver_cfg_extended,perf_cfg,envelope_bounds);

    b737_extended_results.mach_altitude_envelope = envelope;

    valid_env = envelope.valid & isfinite(envelope.max_level_Mach);

    fprintf("\n--- MACH-ALTITUDE ENVELOPE ---\n");
    fprintf("Valid altitude points: %d / %d\n", nnz(valid_env),numel(valid_env));

    figure;
    plot(envelope.stall_Mach,envelope.altitude_ft, "o-","LineWidth",1.5);
    hold on;
    plot(envelope.max_level_Mach(valid_env), envelope.altitude_ft(valid_env), "s-","LineWidth",1.5);
    grid on;
    xlabel("Mach");
    ylabel("Altitude [ft]");
    legend("Stall boundary","Maximum level speed","Location","best");
    title("B737 Mach-altitude envelope");

catch ME
    b737_extended_results.mach_altitude_envelope = struct();
    b737_extended_results.mach_altitude_envelope_error = string(ME.identifier)+": "+string(ME.message);
    warning("B737:MachAltitudeEnvelope", "Mach-altitude envelope failed: %s",ME.message);
end

%% Engine model diagnostic at cruise trim

eng_L.set_rating_mode("continuous");
eng_R.set_rating_mode("continuous");

ac.state.set_full_state(x_trim);
ac.set_controls_from_vector(u_trim);
ac.compute_total_loads(x_trim,u_trim);

fprintf("\n--- CFM56-7B26/3 CRUISE-TRIM DIAGNOSTIC ---\n");
fprintf("Left engine thrust  : %.3f kN\n",eng_L.last_debug.thrust_N/1000);
fprintf("Right engine thrust : %.3f kN\n",eng_R.last_debug.thrust_N/1000);
fprintf("Left TSFC           : %.4f lb/(lbf hr)\n", eng_L.last_debug.tsfc_lb_lbf_hr);
fprintf("Right TSFC          : %.4f lb/(lbf hr)\n", eng_R.last_debug.tsfc_lb_lbf_hr);
fprintf("Total fuel flow     : %.4f kg/s  %.1f kg/hr\n", eng_L.last_debug.fuel_flow_kgps+eng_R.last_debug.fuel_flow_kgps, 3600*(eng_L.last_debug.fuel_flow_kgps+ eng_R.last_debug.fuel_flow_kgps));

%% Save

b737_results = struct();

b737_results.configuration = struct( ...
    'datcom_directory',datcom_dir, ...
    'datcom_files',b737_datcom_data.source_files, ...
    'datcom_loaded_mach',loaded_mach, ...
    'datcom_missing_mach',missing_mach, ...
    'datcom_control_model_status', ...
        b737_datcom_data.control_model_status, ...
    'engine_variant',eng_L.variant, ...
    'engine_deck_status',eng_L.performance_deck_status, ...
    'performance_altitude_m',perf_alt_m, ...
    'trim_mach',mach_trim, ...
    'aero_mach_limit',aero_mach_limit, ...
    'M_sweep',M_sweep);

b737_results.baseline = struct('x_trim',x_trim, 'u_trim',u_trim, 'converged',converged, 'fuel_flow_kgps',fuel_flow_trim, 'coeff',coeff_trim);

b737_results.takeoff = res_to;
b737_results.landing = res_ld;
b737_results.stability = s;
b737_results.debug_table = b737_debug_table;
b737_results.performance_table = b737_perf_table;
b737_results.key_table = b737_key_table;
b737_results.extended = b737_extended_results;

assignin("base","ac",ac);
assignin("base","x_trim",x_trim);
assignin("base","u_trim",u_trim);
assignin("base","b737_results",b737_results);
assignin("base","b737_perf_table",b737_perf_table);
assignin("base","b737_key_table",b737_key_table);
assignin("base","b737_debug_table",b737_debug_table);
assignin("base","b737_extended_results",b737_extended_results);

script_path = mfilename("fullpath");
script_dir = fileparts(script_path);

if strlength(string(script_dir)) == 0 || ~isfolder(script_dir)
    script_dir = pwd;
end

results_dir = fullfile(script_dir,"results");

if ~isfolder(results_dir)
    [ok,msg] = mkdir(results_dir);

    if ~ok
        warning("B737:ResultsDirectory", "Could not create results directory: %s",msg);
        results_dir = fullfile(tempdir,"B737_Performance_Results");
        if ~isfolder(results_dir)
            mkdir(results_dir);
        end
    end
end

results_file = fullfile(results_dir, "b737_800_cfm56_7b26_performance_results.mat");
b737_results.results_file = string(results_file);

try
    save(char(results_file),"b737_results","-v7.3");
catch ME
    warning("B737:PrimarySaveFailed", "Primary save failed: %s",ME.message);

    results_dir = fullfile(tempdir,"B737_Performance_Results");

    if ~isfolder(results_dir)
        mkdir(results_dir);
    end

    results_file = fullfile(results_dir, "b737_800_cfm56_7b26_performance_results.mat");

    b737_results.results_file = string(results_file);
    save(char(results_file),"b737_results","-v7.3");
end

assignin("base","b737_results_file",string(results_file));

fprintf("\n=== B737-800 CFM56-7B26/3 PERFORMANCE COMPLETE ===\n");
fprintf("Saved: %s\n",char(results_file));


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

function out = tf(cond,a,b)
    if cond
        out = a;
    else
        out = b;
    end
end

function opt = robust_level_metric_solution(perf,condition,solver_cfg,perf_cfg,bounds,objective, perf_table,preferred_index,direct_opt)

    if isstruct(direct_opt) && isfield(direct_opt,'converged') && direct_opt.converged
        opt = direct_opt;
        opt.solution_source = "direct_optimization";
        return;
    end

    V_lb = bounds.V_lb;
    V_ub = bounds.V_ub;
    V_pref = perf_table.V_mps(preferred_index);

    V_seeds = unique(max(V_lb,min(V_ub,[V_pref; 0.92*V_pref; 1.08*V_pref; V_lb+0.35*(V_ub-V_lb); V_lb+0.65*(V_ub-V_lb)])));

    best_opt = [];
    best_J = Inf;

    perf_local = perf_cfg;
    perf_local.range_endurance_metric = "classical";

    condition_free = condition;
    if isfield(condition_free,'mach')
        condition_free.mach = [];
    end
    if isfield(condition_free,'velocity_mps')
        condition_free.velocity_mps = [];
    end

    for k = 1:numel(V_seeds)
        [~,idx_near] = min(abs(perf_table.V_mps-V_seeds(k)));

        problem = struct();
        problem.mode = "level";
        problem.objective = string(objective);
        problem.variables = ["velocity","alpha","control_pitch","throttle"];

        problem.initial_guess = [V_seeds(k); deg2rad(perf_table.Alpha_deg(idx_near)); deg2rad(perf_table.Elevator_deg(idx_near)); perf_table.ThrottleMean(idx_near)];

        problem.lb = [V_lb; deg2rad(bounds.alpha_min_deg); deg2rad(bounds.elevator_min_deg); bounds.throttle_min];

        problem.ub = [V_ub; deg2rad(bounds.alpha_max_deg); deg2rad(bounds.elevator_max_deg); bounds.throttle_max];

        problem.fixed = struct('beta',0, 'phi',0, 'psi',0, 'gamma',0, 'control_roll',0, 'control_yaw',0);

        try
            trial = perf.solve_performance_point(condition_free,problem,solver_cfg,perf_local);

            if trial.converged && isfinite(trial.objective_value) && trial.objective_value < best_J
                best_opt = trial;
                best_J = trial.objective_value;
            end
        catch
        end
    end

    if ~isempty(best_opt)
        opt = best_opt;
        opt.solution_source = "multistart_exact_optimization";
        return;
    end

    % Gradient fallback for piecewise-linear aerodynamic tables:
    % choose the sweep optimum, then solve the exact trim equations at that
    % fixed speed. This preserves exact equilibrium while selecting the
    % performance speed from the transparent sweep.
    V_target = perf_table.V_mps(preferred_index);

    condition_fixed = condition;
    condition_fixed.velocity_mps = V_target;
    if isfield(condition_fixed,'mach')
        condition_fixed.mach = [];
    end

    spec = struct();
    spec.mode = "level";
    spec.variables = ["alpha","control_pitch","throttle"];
    spec.initial_guess = [deg2rad(perf_table.Alpha_deg(preferred_index)); deg2rad(perf_table.Elevator_deg(preferred_index)); perf_table.ThrottleMean(preferred_index)];

    spec.lb = [deg2rad(bounds.alpha_min_deg); deg2rad(bounds.elevator_min_deg); bounds.throttle_min];

    spec.ub = [deg2rad(bounds.alpha_max_deg); deg2rad(bounds.elevator_max_deg); bounds.throttle_max];

    spec.fixed = struct('beta',0, 'phi',0, 'psi',0, 'gamma',0, 'control_roll',0, 'control_yaw',0);

    spec.reference_frame_name = "cg";
    spec.residual_scale = [perf.aircraft.g;perf.aircraft.g;1];

    [x,u,conv,info] = perf.solve_trim(condition_fixed,spec,solver_cfg);

    point = perf.evaluate_trim_point(x,u,condition_fixed,perf_local,info.extras);

    point.trim_converged = conv;
    point.performance_converged = conv;
    point.valid = conv;
    point.residual = info.residual;
    point.residual_scaled = info.residual_scaled;
    point.residual_norm = info.residual_norm;
    point.residual_inf = info.residual_inf;
    point.objective = string(objective);
    point.objective_value = perf.compute_objective_value(point,objective,perf_local);

    opt = struct();
    opt.objective = string(objective);
    opt.objective_value = point.objective_value;
    opt.V_opt_mps = point.V_mps;
    opt.point_opt = point;
    opt.converged = conv;
    opt.z_star = info.z_star;
    opt.variables = info.variables;
    opt.exitflag = info.exitflag;
    opt.output = info.output;
    opt.constraint_violation = info.residual_inf;
    opt.solution_source = "sweep_selected_exact_trim_fallback";
end

function opt = robust_glide_solution(perf,condition,solver_cfg,perf_cfg,bounds, perf_table,preferred_index,direct_opt)

    if isstruct(direct_opt) && isfield(direct_opt,'converged') && direct_opt.converged
        opt = direct_opt;
        opt.solution_source = "direct_optimization";
        return;
    end

    V_lb = bounds.V_lb;
    V_ub = bounds.V_ub;
    V_pref = perf_table.V_mps(preferred_index);
    LD_pref = perf_table.L_D(preferred_index);

    gamma_pref_deg = -rad2deg(atan(1/max(LD_pref,1e-12)));

    V_seeds = unique(max(V_lb,min(V_ub,[V_pref;0.92*V_pref;1.08*V_pref; V_lb+0.35*(V_ub-V_lb); V_lb+0.65*(V_ub-V_lb)])));

    best_opt = [];
    best_J = Inf;

    perf_off = perf_cfg;
    perf_off.engine_off = true;

    condition_free = condition;
    if isfield(condition_free,'mach')
        condition_free.mach = [];
    end
    if isfield(condition_free,'velocity_mps')
        condition_free.velocity_mps = [];
    end

    for k = 1:numel(V_seeds)
        [~,idx_near] = min(abs(perf_table.V_mps-V_seeds(k)));

        problem = struct();
        problem.mode = "glide";
        problem.objective = "max_l_over_d";
        problem.engine_off = true;
        problem.variables = ["velocity","alpha","control_pitch","gamma"];

        problem.initial_guess = [V_seeds(k); deg2rad(perf_table.Alpha_deg(idx_near)); deg2rad(perf_table.Elevator_deg(idx_near)); deg2rad(gamma_pref_deg)];

        problem.lb = [V_lb; deg2rad(bounds.alpha_min_deg); deg2rad(bounds.elevator_min_deg); deg2rad(bounds.gamma_min_deg)];

        problem.ub = [V_ub; deg2rad(bounds.alpha_max_deg); deg2rad(bounds.elevator_max_deg); deg2rad(bounds.gamma_max_deg)];

        problem.fixed = struct('beta',0, 'phi',0, 'psi',0, 'control_roll',0, 'control_yaw',0, 'throttle',0);

        try
            trial = perf.solve_performance_point(condition_free,problem,solver_cfg,perf_off);

            if trial.converged && isfinite(trial.objective_value) && trial.objective_value < best_J
                best_opt = trial;
                best_J = trial.objective_value;
            end
        catch
        end
    end

    if ~isempty(best_opt)
        opt = best_opt;
        opt.solution_source = "multistart_exact_optimization";
        return;
    end

    % Fixed-speed exact engine-off trim fallback at the sweep max-L/D speed.
    condition_fixed = condition;
    condition_fixed.velocity_mps = V_pref;
    if isfield(condition_fixed,'mach')
        condition_fixed.mach = [];
    end

    spec = struct();
    spec.mode = "glide";
    spec.engine_off = true;
    spec.variables = ["alpha","control_pitch","gamma"];
    spec.initial_guess = [deg2rad(perf_table.Alpha_deg(preferred_index)); deg2rad(perf_table.Elevator_deg(preferred_index)); deg2rad(gamma_pref_deg)];

    spec.lb = [deg2rad(bounds.alpha_min_deg); deg2rad(bounds.elevator_min_deg); deg2rad(bounds.gamma_min_deg)];

    spec.ub = [deg2rad(bounds.alpha_max_deg); deg2rad(bounds.elevator_max_deg); deg2rad(bounds.gamma_max_deg)];

    spec.fixed = struct('beta',0, 'phi',0, 'psi',0, 'control_roll',0, 'control_yaw',0, 'throttle',0);

    spec.reference_frame_name = "cg";
    spec.residual_scale = [perf.aircraft.g;perf.aircraft.g;1];

    [x,u,conv,info] = perf.solve_trim(condition_fixed,spec,solver_cfg);

    point = perf.evaluate_trim_point(x,u,condition_fixed,perf_off,info.extras);

    point.trim_converged = conv;
    point.performance_converged = conv;
    point.valid = conv;
    point.residual = info.residual;
    point.residual_scaled = info.residual_scaled;
    point.residual_norm = info.residual_norm;
    point.residual_inf = info.residual_inf;
    point.objective = "max_l_over_d";
    point.objective_value = -point.L_over_D;

    opt = struct();
    opt.objective = "max_l_over_d";
    opt.objective_value = point.objective_value;
    opt.V_opt_mps = point.V_mps;
    opt.point_opt = point;
    opt.converged = conv;
    opt.z_star = info.z_star;
    opt.variables = info.variables;
    opt.exitflag = info.exitflag;
    opt.output = info.output;
    opt.constraint_violation = info.residual_inf;
    opt.solution_source = "sweep_selected_exact_glide_trim_fallback";
end



function sweep = run_b737_level_sweep(ac,lookup,solver,trim_cfg,solver_cfg,alt_m,M_grid, M_seed,z_seed,n_cs,n_pe,S_ref,eng_L,eng_R)

    M_grid = M_grid(:);
    N = numel(M_grid);

    [~,a_atm,~,~] = ac.get_atmosphere(alt_m);
    V_grid = M_grid*a_atm;

    valid = false(N,1);
    converged_vec = false(N,1);
    residual_norm = nan(N,1);
    reason = strings(N,1);

    alpha_deg = nan(N,1);
    elevator_deg = nan(N,1);
    throttle_L = nan(N,1);
    throttle_R = nan(N,1);
    throttle_mean = nan(N,1);

    CL = nan(N,1);
    CD = nan(N,1);
    LD = nan(N,1);

    Drag_N = nan(N,1);
    T_trim_N = nan(N,1);
    T_available_N = nan(N,1);
    T_excess_N = nan(N,1);

    P_required_kW = nan(N,1);
    P_available_kW = nan(N,1);
    P_excess_kW = nan(N,1);

    Ps_mps = nan(N,1);
    ROC_fpm = nan(N,1);

    fuel_flow_kgps = nan(N,1);
    fuel_flow_kgph = nan(N,1);
    fuel_range_nmpkg = nan(N,1);
    fuel_endurance_hrpkg = nan(N,1);

    z_solution = nan(N,3);

    eng_L.set_rating_mode("continuous");
    eng_R.set_rating_mode("continuous");

    [~,seed_index] = min(abs(M_grid-M_seed));
    branch_indices = {seed_index:-1:1,seed_index+1:N};

    for branch = 1:numel(branch_indices)

        z_guess = z_seed(:);
        indices = branch_indices{branch};

        for kk = 1:numel(indices)

            i = indices(kk);
            M_i = M_grid(i);
            V_i = V_grid(i);

            condition_i = struct();
            condition_i.altitude_m = alt_m;
            condition_i.mach = M_i;
            condition_i.velocity_mps = V_i;

            trim_cfg_i = trim_cfg;
            trim_cfg_i.initial_guess = z_guess;

            try
                [x_i,u_i,conv_i,info_i] = solver.solve_trim(condition_i,trim_cfg_i,solver_cfg);

                converged_vec(i) = conv_i;
                residual_norm(i) = get_residual_norm(info_i);

                ac.state.set_full_state(x_i);
                ac.set_controls_from_vector(u_i);

                [m_i,cg_i,~] = ac.compute_total_mass_properties(x_i);
                ac.update_frame_position("cg",cg_i);
                ac.ensure_gravity_source(x_i);
                ac.set_reference_frame("cg");

                [~,~,ff_i] = ac.compute_total_loads_about_cg(x_i,u_i);

                air_i = ac.get_air_data(x_i);
                x_air_i = x_i;
                x_air_i(4:6) = air_i.air_velocity_body_mps;
                c_i = lookup(x_air_i,u_i,ac.geometry);

                V_body = air_i.air_velocity_body_mps;
                V_actual = norm(V_body);

                if V_actual < 1e-9
                    reason(i) = "zero velocity";
                    continue;
                end

                Vhat = V_body/V_actual;
                alpha_i = rad2deg(atan2(V_body(3),V_body(1)));

                throttle_values = u_i(n_cs+1:n_cs+n_pe);
                throttle_L_i = throttle_values(1);
                throttle_R_i = throttle_values(min(2,numel(throttle_values)));
                throttle_mean_i = mean(throttle_values);

                qbar = air_i.qbar_Pa;
                D_i = qbar*S_ref*c_i.CD;

                [F_prop_trim,~,~,~] = ac.compute_loads_by_solver(x_i,u_i,"PropulsionLoadSolver", ac.get_body_frame(),ac.get_reference_frame());

                T_trim_i = dot(F_prop_trim(:),Vhat(:));

                u_available = u_i;
                u_available(n_cs+1:n_cs+n_pe) = 1;

                ac.set_controls_from_vector(u_available);

                [F_prop_available,~,~,~] = ac.compute_loads_by_solver(x_i,u_available,"PropulsionLoadSolver", ac.get_body_frame(),ac.get_reference_frame());

                ac.set_controls_from_vector(u_i);

                T_available_i = dot(F_prop_available(:),Vhat(:));
                T_excess_i = T_available_i-D_i;

                P_required_i = D_i*V_actual;
                P_available_i = T_available_i*V_actual;
                P_excess_i = T_excess_i*V_actual;

                W_i = m_i*ac.g;
                Ps_i = P_excess_i/max(W_i,1e-12);

                trim_ok = conv_i || (isfinite(residual_norm(i)) && residual_norm(i) < 1e-3);

                model_ok = ...
                    isfield(c_i,'datcom_valid') && ...
                    c_i.datcom_valid && ...
                    isfinite(alpha_i) && ...
                    alpha_i >= -5-1e-6 && ...
                    alpha_i <= 15+1e-6 && ...
                    isfinite(c_i.CL) && ...
                    c_i.CL > 0 && ...
                    c_i.CL < 2.5 && ...
                    isfinite(c_i.CD) && ...
                    c_i.CD > 0 && ...
                    c_i.CD < 0.3 && ...
                    all(throttle_values >= -1e-6) && ...
                    all(throttle_values <= 1.0001) && ...
                    isfinite(ff_i) && ...
                    ff_i > 0;

                if ~(trim_ok && model_ok)
                    if ~trim_ok
                        reason(i) = "trim not converged";
                    else
                        reason(i) = "model validity failure";
                    end
                    continue;
                end

                valid(i) = true;
                reason(i) = "valid";

                alpha_deg(i) = alpha_i;
                elevator_deg(i) = rad2deg(u_i(2));
                throttle_L(i) = throttle_L_i;
                throttle_R(i) = throttle_R_i;
                throttle_mean(i) = throttle_mean_i;

                CL(i) = c_i.CL;
                CD(i) = c_i.CD;
                LD(i) = c_i.CL/max(c_i.CD,1e-12);

                Drag_N(i) = D_i;
                T_trim_N(i) = T_trim_i;
                T_available_N(i) = T_available_i;
                T_excess_N(i) = T_excess_i;

                P_required_kW(i) = P_required_i/1000;
                P_available_kW(i) = P_available_i/1000;
                P_excess_kW(i) = P_excess_i/1000;

                Ps_mps(i) = Ps_i;
                ROC_fpm(i) = Ps_i*196.8503937;

                fuel_flow_kgps(i) = ff_i;
                fuel_flow_kgph(i) = ff_i*3600;
                fuel_range_nmpkg(i) = V_actual/ff_i/1852;
                fuel_endurance_hrpkg(i) = 1/ff_i/3600;

                z_solution(i,:) = [atan2(V_body(3),V_body(1)), u_i(2), throttle_mean_i];

                z_guess = z_solution(i,:).';

            catch ME
                reason(i) = string(ME.identifier)+": "+string(ME.message);
            end
        end
    end

    debug_table = table( ...
        M_grid,V_grid,V_grid*1.943844492, ...
        valid,converged_vec,residual_norm, ...
        alpha_deg,CL,CD,throttle_mean, ...
        fuel_flow_kgps,reason, ...
        'VariableNames',{ ...
        'Mach','V_mps','V_kts','Valid','TrimConverged', ...
        'ResidualNorm','Alpha_deg','CL','CD','ThrottleMean', ...
        'FuelFlow_kgps','Reason'});

    idx = find(valid);

    if isempty(idx)
        perf_table = table();
        seed_z_out = z_seed(:);
    else
        perf_table = table( ...
            idx,M_grid(idx),V_grid(idx),V_grid(idx)*1.943844492, ...
            alpha_deg(idx),elevator_deg(idx), ...
            throttle_L(idx),throttle_R(idx),throttle_mean(idx), ...
            CL(idx),CD(idx),LD(idx), ...
            Drag_N(idx),T_trim_N(idx),T_available_N(idx), ...
            T_excess_N(idx), ...
            P_required_kW(idx),P_available_kW(idx),P_excess_kW(idx), ...
            Ps_mps(idx),ROC_fpm(idx), ...
            fuel_flow_kgps(idx),fuel_flow_kgph(idx), ...
            fuel_range_nmpkg(idx),fuel_endurance_hrpkg(idx), ...
            residual_norm(idx), ...
            'VariableNames',{ ...
            'SourceIndex','Mach','V_mps','V_kts','Alpha_deg', ...
            'Elevator_deg','Throttle_L','Throttle_R','ThrottleMean', ...
            'CL','CD','L_D','Drag_N','T_trim_N','T_available_N', ...
            'T_excess_N','P_required_kW','P_available_kW', ...
            'P_excess_kW','Ps_mps','ROC_fpm','FuelFlow_kgps', ...
            'FuelFlow_kgph','FuelRange_nmpkg', ...
            'FuelEndurance_hrpkg','ResidualNorm'});

        [~,nearest] = min(abs(perf_table.Mach-M_seed));
        source_index = perf_table.SourceIndex(nearest);
        seed_z_out = z_solution(source_index,:).';
    end

    sweep = struct();
    sweep.altitude_m = alt_m;
    sweep.M_grid = M_grid;
    sweep.V_grid_mps = V_grid;
    sweep.debug_table = debug_table;
    sweep.perf_table = perf_table;
    sweep.seed_z_out = seed_z_out;
end

function r = get_residual_norm(info)
    r = NaN;

    if isstruct(info)
        if isfield(info,"residual_norm")
            r = info.residual_norm;
        elseif isfield(info,"final_residual_norm")
            r = info.final_residual_norm;
        elseif isfield(info,"residual")
            r = norm(info.residual(:));
        elseif isfield(info,"objective_value")
            r = sqrt(max(info.objective_value,0));
        end
    end

    if isempty(r) || ~isfinite(r)
        r = NaN;
    end
end

