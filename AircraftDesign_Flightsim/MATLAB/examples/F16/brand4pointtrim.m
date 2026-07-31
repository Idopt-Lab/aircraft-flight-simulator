%% F-16 TRIM-CONSTRAINED PERFORMANCE ANALYSIS
% Aircraft definition below; PerformanceAnalysis returns all trim and
% performance quantities. This driver runs: single-point performance
% metrics (stall, best range/endurance, best glide, best angle/rate
% climb, max level speed, sustained/coordinated turn), a baseline trim
% with stability analysis, takeoff and landing, and a velocity sweep
% (dry and afterburner) that feeds the aerodynamics/thrust-power plots.

clear;
clc;
close all;

%% PATH BOOTSTRAP
% Ensures the shared class library and example helper functions are on the
% MATLAB path regardless of the current working directory or session state.
this_script_dir = fileparts(mfilename('fullpath'));
matlab_root_dir = fileparts(fileparts(this_script_dir));
addpath(genpath(fullfile(matlab_root_dir,'classes')));
addpath(genpath(fullfile(matlab_root_dir,'examples')));

%% AIRCRAFT DEFINITION
% Geometry, mass properties, frames, controls, propulsion and load sources.

%% Reference geometry and mass properties

ft_to_m = 0.3048;
ft2_to_m2 = 0.09290304;
lbf_to_N = 4.4482216152605;
g = 9.80665;

aircraft = Aircraft();
aircraft.geometry.set_reference_geometry(300*ft2_to_m2,30*ft_to_m,11.3201786951274*ft_to_m);
aircraft.geometry.set_reference_point([0 0 0]);
aircraft.geometry.aero_mach_max = 2;

analysis_weight_lbf = 0.899666962714079*31377;
empty_weight_lbf = 19980.7005781593;
permanent_payload_lbf = 700;
expendable_payload_lbf = 4400;
engine_weight_lbf = 0.199*23770;
fuel_weight_lbf = max(analysis_weight_lbf-empty_weight_lbf- permanent_payload_lbf-expendable_payload_lbf,0);
airframe_weight_lbf = empty_weight_lbf-engine_weight_lbf;

airframe_mass_kg = airframe_weight_lbf*lbf_to_N/g;
fuel_mass_kg = fuel_weight_lbf*lbf_to_N/g;
engine_mass_kg = engine_weight_lbf*lbf_to_N/g;
permanent_payload_mass_kg = permanent_payload_lbf*lbf_to_N/g;
expendable_payload_mass_kg = expendable_payload_lbf*lbf_to_N/g;

inertia_scale = analysis_weight_lbf/20500;
inertia_correction = 1.3558179483314;
airframe_inertia_kgm2 = inertia_scale*[9496*inertia_correction,0,982*inertia_correction; 0,55814*inertia_correction,0; 982*inertia_correction,0,63100*inertia_correction];

%% Reference frames and component hierarchy

aircraft.add_frame("fuel_tank","body",[0;0;0],@(x) eye(3));
aircraft.add_frame("engine","body",[0;0;0],@(x) eye(3));
aircraft.add_frame("cg","body",[0;0;0],@(x) eye(3));
aircraft.add_frame("gravity_cg","body",[0;0;0], @(x) ReferenceFrame.ned_to_body_dcm(x));

airframe = Component("airframe",airframe_mass_kg,[0;0;0], airframe_inertia_kgm2,aircraft.get_frame("body"),"airframe");
fuel = Component("fuel",fuel_mass_kg,[0;0;0],zeros(3), aircraft.get_frame("fuel_tank"),"fuel");
engine_component = Component("engine",engine_mass_kg,[0;0;0],zeros(3), aircraft.get_frame("engine"),"engine");
permanent_payload = Component("permanent_payload", permanent_payload_mass_kg,[0;0;0],zeros(3), aircraft.get_frame("body"),"payload");
expendable_payload = Component("expendable_payload", expendable_payload_mass_kg,[0;0;0],zeros(3), aircraft.get_frame("body"),"payload");

aircraft.add_component(airframe);
aircraft.add_component(fuel);
aircraft.add_component(engine_component);
aircraft.add_component(permanent_payload);
aircraft.add_component(expendable_payload);

%% Control surfaces

aircraft.add_control_surface(ControlSurface("aileron","aileron","primary",[1 0 0], deg2rad(21.5),deg2rad(-21.5),0,0,0));
aircraft.add_control_surface(ControlSurface("stabilator","elevator","primary",[0 1 0], deg2rad(25),deg2rad(-25),0,0,0));
aircraft.add_control_surface(ControlSurface("rudder","rudder","primary",[0 0 1], deg2rad(30),deg2rad(-30),0,0,0));

%% Propulsion model

engine = BrandtAfterburningEngine("F16A_engine",aircraft.get_frame("engine"),[1;0;0]);
engine.set_rating_mode("dry");
aircraft.add_propulsive_element(engine);
engine_component.add_load_source(PropulsionLoadSolver(engine,aircraft.get_frame("engine")));

%% Aerodynamic and gravity load sources

aerodynamics = CoefficientAerodynamics(@F16ABrandtLookup);
aircraft.add_load_source(AeroLoadSolver(aerodynamics,aircraft.geometry,aircraft,aircraft.get_frame("body")));
aircraft.add_load_source(GravityLoadSolver(aircraft,aircraft.get_frame("gravity_cg")));

%% Final mass-property reference

[~,cg_body_m,~] = aircraft.compute_total_mass_properties();
aircraft.update_frame_position("cg",cg_body_m);
aircraft.update_frame_position("gravity_cg",cg_body_m);
aircraft.set_reference_frame("cg");

%% PERFORMANCE ANALYSIS CONFIGURATION

setup = F16PerformanceConfig(aircraft);
performance = aircraft.get_performance();

condition = setup.condition;
cfg = setup.performance;
solver = setup.solver;
fast_solver = setup.fast_solver;
bounds = setup.bounds;

%% DRY-POWER PERFORMANCE

engine.set_rating_mode("dry");

results.dry.stall = performance.optimize_stall_speed(condition,solver,cfg,bounds.stall);
results.dry.best_range = performance.optimize_paper_best_range(condition,solver,cfg,bounds.range);
results.dry.best_endurance = performance.optimize_paper_best_endurance(condition,solver,cfg,bounds.endurance);
results.glide = performance.analyze_glide_performance(condition,solver,cfg,bounds.glide);
results.dry.best_angle_climb = performance.optimize_best_angle_climb(condition,solver,cfg,bounds.climb);
results.dry.best_rate_climb = performance.optimize_best_rate_climb(condition,solver,cfg,bounds.climb);
results.dry.velocity_sweep = performance.evaluate_velocity_sweep(condition,setup.velocity_grid_mps,setup.level_spec,fast_solver,cfg);

%% BASELINE TRIM AND STABILITY
% Dry rating, matching a subsonic cruise condition (no afterburner needed
% at Ma 0.6). Same Brandt-lookup / BrandtAfterburningEngine aircraft
% object used for performance above; no separate F16Lookup/
% TurbofanPropulsion model.

[~,a_baseline,~,~] = aircraft.get_atmosphere(setup.baseline_altitude_m);
V_baseline_mps = setup.baseline_mach*a_baseline;
baseline_condition = struct('altitude_m',setup.baseline_altitude_m,'velocity_mps',V_baseline_mps);

fprintf("\n=== F-16 BASELINE CRUISE TRIM ===\n");
fprintf("Altitude : %.1f m\n",setup.baseline_altitude_m);
fprintf("Mach     : %.3f\n",setup.baseline_mach);
fprintf("V        : %.3f m/s\n",V_baseline_mps);

[baseline_x,baseline_u,baseline_converged,baseline_info] = performance.solve_trim(baseline_condition,setup.level_spec,solver);
baseline = performance.evaluate_trim_point(baseline_x,baseline_u,baseline_condition,cfg,baseline_info.extras);
baseline.trim_converged = baseline_converged;
baseline.trim_info = baseline_info;

performance.print_point(baseline);

fprintf("\n=== F-16 STABILITY ANALYSIS ===\n");

stab = aircraft.get_stability();
stab.set_trim(baseline_x,baseline_u);

[~,~,Cm_alpha] = stab.compute_pitch_static_stability(deg2rad(-5:0.5:15));
[~,~,Cn_beta] = stab.compute_yaw_static_stability(deg2rad(-10:0.5:10));

stab.linearize(1e-5,1e-5);
stab.analyze_modes();

try
    stab.print_modes();
catch
end

fprintf("\n=== F-16 STATIC STABILITY ===\n");
fprintf("Cm_alpha : %.6f 1/rad\n",Cm_alpha);
fprintf("Cn_beta  : %.6f 1/rad\n",Cn_beta);

results.baseline = baseline;
results.stability = struct('stab',stab,'Cm_alpha',Cm_alpha,'Cn_beta',Cn_beta);

%% AFTERBURNER PERFORMANCE

engine.set_rating_mode("afterburner");

results.afterburner.best_angle_climb = performance.optimize_best_angle_climb(condition,solver,cfg,bounds.climb);
results.afterburner.best_rate_climb = performance.optimize_best_rate_climb(condition,solver,cfg,bounds.climb);
results.afterburner.maximum_level_speed = performance.optimize_max_level_speed(condition,solver,cfg,bounds.maximum_speed);
results.afterburner.velocity_sweep = performance.evaluate_velocity_sweep(condition,setup.velocity_grid_mps,setup.level_spec,fast_solver,cfg);
results.turn.sustained_optimum = performance.optimize_sustained_turn(condition,solver,cfg,bounds.turn);
results.turn.coordinated = performance.optimize_coordinated_turn(setup.turn_condition,setup.coordinated_turn_bank_deg, solver,cfg,bounds.turn);

%% TAKEOFF AND LANDING
% Same aircraft object; takeoff uses dry rating, landing uses idle.

fprintf("\n=== F-16 TAKEOFF ANALYSIS ===\n");

engine.set_rating_mode("dry");

to = aircraft.get_takeoff();
[TO_m,to_res] = to.calculate_takeoff(0,0,3000); %#ok<NASGU>
to.print_takeoff_summary(to_res);

fprintf("\n=== F-16 LANDING ANALYSIS ===\n");

engine.set_rating_mode("idle");

ld = aircraft.get_landing();
[LD_m,ld_res] = ld.calculate_landing(0,0,3000); %#ok<NASGU>
ld.print_landing_summary(ld_res);

engine.set_rating_mode("afterburner");

results.takeoff = to_res;
results.landing = ld_res;

%% COLLECT COMPLETE FRAMEWORK OUTPUTS

results.aircraft = aircraft;
results.engine = engine;
results.setup = setup;

%% SUMMARY AND PLOTS

printF16PerformanceSummary(results);
plotF16VelocitySweep(results);

F16Performance = results;

%% LOCAL PERFORMANCE-CONFIGURATION FUNCTION
function setup = F16PerformanceConfig(aircraft)

setup.condition = struct('altitude_m',9144);
setup.turn_condition = struct('altitude_m',9144,'velocity_mps',250);
setup.coordinated_turn_bank_deg = 60;

setup.baseline_altitude_m = 9144;
setup.baseline_mach = 0.600;

setup.velocity_grid_mps = (140:10:550).';

setup.performance = struct( ...
    'available_throttle',1, ...
    'continuous_throttle',1, ...
    'takeoff_throttle',1, ...
    'CL_max',0.9840156989811455, ...
    'CL_min',-0.8, ...
    'max_Mach',2, ...
    'propulsion_type',"thrust", ...
    'range_endurance_metric',"actual_flow", ...
    'fuel_LHV_Jpkg',43e6, ...
    'n_limit_pos',9, ...
    'n_limit_neg',-3, ...
    'max_dynamic_pressure_Pa',2133*47.88025898033584, ...
    'service_ceiling_threshold_fpm',100, ...
    'enforce_propulsion_constraints',true);

setup.solver = struct( ...
    'residual_tolerance',1e-5, ...
    'inequality_tolerance',1e-6, ...
    'equality_tolerance',1e-5, ...
    'optimality_tolerance',1e-5, ...
    'fmincon_options',optimoptions("fmincon", ...
        "Algorithm","sqp", ...
        "Display","none", ...
        "MaxIterations",1200, ...
        "MaxFunctionEvaluations",20000, ...
        "OptimalityTolerance",1e-8, ...
        "ConstraintTolerance",1e-8, ...
        "StepTolerance",1e-10, ...
        "FiniteDifferenceType","central", ...
        "FiniteDifferenceStepSize",1e-4));

setup.fast_solver = struct( ...
    'residual_tolerance',1e-5, ...
    'inequality_tolerance',1e-6, ...
    'equality_tolerance',1e-5, ...
    'optimality_tolerance',1e-5, ...
    'fmincon_options',optimoptions("fmincon", ...
        "Algorithm","sqp", ...
        "Display","none", ...
        "MaxIterations",350, ...
        "MaxFunctionEvaluations",6000, ...
        "OptimalityTolerance",1e-6, ...
        "ConstraintTolerance",1e-6, ...
        "StepTolerance",1e-8, ...
        "FiniteDifferenceType","forward", ...
        "FiniteDifferenceStepSize",1e-4));

setup.level_spec = struct( ...
    'mode',"level", ...
    'variables',["alpha";"control_pitch";"throttle"], ...
    'initial_guess',[deg2rad(8);deg2rad(1.5);0.4], ...
    'lb',[deg2rad(-5);deg2rad(-25);0], ...
    'ub',[deg2rad(20);deg2rad(25);1], ...
    'fixed',struct('beta',0,'phi',0,'psi',0,'gamma',0, ...
        'control_roll',0,'control_yaw',0), ...
    'reference_frame_name',"cg", ...
    'residual_scale',[aircraft.g;aircraft.g;1]);

%% Optimization bounds

setup.bounds.stall = struct( ...
    'V_lb',40,'V_ub',250,'V0',120, ...
    'alpha0_deg',15,'alpha_min_deg',-5,'alpha_max_deg',20, ...
    'elevator0_deg',1,'elevator_min_deg',-25,'elevator_max_deg',25, ...
    'throttle0',0.4,'throttle_min',0,'throttle_max',1);

setup.bounds.climb = struct( ...
    'V_lb',60,'V_ub',600,'V0',350, ...
    'alpha0_deg',1.5,'alpha_min_deg',-5,'alpha_max_deg',20, ...
    'elevator0_deg',1.7,'elevator_min_deg',-25,'elevator_max_deg',25, ...
    'gamma0_deg',12,'gamma_min_deg',0,'gamma_max_deg',45, ...
    'throttle0',1);

setup.bounds.range = struct( ...
    'V_lb',60,'V_ub',400,'V0',230, ...
    'alpha0_deg',5,'alpha_min_deg',-5,'alpha_max_deg',20, ...
    'elevator0_deg',1,'elevator_min_deg',-25,'elevator_max_deg',25, ...
    'throttle0',0.4,'throttle_min',0,'throttle_max',1);

setup.bounds.endurance = struct( ...
    'V_lb',60,'V_ub',400,'V0',190, ...
    'alpha0_deg',6,'alpha_min_deg',-5,'alpha_max_deg',20, ...
    'elevator0_deg',1,'elevator_min_deg',-25,'elevator_max_deg',25, ...
    'throttle0',0.35,'throttle_min',0,'throttle_max',1);

setup.bounds.glide = struct( ...
    'V_lb',50,'V_ub',350,'V0',190, ...
    'alpha0_deg',6,'alpha_min_deg',-5,'alpha_max_deg',20, ...
    'elevator0_deg',1,'elevator_min_deg',-25,'elevator_max_deg',25, ...
    'gamma0_deg',-5,'gamma_min_deg',-30,'gamma_max_deg',-0.01);

setup.bounds.maximum_speed = struct( ...
    'V_lb',100,'V_ub',600,'V0',450, ...
    'alpha0_deg',1,'alpha_min_deg',-5,'alpha_max_deg',20, ...
    'elevator0_deg',1,'elevator_min_deg',-25,'elevator_max_deg',25, ...
    'throttle0',1,'throttle_min',0,'throttle_max',1);

setup.bounds.turn = struct( ...
    'V_lb',100,'V_ub',350,'V0',240, ...
    'bank0_deg',45,'bank_min_deg',1,'bank_max_deg',80, ...
    'throttle0',1,'throttle_min',0,'throttle_max',1, ...
    'turn_rate_ub_radps',1.2, ...
    'alpha0_deg',8,'alpha_min_deg',-5,'alpha_max_deg',20, ...
    'elevator0_deg',0,'elevator_min_deg',-25,'elevator_max_deg',25, ...
    'aileron0_deg',0,'aileron_min_deg',-21.5,'aileron_max_deg',21.5, ...
    'rudder0_deg',0,'rudder_min_deg',-30,'rudder_max_deg',30, ...
    'Ps_trim_tolerance_mps',0.05);

end

%% LOCAL PLOTTING FUNCTION
function plotF16VelocitySweep(results)

dry = sweepSeries(results.dry.velocity_sweep);
afterburner = sweepSeries(results.afterburner.velocity_sweep);

%% Aerodynamic plots

figure("Name","Aerodynamics");
tiledlayout(2,2,"TileSpacing","compact");
nexttile;
plot(dry.V_mps,dry.L_N/1000,"LineWidth",1.5);
hold on;
yline(dry.weight_N(1)/1000,"--");
grid on;
xlabel("V [m/s]");
ylabel("Force [kN]");
legend("Lift","Weight","Location","best");
nexttile;
plot(dry.V_mps,dry.D_N/1000,"LineWidth",1.5);
grid on;
xlabel("V [m/s]");
ylabel("Drag [kN]");
nexttile;
plot(dry.V_mps,dry.CL,"LineWidth",1.5);
hold on;
plot(dry.V_mps,dry.CD,"LineWidth",1.5);
grid on;
xlabel("V [m/s]");
ylabel("Coefficient");
legend("C_L","C_D","Location","best");
nexttile;
plot(dry.V_mps,dry.L_over_D,"LineWidth",1.5);
grid on;
xlabel("V [m/s]");
ylabel("L/D");

%% Propulsion and power plots

figure("Name","Thrust and Power");
tiledlayout(2,2,"TileSpacing","compact");
nexttile;
plot(dry.V_mps,dry.T_required_N/1000,"LineWidth",1.5);
hold on;
plot(dry.V_mps,dry.T_available_N/1000,"LineWidth",1.5);
plot(afterburner.V_mps,afterburner.T_available_N/1000,"LineWidth",1.5);
grid on;
xlabel("V [m/s]");
ylabel("Thrust [kN]");
legend("Required","Dry available","AB available","Location","best");
nexttile;
plot(dry.V_mps,dry.T_excess_N/1000,"LineWidth",1.5);
hold on;
plot(afterburner.V_mps,afterburner.T_excess_N/1000,"LineWidth",1.5);
yline(0,"--");
grid on;
xlabel("V [m/s]");
ylabel("Excess thrust [kN]");
legend("Dry","Afterburner","Location","best");
nexttile;
plot(dry.V_mps,dry.P_thrust_required_W/1e6,"LineWidth",1.5);
hold on;
plot(dry.V_mps,dry.P_thrust_available_W/1e6,"LineWidth",1.5);
plot(afterburner.V_mps,afterburner.P_thrust_available_W/1e6, "LineWidth",1.5);
grid on;
xlabel("V [m/s]");
ylabel("Power [MW]");
legend("Required","Dry available","AB available","Location","best");
nexttile;
plot(dry.V_mps,dry.P_excess_W/1e6,"LineWidth",1.5);
hold on;
plot(afterburner.V_mps,afterburner.P_excess_W/1e6,"LineWidth",1.5);
yline(0,"--");
grid on;
xlabel("V [m/s]");
ylabel("Excess power [MW]");
legend("Dry","Afterburner","Location","best");

end

function data = sweepSeries(sweep)

indices = find(sweep.valid);
if isempty(indices)
    error('plotF16VelocitySweep:NoValidSweepPoints', 'The velocity sweep contains no valid performance points.');
end
points = [sweep.point{indices}];

fields = [ ...
    "V_mps";"weight_N";"L_N";"D_N";"CL";"CD";"L_over_D"; ...
    "alpha_deg";"theta_deg";"elevator_deg";"throttle"; ...
    "fuel_flow_trim";"T_required_N"; ...
    "T_available_N";"T_excess_N";"P_thrust_required_W"; ...
    "P_thrust_available_W";"P_excess_W";"Ps_mps"];

for k = 1:numel(fields)
    name = fields(k);
    data.(name) = reshape([points.(name)],[],1);
end
data.residual_inf = sweep.residual_norm(indices);

end

%% LOCAL RESULTS-SUMMARY FUNCTION
function printF16PerformanceSummary(results)

rows = [ ...
    pointTableRow("Stall",results.dry.stall); ...
    pointTableRow("Best range",results.dry.best_range); ...
    pointTableRow("Best endurance",results.dry.best_endurance); ...
    pointTableRow("Best glide",results.glide.best_glide); ...
    pointTableRow("Minimum sink",results.glide.min_sink); ...
    pointTableRow("Best angle climb (dry)", ...
        results.dry.best_angle_climb); ...
    pointTableRow("Best rate climb (dry)", ...
        results.dry.best_rate_climb); ...
    pointTableRow("Best angle climb (afterburner)", ...
        results.afterburner.best_angle_climb); ...
    pointTableRow("Best rate climb (afterburner)", ...
        results.afterburner.best_rate_climb); ...
    pointTableRow("Maximum level speed", ...
        results.afterburner.maximum_level_speed); ...
    pointTableRow("Maximum sustained turn", ...
        results.turn.sustained_optimum); ...
    pointTableRow("Coordinated turn",results.turn.coordinated)];

fprintf('\n=== F-16 PERFORMANCE POINTS AND TRIM CONTROLS ===\n');
disp(rows);

h = results.setup.condition.altitude_m;
print_one_sweep("DRY",results.dry.velocity_sweep,h);
print_one_sweep("AFTERBURNER",results.afterburner.velocity_sweep,h);

end

function print_one_sweep(label,sweep,altitude_m)

fprintf('\n=== F-16 VELOCITY SWEEP (%s), h = %.0f m ===\n',label,altitude_m);

for i = 1:numel(sweep.velocity_mps)
    if sweep.valid(i)
        p = sweep.point{i};
        fprintf('V=%6.1f m/s | Mach=%.3f | alpha=%6.2f deg | throttle=%.4f | CL=%.4f | CD=%.4f | L/D=%6.2f | T_excess=%8.1f N | Ps=%6.2f m/s | residual=%.2e\n', ...
            sweep.velocity_mps(i),p.Mach,p.alpha_deg,p.throttle,p.CL,p.CD,p.L_over_D,p.T_excess_N,p.Ps_mps,sweep.residual_norm(i));
    else
        fprintf('V=%6.1f m/s | INVALID | %s\n',sweep.velocity_mps(i),sweep.error_message(i));
    end
end

end

function row = pointTableRow(metric,opt)

p = opt.point_opt;
row = table( ...
    string(metric),logical(opt.converged),p.V_mps,p.Mach,p.alpha_deg, ...
    p.beta_deg,p.gamma_trim_deg,p.phi_deg,p.theta_deg, ...
    p.aileron_deg,p.elevator_deg,p.rudder_deg,p.throttle, ...
    p.ROC_trim_mps,p.fuel_flow_trim,{p.u_trim(:).'}, ...
    'VariableNames',{'Metric','Converged','V_mps','Mach','Alpha_deg', 'Beta_deg','Gamma_deg','Phi_deg','Theta_deg','Aileron_deg', 'Stabilator_deg','Rudder_deg','Throttle','ROC_mps', 'FuelFlow_kgps','ControlVector'});

end
