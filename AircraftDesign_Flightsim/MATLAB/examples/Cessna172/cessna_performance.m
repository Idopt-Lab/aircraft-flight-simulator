
clear
clc
close all

%% PATH BOOTSTRAP
% Ensures the shared class library and example helper functions are on the
% MATLAB path regardless of the current working directory or session state.
this_script_dir = fileparts(mfilename('fullpath'));
matlab_root_dir = fileparts(fileparts(this_script_dir));
addpath(genpath(fullfile(matlab_root_dir,'classes')));
addpath(genpath(fullfile(matlab_root_dir,'examples')));

datcom_dir = "D:\aircraft-flight-simulator - Copy\AircraftDesign_Flightsim\MATLAB\examples\Cessna172";
datcom_files = fullfile(datcom_dir,["cessna010.out"; "cessna012.out"; "cessna014.out"; "cessna016.out"; "cessna018.out"; "cessna020.out"; "cessna022.out"; "cessna024.out"]);

target_gross_lb = 2550;

performance_altitude_ft = 5000;
baseline_altitude_ft = 6000;
baseline_speed_ktas = 112;
baseline_reference_rpm = 2500;
baseline_reference_bhp_fraction = 0.65;
baseline_reference_fuel_gph = 8.8;

baseline_CD_increment = 0.01681;
propeller_coefficient_scale = 0.84091;

max_engine_power_hp = 180;
avgas_density_kg_per_usgal = 2.72;

moment_reference_frame_name = "cg";


ac = Aircraft();
ac.set_body_frame("body");

S_ref = 16.17;
b_ref = 10.98;
c_ref = 1.50;

ac.geometry.set_reference_geometry(S_ref,b_ref,c_ref);

NOSE_REF_CG = [1.95;0.00;1.26];
airframe_nose_ref = [1.95;0.00;1.26];
payload_nose_ref = [2.10;0.00;1.20];
datcom_moment_nose_ref = [2.11;0.00;1.26];
engine_nose_ref = [1.68;0.00;1.26];

cg_to_airframe = airframe_nose_ref-NOSE_REF_CG;
cg_to_payload = payload_nose_ref-NOSE_REF_CG;
cg_to_aero_ref = datcom_moment_nose_ref-NOSE_REF_CG;
cg_to_engine = engine_nose_ref-NOSE_REF_CG;

ac.geometry.set_reference_point(cg_to_aero_ref.');

add_or_update_frame(ac,"cg","body",[0;0;0],@(x) eye(3));

cfg = ac.get_configurator();


target_gross_kg = target_gross_lb*0.45359237;
empty_mass_kg = 767;
fuel_mass_kg = 144;
payload_kg = target_gross_kg-empty_mass_kg-fuel_mass_kg;

if payload_kg < 0
    error("C172:InvalidMass", "Target gross mass is lower than empty plus fuel mass.");
end

I_body = [1285 0 0;0 1824 0;0 0 2666];

cfg.add_component('name','empty_airframe', 'type','airframe', 'mass',empty_mass_kg, 'position',cg_to_airframe, 'cg_local',[0;0;0], 'inertia',I_body, 'parent_frame','body', 'frame_name','airframe_frame');

cfg.add_component('name','payload', 'type','payload', 'mass',payload_kg, 'position',cg_to_payload, 'cg_local',[0;0;0], 'inertia',zeros(3,3), 'parent_frame','body', 'frame_name','payload_frame');

cfg.add_component('name','fuel_fixed', 'type','fuel', 'mass',fuel_mass_kg, 'position',[0;0;0], 'cg_local',[0;0;0], 'inertia',zeros(3,3), 'parent_frame','body', 'frame_name','fuel_frame');


c172_model_options = struct('CD0_add',baseline_CD_increment, 'CD_alpha_extra',0, 'CL_max_abs',1.9);
[c172_lookup,c172_data] = c172datcom(datcom_files,false,c172_model_options);
expected_mach_grid = (0.10:0.02:0.24).';
if numel(c172_data.mach) ~= numel(expected_mach_grid) || max(abs(c172_data.mach(:)-expected_mach_grid)) > 1e-12
    error("C172:UnexpectedMachGrid", "Expected DATCOM Mach cases 0.10:0.02:0.24.");
end
aero_mach_limit = max(c172_data.mach);

datcom_reference = [c172_data.reference_area_m2, c172_data.reference_chord_m, c172_data.reference_span_m, c172_data.moment_reference_horizontal_m, c172_data.moment_reference_vertical_m];
expected_reference = [S_ref,c_ref,b_ref, datcom_moment_nose_ref(1),datcom_moment_nose_ref(3)];
if max(abs(datcom_reference-expected_reference)) > 1e-6
    error("C172:DATCOMReferenceMismatch", "Aircraft geometry does not match the DATCOM reference definition.");
end
aero_model = CoefficientAerodynamics(c172_lookup);

cfg.add_aero_solver(aero_model, 'name','c172_aero', 'position',cg_to_aero_ref, 'dcm_fn',@(x) eye(3), 'geom',ac.geometry, 'parent_frame','body', 'frame_name','aero_ref');


cfg.add_control_surface('name','aileron', 'surface_type','aileron', 'classification','primary', 'axis',[1 0 0], 'max_deflection',deg2rad(20), 'min_deflection',deg2rad(-20), 'dCl',0,'dCm',0,'dCn',0);

cfg.add_control_surface('name','elevator', 'surface_type','elevator', 'classification','primary', 'axis',[0 1 0], 'max_deflection',deg2rad(28), 'min_deflection',deg2rad(-26), 'dCl',0,'dCm',0,'dCn',0);

cfg.add_control_surface('name','rudder', 'surface_type','rudder', 'classification','primary', 'axis',[0 0 1], 'max_deflection',deg2rad(30), 'min_deflection',deg2rad(-30), 'dCl',0,'dCm',0,'dCn',0);


o360 = PistonEngineModel(max_engine_power_hp*745.7,2700,"Lycoming_O_360_A4M");
o360.minimum_rpm = 600;
o360.maximum_rpm = 2700;
o360.density_lapse_exponent = 1.0;
o360.installation_power_factor = 1.0;
o360.set_rpm_power_map([0.00,0.25,0.40,0.55,0.70,0.85,1.00], [0.00,0.18,0.38,0.60,0.78,0.91,1.00]);
o360.set_throttle_map([0.00,0.25,0.50,0.75,1.00], [0.00,0.25,0.50,0.75,1.00]);
o360.set_bsfc_map([0.10,0.40,0.47,0.53,0.59,0.65,0.72,0.76,1.00], [0.42,0.32,0.289,0.279,0.275,0.27435,0.270,0.269,0.290]);
o360.idle_fuel_flow_kgps = 1.2*avgas_density_kg_per_usgal/3600;

sensenich_76em8s14 = PropellerPerformanceModel( ...
    1.9304, ...
    [0.00,0.20,0.40,0.60,0.80,1.00,1.20,1.40], ...
    propeller_coefficient_scale* ...
        [0.115,0.105,0.088,0.065,0.050,0.033,0.015,0.000], ...
    propeller_coefficient_scale* ...
        [0.060,0.059,0.057,0.054,0.051,0.048,0.045,0.043], ...
    2, ...
    "Sensenich_76EM8S14_0_60");
sensenich_76em8s14.minimum_rpm = 600;
sensenich_76em8s14.maximum_rpm = 2700;
sensenich_76em8s14.set_map_metadata("calibrated provisional coefficient shape", "calibrated_provisional",false, "Uniform Ct/Cp scale fitted to the 6000-ft, 112-KTAS cruise point; not manufacturer-sourced.");

propulsion_arguments = { ...
    'name','O360', ...
    'element_type','propeller', ...
    'max_output',max_engine_power_hp*745.7, ...
    'position',cg_to_engine, ...
    'mount_euler',[0 0 0], ...
    'direction',[1 0 0], ...
    'fuel_rate',0, ...
    'mass',0, ...
    'parent_frame','body', ...
    'frame_name','engine'};

configurator_path = which('AircraftConfigurator');
if isempty(configurator_path)
    error("C172:MissingAircraftConfigurator", "AircraftConfigurator is not on the MATLAB path.");
end
configurator_supports_models = contains(fileread(configurator_path),"engine_model");

if configurator_supports_models
    [engine_pe,~] = cfg.add_propulsive_element(propulsion_arguments{:}, 'engine_model',o360, 'propeller_model',sensenich_76em8s14, 'rotation_direction',1);
else
    [engine_pe,~] = cfg.add_propulsive_element(propulsion_arguments{:});
    engine_pe.set_models(o360,sensenich_76em8s14);
end
engine_pe.set_rotation_direction(1);
engine_pe.set_air_velocity_frame(ac.get_frame(ac.body_frame_name));
engine_pe.rpm_grid_size = 160;


[m_total,cg_total,I_total] = ac.compute_total_mass_properties();
ac.update_frame_position("cg",cg_total);
ac.ensure_gravity_source();
ac.set_reference_frame(moment_reference_frame_name);

W = m_total*ac.g;
perf = PerformanceAnalysis(ac);

fprintf("\n=== C172 READY ===\n");
fprintf("Gross mass       : %.3f kg  %.1f lb\n",m_total,m_total/0.45359237);
fprintf("Weight           : %.3f N\n",W);
fprintf("CG               : [% .4f % .4f % .4f] m\n",cg_total);
fprintf("Sref/bref/cref   : %.3f / %.3f / %.3f m\n",S_ref,b_ref,c_ref);
fprintf("Alpha data range : %.3f to %.3f deg\n", rad2deg(min(c172_data.alpha)),rad2deg(max(c172_data.alpha)));
fprintf("CL table range   : %.4f to %.4f\n", min(c172_data.CL_grid(:),[],'omitnan'), max(c172_data.CL_grid(:),[],'omitnan'));
fprintf("Mach data range  : %.3f to %.3f\n", min(c172_data.mach),max(c172_data.mach));
fprintf("DATCOM grid      : %d Mach x %d altitude x %d alpha\n", c172_data.database_size);
fprintf("Alpha intersect. : %d\n",c172_data.alpha_intersection_applied);
fprintf("Ragged alpha grid: %d\n",c172_data.ragged_alpha_grid);
fprintf("Cruise CD add    : %.5f\n",baseline_CD_increment);
fprintf("Prop. map scale  : %.5f\n",propeller_coefficient_scale);
fprintf("Prop. map status : %s | manufacturer validated: %d\n", sensenich_76em8s14.map_status, sensenich_76em8s14.manufacturer_validated);
fprintf("Aero moment ref. : [% .4f % .4f % .4f] m body\n",cg_to_aero_ref);
fprintf("Controls/engines : %d / %d\n", numel(ac.control_surfaces),numel(ac.propulsive_elements));


solver_cfg = struct();
solver_cfg.residual_tolerance = 1e-5;
solver_cfg.inequality_tolerance = 1e-6;
solver_cfg.equality_tolerance = 1e-5;
solver_cfg.fmincon_options = optimoptions( ...
    "fmincon", ...
    "Algorithm","sqp", ...
    "Display","none", ...
    "MaxIterations",2500, ...
    "MaxFunctionEvaluations",35000, ...
    "OptimalityTolerance",1e-9, ...
    "StepTolerance",1e-10, ...
    "ConstraintTolerance",1e-8, ...
    "FunctionTolerance",1e-10, ...
    "FiniteDifferenceType","central", ...
    "FiniteDifferenceStepSize",1e-4);

alpha_data_min_deg = rad2deg(min(c172_data.alpha));
alpha_data_max_deg = rad2deg(max(c172_data.alpha));
alpha_min_deg = max(-5,alpha_data_min_deg);
alpha_max_deg = min(18,alpha_data_max_deg);

if alpha_max_deg <= alpha_min_deg
    error("C172:InvalidAlphaDatabase", "DATCOM alpha range is insufficient for trim optimization.");
end

perf_cfg = struct();
perf_cfg.available_throttle = 1;
perf_cfg.takeoff_throttle = 1;
perf_cfg.continuous_throttle = 1;
perf_cfg.stall_throttle = 0.5;
perf_cfg.CL_max = max(c172_data.CL_grid(:),[],'omitnan');
perf_cfg.CL_min = min(c172_data.CL_grid(:),[],'omitnan');
perf_cfg.propulsion_type = "power";
perf_cfg.range_endurance_metric = "classical";
perf_cfg.max_Mach = aero_mach_limit;
perf_cfg.reference_frame_name = moment_reference_frame_name;
perf_cfg.enforce_propulsion_constraints = true;


baseline_condition = struct('altitude_m',baseline_altitude_ft*0.3048, 'velocity_mps',baseline_speed_ktas*0.514444);

baseline_spec = make_level_spec([deg2rad(3);deg2rad(0);0.6], alpha_min_deg,alpha_max_deg,-26,28,ac.g, moment_reference_frame_name);

[x_base,u_base,conv_base,info_base] = perf.solve_trim(baseline_condition,baseline_spec,solver_cfg);
baseline_point = perf.evaluate_trim_point(x_base,u_base,baseline_condition,perf_cfg,info_base.extras);
baseline_point.objective = "baseline_cruise_trim";
baseline_point.trim_converged = conv_base;
baseline_point.valid = conv_base && point_in_c172_database(baseline_point,alpha_data_min_deg,alpha_data_max_deg,aero_mach_limit);
baseline_point.residual = info_base.residual;
baseline_point.residual_norm = info_base.residual_norm;
baseline_point.residual_inf = info_base.residual_inf;

perf.print_point(baseline_point);

baseline_engine = evaluate_engine_point(engine_pe,baseline_point);
baseline_model_bhp_fraction = baseline_engine.brake_power_W/ max(o360.rated_power_W,1e-9);
baseline_model_fuel_gph = baseline_point.fuel_flow_trim*3600/ avgas_density_kg_per_usgal;

validation_quantity = ["KTAS";"RPM";"BHP_percent";"Fuel_GPH"];
sheet_target = [baseline_speed_ktas; baseline_reference_rpm; 100*baseline_reference_bhp_fraction; baseline_reference_fuel_gph];
model_value = [baseline_point.V_kts; baseline_engine.rpm; 100*baseline_model_bhp_fraction; baseline_model_fuel_gph];
model_error = model_value-sheet_target;

baseline_sheet_validation = table(validation_quantity,sheet_target,model_value,model_error, 'VariableNames',{'Quantity','SheetTarget','ModelValue','ModelError'});

fprintf("\n=== 6000 FT STANDARD-TEMPERATURE CRUISE CHECK ===\n");
fprintf("Condition: 2550 lb, 112 KTAS, C172 180 hp performance sheet.\n");
disp(baseline_sheet_validation);
fprintf("Propeller efficiency: %.6f\n", baseline_engine.propeller_efficiency);


required_methods = ["paper_formulation_catalog"; "analyze_paper_trim_suite"; "optimize_paper_stall_speed"; "optimize_paper_best_range"; "optimize_paper_best_endurance"; "compute_paper_required_available_curve"];

method_list = string(methods("PerformanceAnalysis"));
missing_methods = required_methods(~ismember(required_methods,method_list));

if ~isempty(missing_methods)
    error("C172:OldPerformanceAnalysis", "PerformanceAnalysis is missing: %s",strjoin(missing_methods,", "));
end

performance_condition = struct();
performance_condition.altitude_m = performance_altitude_ft*0.3048;

stall_reference = perf.compute_stall_speed(performance_condition,perf_cfg);

V_min_ktas = 45;
V_seed_ktas = 80;
[~,speed_of_sound_performance,~,~] = ac.get_atmosphere(performance_condition.altitude_m);
V_datcom_max_mps = aero_mach_limit*speed_of_sound_performance;
V_max_ktas = min(140,V_datcom_max_mps/0.514444);

V_min = V_min_ktas*0.514444;
V_seed = V_seed_ktas*0.514444;
V_max = V_max_ktas*0.514444;

performance_bounds = struct();

performance_bounds.stall = struct( ...
    'V_lb',V_min, ...
    'V_ub',90*0.514444, ...
    'V0',max(55*0.514444,stall_reference.Vs_mps), ...
    'alpha0_deg',min(12,alpha_max_deg-1), ...
    'alpha_min_deg',alpha_min_deg, ...
    'alpha_max_deg',alpha_max_deg, ...
    'elevator0_deg',0, ...
    'elevator_min_deg',-26, ...
    'elevator_max_deg',28, ...
    'throttle0',0.55, ...
    'throttle_min',0, ...
    'throttle_max',1);

performance_bounds.climb = struct( ...
    'V_lb',50*0.514444, ...
    'V_ub',min(120,V_max_ktas)*0.514444, ...
    'V0',75*0.514444, ...
    'alpha0_deg',5, ...
    'alpha_min_deg',alpha_min_deg, ...
    'alpha_max_deg',alpha_max_deg, ...
    'elevator0_deg',0, ...
    'elevator_min_deg',-26, ...
    'elevator_max_deg',28, ...
    'gamma0_deg',5, ...
    'gamma_min_deg',0, ...
    'gamma_max_deg',20, ...
    'throttle0',1);

performance_bounds.range = struct( ...
    'V_lb',50*0.514444, ...
    'V_ub',min(135,V_max_ktas)*0.514444, ...
    'V0',80*0.514444, ...
    'alpha0_deg',4, ...
    'alpha_min_deg',alpha_min_deg, ...
    'alpha_max_deg',alpha_max_deg, ...
    'elevator0_deg',0, ...
    'elevator_min_deg',-26, ...
    'elevator_max_deg',28, ...
    'throttle0',0.55, ...
    'throttle_min',0, ...
    'throttle_max',1);

performance_bounds.endurance = performance_bounds.range;
performance_bounds.endurance.V0 = 65*0.514444;

performance_bounds.glide = struct( ...
    'V_lb',50*0.514444, ...
    'V_ub',min(120,V_max_ktas)*0.514444, ...
    'V0',75*0.514444, ...
    'alpha0_deg',5, ...
    'alpha_min_deg',alpha_min_deg, ...
    'alpha_max_deg',alpha_max_deg, ...
    'elevator0_deg',0, ...
    'elevator_min_deg',-26, ...
    'elevator_max_deg',28, ...
    'gamma0_deg',-5, ...
    'gamma_min_deg',-20, ...
    'gamma_max_deg',-0.01);

performance_bounds.max_speed = struct( ...
    'V_lb',80*0.514444, ...
    'V_ub',V_max, ...
    'V0',min(125,V_max_ktas)*0.514444, ...
    'alpha0_deg',1, ...
    'alpha_min_deg',alpha_min_deg, ...
    'alpha_max_deg',alpha_max_deg, ...
    'elevator0_deg',0, ...
    'elevator_min_deg',-26, ...
    'elevator_max_deg',28, ...
    'throttle0',0.95, ...
    'throttle_min',0, ...
    'throttle_max',1);

performance_bounds.alpha0_deg = 5;
performance_bounds.alpha_min_deg = alpha_min_deg;
performance_bounds.alpha_max_deg = alpha_max_deg;
performance_bounds.elevator0_deg = 0;
performance_bounds.elevator_min_deg = -26;
performance_bounds.elevator_max_deg = 28;
performance_bounds.throttle0 = 0.55;
performance_bounds.throttle_min = 0;
performance_bounds.throttle_max = 1;
performance_bounds.gamma0_deg = 5;
performance_bounds.gamma_min_deg = -10;
performance_bounds.gamma_max_deg = 20;

performance_bounds.V_sweep_mps = [];

performance_suite = perf.analyze_paper_trim_suite(performance_condition,solver_cfg,perf_cfg,performance_bounds);

V_down = (V_seed_ktas:-2:V_min_ktas).'*0.514444;
V_up = (V_seed_ktas:2:V_max_ktas).'*0.514444;

curve_down = perf.compute_paper_required_available_curve(performance_condition,V_down,solver_cfg,perf_cfg,performance_bounds);
curve_up = perf.compute_paper_required_available_curve(performance_condition,V_up,solver_cfg,perf_cfg,performance_bounds);

performance_suite.required_available_curve = merge_bidirectional_curves(curve_down,curve_up);

fprintf("\n=== C172 AIRBORNE PERFORMANCE FORMULATIONS ===\n");
disp(performance_suite.catalog);
fprintf("Excluded deliberately: VR and VMU.\n");


metric_name = ["Minimum powered trimmed speed"; "Best angle of climb"; "Best rate of climb"; "Best range"; "Best endurance"; "Best engine-off glide"; "Maximum level speed"];

symbol = ["VS";"Vx";"Vy";"Vbr";"Vbe";"Vbg";"Vh"];

opt_list = { ...
    performance_suite.stall_speed; ...
    performance_suite.best_angle_climb; ...
    performance_suite.best_rate_climb; ...
    performance_suite.best_range; ...
    performance_suite.best_endurance; ...
    performance_suite.best_glide; ...
    performance_suite.maximum_level_speed};

n_metric = numel(opt_list);
converged = false(n_metric,1);
V_mps = nan(n_metric,1);
V_KTAS = nan(n_metric,1);
Mach = nan(n_metric,1);
alpha_deg = nan(n_metric,1);
gamma_deg = nan(n_metric,1);
elevator_deg = nan(n_metric,1);
throttle = nan(n_metric,1);
CL = nan(n_metric,1);
CD = nan(n_metric,1);
L_over_D = nan(n_metric,1);
CL32_over_CD = nan(n_metric,1);
ROC_fpm = nan(n_metric,1);
shaft_power_hp = nan(n_metric,1);
engine_rpm = nan(n_metric,1);
propeller_efficiency = nan(n_metric,1);
fuel_flow_gph = nan(n_metric,1);
power_balance_error_W = nan(n_metric,1);
propulsion_valid = false(n_metric,1);
active_limit = strings(n_metric,1);

for i = 1:n_metric
    opt_i = opt_list{i};
    converged(i) = isstruct(opt_i) && isfield(opt_i,'converged') && logical(opt_i.converged);

    if ~converged(i) || ~isfield(opt_i,'point_opt') || isempty(opt_i.point_opt)
        active_limit(i) = "failed";
        continue
    end

    p = opt_i.point_opt;
    V_mps(i) = p.V_mps;
    V_KTAS(i) = p.V_kts;
    Mach(i) = p.Mach;
    alpha_deg(i) = p.alpha_deg;
    gamma_deg(i) = p.gamma_deg;
    elevator_deg(i) = p.elevator_deg;
    throttle(i) = p.throttle;
    CL(i) = p.CL;
    CD(i) = p.CD;
    L_over_D(i) = p.L_over_D;
    CL32_over_CD(i) = max(p.CL,0)^(3/2)/max(p.CD,1e-12);
    ROC_fpm(i) = p.ROC_trim_fpm;
    engine_i = evaluate_engine_point(engine_pe,p);
    shaft_power_hp(i) = engine_i.shaft_power_W/745.7;
    engine_rpm(i) = engine_i.rpm;
    propeller_efficiency(i) = engine_i.propeller_efficiency;
    power_balance_error_W(i) = engine_i.power_balance_error_W;
    propulsion_valid(i) = engine_i.operating_valid;
    if isfinite(p.fuel_flow_trim) && p.fuel_flow_trim >= 0
        fuel_flow_gph(i) = p.fuel_flow_trim*3600/avgas_density_kg_per_usgal;
    end
    active_limit(i) = active_constraint_status(opt_i,p,alpha_min_deg,alpha_max_deg);
end

performance_summary = table( ...
    metric_name,symbol,converged,V_mps,V_KTAS,Mach, ...
    alpha_deg,gamma_deg,elevator_deg,throttle,CL,CD,L_over_D, ...
    CL32_over_CD,ROC_fpm,shaft_power_hp,engine_rpm, ...
    propeller_efficiency,fuel_flow_gph,power_balance_error_W, ...
    propulsion_valid,active_limit, ...
    'VariableNames',{ ...
    'Metric','Symbol','Converged','V_mps','V_KTAS','Mach', ...
    'Alpha_deg','Gamma_deg','Elevator_deg','Throttle','CL','CD', ...
    'L_over_D','CL32_over_CD','ROC_fpm','ShaftPower_hp', ...
    'Engine_RPM','PropellerEfficiency', ...
    'FuelFlow_GPH','PowerBalanceError_W','PropulsionValid', ...
    'ActiveLimit'});

fprintf("\n=== C172 PERFORMANCE SUMMARY ===\n");
fprintf("Analysis altitude: %.0f ft | gross weight: %.0f lb\n", performance_altitude_ft,target_gross_lb);
disp(performance_summary);

print_result("VS  minimum powered trimmed speed", performance_suite.stall_speed,engine_pe,avgas_density_kg_per_usgal);
print_result("Vx  best angle of climb", performance_suite.best_angle_climb,engine_pe,avgas_density_kg_per_usgal);
print_result("Vy  best rate of climb", performance_suite.best_rate_climb,engine_pe,avgas_density_kg_per_usgal);
print_result("Vbr best range", performance_suite.best_range,engine_pe,avgas_density_kg_per_usgal);
print_result("Vbe best endurance", performance_suite.best_endurance,engine_pe,avgas_density_kg_per_usgal);
print_result("Vbg best engine-off glide", performance_suite.best_glide,engine_pe,avgas_density_kg_per_usgal);
print_result("Vh  maximum level speed", performance_suite.maximum_level_speed,engine_pe,avgas_density_kg_per_usgal);

if performance_suite.stall_speed.converged && isfield(performance_suite.stall_speed,'CLmax_reference') && isfield(performance_suite.stall_speed.CLmax_reference,'Vs_mps')
    fprintf("\nMinimum powered trim V : %.3f KTAS\n", performance_suite.stall_speed.point_opt.V_kts);
    fprintf("CLmax reference VS     : %.3f KTAS\n", performance_suite.stall_speed.CLmax_reference.Vs_kts);
end


curve = performance_suite.required_available_curve;
req = curve.required;
av = curve.available;

valid_req = req.valid(:) & isfinite(req.power_W(:));
valid_av = av.valid(:) & isfinite(av.power_W(:));
valid_both = valid_req & valid_av;

required_thrust_hp = req.power_W(:)/745.7;
available_thrust_hp = av.power_W(:)/745.7;
excess_thrust_hp = curve.excess_power_W(:)/745.7;

CL_curve = nan(size(curve.V_mps));
CD_curve = nan(size(curve.V_mps));
for i = 1:numel(curve.V_mps)
    if req.valid(i) && ~isempty(req.point{i})
        CL_curve(i) = req.point{i}.CL;
        CD_curve(i) = req.point{i}.CD;
    end
end

LD_curve = CL_curve./max(CD_curve,1e-12);
CL32CD_curve = max(CL_curve,0).^(3/2)./max(CD_curve,1e-12);

figure;
plot(curve.V_kts(valid_req),required_thrust_hp(valid_req), "o-","LineWidth",1.5);
hold on;
plot(curve.V_kts(valid_av),available_thrust_hp(valid_av), "s-","LineWidth",1.5);
grid on;
xlabel("True airspeed (kt)");
ylabel("Thrust power (hp)");
title("C172 required and available power");
legend("Required","Available","Location","best");

figure;
plot(curve.V_kts(valid_req),LD_curve(valid_req), "o-","LineWidth",1.5);
hold on;
plot(curve.V_kts(valid_req),CL32CD_curve(valid_req), "s-","LineWidth",1.5);
grid on;
xlabel("True airspeed (kt)");
ylabel("Aerodynamic performance metric");
title("C172 best-range and best-endurance metrics");
legend("C_L/C_D","C_L^{3/2}/C_D","Location","best");

figure;
plot(curve.V_kts(valid_av),av.gamma_deg(valid_av), "o-","LineWidth",1.5);
grid on;
xlabel("True airspeed (kt)");
ylabel("Maximum steady climb angle (deg)");
title("C172 climb-angle curve");

figure;
plot(curve.V_kts(valid_av),av.ROC_fpm(valid_av), "o-","LineWidth",1.5);
hold on;
yline(0,"--");
grid on;
xlabel("True airspeed (kt)");
ylabel("Maximum steady rate of climb (ft/min)");
title("C172 climb-rate curve");

figure;
plot(curve.V_kts(valid_both),excess_thrust_hp(valid_both), "o-","LineWidth",1.5);
hold on;
yline(0,"--");
grid on;
xlabel("True airspeed (kt)");
ylabel("Excess thrust power (hp)");
title("C172 excess power");

if any(valid_req)
    idx_req = find(valid_req);

    [min_power_hp,k_power] = min(required_thrust_hp(valid_req));
    idx_min_power = idx_req(k_power);

    [max_ld,k_ld] = max(LD_curve(valid_req));
    idx_max_ld = idx_req(k_ld);

    [max_cl32cd,k_end] = max(CL32CD_curve(valid_req));
    idx_max_end = idx_req(k_end);

    fprintf("\n=== REQUIRED POWER / AERODYNAMIC CURVE ===\n");
    fprintf("Valid required points : %d / %d\n", nnz(valid_req),numel(valid_req));
    fprintf("First/last valid speed: %.3f / %.3f KTAS\n", curve.V_kts(find(valid_req,1,'first')), curve.V_kts(find(valid_req,1,'last')));
    fprintf("Minimum required thrust power: %.3f hp at %.3f KTAS\n", min_power_hp,curve.V_kts(idx_min_power));
    fprintf("Max C_L/C_D           : %.5f at %.3f KTAS\n", max_ld,curve.V_kts(idx_max_ld));
    fprintf("Max C_L^(3/2)/C_D     : %.5f at %.3f KTAS\n", max_cl32cd,curve.V_kts(idx_max_end));
end

if any(valid_av)
    idx_av = find(valid_av);
    [max_roc,k_roc] = max(av.ROC_fpm(valid_av));
    idx_max_roc = idx_av(k_roc);
    [max_gamma,k_gamma] = max(av.gamma_deg(valid_av));
    idx_max_gamma = idx_av(k_gamma);

    fprintf("\n=== AVAILABLE CLIMB CURVE ===\n");
    fprintf("Maximum climb angle   : %.4f deg at %.3f KTAS\n", max_gamma,curve.V_kts(idx_max_gamma));
    fprintf("Maximum rate of climb : %.3f ft/min at %.3f KTAS\n", max_roc,curve.V_kts(idx_max_roc));
end


performance_results = struct();
performance_results.aircraft = ac;
performance_results.mass_kg = m_total;
performance_results.cg_body_m = cg_total;
performance_results.inertia_kgm2 = I_total;
performance_results.baseline_cruise = baseline_point;
performance_results.baseline_sheet_validation = baseline_sheet_validation;
performance_results.performance_suite = performance_suite;
performance_results.performance_summary = performance_summary;
performance_results.assumptions = struct( ...
    'gross_weight_lb',target_gross_lb, ...
    'performance_altitude_ft',performance_altitude_ft, ...
    'avgas_density_kg_per_usgal',avgas_density_kg_per_usgal, ...
    'baseline_reference_rpm',baseline_reference_rpm, ...
    'baseline_reference_bhp_fraction',baseline_reference_bhp_fraction, ...
    'baseline_reference_fuel_gph',baseline_reference_fuel_gph, ...
    'baseline_CD_increment',baseline_CD_increment, ...
    'propeller_coefficient_scale',propeller_coefficient_scale, ...
    'aero_mach_limit',aero_mach_limit, ...
    'engine_model',string(o360.name), ...
    'propeller_model',string(sensenich_76em8s14.name), ...
    'propeller_coefficient_map_source',sensenich_76em8s14.map_source, ...
    'propeller_coefficient_map_status',sensenich_76em8s14.map_status, ...
    'propeller_manufacturer_validated', ...
        sensenich_76em8s14.manufacturer_validated, ...
    'propeller_calibration_description', ...
        sensenich_76em8s14.calibration_description, ...
    'aerodynamic_drag_calibration_source', ...
        "constant CD increment from the 6000-ft cruise sheet point", ...
    'moment_reference_frame_name',moment_reference_frame_name, ...
    'alpha_data_min_deg',alpha_data_min_deg, ...
    'alpha_data_max_deg',alpha_data_max_deg);

fprintf("\n=== C172 AIRBORNE PERFORMANCE COMPLETE ===\n");
fprintf("Included : VS, Vx, Vy, Vbr, Vbe, Vbg, Vh, power/climb curves\n");
fprintf("Excluded : VR and VMU\n");
fprintf("Primary prop-aircraft metrics: C_L/C_D and C_L^(3/2)/C_D\n");

function spec = make_level_spec(z0,alpha_min_deg,alpha_max_deg, elevator_min_deg,elevator_max_deg,g,reference_frame_name)
    spec = struct();
    spec.mode = "level";
    spec.variables = ["alpha";"control_pitch";"throttle"];
    spec.initial_guess = z0(:);
    spec.lb = [deg2rad(alpha_min_deg); deg2rad(elevator_min_deg); 0];
    spec.ub = [deg2rad(alpha_max_deg); deg2rad(elevator_max_deg); 1];
    spec.fixed = struct('beta',0,'phi',0,'psi',0,'gamma',0, 'p',0,'q',0,'r',0, 'control_roll',0,'control_yaw',0);
    spec.reference_frame_name = string(reference_frame_name);
    spec.residual_scale = [g;g;1];
end

function tf = point_in_c172_database(p,alpha_min_deg,alpha_max_deg,mach_limit)
    tf = isfinite(p.Mach) && p.Mach <= mach_limit+1e-8 && isfinite(p.alpha_deg) && p.alpha_deg >= alpha_min_deg-1e-6 && p.alpha_deg <= alpha_max_deg+1e-6;
end

function status = active_constraint_status(opt,p,alpha_min_deg,alpha_max_deg)
    status = "interior";

    if isfield(p,'propulsion_operating_report')
        report = p.propulsion_operating_report;
        active_propulsion = report.constraint_values >= -1e-4;
        if any(active_propulsion)
            status = strjoin(report.constraint_names(active_propulsion),"+");
            return
        end
    end

    if isfield(opt,'problem') && isstruct(opt.problem) && isfield(opt.problem,'variables') && isfield(opt.problem,'lb') && isfield(opt.problem,'ub')
        names = string(opt.problem.variables(:));
        z = opt.z_star(:);
        lb = opt.problem.lb(:);
        ub = opt.problem.ub(:);

        distance = min(abs(z-lb),abs(ub-z));
        scale = max(abs(ub-lb),1);
        active = distance <= 1e-4.*scale;

        if any(active)
            status = strjoin(names(active),"+");
            return
        end
    end

    if abs(p.alpha_deg-alpha_min_deg) <= 1e-3
        status = "alpha_lower";
    elseif abs(p.alpha_deg-alpha_max_deg) <= 1e-3
        status = "alpha_upper";
    elseif abs(p.throttle) <= 1e-5
        status = "throttle_lower";
    elseif abs(p.throttle-1) <= 1e-5
        status = "throttle_upper";
    end
end

function print_result(label,opt,engine_pe,avgas_density)
    fprintf("\n=== %s ===\n",label);

    if ~isstruct(opt) || ~isfield(opt,'converged') || ~opt.converged || ~isfield(opt,'point_opt') || isempty(opt.point_opt)
        fprintf("Converged: 0\n");
        return
    end

    p = opt.point_opt;
    engine = evaluate_engine_point(engine_pe,p);
    shaft_hp = engine.shaft_power_W/745.7;
    fuel_gph = NaN;
    if isfinite(p.fuel_flow_trim) && p.fuel_flow_trim >= 0
        fuel_gph = p.fuel_flow_trim*3600/avgas_density;
    end

    fprintf("Converged       : 1\n");
    fprintf("V               : %.6f m/s  %.3f KTAS\n",p.V_mps,p.V_kts);
    fprintf("Mach            : %.6f\n",p.Mach);
    fprintf("Alpha / gamma   : %.6f / %.6f deg\n",p.alpha_deg,p.gamma_deg);
    fprintf("Elevator        : %.6f deg\n",p.elevator_deg);
    fprintf("Throttle        : %.8f\n",p.throttle);
    fprintf("CL / CD / L-D  : %.8f / %.8f / %.8f\n", p.CL,p.CD,p.L_over_D);
    fprintf("CL^1.5 / CD     : %.8f\n", max(p.CL,0)^(3/2)/max(p.CD,1e-12));
    fprintf("ROC             : %.6f ft/min\n",p.ROC_trim_fpm);
    fprintf("Shaft power est.: %.6f hp\n",shaft_hp);
    fprintf("Engine RPM      : %.3f\n",engine.rpm);
    fprintf("Prop efficiency : %.6f\n",engine.propeller_efficiency);
    fprintf("Power residual  : %.3f W\n",engine.power_balance_error_W);
    fprintf("Propulsion valid: %d  %s\n", engine.operating_valid,engine.limit_state);
    fprintf("Fuel flow       : %.6f GPH\n",fuel_gph);
    fprintf("Residual inf.   : %.3e\n",p.residual_inf);
end

function out = evaluate_engine_point(engine_pe,point)
    out = struct('rpm',NaN, 'brake_power_W',NaN, 'shaft_power_W',NaN, 'propeller_efficiency',NaN, 'power_balance_error_W',NaN, 'operating_valid',false, 'limit_state',"not_evaluated");

    if ~isstruct(point) || ~isfield(point,'x_trim') || ~isfield(point,'u_trim') || ~isfield(point,'throttle') || ~isfinite(point.throttle) || point.throttle <= 0
        if isstruct(point) && isfield(point,'throttle') && isfinite(point.throttle) && point.throttle <= 0
            out.rpm = 0;
            out.brake_power_W = 0;
            out.shaft_power_W = 0;
            out.propeller_efficiency = 0;
            out.power_balance_error_W = 0;
            out.operating_valid = true;
            out.limit_state = "engine_off";
        end
        return
    end

    engine_pe.set_throttle(point.throttle);
    engine_pe.get_FM(point.x_trim,point.u_trim);
    operating_point = engine_pe.last_operating_point;
    operating_status = engine_pe.get_operating_status();
    out.operating_valid = operating_status.valid;
    out.limit_state = operating_status.limit_state;

    if isfield(operating_point,'rpm')
        out.rpm = operating_point.rpm;
    end
    if isfield(operating_point,'power_balance_error_W')
        out.power_balance_error_W = operating_point.power_balance_error_W;
    end
    if isfield(operating_point,'engine')
        if isfield(operating_point.engine,'brake_power_W')
            out.brake_power_W = operating_point.engine.brake_power_W;
        end
        if isfield(operating_point.engine,'shaft_power_W')
            out.shaft_power_W = operating_point.engine.shaft_power_W;
        end
    end
    if isfield(operating_point,'propeller') && isfield(operating_point.propeller,'efficiency')
        out.propeller_efficiency = operating_point.propeller.efficiency;
    end
end

function merged = merge_bidirectional_curves(down,up)
    down = reverse_curve(down);

    start_up = 1;
    if ~isempty(down.V_mps) && ~isempty(up.V_mps) && abs(down.V_mps(end)-up.V_mps(1)) <= 1e-9
        start_up = 2;
    end

    idx_up = start_up:numel(up.V_mps);

    merged = struct();
    merged.formulation = "required_available_bidirectional";
    merged.V_mps = [down.V_mps(:);up.V_mps(idx_up)];
    merged.V_kts = [down.V_kts(:);up.V_kts(idx_up)];

    merged.required = merge_branch(down.required,up.required,idx_up);
    merged.available = merge_branch(down.available,up.available,idx_up);
    merged.excess_power_W = [down.excess_power_W(:);up.excess_power_W(idx_up)];
    merged.excess_thrust_N = [down.excess_thrust_N(:);up.excess_thrust_N(idx_up)];
end

function curve = reverse_curve(curve)
    curve.V_mps = flipud(curve.V_mps(:));
    curve.V_kts = flipud(curve.V_kts(:));
    curve.required = reverse_branch(curve.required);
    curve.available = reverse_branch(curve.available);
    curve.excess_power_W = flipud(curve.excess_power_W(:));
    curve.excess_thrust_N = flipud(curve.excess_thrust_N(:));
end

function branch = reverse_branch(branch)
    fields = fieldnames(branch);
    for k = 1:numel(fields)
        name = fields{k};
        value = branch.(name);
        if iscell(value) || (isnumeric(value) && ~isscalar(value)) || (islogical(value) && ~isscalar(value)) || (isstring(value) && ~isscalar(value))
            branch.(name) = flipud(value(:));
        end
    end
end

function out = merge_branch(a,b,idx_b)
    out = struct();
    fields = union(fieldnames(a),fieldnames(b),'stable');

    for k = 1:numel(fields)
        name = fields{k};

        if ~isfield(a,name)
            out.(name) = b.(name);
            continue
        elseif ~isfield(b,name)
            out.(name) = a.(name);
            continue
        end

        va = a.(name);
        vb = b.(name);

        if iscell(va)
            out.(name) = [va(:);vb(idx_b)];
        elseif (isnumeric(va) || islogical(va) || isstring(va)) && ~isscalar(va)
            out.(name) = [va(:);vb(idx_b)];
        elseif isscalar(va) && isscalar(vb)
            out.(name) = va;
        else
            out.(name) = va;
        end
    end
end

function add_or_update_frame(ac,name,parent,position,dcm_fn)
    if ac.has_frame(name)
        ac.update_frame_position(name,position(:));
        ac.update_frame_orientation(name,dcm_fn);
    else
        ac.add_frame(name,parent,position(:),dcm_fn);
    end
end

function C = ned_to_body_dcm(x)
    phi = x(7);
    theta = x(8);
    psi = x(9);

    cphi = cos(phi);
    sphi = sin(phi);
    ctheta = cos(theta);
    stheta = sin(theta);
    cpsi = cos(psi);
    spsi = sin(psi);

    C = [ctheta*cpsi, ctheta*spsi, -stheta; sphi*stheta*cpsi-cphi*spsi, sphi*stheta*spsi+cphi*cpsi, sphi*ctheta; cphi*stheta*cpsi+sphi*spsi, cphi*stheta*spsi-sphi*cpsi, cphi*ctheta];
end
