
%% PATH BOOTSTRAP
% Ensures the shared class library and example helper functions are on the
% MATLAB path regardless of the current working directory or session state.
this_script_dir = fileparts(mfilename('fullpath'));
matlab_root_dir = fileparts(fileparts(this_script_dir));
addpath(genpath(fullfile(matlab_root_dir,'classes')));
addpath(genpath(fullfile(matlab_root_dir,'examples')));

alt_trim_m = 9144;
mach_trim = 0.60;
sim_time = 100.0;
dt_fc = 0.01;
disable_rate_limiting = true;

%% Aircraft object

ac = Aircraft();

%% Reference geometry

S_ref = 27.87;
b_ref = 9.14;
c_ref = 3.45;
ac.geometry.set_reference_geometry(S_ref,b_ref,c_ref);

nozzle_cg = 193.0;
nozzle_aerorp = 189.5;
nozzle_fueltank = 174.4;
nozzle_engine = 0.0;

cg_to_aerorp = (nozzle_aerorp-nozzle_cg)*0.0254*[1;0;0];
cg_to_tank = (nozzle_fueltank-nozzle_cg)*0.0254*[1;0;0];
cg_to_engine = (nozzle_engine-nozzle_cg)*0.0254*[1;0;0];

ac.geometry.set_reference_point(cg_to_aerorp.');

%% Frames

ac.set_body_frame("body");
ac.set_reference_frame("body");

add_or_update_frame(ac,"aero_ref","body",cg_to_aerorp,@(x) eye(3));
add_or_update_frame(ac,"fuel_tank","body",cg_to_tank,@(x) eye(3));
add_or_update_frame(ac,"engine","body",cg_to_engine,@(x) eye(3));
add_or_update_frame(ac,"cg","body",[0;0;0],@(x) eye(3));
add_or_update_frame(ac,"gravity_cg","body",[0;0;0],@(x) ned_to_body_dcm(x));

%% Mass components

slug_ft2_to_kg_m2 = 1.35582;
Ixx = 9496*slug_ft2_to_kg_m2;
Iyy = 55814*slug_ft2_to_kg_m2;
Izz = 63100*slug_ft2_to_kg_m2;
Ixz = -982*slug_ft2_to_kg_m2;
I_body = [Ixx 0 -Ixz; 0 Iyy 0; -Ixz 0 Izz];

empty_mass = 9300;
fuel_mass = 2400;
tank_mass = fuel_mass/2;

airframe = Component("airframe",empty_mass,[0;0;0],I_body,ac.get_frame("body"),"airframe");
tank1 = Component("tank1",tank_mass,[0;0;0],zeros(3,3),ac.get_frame("fuel_tank"),"fuel");
tank2 = Component("tank2",tank_mass,[0;0;0],zeros(3,3),ac.get_frame("fuel_tank"),"fuel");
engine_comp = Component("engine",0,[0;0;0],zeros(3,3),ac.get_frame("engine"),"engine");

airframe.set_metadata("description","F16 airframe");
tank1.set_metadata("fuel_mass",tank_mass);
tank2.set_metadata("fuel_mass",tank_mass);
engine_comp.set_metadata("max_thrust",128992);
engine_comp.set_metadata("engine_type","F100-PW-229");

ac.add_component(airframe);
ac.add_component(tank1);
ac.add_component(tank2);
ac.add_component(engine_comp);

%% Controls

aileron = ControlSurface("aileron","aileron","primary",[1 0 0],deg2rad(21.5),deg2rad(-21.5),0.05,0,0);
elevator = ControlSurface("elevator","elevator","primary",[0 1 0],deg2rad(25),deg2rad(-25),0,-0.10,0);
rudder = ControlSurface("rudder","rudder","primary",[0 0 1],deg2rad(30),deg2rad(-30),0,0,-0.08);

ac.add_control_surface(aileron);
ac.add_control_surface(elevator);
ac.add_control_surface(rudder);

%% Propulsion

engine_pe = TurbofanPropulsion("F100_PW_229",ac.get_frame("engine"),[1;0;0],128992,2.5,@(thr,M,alt,V) F100ThrustModel(thr,M,alt,V));
ac.add_propulsive_element(engine_pe);

prop_solver = PropulsionLoadSolver(engine_pe,ac.get_frame("engine"));
engine_comp.add_load_source(prop_solver);

%% Aerodynamics and gravity loads

aero_model = CoefficientAerodynamics(@F16Lookup);
aero_solver = AeroLoadSolver(aero_model,ac.geometry,ac,ac.get_frame("aero_ref"));
gravity_solver = GravityLoadSolver(ac,ac.get_frame("gravity_cg"));

ac.add_load_source(aero_solver);
ac.add_load_source(gravity_solver);

%% Sync mass and CG frames

[m_total,cg_total,I_total] = ac.compute_total_mass_properties();
ac.update_frame_position("cg",cg_total);
ac.update_frame_position("gravity_cg",cg_total);

W = m_total*ac.g;

fprintf("\n=== AIRCRAFT READY ===\n");
fprintf("Mass    : %.2f kg\n",m_total);
fprintf("Weight  : %.2f N\n",W);
fprintf("CG body : [% .4f % .4f % .4f] m\n",cg_total);

%% Trim setup

[~,a_trim,~,~] = atmosisa(alt_trim_m);
V_trim = mach_trim*a_trim;

alpha0 = deg2rad(7.0);
elev0 = deg2rad(-0.5);
thr0 = 0.05;

condition = struct();
condition.altitude_m = alt_trim_m;
condition.mach = mach_trim;

trim_cfg = struct();
trim_cfg.variables = ["alpha","elevator","throttle"];
trim_cfg.residuals = ["Fx","Fz","My"];
trim_cfg.initial_guess = [alpha0; elev0; thr0];
trim_cfg.lb = [deg2rad(-5); deg2rad(-25); 0.0];
trim_cfg.ub = [deg2rad(20); deg2rad(25); 1.0];
trim_cfg.weights = [1; 1; 1];

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
solver_cfg.fmincon_options = optimoptions("fmincon","Algorithm","sqp","Display","iter","MaxIterations",15000,"MaxFunctionEvaluations",100000,"OptimalityTolerance",1e-10,"StepTolerance",1e-10,"ConstraintTolerance",1e-10,"FunctionTolerance",1e-10,"FiniteDifferenceStepSize",1e-8);

%% Run trim

fprintf("\n=== RUNNING CONFIG-DRIVEN CRUISE TRIM ===\n");
fprintf("Altitude : %.1f m\n",alt_trim_m);
fprintf("Mach     : %.3f\n",mach_trim);
fprintf("V        : %.3f m/s\n",V_trim);

ac.set_reference_frame("cg");

solver = ac.get_trim_solver();
[x_trim,u_trim,converged,info] = solver.solve_trim(condition,trim_cfg,solver_cfg);

if converged
    fprintf("\n=== TRIM CONVERGED ===\n");
else
    fprintf("\n=== TRIM NOT CONVERGED: using returned best solution ===\n");
end

solver.print_summary();

%% Final verification about CG

ac.state.set_full_state(x_trim);
ac.set_controls_from_vector(u_trim);

[m_trim,cg_trim,I_trim] = ac.compute_total_mass_properties(x_trim);
ac.update_frame_position("cg",cg_trim);
ac.update_frame_position("gravity_cg",cg_trim);
ac.set_reference_frame("cg");

[F_trim,M_trim,fuel_flow_trim] = ac.compute_total_loads(x_trim,u_trim);

W_trim = m_trim*ac.g;
cbar_trim = ac.geometry.mean_aerodynamic_chord;

fprintf("\n=== TRIM VERIFICATION ABOUT CG ===\n");
fprintf("Mass       : %.2f kg\n",m_trim);
fprintf("CG         : [% .4f % .4f % .4f] m\n",cg_trim);
fprintf("Alpha      : %.4f deg\n",rad2deg(atan2(x_trim(6),x_trim(4))));
fprintf("Theta      : %.4f deg\n",rad2deg(x_trim(8)));
fprintf("Aileron    : %.4f deg\n",rad2deg(u_trim(1)));
fprintf("Elevator   : %.4f deg\n",rad2deg(u_trim(2)));
fprintf("Rudder     : %.4f deg\n",rad2deg(u_trim(3)));
fprintf("Throttle   : %.6f\n",u_trim(4));
fprintf("F [N]      : [% .4e % .4e % .4e]\n",F_trim);
fprintf("M [Nm]     : [% .4e % .4e % .4e]\n",M_trim);
fprintf("Fx/W       : %.6e\n",F_trim(1)/W_trim);
fprintf("Fz/W       : %.6e\n",F_trim(3)/W_trim);
fprintf("My/Wc      : %.6e\n",M_trim(2)/(W_trim*cbar_trim));
fprintf("Fuel flow  : %.6f kg/s\n",fuel_flow_trim);

%% Performance analysis

fprintf("\n=== RUNNING PERFORMANCE ANALYSIS ===\n");

perf = ac.get_performance();

perf_cfg = struct('available_throttle',1, 'propulsion_type',"thrust", 'reference_frame_name',"cg");

perf_point = perf.evaluate_trim_point(x_trim,u_trim,condition,perf_cfg);
perf_point.trim_converged = converged;

perf.print_point(perf_point);

%% Stability analysis

fprintf("\n=== RUNNING STABILITY ANALYSIS ===\n");

% Make sure all reference frames are synchronized about the actual trim CG
[m_trim,cg_trim,I_trim] = ac.compute_total_mass_properties(x_trim);
ac.update_frame_position("cg",cg_trim);
ac.update_frame_position("gravity_cg",cg_trim);
ac.set_reference_frame("cg");

% Restore trim condition on aircraft object
ac.state.set_full_state(x_trim);
ac.set_controls_from_vector(u_trim);

% Get stability object from Aircraft lazy object system
stab = ac.get_stability();
stab.set_trim(x_trim,u_trim);

% Numerical linearization step sizes
dx_lin = 1e-5;
du_lin = 1e-5;

[A,B] = stab.linearize(dx_lin,du_lin);

stab.analyze_modes();
stab.print_modes();

% Static stability derivatives
alpha_range = deg2rad(-5:0.5:15);
beta_range  = deg2rad(-10:0.5:10);

[alpha_vec,Cm_vec,Cm_alpha] = stab.compute_pitch_static_stability(alpha_range);
[beta_vec,Cn_vec,Cn_beta]   = stab.compute_yaw_static_stability(beta_range);

fprintf("\n=== STATIC STABILITY ===\n");
fprintf("Cm_alpha : %.6f 1/rad\n",Cm_alpha);
fprintf("Cn_beta  : %.6f 1/rad\n",Cn_beta);

if Cm_alpha < 0
    fprintf("Pitch static stability : stable tendency\n");
else
    fprintf("Pitch static stability : unstable/neutral tendency\n");
end

if Cn_beta > 0
    fprintf("Directional static stability : stable tendency\n");
else
    fprintf("Directional static stability : unstable/neutral tendency\n");
end

% Participation factors
stab.print_participation();

% Optional eigenvalue plot
stab.plot_eigenvalues();

%% Takeoff analysis

fprintf("\n=== RUNNING TAKEOFF ANALYSIS ===\n");

to = ac.get_takeoff();

takeoff_altitude_m = 0;
runway_slope_deg = 0;
runway_length_m = 3000;

[TO_m, to_res] = to.calculate_takeoff(takeoff_altitude_m, runway_slope_deg, runway_length_m);

to.print_takeoff_summary(to_res);

%% Landing analysis

fprintf("\n=== RUNNING LANDING ANALYSIS ===\n");

ld = ac.get_landing();

landing_altitude_m = 0;
landing_temp_offset_K = 0;
landing_runway_length_m = 3000;

[LD_m, ld_res] = ld.calculate_landing(landing_altitude_m, landing_temp_offset_K, landing_runway_length_m);

ld.print_landing_summary(ld_res);

%% Simulink / base workspace export

n_cs = numel(ac.control_surfaces);
n_pe = numel(ac.propulsive_elements);
n_total = n_cs+n_pe;

t_vec = (0:dt_fc:sim_time).';

x0 = x_trim(:);
u0 = u_trim(:);

Initialpos = x0(1:3);
InitialVel = x0(4:6);
InitialOri = x0(7:9);
InitialRot = x0(10:12);

u_trim_mat = repmat(u0.',numel(t_vec),1);

control_input_data = struct();
control_input_data.time = t_vec;
control_input_data.signals.values = u_trim_mat;
control_input_data.signals.dimensions = n_total;

ground_k = ac.ground_k;
ground_c = ac.ground_c;

assignin("base","ac",ac);
assignin("base","n_cs",n_cs);
assignin("base","n_pe",n_pe);
assignin("base","n_total",n_total);
assignin("base","initial_state",x0);
assignin("base","Initialpos",Initialpos);
assignin("base","InitialVel",InitialVel);
assignin("base","InitialOri",InitialOri);
assignin("base","InitialRot",InitialRot);
assignin("base","u_trim",u0);
assignin("base","control_input_data",control_input_data);
assignin("base","sim_stop_time",sim_time);
assignin("base","dt_fc",dt_fc);
assignin("base","disable_rate_limiting",disable_rate_limiting);
assignin("base","ground_k",ground_k);
assignin("base","ground_c",ground_c);

assignin("base","perf",perf);
assignin("base","perf_point",perf_point);

assignin("base","stab",stab);
assignin("base","A_matrix",A);
assignin("base","B_matrix",B);
assignin("base","eigenvalues",stab.eigenvalues);
assignin("base","eigenvectors",stab.eigenvectors);
assignin("base","alpha_static_vec",alpha_vec);
assignin("base","Cm_static_vec",Cm_vec);
assignin("base","Cm_alpha",Cm_alpha);
assignin("base","beta_static_vec",beta_vec);
assignin("base","Cn_static_vec",Cn_vec);
assignin("base","Cn_beta",Cn_beta);

assignin("base","to",to);
assignin("base","TO_m",TO_m);
assignin("base","takeoff_results",to_res);

assignin("base","ld",ld);
assignin("base","LD_m",LD_m);
assignin("base","landing_results",ld_res);

fprintf("\n=== SIMULINK VARIABLES READY ===\n");
fprintf("n_cs             = %d\n",n_cs);
fprintf("n_pe             = %d\n",n_pe);
fprintf("n_total          = %d\n",n_total);
fprintf("sim_stop_time    = %.2f s\n",sim_time);
fprintf("dt_fc            = %.4f s\n",dt_fc);
fprintf("rate limiting off= %d\n",disable_rate_limiting);
fprintf("x0 alpha         = %.4f deg\n",rad2deg(atan2(x0(6),x0(4))));
fprintf("A_matrix size    : %d x %d\n",size(A,1),size(A,2));
fprintf("B_matrix size    : %d x %d\n",size(B,1),size(B,2));
fprintf("Eigenvalues      : %d\n",numel(stab.eigenvalues));

%% ============================================================
% Local functions
% ============================================================

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
    C = [ct*cs ct*ss -st; sp*st*cs-cp*ss sp*st*ss+cp*cs sp*ct; cp*st*cs+sp*ss cp*st*ss-sp*cs cp*ct];
end

function add_or_update_frame(ac,name,parent_name,r_parent,dcm_fn)
    if ac.has_frame(name)
        ac.update_frame_position(name,r_parent);
        ac.update_frame_orientation(name,dcm_fn);
    else
        ac.add_frame(name,parent_name,r_parent,dcm_fn);
    end
end
