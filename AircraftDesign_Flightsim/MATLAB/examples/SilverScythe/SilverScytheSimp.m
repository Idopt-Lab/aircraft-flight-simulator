function out = SilverScytheSimp(fidelity,options)
% SILVERSCYTHESIMP  Silver Scythe trim, analysis and Simulink export.
%
% Same role and structure as examples/F16/F16Simp.m: build the aircraft,
% trim it at a cruise point, run the framework's performance, stability,
% takeoff and landing analyses, then publish the base-workspace variables
% that classes/FlightSim_Aircraft.slx and aircraft_sfunc_generalized expect.
%
%   SilverScytheSimp                 % high fidelity, cruise at 30 kft M0.80
%   SilverScytheSimp("low")          % textbook aerodynamic and engine model
%   out = SilverScytheSimp("high",struct('mach_trim',0.90));
%
% After it returns:
%
%   open_system(fullfile(fileparts(which('Aircraft')),'FlightSim_Aircraft.slx'))
%   sim('FlightSim_Aircraft')
%
% and, for the CesiumJS view, start the web server in
% AircraftDesign_Flightsim/CesiumJS (python main.py), open
% http://localhost:8000, then run the model.
%
% Options
%   alt_trim_m, mach_trim      trim condition (default 9144 m, Mach 0.80)
%   sim_time, dt_fc            Simulink stop time and control sample time
%   disable_rate_limiting      passed through to the S-function
%   run_performance            framework performance point   (default true)
%   run_stability              linearisation and modes       (default true)
%   run_takeoff, run_landing   field-performance analyses    (default true)
%   plots                      draw the eigenvalue plot      (default false)
%   longitudinal_source        see SilverScytheContext
%   engine_mode                see SilverScytheContext

%% PATH BOOTSTRAP
% Ensures the shared class library and example helper functions are on the
% MATLAB path regardless of the current working directory or session state.
this_script_dir = fileparts(mfilename('fullpath'));
matlab_root_dir = fileparts(fileparts(this_script_dir));
addpath(genpath(fullfile(matlab_root_dir,'classes')));
addpath(genpath(fullfile(matlab_root_dir,'examples')));

if nargin < 1 || isempty(fidelity), fidelity = "high"; end
if nargin < 2 || isempty(options), options = struct(); end

defaults = struct( ...
    'alt_trim_m',9144, ...            % 30,000 ft, the report cruise altitude
    'mach_trim',0.80, ...
    'sim_time',100.0, ...
    'dt_fc',0.01, ...
    'disable_rate_limiting',true, ...
    'run_performance',true, ...
    'run_stability',true, ...
    'run_takeoff',true, ...
    'run_landing',true, ...
    'plots',false, ...
    'verbose',true);
f = fieldnames(defaults);
for i = 1:numel(f)
    if ~isfield(options,f{i}) || isempty(options.(f{i}))
        options.(f{i}) = defaults.(f{i});
    end
end

clear aircraft_sfunc_generalized

%% ------------------------------------------------------------------
%  Aircraft
%% ------------------------------------------------------------------
[ac,ctx,info] = SilverScytheBuild(fidelity,options);
ss = ctx.ss;

W = info.mass_kg*ac.g;

fprintf("\n=== SILVER SCYTHE READY (%s fidelity) ===\n",upper(char(ctx.fidelity)));
fprintf("Mass       : %.2f kg (%.0f lb, report %.0f lb)\n", ...
    info.mass_kg,info.mass_lb,info.report_mass_lb);
fprintf("Weight     : %.2f N\n",W);
fprintf("CG body    : [% .4f % .4f % .4f] m\n",info.cg_body_m);
fprintf("CG station : %.3f ft aft of nose (report %.2f ft)\n", ...
    info.cg_station_ft,info.report_cg_station_ft);
fprintf("S, b, cbar : %.3f m^2, %.3f m, %.3f m\n",ss.S_ref,ss.b_ref,ss.cbar);
fprintf("Aero ref   : station %.4f ft (%s longitudinal layer)\n", ...
    info.aero_reference_station_ft,char(info.longitudinal_source));
fprintf("Static marg: %+.4f MAC about the aero reference, %+.4f about the CG\n", ...
    info.static_margin_at_aero_ref,info.static_margin_at_cg);
fprintf("Controls   : %s\n",strjoin(cellstr(info.control_order),", "));

%% ------------------------------------------------------------------
%  Trim
%% ------------------------------------------------------------------
alt_trim_m = options.alt_trim_m;
mach_trim = options.mach_trim;
[~,a_trim,~,~] = atmosisa(alt_trim_m);
V_trim = mach_trim*a_trim;

condition = struct('altitude_m',alt_trim_m,'mach',mach_trim);

trim_cfg = struct();
trim_cfg.variables = ["alpha","elevator","throttle"];
trim_cfg.residuals = ["Fx","Fz","My"];
trim_cfg.initial_guess = [deg2rad(3.0); deg2rad(-5.0); 0.20];
trim_cfg.lb = [deg2rad(-5); -deg2rad(ss.controls.ruddervator.max_deg); 0.0];
trim_cfg.ub = [deg2rad(20);  deg2rad(ss.controls.ruddervator.max_deg); 1.0];
trim_cfg.weights = [1; 1; 1];

trim_cfg.fixed = struct();
trim_cfg.fixed.beta = 0;
trim_cfg.fixed.phi = 0;
trim_cfg.fixed.psi = 0;
trim_cfg.fixed.gamma = 0;
trim_cfg.fixed.aileron = 0;
trim_cfg.fixed.rudder = 0;
trim_cfg.fixed.nozzle_pitch = 0;
trim_cfg.fixed.nozzle_yaw = 0;
trim_cfg.fixed.p = 0;
trim_cfg.fixed.q = 0;
trim_cfg.fixed.r = 0;

trim_cfg.reference_frame_name = "cg";

solver_cfg = struct();
solver_cfg.residual_tolerance = 1e-5;
solver_cfg.fmincon_options = optimoptions("fmincon", ...
    "Algorithm","sqp","Display","iter", ...
    "MaxIterations",15000,"MaxFunctionEvaluations",100000, ...
    "OptimalityTolerance",1e-10,"StepTolerance",1e-12, ...
    "ConstraintTolerance",1e-10,"FunctionTolerance",1e-12, ...
    "FiniteDifferenceStepSize",1e-8);

fprintf("\n=== RUNNING CRUISE TRIM ===\n");
fprintf("Altitude : %.1f m (%.0f ft)\n",alt_trim_m,alt_trim_m/ss.units.ft_to_m);
fprintf("Mach     : %.3f\n",mach_trim);
fprintf("V        : %.3f m/s\n",V_trim);

ac.set_reference_frame("cg");
solver = ac.get_trim_solver();
[x_trim,u_trim,converged,trim_info] = solver.solve_trim(condition,trim_cfg,solver_cfg);

if converged
    fprintf("\n=== TRIM CONVERGED ===\n");
else
    fprintf("\n=== TRIM NOT CONVERGED: using the returned best solution ===\n");
    if abs(rad2deg(u_trim(2))) > 0.99*ss.controls.ruddervator.max_deg
        fprintf(['The ruddervator is on its %.0f deg limit.  With the ', ...
                 'textbook\nlongitudinal layer the report''s own static ', ...
                 'margin of %.2f MAC exceeds\nthe pitch authority of the ', ...
                 'report''s own control-surface sizing.  Re-run\nwith ', ...
                 'options.longitudinal_source = "datcom" for the ', ...
                 'flyable\nconfiguration; see SilverScytheValidation.\n'], ...
                 ss.controls.ruddervator.max_deg,ss.static.static_margin_MAC);
    end
end

solver.print_summary();

%% ------------------------------------------------------------------
%  Verification about the trim CG
%% ------------------------------------------------------------------
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
fprintf("Ruddervator: %.4f deg (symmetric)\n",rad2deg(u_trim(2)));
fprintf("Rudder     : %.4f deg (differential)\n",rad2deg(u_trim(3)));
fprintf("Nozzle     : %.4f deg pitch, %.4f deg yaw\n", ...
    rad2deg(u_trim(4)),rad2deg(u_trim(5)));
fprintf("Throttle   : %.6f / %.6f\n",u_trim(6),u_trim(7));
fprintf("F [N]      : [% .4e % .4e % .4e]\n",F_trim);
fprintf("M [Nm]     : [% .4e % .4e % .4e]\n",M_trim);
fprintf("Fx/W       : %.6e\n",F_trim(1)/W_trim);
fprintf("Fz/W       : %.6e\n",F_trim(3)/W_trim);
fprintf("My/Wc      : %.6e\n",M_trim(2)/(W_trim*cbar_trim));
fprintf("Fuel flow  : %.6f kg/s (%.0f lb/h)\n",fuel_flow_trim, ...
    fuel_flow_trim*3600/ss.units.lbm_to_kg);

coeff = report_trim_coefficients(ac,ctx,x_trim,u_trim);
fprintf("CL / CD    : %.5f / %.6f   L/D = %.3f\n", ...
    coeff.CL,coeff.CD,coeff.CL/coeff.CD);
fprintf("CD split   : profile %.6f, induced %.6f, wave %.6f\n", ...
    coeff.CD_profile,coeff.CD_induced,coeff.CD_wave);

out = struct();
out.fidelity = ctx.fidelity;
out.ctx = ctx;
out.info = info;
out.trim = struct('state',x_trim,'controls',u_trim,'converged',converged, ...
    'info',trim_info,'condition',condition,'coefficients',coeff, ...
    'force_N',F_trim,'moment_Nm',M_trim,'fuel_flow_kgps',fuel_flow_trim);

%% ------------------------------------------------------------------
%  Performance point
%% ------------------------------------------------------------------
perf = [];
perf_point = [];
if options.run_performance
    fprintf("\n=== RUNNING PERFORMANCE ANALYSIS ===\n");
    perf = ac.get_performance();
    perf_cfg = struct('available_throttle',1,'propulsion_type',"thrust", ...
        'reference_frame_name',"cg");
    perf_point = perf.evaluate_trim_point(x_trim,u_trim,condition,perf_cfg);
    perf_point.trim_converged = converged;
    perf.print_point(perf_point);
    out.performance_point = perf_point;
end

%% ------------------------------------------------------------------
%  Stability
%% ------------------------------------------------------------------
stab = [];
A = [];
B = [];
alpha_vec = [];
Cm_vec = [];
Cm_alpha = NaN;
beta_vec = [];
Cn_vec = [];
Cn_beta = NaN;
if options.run_stability
    fprintf("\n=== RUNNING STABILITY ANALYSIS ===\n");

    [m_trim,cg_trim,I_trim] = ac.compute_total_mass_properties(x_trim); %#ok<ASGLU>
    ac.update_frame_position("cg",cg_trim);
    ac.update_frame_position("gravity_cg",cg_trim);
    ac.set_reference_frame("cg");
    ac.state.set_full_state(x_trim);
    ac.set_controls_from_vector(u_trim);

    stab = ac.get_stability();
    stab.set_trim(x_trim,u_trim);

    [A,B] = stab.linearize(1e-5,1e-5);
    stab.analyze_modes();
    stab.print_modes();

    alpha_range = deg2rad(-5:0.5:15);
    beta_range  = deg2rad(-10:0.5:10);
    [alpha_vec,Cm_vec,Cm_alpha] = stab.compute_pitch_static_stability(alpha_range);
    [beta_vec,Cn_vec,Cn_beta]   = stab.compute_yaw_static_stability(beta_range);

    fprintf("\n=== STATIC STABILITY ===\n");
    fprintf("Cm_alpha : %+.6f 1/rad (model derivative layer %+.6f)\n", ...
        Cm_alpha,ctx.sd.Cm_alpha);
    fprintf("Cn_beta  : %+.6f 1/rad (model derivative layer %+.6f)\n", ...
        Cn_beta,ctx.sd.Cn_beta);
    if Cm_alpha < 0
        fprintf("Pitch static stability       : stable tendency\n");
    else
        fprintf(['Pitch static stability       : unstable tendency, as ', ...
                 'designed (report 12.3:\n  relaxed static stability with ', ...
                 'full-authority digital fly-by-wire)\n']);
    end
    if Cn_beta > 0
        fprintf("Directional static stability : stable tendency\n");
    else
        fprintf("Directional static stability : unstable/neutral tendency\n");
    end

    stab.print_participation();
    if options.plots
        stab.plot_eigenvalues();
    end

    out.stability = struct('A',A,'B',B,'Cm_alpha',Cm_alpha,'Cn_beta',Cn_beta, ...
        'eigenvalues',stab.eigenvalues);
end

%% ------------------------------------------------------------------
%  Takeoff and landing
%% ------------------------------------------------------------------
to = [];
TO_m = NaN;
to_res = struct();
if options.run_takeoff
    fprintf("\n=== RUNNING TAKEOFF ANALYSIS ===\n");
    to = ac.get_takeoff();
    [TO_m,to_res] = to.calculate_takeoff(0,0,3000);
    to.print_takeoff_summary(to_res);
    fprintf(['Carrier note: the RFP case is a catapult launch, not a ', ...
             'runway takeoff.\n  Catapult end speed and the launch SEROC ', ...
             'are in SilverScythePerformance.\n']);
    out.takeoff = struct('distance_m',TO_m,'results',to_res);
end

ld = [];
LD_m = NaN;
ld_res = struct();
if options.run_landing
    fprintf("\n=== RUNNING LANDING ANALYSIS ===\n");
    ld = ac.get_landing();
    [LD_m,ld_res] = ld.calculate_landing(0,0,3000);
    ld.print_landing_summary(ld_res);
    fprintf(['Carrier note: recovery is arrested, not a runway rollout.\n', ...
             '  Approach and engagement speeds are in ', ...
             'SilverScythePerformance.\n']);
    out.landing = struct('distance_m',LD_m,'results',ld_res);
end

%% ------------------------------------------------------------------
%  Simulink / base workspace export
%% ------------------------------------------------------------------
n_cs = numel(ac.control_surfaces);
n_pe = numel(ac.propulsive_elements);
n_total = n_cs + n_pe;

sim_time = options.sim_time;
dt_fc = options.dt_fc;
disable_rate_limiting = options.disable_rate_limiting;

t_vec = (0:dt_fc:sim_time).';
x0 = x_trim(:);
u0 = u_trim(:);

Initialpos = x0(1:3);
InitialVel = x0(4:6);
InitialOri = x0(7:9);
InitialRot = x0(10:12);

control_input_data = struct();
control_input_data.time = t_vec;
control_input_data.signals.values = repmat(u0.',numel(t_vec),1);
control_input_data.signals.dimensions = n_total;

ground_k = ac.ground_k;
ground_c = ac.ground_c;

assignin("base","ac",ac);
assignin("base","silver_scythe_fidelity",string(ctx.fidelity));
assignin("base","silver_scythe_ctx",ctx);
assignin("base","silver_scythe_out",out);
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
assignin("base","control_order",info.control_order);

if ~isempty(perf)
    assignin("base","perf",perf);
    assignin("base","perf_point",perf_point);
end
if ~isempty(stab)
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
end
if ~isempty(to)
    assignin("base","to",to);
    assignin("base","TO_m",TO_m);
    assignin("base","takeoff_results",to_res);
end
if ~isempty(ld)
    assignin("base","ld",ld);
    assignin("base","LD_m",LD_m);
    assignin("base","landing_results",ld_res);
end

fprintf("\n=== SIMULINK VARIABLES READY ===\n");
fprintf("fidelity         = %s\n",char(ctx.fidelity));
fprintf("engine model     = %s\n",char(ctx.thrust_mode));
fprintf("n_cs             = %d\n",n_cs);
fprintf("n_pe             = %d\n",n_pe);
fprintf("n_total          = %d\n",n_total);
fprintf("sim_stop_time    = %.2f s\n",sim_time);
fprintf("dt_fc            = %.4f s\n",dt_fc);
fprintf("rate limiting off= %d\n",disable_rate_limiting);
fprintf("x0 alpha         = %.4f deg\n",rad2deg(atan2(x0(6),x0(4))));
fprintf("control order    = %s\n",strjoin(cellstr(info.control_order),", "));
if options.run_stability
    fprintf("A_matrix size    : %d x %d\n",size(A,1),size(A,2));
    fprintf("B_matrix size    : %d x %d\n",size(B,1),size(B,2));
end
fprintf(['\nNext: open classes/FlightSim_Aircraft.slx and run it, or start\n', ...
         'AircraftDesign_Flightsim/CesiumJS (python main.py) and open\n', ...
         'http://localhost:8000 for the live world view.\n\n']);

end

% =======================================================================
function coeff = report_trim_coefficients(ac,ctx,x,u)
% Re-evaluate the coefficient lookup at the trim point for reporting.
C = SilverScytheLookup(x,u,ac.geometry,ctx);
coeff = struct('CL',C.CL,'CD',C.CD,'CY',C.CY, ...
    'Cl',C.Cl,'Cm',C.Cm,'Cn',C.Cn, ...
    'CD_profile',C.debug_CD_profile,'CD_induced',C.debug_CD_induced, ...
    'CD_wave',C.debug_CD_wave,'mach',C.debug_mach, ...
    'alpha_deg',rad2deg(C.debug_alpha_rad),'e_span',C.debug_e_span, ...
    'k_induced',C.debug_k_induced);
end
