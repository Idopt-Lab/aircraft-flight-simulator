%% F16_Paper_Aligned_Performance_Clean.m
% Full airborne trim-performance analysis:
% VS, Vx, Vy, best range, best endurance, best glide, Vh and performance curves.
% VR and VMU are excluded because they require ground-contact models.

clear; clc; close all; clear classes; clear functions; rehash toolboxcache;

%% PATH BOOTSTRAP
% Ensures the shared class library and example helper functions are on the
% MATLAB path regardless of the current working directory or session state.
this_script_dir = fileparts(mfilename('fullpath'));
matlab_root_dir = fileparts(fileparts(this_script_dir));
addpath(genpath(fullfile(matlab_root_dir,'classes')));
addpath(genpath(fullfile(matlab_root_dir,'examples')));

%% Aircraft definition

ac=Aircraft();
S=27.87; b=9.14; c=3.45;
ac.geometry.set_reference_geometry(S,b,c);

nozzle_cg=193; nozzle_aero=189.5; nozzle_tank=174.4; nozzle_engine=0;
r_aero=(nozzle_aero-nozzle_cg)*0.0254*[1;0;0];
r_tank=(nozzle_tank-nozzle_cg)*0.0254*[1;0;0];
r_engine=(nozzle_engine-nozzle_cg)*0.0254*[1;0;0];
ac.geometry.set_reference_point(r_aero.');

ac.set_body_frame("body");
ac.set_reference_frame("body");
ac.add_frame("aero_ref","body",r_aero,@(x) eye(3));
ac.add_frame("fuel_tank","body",r_tank,@(x) eye(3));
ac.add_frame("engine","body",r_engine,@(x) eye(3));
ac.add_frame("cg","body",[0;0;0],@(x) eye(3));
ac.add_frame("gravity_cg","body",[0;0;0],@(x) ReferenceFrame.ned_to_body_dcm(x));

kI=1.35582;
I=[9496*kI,0,982*kI;0,55814*kI,0;982*kI,0,63100*kI];
empty_mass=9300; fuel_mass=2400; tank_mass=fuel_mass/2;

airframe=Component("airframe",empty_mass,[0;0;0],I,ac.get_frame("body"),"airframe");
tank1=Component("tank1",tank_mass,[0;0;0],zeros(3),ac.get_frame("fuel_tank"),"fuel");
tank2=Component("tank2",tank_mass,[0;0;0],zeros(3),ac.get_frame("fuel_tank"),"fuel");
engine_component=Component("engine",0,[0;0;0],zeros(3),ac.get_frame("engine"),"engine");

airframe.set_metadata("description","F-16 airframe");
tank1.set_metadata("fuel_mass",tank_mass);
tank2.set_metadata("fuel_mass",tank_mass);
engine_component.set_metadata("max_thrust",128992);
engine_component.set_metadata("engine_type","F100-PW-229");

ac.add_component(airframe);
ac.add_component(tank1);
ac.add_component(tank2);
ac.add_component(engine_component);

ac.add_control_surface(ControlSurface("aileron","aileron","primary",[1 0 0],deg2rad(21.5),deg2rad(-21.5),0,0,0));
ac.add_control_surface(ControlSurface("elevator","elevator","primary",[0 1 0],deg2rad(25),deg2rad(-25),0,0,0));
ac.add_control_surface(ControlSurface("rudder","rudder","primary",[0 0 1],deg2rad(30),deg2rad(-30),0,0,0));

engine=TurbofanPropulsion("F100_PW_229",ac.get_frame("engine"),[1;0;0],128992,2.5,@(thr,M,alt,V) F100ThrustModel(thr,M,alt,V));
ac.add_propulsive_element(engine);
engine_component.add_load_source(PropulsionLoadSolver(engine,ac.get_frame("engine")));

aero=CoefficientAerodynamics(@F16Lookup);
ac.add_load_source(AeroLoadSolver(aero,ac.geometry,ac,ac.get_frame("aero_ref")));
ac.add_load_source(GravityLoadSolver(ac,ac.get_frame("gravity_cg")));

[mass,cg,Icg]=ac.compute_total_mass_properties();
ac.update_frame_position("cg",cg);
ac.update_frame_position("gravity_cg",cg);
ac.set_reference_frame("cg");

fprintf("\n=== F-16 MODEL ===\n");
fprintf("Mass %.2f kg | Weight %.2f N | CG [% .4f % .4f % .4f] m\n",mass,mass*ac.g,cg);
fprintf("S %.3f m^2 | b %.3f m | c %.3f m\n",S,b,c);

%% Analysis configuration

perf=PerformanceAnalysis(ac);
condition=struct('altitude_m',9144,'mach',0.60);
aero_mach_limit=1.20;
[~,a,~,~]=isa1976(condition.altitude_m);
V_database_max=aero_mach_limit*a;

solver=struct('residual_tolerance',1e-5,'inequality_tolerance',1e-6, 'equality_tolerance',1e-5,'optimality_tolerance',1e-5);

solver.fmincon_options=optimoptions("fmincon","Algorithm","sqp","Display","none", ...
    "MaxIterations",1000,"MaxFunctionEvaluations",15000, ...
    "OptimalityTolerance",1e-8,"ConstraintTolerance",1e-8, ...
    "StepTolerance",1e-9,"FiniteDifferenceType","central", ...
    "FiniteDifferenceStepSize",1e-4);

perf_cfg=struct('available_throttle',1,'stall_throttle',0.25, ...
    'takeoff_throttle',1,'continuous_throttle',1, ...
    'range_endurance_metric',"actual_flow",'propulsion_type',"thrust", ...
    'CL_max',1.6,'CL_min',-1,'max_Mach',aero_mach_limit, ...
    'max_dynamic_pressure_Pa',[]);

%% Baseline cruise trim

spec=struct();
spec.mode="level";
spec.variables=["alpha";"control_pitch";"throttle"];
spec.initial_guess=[deg2rad(7);deg2rad(-0.5);0.25];
spec.lb=[deg2rad(-5);deg2rad(-25);0];
spec.ub=[deg2rad(20);deg2rad(25);1];
spec.fixed=struct('beta',0,'phi',0,'psi',0,'gamma',0,'p',0,'q',0,'r',0, 'control_roll',0,'control_yaw',0);
spec.reference_frame_name="cg";
spec.residual_scale=[ac.g;ac.g;1];

fprintf("\nRunning baseline cruise trim...\n");
[x,u,trim_ok,trim_info]=perf.solve_trim(condition,spec,solver);
cruise=perf.evaluate_trim_point(x,u,condition,perf_cfg,trim_info.extras);
cruise.trim_converged=trim_ok;
cruise.residual_inf=trim_info.residual_inf;

fprintf("\n=== BASELINE CRUISE TRIM ===\n");
fprintf("Converged %d | Residual %.3e\n",trim_ok,trim_info.residual_inf);
fprintf("Altitude %.1f m | Mach %.4f | V %.3f m/s\n",cruise.altitude_m,cruise.Mach,cruise.V_mps);
fprintf("Alpha %.4f deg | Theta %.4f deg | Elevator %.4f deg | Throttle %.6f\n", cruise.alpha_deg,cruise.theta_deg,cruise.elevator_deg,cruise.throttle);
fprintf("CL %.5f | CD %.5f | L/D %.4f | Fuel %.5f kg/s\n", cruise.CL,cruise.CD,cruise.L_over_D,cruise.fuel_flow_trim);

%% Paper-aligned bounds

stall_ref=perf.compute_stall_speed(condition,perf_cfg);
Vmin=max(1.02*stall_ref.Vs_mps,80);
Vmax=min(350,V_database_max);

bounds=struct();

bounds.stall=struct('V_lb',60,'V_ub',min(170,V_database_max), ...
    'V0',min(max(110,stall_ref.Vs_mps),V_database_max), ...
    'alpha0_deg',15,'alpha_min_deg',-5,'alpha_max_deg',25, ...
    'elevator0_deg',0,'elevator_min_deg',-25,'elevator_max_deg',25, ...
    'throttle0',perf_cfg.stall_throttle,'throttle_min',0,'throttle_max',1);

bounds.climb=struct('V_lb',Vmin,'V_ub',Vmax,'V0',min(max(cruise.V_mps,170),Vmax), ...
    'alpha0_deg',8,'alpha_min_deg',-5,'alpha_max_deg',20, ...
    'elevator0_deg',-1,'elevator_min_deg',-25,'elevator_max_deg',25, ...
    'gamma0_deg',15,'gamma_min_deg',0,'gamma_max_deg',45,'throttle0',1);

bounds.range=struct('V_lb',Vmin,'V_ub',Vmax,'V0',min(max(cruise.V_mps,200),Vmax), ...
    'alpha0_deg',5,'alpha_min_deg',-5,'alpha_max_deg',20, ...
    'elevator0_deg',0,'elevator_min_deg',-25,'elevator_max_deg',25, ...
    'throttle0',cruise.throttle,'throttle_min',0,'throttle_max',1);

bounds.endurance=bounds.range;

bounds.glide=struct('V_lb',Vmin,'V_ub',min(320,V_database_max), ...
    'V0',min(max(cruise.V_mps,190),min(320,V_database_max)), ...
    'alpha0_deg',5,'alpha_min_deg',-5,'alpha_max_deg',20, ...
    'elevator0_deg',0,'elevator_min_deg',-25,'elevator_max_deg',25, ...
    'gamma0_deg',-7,'gamma_min_deg',-30,'gamma_max_deg',-0.01);

bounds.max_speed=struct('V_lb',Vmin,'V_ub',V_database_max,'V0',0.90*V_database_max, ...
    'alpha0_deg',1,'alpha_min_deg',-5,'alpha_max_deg',20, ...
    'elevator0_deg',0,'elevator_min_deg',-25,'elevator_max_deg',25, ...
    'throttle0',1,'throttle_min',0,'throttle_max',1);

bounds.alpha0_deg=5;
bounds.alpha_min_deg=-5;
bounds.alpha_max_deg=25;
bounds.elevator0_deg=0;
bounds.elevator_min_deg=-25;
bounds.elevator_max_deg=25;
bounds.throttle0=cruise.throttle;
bounds.throttle_min=0;
bounds.throttle_max=1;
bounds.gamma0_deg=5;
bounds.gamma_min_deg=-20;
bounds.gamma_max_deg=45;
bounds.V_sweep_mps=linspace(Vmin,V_database_max,31)';

%% Run performance suite

fprintf("\nRunning paper-aligned performance suite...\n");
tic;
suite=perf.analyze_paper_trim_suite(condition,solver,perf_cfg,bounds);
fprintf("Performance suite completed in %.1f s\n",toc);

names=["Minimum powered trim speed";"Best angle climb";"Best rate climb"; "Best range";"Best endurance";"Best engine-off glide";"Maximum level speed"];
symbols=["VS";"Vx";"Vy";"Vbr";"Vbe";"Vbg";"Vh"];
opts={suite.stall_speed;suite.best_angle_climb;suite.best_rate_climb; suite.best_range;suite.best_endurance;suite.best_glide;suite.maximum_level_speed};

n=numel(opts);
Converged=false(n,1); InDatabase=false(n,1); V_mps=nan(n,1); V_KTAS=nan(n,1);
Mach=nan(n,1); Alpha_deg=nan(n,1); Gamma_deg=nan(n,1); Elevator_deg=nan(n,1);
Throttle=nan(n,1); CL=nan(n,1); CD=nan(n,1); L_D=nan(n,1);
ROC_mps=nan(n,1); FuelFlow_kgps=nan(n,1); VelocityBound=strings(n,1);

for i=1:n
    opt=opts{i};
    Converged(i)=isstruct(opt)&&isfield(opt,"converged")&&logical(opt.converged);
    if ~Converged(i)||~isfield(opt,"point_opt")||isempty(opt.point_opt)
        VelocityBound(i)="failed"; continue;
    end

    p=opt.point_opt;
    V_mps(i)=p.V_mps; V_KTAS(i)=p.V_kts; Mach(i)=p.Mach;
    Alpha_deg(i)=p.alpha_deg; Gamma_deg(i)=p.gamma_deg; Elevator_deg(i)=p.elevator_deg;
    Throttle(i)=p.throttle; CL(i)=p.CL; CD(i)=p.CD; L_D(i)=p.L_over_D;
    ROC_mps(i)=p.ROC_trim_mps; FuelFlow_kgps(i)=p.fuel_flow_trim;

    InDatabase(i)=p.Mach<=aero_mach_limit+1e-8 && ...
        p.alpha_deg>=rad2deg(-0.175)-1e-8 && p.alpha_deg<=rad2deg(0.785)+1e-8 && ...
        abs(p.beta_deg)<=15+1e-8 && abs(p.elevator_deg)<=25+1e-8 && ...
        abs(p.aileron_deg)<=21.5+1e-8 && abs(p.rudder_deg)<=30+1e-8;

    VelocityBound(i)="interior";
    if isfield(opt,"problem")&&isfield(opt,"z_star")
        vars=string(opt.problem.variables(:));
        j=find(strcmpi(vars,"velocity"),1);
        if isempty(j)
            VelocityBound(i)="not_velocity_optimization";
        else
            z=opt.z_star(j); lb=opt.problem.lb(j); ub=opt.problem.ub(j);
            tol=max(1e-5,1e-4*max(1,ub-lb));
            if abs(z-lb)<=tol, VelocityBound(i)="lower_bound";
            elseif abs(z-ub)<=tol, VelocityBound(i)="upper_bound";
            end
        end
    end
end

Summary=table(names,symbols,Converged,InDatabase,V_mps,V_KTAS,Mach,Alpha_deg, ...
    Gamma_deg,Elevator_deg,Throttle,CL,CD,L_D,ROC_mps,FuelFlow_kgps,VelocityBound, ...
    'VariableNames',{'Metric','Symbol','Converged','InAeroDatabase','V_mps', ...
    'V_KTAS','Mach','Alpha_deg','Gamma_deg','Elevator_deg','Throttle','CL', ...
    'CD','L_D','ROC_mps','FuelFlow_kgps','VelocityBound'});

fprintf("\n=== PAPER-ALIGNED PERFORMANCE SUMMARY ===\n");
disp(Summary);

%% Required and available curves

curve=suite.required_available_curve;

if isstruct(curve)&&isfield(curve,"V_mps")&&~isempty(curve.V_mps)
    req=curve.required; av=curve.available;
    mr=req.valid(:)&isfinite(req.thrust_N(:))&isfinite(req.power_W(:));
    ma=av.valid(:)&isfinite(av.thrust_N(:))&isfinite(av.power_W(:));
    me=mr&ma&isfinite(curve.excess_thrust_N(:))&isfinite(curve.excess_power_W(:));

    figure("Name","Thrust Performance");
    tiledlayout(2,1,"TileSpacing","compact");
    nexttile;
    plot(curve.V_kts(mr),req.thrust_N(mr)/1000,"LineWidth",1.5); hold on;
    plot(curve.V_kts(ma),av.thrust_N(ma)/1000,"LineWidth",1.5);
    grid on; xlabel("TAS [kt]"); ylabel("Thrust [kN]");
    legend("Required","Available","Location","best"); title("Thrust Required and Available");
    nexttile;
    plot(curve.V_kts(me),curve.excess_thrust_N(me)/1000,"LineWidth",1.5); yline(0,"--");
    grid on; xlabel("TAS [kt]"); ylabel("Excess thrust [kN]"); title("Excess Thrust");

    figure("Name","Power and Climb Performance");
    tiledlayout(3,1,"TileSpacing","compact");
    nexttile;
    plot(curve.V_kts(mr),req.power_W(mr)/1e6,"LineWidth",1.5); hold on;
    plot(curve.V_kts(ma),av.power_W(ma)/1e6,"LineWidth",1.5);
    grid on; xlabel("TAS [kt]"); ylabel("Power [MW]");
    legend("Required","Available","Location","best"); title("Power Required and Available");
    nexttile;
    plot(curve.V_kts(me),curve.excess_power_W(me)/1e6,"LineWidth",1.5); yline(0,"--");
    grid on; xlabel("TAS [kt]"); ylabel("Excess power [MW]"); title("Excess Power");
    nexttile;
    plot(curve.V_kts(ma),av.ROC_fpm(ma),"LineWidth",1.5); yline(0,"--");
    grid on; xlabel("TAS [kt]"); ylabel("ROC [ft/min]"); title("Maximum Steady Climb Rate");
end

%% Save results

Results=struct();
Results.formulation="Gould_2025_airborne_trim_performance";
Results.excluded=["rotation_speed_VR";"minimum_unstick_speed_VMU"];
Results.aero_mach_limit=aero_mach_limit;
Results.baseline_cruise=cruise;
Results.suite=suite;
Results.summary=Summary;
Results.required_available_curve=curve;

assignin("base","Results",Results);
assignin("base","PerformanceSummary",Summary);

save_results=false;
if save_results
    save("F16_performance_results.mat","Results","-v7.3");
end

fprintf("\nAnalysis complete. Included VS, Vx, Vy, Vbr, Vbe, Vbg and Vh.\n");
fprintf("VR and VMU remain excluded. Aerodynamic database limit: Mach %.2f.\n",aero_mach_limit);
