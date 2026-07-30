clear; clc; close all; clear classes; clear functions; rehash toolboxcache;

%% PATH BOOTSTRAP
% Ensures the shared class library and example helper functions are on the
% MATLAB path regardless of the current working directory or session state.
this_script_dir = fileparts(mfilename('fullpath'));
matlab_root_dir = fileparts(fileparts(this_script_dir));
addpath(genpath(fullfile(matlab_root_dir,'classes')));
addpath(genpath(fullfile(matlab_root_dir,'examples')));

lbfN=4.4482216152605; g=9.80665; S=300*0.09290304; b=30*0.3048; c=11.3201786951274*0.3048;
Wto=31377; W=0.899666962714079*Wto; Wempty=19980.7005781593; Wperm=700; Wexp=4400; Weng=0.199*23770; Wfuel=W-Wempty-Wperm-Wexp; Wairframe=Wempty-Weng;
scale=W/20500; k=1.3558179483314; I=scale*[9496*k 0 982*k;0 55814*k 0;982*k 0 63100*k];

ac=Aircraft(); ac.geometry.set_reference_geometry(S,b,c); ac.geometry.set_reference_point([0 0 0]);
ac.add_frame("fuel_tank","body",[0;0;0],@(x) eye(3)); ac.add_frame("engine","body",[0;0;0],@(x) eye(3)); ac.add_frame("cg","body",[0;0;0],@(x) eye(3)); ac.add_frame("gravity_cg","body",[0;0;0],@(x) ReferenceFrame.ned_to_body_dcm(x));

airframe=Component("airframe",Wairframe*lbfN/g,[0;0;0],I,ac.get_frame("body"),"airframe");
fuel=Component("fuel",Wfuel*lbfN/g,[0;0;0],zeros(3),ac.get_frame("fuel_tank"),"fuel");
engineComponent=Component("engine",Weng*lbfN/g,[0;0;0],zeros(3),ac.get_frame("engine"),"engine");
permanentPayload=Component("permanent_payload",Wperm*lbfN/g,[0;0;0],zeros(3),ac.get_frame("body"),"payload");
expendablePayload=Component("expendable_payload",Wexp*lbfN/g,[0;0;0],zeros(3),ac.get_frame("body"),"payload");
ac.add_component(airframe); ac.add_component(fuel); ac.add_component(engineComponent); ac.add_component(permanentPayload); ac.add_component(expendablePayload);

ac.add_control_surface(ControlSurface("aileron","aileron","primary",[1 0 0],deg2rad(21.5),deg2rad(-21.5),0,0,0));
ac.add_control_surface(ControlSurface("stabilator","elevator","primary",[0 1 0],deg2rad(25),deg2rad(-25),0,0,0));
ac.add_control_surface(ControlSurface("rudder","rudder","primary",[0 0 1],deg2rad(30),deg2rad(-30),0,0,0));

tsfcConv=0.45359237/(lbfN*3600); p=struct();
p.number_of_engines=1; p.rated_dry_thrust_N=15000*lbfN; p.rated_afterburner_thrust_N=23770*lbfN;
p.tsfc_sl_dry_kg_N_s=0.7*tsfcConv; p.tsfc_sl_afterburner_kg_N_s=2.2*tsfcConv; p.temperature_ratio_limit=1;
p.dry_lapse_coefficients=[0.3 1 1.7]; p.afterburner_lapse_coefficients=[0.1 0.5 2.2];
p.tsfc_mach_factor=0.35; p.tsfc_dry_mach_reference=0; p.tsfc_afterburner_mach_reference=0.4;
p.tsfc_mach_exponent=1; p.tsfc_theta_exponent=0.5; p.maximum_mach=2; p.throttle_exponent=1;
p.installation_loss_factor=1; p.engine_health_factor=1; p.idle_thrust_fraction=0.025; p.idle_fuel_flow_kgps=0.035;

engine=AfterburningTurbojet("F16A_engine",ac.get_frame("engine"),[1;0;0],p); ac.add_propulsive_element(engine);
engineSolver=PropulsionLoadSolver(engine,ac.get_frame("engine")); engineComponent.add_load_source(engineSolver);
aero=CoefficientAerodynamics(@F16ABrandtLookup); aeroSolver=AeroLoadSolver(aero,ac.geometry,ac,ac.get_frame("body")); gravitySolver=GravityLoadSolver(ac,ac.get_frame("gravity_cg"));
ac.add_load_source(aeroSolver); ac.add_load_source(gravitySolver);

[mass,cg,Icg]=ac.compute_total_mass_properties(); ac.update_frame_position("cg",cg); ac.update_frame_position("gravity_cg",cg); ac.set_reference_frame("cg");

altitude_m=9000; V=250; alpha=deg2rad(5); beta=0; throttle=1;
x=zeros(12,1); x(3)=-altitude_m; x(4)=V*cos(alpha)*cos(beta); x(5)=V*sin(beta); x(6)=V*sin(alpha)*cos(beta); x(8)=alpha;
u=[0;0;0;throttle]; ac.set_controls_from_vector(u);

[FaDirect,MaDirect,coeff]=aero.get_FM(x,u,ac.geometry,ac);
assert(all(isfinite([FaDirect;MaDirect])) && coeff.CD>0,"Aerodynamic model failed.");
fprintf("\nAERODYNAMIC CHECK\n");
fprintf("Mach %.6f  alpha %.3f deg  CL %.6f  CD %.6f  valid %d\n",coeff.mach,rad2deg(coeff.alpha_rad),coeff.CL,coeff.CD,coeff.valid);
fprintf("F_aero body [N]   = [% .6e % .6e % .6e]\n",FaDirect);
fprintf("M_aero body [Nm]  = [% .6e % .6e % .6e]\n",MaDirect);

engine.set_throttle(1); engine.set_rating_mode("dry"); [Fd,Md,mdotD]=engine.get_FM(x,u); dryDebug=engine.last_debug;
engine.set_rating_mode("afterburner"); [Fab,Mab,mdotAB]=engine.get_FM(x,u); abDebug=engine.last_debug;
assert(all(isfinite([Fd;Md;Fab;Mab;mdotD;mdotAB])) && Fd(1)>0 && Fab(1)>Fd(1),"Engine model failed.");
fprintf("\nENGINE CHECK\n");
fprintf("Dry:         thrust %.6f kN  fuel %.6f kg/s  lapse %.8f\n",Fd(1)/1000,mdotD,dryDebug.lapse_raw);
fprintf("Afterburner: thrust %.6f kN  fuel %.6f kg/s  lapse %.8f\n",Fab(1)/1000,mdotAB,abDebug.lapse_raw);

engine.set_rating_mode("dry"); ac.set_controls_from_vector(u);
[FaD,MaD]=ac.compute_loads_by_solver(x,u,"AeroLoadSolver","body","cg");
[FpD,MpD,ffD]=ac.compute_loads_by_solver(x,u,"PropulsionLoadSolver","body","cg");
[FgD,MgD]=ac.compute_loads_by_solver(x,u,"GravityLoadSolver","body","cg");
[FtD,MtD,ffTotalD]=ac.compute_total_loads(x,u);
errD=norm([FtD-(FaD+FpD+FgD);MtD-(MaD+MpD+MgD)]);

engine.set_rating_mode("afterburner"); ac.set_controls_from_vector(u);
[FaA,MaA]=ac.compute_loads_by_solver(x,u,"AeroLoadSolver","body","cg");
[FpA,MpA,ffA]=ac.compute_loads_by_solver(x,u,"PropulsionLoadSolver","body","cg");
[FgA,MgA]=ac.compute_loads_by_solver(x,u,"GravityLoadSolver","body","cg");
[FtA,MtA,ffTotalA]=ac.compute_total_loads(x,u);
errA=norm([FtA-(FaA+FpA+FgA);MtA-(MaA+MpA+MgA)]);

Name=["Aerodynamics";"Gravity";"Engine dry";"Total dry";"Engine afterburner";"Total afterburner"];
FM=[FaD.' MaD.' 0;FgD.' MgD.' 0;FpD.' MpD.' ffD;FtD.' MtD.' ffTotalD;FpA.' MpA.' ffA;FtA.' MtA.' ffTotalA];
LoadCheck=table(Name,FM(:,1),FM(:,2),FM(:,3),FM(:,4),FM(:,5),FM(:,6),FM(:,7),'VariableNames',{'Source','Fx_N','Fy_N','Fz_N','Mx_Nm','My_Nm','Mz_Nm','FuelFlow_kgps'});
fprintf("\nTOTAL AIRCRAFT LOAD CHECK ABOUT CG\n"); disp(LoadCheck);
fprintf("Dry source-sum error         %.6e\n",errD);
fprintf("Afterburner source-sum error %.6e\n",errA);
fprintf("Mass %.6f kg  weight %.6f kN\n",mass,mass*g/1000);
assert(errD<1e-8*max(1,norm([FtD;MtD])) && errA<1e-8*max(1,norm([FtA;MtA])),"Aircraft load-source summation failed.");

figure("Name","F-16A Load Verification");
bar(categorical(["Aero","Gravity","Engine dry","Total dry","Engine AB","Total AB"]),FM(:,1:3)/1000); grid on; ylabel("Body force [kN]"); legend("F_x","F_y","F_z","Location","best");

F16LoadVerification=struct('condition',struct('altitude_m',altitude_m,'V_mps',V,'alpha_deg',rad2deg(alpha)),'mass_kg',mass,'cg_m',cg,'inertia_kgm2',Icg,'aero_coefficients',coeff,'loads',LoadCheck,'dry_reconciliation_error',errD,'afterburner_reconciliation_error',errA);
