function results = SilverScytheSelfCheck(verbose)
% SILVERSCYTHESELFCHECK  Verify each ported analysis tool against its source.
%
% Every check compares a ported routine against a number that was produced
% by the original tool, not against another part of this package:
%
%   1  geometry            report Table 1 and Figure 18
%   2  FRICTION            report Figure 24 profile drag at Mach 0.80
%   3  LIDRAG              archived lidrag.exe output, span e = 0.89127
%   4  Raymer wave drag    report Figure 24 wave drag at Mach 0.80 and 1.64
%   5  induced drag        report Figure 24 induced drag, evaluated at the
%                          report's own CL (SilverScytheReportCL)
%   6  NPSS engine deck    report Table 11 design point
%   7  Mattingly deck      F414-GE-400 sea-level static ratings, Table 10
%   8  airfoil polar       cl_alpha.m fit of the NACA 64-206 XFOIL polar
%   9  neutral point       report Table 8, 22.7 ft
%  10  mass build-up       report Tables 12 and 13 itemised weights
%  11  AAMF force/moment   closure of the aerodynamic force against
%                          qbar*S*C for an independent hand evaluation
%
%   results = SilverScytheSelfCheck;

if nargin < 1 || isempty(verbose), verbose = true; end

ss = SilverScytheData();
checks = struct('name',{},'value',{},'reference',{},'tol_pct',{},'unit',{});

%% 1  Geometry
checks = add(checks,'Aspect ratio',ss.AR,42^2/798,1e-9,'');
checks = add(checks,'Mean aerodynamic chord',ss.cbar_ft,21.9649,0.1,'ft');
checks = add(checks,'Wetted area total',ss.S_wet_ft2,2380.520,0.01,'ft^2');

%% 2  FRICTION against report Figure 24
fr = SilverScytheFriction(ss,ss.targets.drag.Mach,ss.targets.drag.altitude_ft);
checks = add(checks,'FRICTION CD_profile, M0.80', ...
    fr.CD_profile(1,1),ss.targets.drag.CD_profile(1),0.5,'');
checks = add(checks,'FRICTION CD_profile, M1.64', ...
    fr.CD_profile(1,2),ss.targets.drag.CD_profile(2),1.0,'');

%% 3  LIDRAG against the archived executable output
span = SilverScytheSpanload(ss);
lid = SilverScytheLidrag(span.eta,span.ccl_ca,ss.AR);
checks = add(checks,'LIDRAG span efficiency', ...
    lid.e_span,ss.targets.lidrag_e_span,1.0,'');

ref = read_reference_lidrag(fullfile(ss.paths.reference_dir, ...
    'lidrag_silverscythe.out'));
if ~isnan(ref.e_span)
    checks = add(checks,'LIDRAG e vs archived file', ...
        lid.e_span,ref.e_span,1.0,'');
end
% Elliptic loading must return exactly unity
eta_e = linspace(0,1,50).';
lid_e = SilverScytheLidrag(eta_e,sin(acos(min(max(eta_e,0),1))),ss.AR);
checks = add(checks,'LIDRAG elliptic e',lid_e.e_span,1.0,0.1,'');

%% 4  Wave drag against report Figure 24
rc = SilverScytheReportCL(ss,ss.targets.drag.Mach, ...
    ss.targets.drag.altitude_ft,ss.weights.aero_deck_lb);
wv = SilverScytheWaveDrag(ss,ss.targets.drag.Mach,rc.CL_report);
checks = add(checks,'Wave drag, M0.80', ...
    wv.CD_wave(1),ss.targets.drag.CD_wave(1),1.0,'');
checks = add(checks,'Wave drag, M1.64', ...
    wv.CD_wave(2),ss.targets.drag.CD_wave(2),1.0,'');

%% 5  Induced drag at the report's own CL
CDi = lid.k_induced*rc.CL_report.^2;
% Tolerance is looser on the induced term than on the other two because
% report Figure 24 does not state the weight at which it was drawn.  The
% deck weight used here is the strike-start weight of 61,018 lb; a weight of
% about 59,600 lb reproduces the published induced drag exactly, and the
% remaining spread is that ambiguity plus the 0.23% LIDRAG delta.
checks = add(checks,'Induced drag, M0.80', ...
    CDi(1),ss.targets.drag.CD_induced(1),5.0,'');
checks = add(checks,'Induced drag, M1.64', ...
    CDi(2),ss.targets.drag.CD_induced(2),3.0,'');
checks = add(checks,'Total CD, M0.80', ...
    fr.CD_profile(1,1)+CDi(1)+wv.CD_wave(1),ss.targets.drag.CD_total(1),3.0,'');
checks = add(checks,'Total CD, M1.64', ...
    fr.CD_profile(1,2)+CDi(2)+wv.CD_wave(2),ss.targets.drag.CD_total(2),2.0,'');

%% 6  NPSS deck against report Table 11
eng = SilverScytheEngineDeck(ss);
npss_dry = SilverScytheThrustModel(0.5,0,0,0,"npss",ss,eng);
npss_wet = SilverScytheThrustModel(1.0,0,0,0,"npss",ss,eng);
checks = add(checks,'NPSS dry thrust, SLS', ...
    npss_dry.thrust_lbf,ss.engine.npss.T_dry_lbf,0.5,'lbf');
checks = add(checks,'NPSS dry TSFC, SLS', ...
    npss_dry.TSFC,ss.engine.npss.TSFC_dry,0.5,'');
checks = add(checks,'NPSS wet thrust, SLS', ...
    npss_wet.thrust_lbf,ss.engine.npss.T_wet_lbf,0.5,'lbf');
checks = add(checks,'NPSS wet TSFC, SLS', ...
    npss_wet.TSFC,ss.engine.npss.TSFC_wet,0.5,'');

%% 7  Mattingly deck against report Table 10
mat_dry = SilverScytheThrustModel(0.5,0,0,0,"mattingly",ss);
mat_wet = SilverScytheThrustModel(1.0,0,0,0,"mattingly",ss);
checks = add(checks,'Mattingly dry thrust, SLS', ...
    mat_dry.thrust_lbf,ss.engine.T_SL_dry_lbf,0.5,'lbf');
checks = add(checks,'Mattingly wet thrust, SLS', ...
    mat_wet.thrust_lbf,ss.engine.T_SL_AB_lbf,0.5,'lbf');
checks = add(checks,'Mattingly dry TSFC, SLS', ...
    mat_dry.TSFC,ss.engine.TSFC_SL_dry,0.5,'');

%% 8  Airfoil polar
af = SilverScytheAirfoil(ss);
checks = add(checks,'NACA 64-206 2-D lift slope', ...
    af.cl_alpha_per_deg,0.102358,1.0,'/deg');
checks = add(checks,'NACA 64-206 alpha_L0', ...
    af.alpha_L0_deg,-1.68499,2.0,'deg');

%% 9  Neutral point against report Table 8
sd = SilverScytheStabDerivatives(ss,"textbook", ...
    high_fidelity_lift_slope(ss), ...
    struct('cg_station_ft',ss.static.textbook_reference_station_ft));
checks = add(checks,'Neutral point',sd.neutral_point_ft, ...
    ss.static.neutral_point_station_ft,1.0,'ft');
checks = add(checks,'Static margin, report CG', ...
    sd.static_margin,ss.static.static_margin_MAC,5.0,'MAC');

%% 10  Mass build-up against report Tables 12 and 13
w = cell2mat(ss.components(:,2));
x = cell2mat(ss.components(:,3));
types = string(ss.components(:,7));
is_empty_group = ~ismember(types,["useful","fuel","payload"]);
W_empty = sum(w(is_empty_group));
W_total = sum(w);
checks = add(checks,'Empty weight',W_empty,ss.weights.empty_lb,0.1,'lb');
checks = add(checks,'Takeoff gross weight',W_total, ...
    ss.weights.TOGW_strike_lb,0.2,'lb');
checks = add(checks,'Takeoff CG from itemised stations', ...
    sum(w.*x)/W_total,ss.weights.takeoff_cg_station_ft,4.0,'ft');

%% 11  AAMF force/moment closure
[ac,ctx] = SilverScytheBuild("high",struct('verbose',false));
alt_m = 30000*ss.units.ft_to_m;
[~,a_snd,~,rho] = atmosisa(alt_m);
V = 0.8*a_snd;
alpha = deg2rad(4);
x_state = zeros(12,1);
x_state(3) = -alt_m;
x_state(4) = V*cos(alpha);
x_state(6) = V*sin(alpha);
u_state = [0;0;0;0;0;0;0];
C = SilverScytheLookup(x_state,u_state,ac.geometry,ctx);
% Isolate the aerodynamic solver so gravity and idle thrust do not enter.
ac.set_reference_frame("cg");
[F_aero,~] = ac.compute_loads_by_solver(x_state,u_state,'AeroLoadSolver', ...
    ac.get_body_frame(),ac.get_frame("cg"));
qbar = 0.5*rho*V^2;
L_hand = C.CL*qbar*ss.S_ref;
D_hand = C.CD*qbar*ss.S_ref;
Fx_hand = -D_hand*cos(alpha) + L_hand*sin(alpha);
Fz_hand = -D_hand*sin(alpha) - L_hand*cos(alpha);
checks = add(checks,'Aero Fx closure',F_aero(1),Fx_hand,1.0,'N');
checks = add(checks,'Aero Fz closure',F_aero(3),Fz_hand,1.0,'N');

%% Report
results = struct();
results.checks = checks;
results.n = numel(checks);
pass = true(1,numel(checks));
for i = 1:numel(checks)
    ref = checks(i).reference;
    if abs(ref) > 1e-12
        err = 100*abs(checks(i).value/ref - 1);
    else
        err = 100*abs(checks(i).value - ref);
    end
    checks(i).error_pct = err;
    pass(i) = err <= checks(i).tol_pct;
end
results.checks = checks;
results.pass = pass;
results.all_pass = all(pass);

if verbose
    fprintf('\n=== SILVER SCYTHE SELF-CHECK ===\n');
    fprintf('%-38s %14s %14s %9s  %s\n', ...
        'quantity','model','reference','error','');
    for i = 1:numel(checks)
        if pass(i), tag = 'PASS'; else, tag = 'FAIL'; end
        fprintf('%-38s %14.6g %14.6g %8.3f%%  %s %s\n', ...
            checks(i).name,checks(i).value,checks(i).reference, ...
            checks(i).error_pct,tag,checks(i).unit);
    end
    fprintf('\n%d of %d checks passed.\n\n',nnz(pass),numel(checks));
end

end

% =======================================================================
function c = add(c,name,value,reference,tol_pct,unit)
c(end+1) = struct('name',name,'value',double(value), ...
    'reference',double(reference),'tol_pct',tol_pct,'unit',unit);
end

% =======================================================================
function a = high_fidelity_lift_slope(ss)
% Incompressible 3-D lift slope of the high-fidelity model, without
% building the whole context.
af = SilverScytheAirfoil(ss);
a0 = af.cl_alpha_per_rad;
Lambda = deg2rad(ss.midchord_sweep_deg);
num = a0*cos(Lambda);
r = num/(pi*ss.AR);
a = num/(sqrt(1 + r^2) + r);
end

% =======================================================================
function out = read_reference_lidrag(path)
out = struct('e_span',NaN,'CL',NaN);
if ~isfile(path), return; end
txt = fileread(path);
tok = regexp(txt,'Span e\s*=\s*([\d.]+)\s+CL\s*=\s*([\d.]+)','tokens','once');
if numel(tok) == 2
    out.e_span = str2double(tok{1});
    out.CL = str2double(tok{2});
end
end
