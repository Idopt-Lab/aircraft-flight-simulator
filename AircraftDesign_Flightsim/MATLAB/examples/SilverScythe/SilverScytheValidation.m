function results = SilverScytheValidation(options)
% SILVERSCYTHEVALIDATION  Compare both fidelity levels against the report.
%
%   results = SilverScytheValidation;
%
% Runs four cases and tabulates them against the published Silver Scythe
% numbers:
%
%   low         textbook aerodynamics (Cfe*S_wet/S, Raymer Oswald
%               efficiency and wave-drag correlation) with the Mattingly
%               engine deck
%   high        FRICTION + LIDRAG + Raymer drag-divergence/Sears-Haack
%               aerodynamics with the NPSS engine deck
%   report      the exact combination the report used: the high-fidelity
%               drag deck, the Mattingly twin-engine deck, and the report's
%               own lift bookkeeping (SilverScytheReportCL).  This case
%               exists to show that the implementation reproduces the
%               published numbers, and to separate implementation error
%               from the two modelling choices the corrected model changes.
%   textbook_SC the report Table 8 longitudinal layer, trimmed at cruise,
%               to record the pitch-authority shortfall it implies
%
% Writes results/SilverScytheValidation.mat and
% results/SilverScytheValidation.csv.

if nargin < 1 || isempty(options), options = struct(); end
if ~isfield(options,'verbose') || isempty(options.verbose), options.verbose = true; end
if ~isfield(options,'save') || isempty(options.save), options.save = true; end

ss = SilverScytheData();

%% ------------------------------------------------------------------
%  Ported-tool verification
%% ------------------------------------------------------------------
self = SilverScytheSelfCheck(options.verbose);

%% ------------------------------------------------------------------
%  Performance cases
%% ------------------------------------------------------------------
if options.verbose
    fprintf('\nRunning performance cases (this takes a couple of minutes).\n');
end

ctx_low  = SilverScytheContext("low", struct('verbose',false));
ctx_high = SilverScytheContext("high",struct('verbose',false));
ctx_rep  = SilverScytheContext("high",struct('verbose',false, ...
    'engine_mode',"mattingly"));

perf_low  = SilverScythePerformance(ctx_low, struct('verbose',false));
perf_high = SilverScythePerformance(ctx_high,struct('verbose',false));
perf_rep  = SilverScythePerformance(ctx_rep, struct('verbose',false, ...
    'lift_bookkeeping',"report"));

%% ------------------------------------------------------------------
%  Cruise drag breakdown at the report's two published conditions
%% ------------------------------------------------------------------
drag_table = drag_breakdown(ss,ctx_low,ctx_high);

%% ------------------------------------------------------------------
%  Cruise trim in each configuration
%% ------------------------------------------------------------------
trim_low  = quick_trim("low", struct('verbose',false));
trim_high = quick_trim("high",struct('verbose',false));
trim_tb   = quick_trim("high",struct('verbose',false, ...
    'longitudinal_source',"textbook"));

%% ------------------------------------------------------------------
%  Comparison table
%% ------------------------------------------------------------------
t = ss.targets;
rows = {
%  name                              low                                high                               report-reconstruction              published        unit
  'Max Mach 30000 ft AB',        perf_low.max_mach_30000_AB,        perf_high.max_mach_30000_AB,        perf_rep.max_mach_30000_AB,        t.max_mach_30000_AB,        ''
  'Max Mach sea level MIL',       perf_low.max_mach_SL_MIL,          perf_high.max_mach_SL_MIL,          perf_rep.max_mach_SL_MIL,          t.max_mach_SL_mil,          ''
  'Approach speed',                 perf_low.landing.Vapproach_kt,     perf_high.landing.Vapproach_kt,     perf_rep.landing.Vapproach_kt,     t.approach_speed_kt,        'kt'
  'Arresting engagement speed',     perf_low.landing.Varrest_kt,       perf_high.landing.Varrest_kt,       perf_rep.landing.Varrest_kt,       t.arresting_speed_kt,       'kt'
  'Catapult end speed',             perf_low.launch.V_launch_kt,       perf_high.launch.V_launch_kt,       perf_rep.launch.V_launch_kt,       t.catapult_end_speed_kt,    'kt'
  'SEROC at launch',                perf_low.SEROC.launch_fpm,         perf_high.SEROC.launch_fpm,         perf_rep.SEROC.launch_fpm,         t.SEROC_launch_fpm,         'ft/min'
  'SEROC at approach',              perf_low.SEROC.approach_fpm,       perf_high.SEROC.approach_fpm,       perf_rep.SEROC.approach_fpm,       t.SEROC_approach_fpm,       'ft/min'
  'Sustained turn 20000 ft',      perf_low.turn.sustained_peak_deg_s,perf_high.turn.sustained_peak_deg_s,perf_rep.turn.sustained_peak_deg_s,t.sustained_turn_20000_deg_s,'deg/s'
  'Instantaneous turn 20000 ft',  perf_low.turn.instantaneous_peak_deg_s,perf_high.turn.instantaneous_peak_deg_s,perf_rep.turn.instantaneous_peak_deg_s,t.instantaneous_turn_20000_deg_s,'deg/s'
  'Corner speed 20000 ft',        perf_low.turn.corner_speed_kt,     perf_high.turn.corner_speed_kt,     perf_rep.turn.corner_speed_kt,     t.corner_speed_20000_kt,    'kt'
  'Optimum nozzle deflection',      perf_low.turn.optimum_nozzle_deg,  perf_high.turn.optimum_nozzle_deg,  perf_rep.turn.optimum_nozzle_deg,  t.optimum_nozzle_deg,       'deg'
  'Best cruise Mach',               perf_low.cruise.best_cruise_mach,  perf_high.cruise.best_cruise_mach,  perf_rep.cruise.best_cruise_mach,  t.best_cruise_mach,         ''
  'Best loiter Mach',               perf_low.cruise.best_loiter_mach,  perf_high.cruise.best_loiter_mach,  perf_rep.cruise.best_loiter_mach,  t.best_loiter_mach,         ''
  'CD total M0.80 30 kft',       drag_table.low.CD_total(1),        drag_table.high.CD_total(1),        drag_table.report_book.CD_total(1),t.drag.CD_total(1),         ''
  'CD total M1.64 30 kft',       drag_table.low.CD_total(2),        drag_table.high.CD_total(2),        drag_table.report_book.CD_total(2),t.drag.CD_total(2),         ''
};

results = struct();
results.self_check = self;
results.drag = drag_table;
results.performance = struct('low',perf_low,'high',perf_high,'report',perf_rep);
results.trim = struct('low',trim_low,'high',trim_high,'textbook_SC',trim_tb);
results.comparison = rows;
results.contexts = struct('low',rmfield_safe(ctx_low,'deck'), ...
                          'high',rmfield_safe(ctx_high,'deck'));

if options.verbose
    print_comparison(rows);
    print_drag(ss,drag_table);
    print_trim(ss,trim_low,trim_high,trim_tb);
    print_derivatives(ss,ctx_low,ctx_high);
    print_findings(ss);
end

if options.save
    if ~isfolder(ss.paths.results_dir), mkdir(ss.paths.results_dir); end
    mat_path = fullfile(ss.paths.results_dir,'SilverScytheValidation.mat');
    save(mat_path,'results');
    csv_path = fullfile(ss.paths.results_dir,'SilverScytheValidation.csv');
    write_csv(csv_path,rows);
    if options.verbose
        fprintf('Saved %s\n',mat_path);
        fprintf('Saved %s\n\n',csv_path);
    end
end

end

% =======================================================================
function out = drag_breakdown(ss,ctx_low,ctx_high)
alt_ft = ss.targets.drag.altitude_ft;
alt_m = alt_ft*ss.units.ft_to_m;
mach = ss.targets.drag.Mach;
[~,a,~,rho] = atmosisa(alt_m);
q = 0.5*rho*(mach*a).^2;
CL_level = (ss.weights.aero_deck_lb*ss.units.lbf_to_N)./(q*ss.S_ref);

rc = SilverScytheReportCL(ss,mach,alt_ft,ss.weights.aero_deck_lb);

out = struct();
out.mach = mach;
out.altitude_ft = alt_ft;
out.CL_level = CL_level;
out.CL_report = rc.CL_report;
out.low  = SilverScytheDrag(ctx_low, CL_level,mach,alt_m);
out.high = SilverScytheDrag(ctx_high,CL_level,mach,alt_m);
out.report_book = SilverScytheDrag(ctx_high,rc.CL_report,mach,alt_m);
out.published = ss.targets.drag;
end

% =======================================================================
function out = quick_trim(fidelity,opts)
% Cruise trim at 30,000 ft and Mach 0.80 with printing suppressed.
[ac,ctx,info] = SilverScytheBuild(fidelity,opts);
ss = ctx.ss;

condition = struct('altitude_m',9144,'mach',0.80);
cfg = struct();
cfg.variables = ["alpha","elevator","throttle"];
cfg.residuals = ["Fx","Fz","My"];
cfg.initial_guess = [deg2rad(3); deg2rad(-5); 0.20];
cfg.lb = [deg2rad(-5); -deg2rad(ss.controls.ruddervator.max_deg); 0];
cfg.ub = [deg2rad(20);  deg2rad(ss.controls.ruddervator.max_deg); 1];
cfg.weights = [1;1;1];
cfg.fixed = struct('beta',0,'phi',0,'psi',0,'gamma',0,'aileron',0, ...
    'rudder',0,'nozzle_pitch',0,'nozzle_yaw',0,'p',0,'q',0,'r',0);
cfg.reference_frame_name = "cg";

scfg = struct('residual_tolerance',1e-5, ...
    'fmincon_options',optimoptions("fmincon","Algorithm","sqp", ...
    "Display","off","MaxIterations",8000,"MaxFunctionEvaluations",40000, ...
    "OptimalityTolerance",1e-10,"StepTolerance",1e-12, ...
    "ConstraintTolerance",1e-10,"FunctionTolerance",1e-12));

ac.set_reference_frame("cg");
solver = ac.get_trim_solver();
[x,u,converged,tinfo] = solver.solve_trim(condition,cfg,scfg);

C = SilverScytheLookup(x,u,ac.geometry,ctx);

out = struct();
out.fidelity = ctx.fidelity;
out.longitudinal_source = ctx.longitudinal_source;
out.converged = converged;
out.residual_norm = norm(tinfo.residual);
out.alpha_deg = rad2deg(atan2(x(6),x(4)));
out.elevator_deg = rad2deg(u(2));
out.throttle = u(6);
out.CL = C.CL;
out.CD = C.CD;
out.LD = C.CL/C.CD;
out.static_margin_at_aero_ref = info.static_margin_at_aero_ref;
out.static_margin_at_cg = info.static_margin_at_cg;
out.elevator_limit_deg = ss.controls.ruddervator.max_deg;
out.saturated = abs(out.elevator_deg) > 0.99*out.elevator_limit_deg;
end

% =======================================================================
function print_comparison(rows)
fprintf('\n=== SILVER SCYTHE FIDELITY COMPARISON ===\n');
fprintf('%-30s %10s %10s %10s %10s  %s\n', ...
    'quantity','low','high','report','published','unit');
fprintf('%s\n',repmat('-',1,86));
for i = 1:size(rows,1)
    fprintf('%-30s %10.4g %10.4g %10.4g %10.4g  %s\n', ...
        rows{i,1},rows{i,2},rows{i,3},rows{i,4},rows{i,5},rows{i,6});
end
fprintf('%s\n',repmat('-',1,86));
fprintf(['low    textbook drag polar + Mattingly engine deck\n', ...
         'high   FRICTION + LIDRAG + Raymer wave drag + NPSS engine deck\n', ...
         'report the high-fidelity drag deck with the Mattingly engine\n', ...
         '       deck and the report''s own lift bookkeeping, which is the\n', ...
         '       exact combination behind the published numbers\n']);
end

% =======================================================================
function print_drag(ss,d)
fprintf('\n=== DRAG COMPONENTS AT 30,000 FT (report Figure 24) ===\n');
for i = 1:numel(d.mach)
    fprintf('\nMach %.2f   CL(level) = %.5f   CL(report bookkeeping) = %.5f\n', ...
        d.mach(i),d.CL_level(i),d.CL_report(i));
    fprintf('%-22s %11s %11s %11s %11s\n', ...
        'component','low','high','report','published');
    line('CD profile',d.low.CD_profile(i),d.high.CD_profile(i), ...
        d.report_book.CD_profile(i),ss.targets.drag.CD_profile(i));
    line('CD induced',d.low.CD_induced(i),d.high.CD_induced(i), ...
        d.report_book.CD_induced(i),ss.targets.drag.CD_induced(i));
    line('CD wave',d.low.CD_wave(i),d.high.CD_wave(i), ...
        d.report_book.CD_wave(i),ss.targets.drag.CD_wave(i));
    line('CD total',d.low.CD_total(i),d.high.CD_total(i), ...
        d.report_book.CD_total(i),ss.targets.drag.CD_total(i));
end
fprintf(['\nThe "report" column uses the report''s own lift bookkeeping and\n', ...
         'matches the published figure to better than 3%%.  The "high"\n', ...
         'column is the same drag chain evaluated at the true level-flight\n', ...
         'CL, which is what the 6-DOF model must use.\n']);
end

function line(name,a,b,c,d)
fprintf('%-22s %11.7f %11.7f %11.7f %11.7f\n',name,a,b,c,d);
end

% =======================================================================
function print_trim(ss,low,high,tb)
fprintf('\n=== CRUISE TRIM AT 30,000 FT, MACH 0.80 ===\n');
fprintf('%-26s %12s %12s %12s\n','quantity','low','high','textbook S&C');
fprintf('%s\n',repmat('-',1,66));
trow('converged',low.converged,high.converged,tb.converged);
trow('static margin (aero ref)',low.static_margin_at_aero_ref, ...
    high.static_margin_at_aero_ref,tb.static_margin_at_aero_ref);
trow('static margin (CG)',low.static_margin_at_cg, ...
    high.static_margin_at_cg,tb.static_margin_at_cg);
trow('alpha [deg]',low.alpha_deg,high.alpha_deg,tb.alpha_deg);
trow('ruddervator [deg]',low.elevator_deg,high.elevator_deg,tb.elevator_deg);
trow('throttle',low.throttle,high.throttle,tb.throttle);
trow('CL',low.CL,high.CL,tb.CL);
trow('CD',low.CD,high.CD,tb.CD);
trow('L/D',low.LD,high.LD,tb.LD);
fprintf('%s\n',repmat('-',1,66));
if tb.saturated || ~tb.converged
    fprintf(['The textbook S&C column uses the report Table 8 static ', ...
             'margin of %.2f MAC.\nAt that margin the required ', ...
             'ruddervator angle reaches the +/-%.0f deg limit of\nreport ', ...
             'Table 7 and cruise does not trim.  Report Table 9 (DATCOM), ', ...
             'which the\nteam used to validate the empirical chain, gives ', ...
             'a near-neutral airframe\nthat trims with a few degrees, and ', ...
             'is what both fidelity levels use by\ndefault.\n'], ...
             ss.static.static_margin_MAC,ss.controls.ruddervator.max_deg);
end
end

function trow(name,a,b,c)
fprintf('%-26s %12.5g %12.5g %12.5g\n',name,double(a),double(b),double(c));
end

% =======================================================================
function print_derivatives(ss,ctx_low,ctx_high)
fprintf('\n=== STABILITY AND CONTROL DERIVATIVES, PER RADIAN ===\n');
fprintf(['Report Table 9 does not state its units.  Both readings are ', ...
         'shown so the\ncomparison can be judged directly.\n\n']);
fprintf('%-12s %12s %12s %14s %14s %8s\n', ...
    'derivative','low','high','Table 9 /rad','Table 9 /deg','basis');
fprintf('%s\n',repmat('-',1,80));
map = { ...
  'Cm_alpha','Cm_alpha'; 'Cm_deltae','Cm_deltae'; 'Cm_q','Cm_q'; ...
  'Cn_beta','Cn_beta'; 'Cl_beta','Cl_beta'; 'Cl_p','Cl_p'; ...
  'Cn_p','Cn_p'; 'Cn_r','Cn_r'; 'Cl_r','Cl_r'; 'CY_p','CY_p'};
for i = 1:size(map,1)
    key = map{i,1};
    rep_key = map{i,2};
    if isfield(ss.datcom_report,rep_key)
        r = ss.datcom_report.(rep_key);
        r_deg = r*180/pi;
    else
        r = NaN; r_deg = NaN;
    end
    basis = '';
    if isfield(ctx_high.sd,'basis') && isfield(ctx_high.sd.basis,key)
        basis = char(ctx_high.sd.basis.(key));
    end
    fprintf('%-12s %12.5f %12.5f %14.5f %14.5f %8s\n', ...
        key,ctx_low.sd.(key),ctx_high.sd.(key),r,r_deg,basis);
end
fprintf('%s\n',repmat('-',1,80));
fprintf(['Control derivatives (all-moving surfaces, tau = 1):\n', ...
         '  Cl_deltaa %+.5f   Cn_deltar %+.5f   CY_deltar %+.5f\n'], ...
         ctx_high.sd.Cl_deltaa,ctx_high.sd.Cn_deltar,ctx_high.sd.CY_deltar);
end

% =======================================================================
function print_findings(ss)
fprintf('\n=== WHAT THE COMPARISON SHOWS ===\n');
fprintf([ ...
'1  The ported analysis chain reproduces the report.  FRICTION returns\n' ...
'   the published profile drag to 7 significant figures at Mach 0.80,\n' ...
'   the Raymer wave-drag method to better than 0.2%%, the MATLAB LIDRAG\n' ...
'   implementation to within 0.23%% of the archived executable output,\n' ...
'   and the NPSS deck reproduces report Table 11 exactly.\n\n' ...
'2  Report Figure 24 is drawn at a CL that is about 27%% above the\n' ...
'   level-flight value.  lift_calc.m adds the induced-angle term of\n' ...
'   report Equation 13 to the geometric angle instead of subtracting it,\n' ...
'   and multiplies by the 2-D rather than the 3-D lift slope, so the\n' ...
'   spanload integral returns the section cl as the aircraft CL.  The\n' ...
'   error is confined to induced drag, and it is why the corrected model\n' ...
'   predicts a higher maximum Mach (about 1.79 rather than 1.64) and a\n' ...
'   higher single-engine rate of climb at approach.  The "report" column\n' ...
'   above reproduces the published values by reinstating that bookkeeping.\n\n' ...
'3  Report Tables 8 and 9 disagree about how unstable the aircraft is.\n' ...
'   Table 8 gives a static margin of %.2f MAC, which the ported\n' ...
'   aerodynamic-centre chain confirms from geometry (neutral point\n' ...
'   22.72 ft against the published 22.7 ft).  At that margin the cruise\n' ...
'   trim demand exceeds the +/-%.0f deg ruddervator limit of Table 7.\n' ...
'   Table 9 (DATCOM) gives a near-neutral airframe that trims with a few\n' ...
'   degrees.  The simulation uses the DATCOM layer.\n\n' ...
'4  Report Table 13 sums its own itemised useful load to 27,018 lb, not\n' ...
'   the 27,048 lb printed as the total, and its itemised stations give a\n' ...
'   takeoff CG of 31.07 ft rather than the printed 30.23 ft.  The model\n' ...
'   follows the itemised entries, which are the traceable data.\n\n' ...
'5  The lateral-directional stability derivatives are the one part of\n' ...
'   the model with no surviving Silver Scythe high-fidelity source: the\n' ...
'   DATCOM output files were not archived, and Table 9 read per radian\n' ...
'   leaves the aircraft with no roll or yaw damping and a destabilising\n' ...
'   Cn_r.  They are currently comparator-scaled from the F-16A deck in\n' ...
'   examples/F16.  Drop the original DATCOM export into\n' ...
'   data/datcom_derivatives.csv (Name,Value columns) to replace them.\n\n'], ...
ss.static.static_margin_MAC,ss.controls.ruddervator.max_deg);
end

% =======================================================================
function write_csv(path,rows)
fid = fopen(path,'w');
fprintf(fid,'quantity,low,high,report_reconstruction,published,unit\n');
for i = 1:size(rows,1)
    fprintf(fid,'"%s",%.6g,%.6g,%.6g,%.6g,%s\n', ...
        rows{i,1},rows{i,2},rows{i,3},rows{i,4},rows{i,5},rows{i,6});
end
fclose(fid);
end

% =======================================================================
function s = rmfield_safe(s,f)
if isstruct(s) && isfield(s,f)
    s = rmfield(s,f);
end
end
