function perf = SilverScythePerformance(ctx,options)
% SILVERSCYTHEPERFORMANCE  Report-style performance analysis of Silver Scythe.
%
% Reproduces the analyses of report section 16 using the ported
% aerodynamic and propulsion models, so that each published headline number
% can be compared against the model directly:
%
%   16.2  carrier catapult end speed, MIL-STD-3013B Eq 17, tropical day
%   16.3  approach and arresting-gear engagement speed
%   16.4  single-engine rate of climb at launch and at approach
%   16.5  level-flight envelope and maximum Mach
%   16.6  instantaneous and sustained turn performance at 20,000 ft
%   16.7  range and endurance factors at 30,000 ft
%
% The methods are the senior-design ones (takeoff_landing.m,
% FlightPerformance.m, climb_func.m, turn_performance.m, breguet.m), all
% recomputed here from the ported aerodynamic and engine models rather than
% from a stored deck.
%
% This function works in imperial units throughout, matching the report and
% the senior-design scripts, and converts only at the interfaces.
%
%   perf = SilverScythePerformance(SilverScytheContext("high"));
%
% options.lift_bookkeeping
%   "state"   (default) level-flight CL = W/(q*S_ref), the physically
%             correct value and the one the 6-DOF model uses.
%   "report"  the report's own CL bookkeeping via SilverScytheReportCL,
%             i.e. report Equations 13-14 as coded in lift_calc.m.  This
%             exists to demonstrate that the two published numbers the
%             corrected model does not reproduce - maximum Mach and the
%             single-engine rate of climb - are both consequences of that
%             bookkeeping and not of the drag methods.  See
%             SilverScytheValidation.

if nargin < 1 || isempty(ctx), ctx = SilverScytheContext("high"); end
if nargin < 2 || isempty(options), options = struct(); end
if ~isfield(options,'verbose') || isempty(options.verbose), options.verbose = true; end
if ~isfield(options,'mach_grid') || isempty(options.mach_grid)
    options.mach_grid = 0.10:0.02:2.20;
end
if ~isfield(options,'altitude_ft_grid') || isempty(options.altitude_ft_grid)
    options.altitude_ft_grid = 0:1000:60000;
end
if ~isfield(options,'lift_bookkeeping') || isempty(options.lift_bookkeeping)
    options.lift_bookkeeping = "state";
end
options.lift_bookkeeping = lower(string(options.lift_bookkeeping));
if ~any(options.lift_bookkeeping == ["state","report"])
    error('SilverScythePerformance:UnknownBookkeeping', ...
        'options.lift_bookkeeping must be "state" or "report".');
end
book = options.lift_bookkeeping;

ss = ctx.ss;
mach = options.mach_grid;
alt_grid = options.altitude_ft_grid;
g = 32.174;

perf = struct();
perf.fidelity = ctx.fidelity;
perf.lift_bookkeeping = book;
perf.mach = mach;
perf.altitude_ft = alt_grid;

%% ------------------------------------------------------------------
%  16.5  Level-flight envelope and maximum Mach
%% ------------------------------------------------------------------
n_alt = numel(alt_grid);
max_mach_AB = nan(1,n_alt);
max_mach_MIL = nan(1,n_alt);
drag_lbf = zeros(n_alt,numel(mach));
thrust_AB_lbf = zeros(n_alt,numel(mach));
thrust_MIL_lbf = zeros(n_alt,numel(mach));
LD = zeros(n_alt,numel(mach));

W_dash = ss.weights.aero_deck_lb;

for i = 1:n_alt
    [q,~] = dynamic_pressure(ss,mach,alt_grid(i));
    CL = level_CL(ss,W_dash,q,mach,alt_grid(i),book);
    D = SilverScytheDrag(ctx,CL,mach,alt_grid(i)*ss.units.ft_to_m);
    drag_lbf(i,:) = D.CD_total.*q*ss.S_ref_ft2;
    LD(i,:) = CL./D.CD_total;
    [thrust_AB_lbf(i,:),thrust_MIL_lbf(i,:)] = both_engines(ctx,mach,alt_grid(i));
    max_mach_AB(i)  = last_crossing(mach,thrust_AB_lbf(i,:)  - drag_lbf(i,:));
    max_mach_MIL(i) = last_crossing(mach,thrust_MIL_lbf(i,:) - drag_lbf(i,:));
end

perf.envelope = struct();
perf.envelope.max_mach_AB  = max_mach_AB;
perf.envelope.max_mach_MIL = max_mach_MIL;
perf.envelope.drag_lbf     = drag_lbf;
perf.envelope.thrust_AB_lbf = thrust_AB_lbf;
perf.envelope.thrust_MIL_lbf = thrust_MIL_lbf;
perf.envelope.LD = LD;
perf.envelope.weight_lb = W_dash;
perf.max_mach_30000_AB = interp_at(alt_grid,max_mach_AB,30000);
perf.max_mach_SL_MIL   = interp_at(alt_grid,max_mach_MIL,0);
[perf.max_mach_overall,i_best] = max(max_mach_AB);
perf.max_mach_altitude_ft = alt_grid(i_best);

%% ------------------------------------------------------------------
%  16.3  Approach and arresting-gear speeds
%
%  Senior-design takeoff_landing.m: stall at the 8.1 deg approach attitude
%  with flaps and slats, less the vertical thrust component from the
%  nozzles at 12 deg, then Vapp = 1.10 Vs and Varrest = 1.05 Vapp.
%% ------------------------------------------------------------------
[~,a_sl,rho_sl] = atmosphere_imperial(ss,0);
W_land = ss.weights.landing_analysis_lb;
CLmax_land = ctx.lift.CL_landing_ref;

thrust_lift = 2*ss.carrier.thrust_lift_per_engine_lbf* ...
    sind(ss.carrier.nozzle_launch_deg);
W_eff = W_land - thrust_lift;

Vs = sqrt(2*W_eff/(rho_sl*ss.S_ref_ft2*CLmax_land));
Vapp = Vs*ss.carrier.Vapp_over_Vs;
Varrest = Vapp*ss.carrier.Varrest_over_Vapp;

perf.landing = struct('weight_lb',W_land,'CLmax',CLmax_land, ...
    'thrust_lift_lbf',thrust_lift, ...
    'Vstall_kt',Vs/1.68781,'Vapproach_kt',Vapp/1.68781, ...
    'Varrest_kt',Varrest/1.68781,'mach_approach',Vapp/a_sl);

%% ------------------------------------------------------------------
%  16.2  Carrier catapult end speed, MIL-STD-3013B Eq 17
%
%  Tropical day: 89.8 degF, density scaled by the ISA/tropical temperature
%  ratio, as in takeoff_landing.m.
%% ------------------------------------------------------------------
[T_isa_K,~,~,rho_isa] = atmosisa(0);
T_trop = ss.carrier.tropical_day_K;
rho_trop = rho_isa*(T_isa_K/T_trop)*0.00194032;
a_trop = sqrt(1.4*287*T_trop)/ss.units.ft_to_m;

W_TO = ss.weights.TOGW_strike_lb;
Wsys = W_TO + ss.carrier.catapult_equiv_wt_lb;
Sc = ss.carrier.catapult_stroke_ft;

VDL_fit = polyfit(ss.carrier.deadload_weight_lb, ...
                  ss.carrier.deadload_endspeed_kt,2);
VDL = polyval(VDL_fit,W_TO)*1.68781;

[T_AB_total,~] = both_engines(ctx,0.20,0);
launch = solve_catapult(ctx,ss,T_AB_total,Wsys,Sc,VDL,rho_trop,a_trop,W_TO,book);
launch_SE = solve_catapult(ctx,ss,T_AB_total/2,Wsys,Sc,VDL,rho_trop,a_trop,W_TO,book);

perf.launch = struct('weight_lb',W_TO,'Wsys_lb',Wsys, ...
    'thrust_AB_lbf',T_AB_total,'VDL_kt',VDL/1.68781, ...
    'V_launch_kt',launch.V/1.68781,'mach_launch',launch.V/a_trop, ...
    'V_launch_SE_kt',launch_SE.V/1.68781, ...
    'drag_at_launch_lbf',launch.D,'rho_tropical_slug_ft3',rho_trop);

%% ------------------------------------------------------------------
%  16.4  Single-engine rate of climb, sea level
%
%  climb_func.m: excess power with one engine at afterburner, divided by
%  the relevant weight.  Evaluated at the catapult end speed against TOGW
%  and at the approach speed against the landing weight.
%% ------------------------------------------------------------------
perf.SEROC = struct();
% In "report" mode the drag reference weight is the strike-start weight at
% which the senior-design deck was generated, because climb_func.m read
% drag from that single deck for every segment while dividing excess power
% by the segment weight.  In "state" mode each condition uses its own
% weight, which is the physically correct treatment.
if book == "report"
    W_drag_launch   = ss.weights.aero_deck_lb;
    W_drag_approach = ss.weights.aero_deck_lb;
else
    W_drag_launch   = W_TO;
    W_drag_approach = W_land;
end
perf.SEROC.launch_fpm = single_engine_roc(ctx,ss,launch.V,W_TO,W_drag_launch,book);
perf.SEROC.approach_fpm = single_engine_roc(ctx,ss,Vapp,W_land,W_drag_approach,book);
perf.SEROC.requirement_launch_fpm = ss.carrier.min_SEROC_launch_fpm;
perf.SEROC.requirement_approach_fpm = ss.carrier.min_SEROC_approach_fpm;

%% ------------------------------------------------------------------
%  16.6  Turn performance at 20,000 ft
%
%  turn_performance.m: the g-limit and stall-limit instantaneous curves,
%  and the sustained curve from report Equation 21 with the nozzles at the
%  optimum 11 deg found by the report's trade study.  The report's
%  sustained solver uses zero-lift drag only, which is preserved here.
%% ------------------------------------------------------------------
alt_turn = 20000;
[~,a_turn,rho_turn] = atmosphere_imperial(ss,alt_turn);
W_combat = ss.weights.combat_lb;
n_limit = ss.static.g_limit;
CLmax_turn = ctx.lift.CL_max_clean;
delta_nozzle = ss.controls.nozzle.optimum_turn_deg;
delta_inst = ss.controls.nozzle.max_deg;

V_turn = mach*a_turn;
q_turn = 0.5*rho_turn*V_turn.^2;
[T_AB_turn,~] = both_engines(ctx,mach,alt_turn);

D_turn = SilverScytheDrag(ctx,zeros(size(mach)),mach,alt_turn*ss.units.ft_to_m);
CD0_turn = D_turn.CD_profile;
k_turn = D_turn.k_induced;

n_aero = (CLmax_turn.*q_turn*ss.S_ref_ft2 + T_AB_turn.*sind(delta_inst))/W_combat;
n_aero = max(n_aero,1);

turn_glimit = (g./V_turn).*sqrt(max(n_limit^2 - 1,0))*180/pi;
turn_stall  = (g./V_turn).*sqrt(max(n_aero.^2 - 1,0))*180/pi;

n_sust = ones(size(mach));
turn_sust = zeros(size(mach));
for i = 1:numel(mach)
    Tx = T_AB_turn(i)*cosd(delta_nozzle);
    Tv = T_AB_turn(i)*sind(delta_nozzle);
    f = @(n) Tx - q_turn(i)*ss.S_ref_ft2*(CD0_turn(i) + k_turn(i)* ...
        ((2*(n*W_combat - Tv))/(rho_turn*V_turn(i)^2*ss.S_ref_ft2))^2);
    n_hi = n_aero(i);
    if f(1) < 0
        n_sol = 1;
    elseif f(n_hi) >= 0
        n_sol = n_hi;
    else
        lo = 1; hi = n_hi;
        for it = 1:200
            mid = 0.5*(lo + hi);
            if f(mid) >= 0, lo = mid; else, hi = mid; end
        end
        n_sol = 0.5*(lo + hi);
    end
    CL_req = 2*(n_sol*W_combat - Tv)/(rho_turn*V_turn(i)^2*ss.S_ref_ft2);
    if CL_req > CLmax_turn
        n_sol = (0.5*rho_turn*V_turn(i)^2*ss.S_ref_ft2*CLmax_turn + Tv)/W_combat;
    end
    n_sol = min(max(n_sol,1),n_limit);
    n_sust(i) = n_sol;
    turn_sust(i) = (g/V_turn(i))*sqrt(max(n_sol^2 - 1,0))*180/pi;
end

corner_mach = last_crossing(mach,turn_glimit - turn_stall);
if isnan(corner_mach)
    corner_turn = NaN;
else
    corner_turn = interp1(mach,turn_glimit,corner_mach);
end
[peak_sust,i_sust] = max(turn_sust);

perf.turn = struct();
perf.turn.altitude_ft   = alt_turn;
perf.turn.weight_lb     = W_combat;
perf.turn.mach          = mach;
perf.turn.turn_glimit_deg_s = turn_glimit;
perf.turn.turn_stall_deg_s  = turn_stall;
perf.turn.turn_sustained_deg_s = turn_sust;
perf.turn.n_sustained   = n_sust;
perf.turn.n_instantaneous = min(n_aero,n_limit);
perf.turn.nozzle_deg    = delta_nozzle;
perf.turn.CLmax         = CLmax_turn;
perf.turn.corner_mach   = corner_mach;
perf.turn.corner_speed_kt = corner_mach*a_turn/1.68781;
perf.turn.instantaneous_peak_deg_s = corner_turn;
perf.turn.sustained_peak_deg_s = peak_sust;
perf.turn.sustained_peak_mach  = mach(i_sust);
perf.turn.min_radius_sustained_ft = min(V_turn.^2./ ...
    (g*sqrt(max(n_sust.^2 - 1,1e-9))));
perf.turn.min_radius_instantaneous_ft = min(V_turn.^2./ ...
    (g*sqrt(max(min(n_aero,n_limit).^2 - 1,1e-9))));

% Nozzle deflection trade study, report Figure 43
nozzle_sweep = 0:1:ss.controls.nozzle.max_deg;
sust_vs_nozzle = zeros(size(nozzle_sweep));
for j = 1:numel(nozzle_sweep)
    sust_vs_nozzle(j) = peak_sustained_turn(ctx,ss,mach,alt_turn, ...
        W_combat,nozzle_sweep(j),CLmax_turn,n_limit);
end
[~,i_opt] = max(sust_vs_nozzle);
perf.turn.nozzle_sweep_deg = nozzle_sweep;
perf.turn.sustained_vs_nozzle_deg_s = sust_vs_nozzle;
perf.turn.optimum_nozzle_deg = nozzle_sweep(i_opt);

%% ------------------------------------------------------------------
%  16.7  Range and endurance factors at 30,000 ft
%% ------------------------------------------------------------------
alt_cruise = 30000;
mach_cruise = mach(mach <= 1.20);
[q_c,a_c] = dynamic_pressure(ss,mach_cruise,alt_cruise);
CL_c = level_CL(ss,ss.weights.cruise_lb,q_c,mach_cruise,alt_cruise,book);
D_c = SilverScytheDrag(ctx,CL_c,mach_cruise,alt_cruise*ss.units.ft_to_m);
LD_c = CL_c./D_c.CD_total;
TSFC_dry = zeros(size(mach_cruise));
for i = 1:numel(mach_cruise)
    o = ctx.thrust_handle(0.5,mach_cruise(i),alt_cruise*ss.units.ft_to_m,0);
    TSFC_dry(i) = o.TSFC;
end
V_c = mach_cruise*a_c;
RF = (V_c./TSFC_dry).*LD_c;
EF = (1./TSFC_dry).*LD_c;
[maxRF,i_RF] = max(RF);
[maxEF,i_EF] = max(EF);
[maxLD,i_LD] = max(LD_c);

perf.cruise = struct('altitude_ft',alt_cruise,'weight_lb',ss.weights.cruise_lb, ...
    'mach',mach_cruise,'LD',LD_c,'TSFC_dry',TSFC_dry, ...
    'range_factor',RF,'endurance_factor',EF, ...
    'best_cruise_mach',mach_cruise(i_RF),'max_range_factor',maxRF, ...
    'best_loiter_mach',mach_cruise(i_EF),'max_endurance_factor',maxEF, ...
    'max_LD',maxLD,'max_LD_mach',mach_cruise(i_LD));

if options.verbose
    print_performance(perf,ss);
end

end

% =======================================================================
function CL = level_CL(ss,W_lb,q_psf,mach,alt_ft,book)
% Level-flight CL, optionally through the report's own bookkeeping.
if book == "report"
    r = SilverScytheReportCL(ss,mach,alt_ft,W_lb);
    CL = r.CL_report;
else
    CL = W_lb./(q_psf*ss.S_ref_ft2);
end
end

% =======================================================================
function [q,a] = dynamic_pressure(ss,mach,alt_ft)
[~,a,rho] = atmosphere_imperial(ss,alt_ft);
q = 0.5*rho*(mach*a).^2;
end

% =======================================================================
function [rho_out,a,rho] = atmosphere_imperial(ss,alt_ft)
[~,a_mps,~,rho_si] = atmosisa(alt_ft*ss.units.ft_to_m);
a = a_mps/ss.units.ft_to_m;
rho = rho_si*0.00194032;
rho_out = rho;
end

% =======================================================================
function [T_AB,T_MIL] = both_engines(ctx,mach,alt_ft)
alt_m = alt_ft*ctx.ss.units.ft_to_m;
mach = mach(:).';
T_AB = zeros(size(mach));
T_MIL = zeros(size(mach));
for i = 1:numel(mach)
    a = ctx.thrust_handle(1.0,mach(i),alt_m,0);
    m = ctx.thrust_handle(0.5,mach(i),alt_m,0);
    T_AB(i)  = 2*a.thrust_lbf;
    T_MIL(i) = 2*m.thrust_lbf;
end
end

% =======================================================================
function m = last_crossing(x,y)
% Largest x at which y changes sign, matching FlightPerformance.m.
idx = find(y(1:end-1).*y(2:end) < 0);
if isempty(idx)
    if all(y > 0)
        m = x(end);
    else
        m = NaN;
    end
    return;
end
i = idx(end);
m = interp1(y(i:i+1),x(i:i+1),0);
end

% =======================================================================
function v = interp_at(x,y,xq)
good = isfinite(y);
if nnz(good) < 2
    v = NaN;
else
    v = interp1(x(good),y(good),xq,'linear',NaN);
end
end

% =======================================================================
function out = solve_catapult(ctx,ss,T_lbf,Wsys,Sc,VDL,rho,a_fps,W_ref,book)
% MIL-STD-3013B catapult end speed, report Eq 17, iterated on drag.
V = 0.20*a_fps;
D = 0;
for it = 1:60
    mach = V/a_fps;
    q = 0.5*rho*V^2;
    CL = level_CL(ss,W_ref,q,mach,0,book);
    CL = min(CL,ctx.lift.CL_max_takeoff);
    Dr = SilverScytheDrag(ctx,CL,mach,0);
    D_new = Dr.CD_total*q*ss.S_ref_ft2;
    term = 22.5888*Sc*(T_lbf - D_new)/Wsys;
    V_new = sqrt((term + sqrt(term^2 + 4*VDL^4))/2);
    if abs(V_new - V) < 1e-6*max(V,1), V = V_new; D = D_new; break; end
    V = V_new; D = D_new;
end
out = struct('V',V,'D',D,'mach',V/a_fps);
end

% =======================================================================
function roc_fpm = single_engine_roc(ctx,ss,V_fps,W_lb,W_ref,book)
[~,a,rho] = atmosphere_imperial(ss,0);
mach = V_fps/a;
q = 0.5*rho*V_fps^2;
CL = level_CL(ss,W_ref,q,mach,0,book);
CL = min(CL,ctx.lift.CL_max_takeoff);
D = SilverScytheDrag(ctx,CL,mach,0);
drag_lbf = D.CD_total*q*ss.S_ref_ft2;
o = ctx.thrust_handle(1.0,mach,0,V_fps);
T_single = o.thrust_lbf;
roc_fpm = 60*(T_single - drag_lbf)*V_fps/W_lb;
end

% =======================================================================
function peak = peak_sustained_turn(ctx,ss,mach,alt_ft,W,nozzle_deg,CLmax,n_lim)
g = 32.174;
[~,a,rho] = atmosphere_imperial(ss,alt_ft);
V = mach*a;
q = 0.5*rho*V.^2;
[T_AB,~] = both_engines(ctx,mach,alt_ft);
D = SilverScytheDrag(ctx,zeros(size(mach)),mach,alt_ft*ss.units.ft_to_m);
CD0 = D.CD_profile;
k = D.k_induced;
turn = zeros(size(mach));
for i = 1:numel(mach)
    Tx = T_AB(i)*cosd(nozzle_deg);
    Tv = T_AB(i)*sind(nozzle_deg);
    n_hi = max((CLmax*q(i)*ss.S_ref_ft2 + T_AB(i)*sind(ss.controls.nozzle.max_deg))/W,1);
    f = @(n) Tx - q(i)*ss.S_ref_ft2*(CD0(i) + k(i)* ...
        ((2*(n*W - Tv))/(rho*V(i)^2*ss.S_ref_ft2))^2);
    if f(1) < 0
        n = 1;
    elseif f(n_hi) >= 0
        n = n_hi;
    else
        lo = 1; hi = n_hi;
        for it = 1:200
            mid = 0.5*(lo + hi);
            if f(mid) >= 0, lo = mid; else, hi = mid; end
        end
        n = 0.5*(lo + hi);
    end
    CL_req = 2*(n*W - Tv)/(rho*V(i)^2*ss.S_ref_ft2);
    if CL_req > CLmax
        n = (0.5*rho*V(i)^2*ss.S_ref_ft2*CLmax + Tv)/W;
    end
    n = min(max(n,1),n_lim);
    turn(i) = (g/V(i))*sqrt(max(n^2 - 1,0))*180/pi;
end
peak = max(turn);
end

% =======================================================================
function print_performance(perf,ss)
t = ss.targets;
fprintf('\n=== SILVER SCYTHE PERFORMANCE (%s fidelity, %s lift bookkeeping) ===\n', ...
    upper(char(perf.fidelity)),char(perf.lift_bookkeeping));
row('Max Mach, 30,000 ft, AB',   perf.max_mach_30000_AB, t.max_mach_30000_AB,'');
row('Max Mach, sea level, MIL',  perf.max_mach_SL_MIL,   t.max_mach_SL_mil,'');
row('Approach speed',            perf.landing.Vapproach_kt, t.approach_speed_kt,'kt');
row('Arresting engagement speed',perf.landing.Varrest_kt,   t.arresting_speed_kt,'kt');
row('Catapult end speed',        perf.launch.V_launch_kt,   t.catapult_end_speed_kt,'kt');
row('SEROC at launch',           perf.SEROC.launch_fpm,     t.SEROC_launch_fpm,'ft/min');
row('SEROC at approach',         perf.SEROC.approach_fpm,   t.SEROC_approach_fpm,'ft/min');
row('Sustained turn, 20,000 ft', perf.turn.sustained_peak_deg_s, ...
    t.sustained_turn_20000_deg_s,'deg/s');
row('Instantaneous turn, 20 kft',perf.turn.instantaneous_peak_deg_s, ...
    t.instantaneous_turn_20000_deg_s,'deg/s');
row('Corner speed, 20,000 ft',   perf.turn.corner_speed_kt, t.corner_speed_20000_kt,'kt');
row('Optimum nozzle deflection', perf.turn.optimum_nozzle_deg, t.optimum_nozzle_deg,'deg');
row('Best cruise Mach',          perf.cruise.best_cruise_mach, t.best_cruise_mach,'');
row('Best loiter Mach',          perf.cruise.best_loiter_mach, t.best_loiter_mach,'');
fprintf('Max L/D at 30,000 ft            : %8.3f at Mach %.2f\n', ...
    perf.cruise.max_LD,perf.cruise.max_LD_mach);
fprintf('\n');
end

function row(name,value,target,unit)
if isnan(value)
    fprintf('%-32s: %8s   report %8.2f %-6s\n',name,'n/a',target,unit);
else
    err = 100*(value/target - 1);
    fprintf('%-32s: %8.2f   report %8.2f %-6s  %+6.1f%%\n', ...
        name,value,target,unit,err);
end
end
