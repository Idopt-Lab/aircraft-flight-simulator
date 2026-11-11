function ctx = SilverScytheContext(fidelity,options)
% SILVERSCYTHECONTEXT  Assemble one fidelity level of the Silver Scythe model.
%
% Everything expensive - the XFOIL polar, the FRICTION/LIDRAG/wave-drag
% deck, the NPSS engine deck, the derivative build-up - is prepared once
% here and handed to the coefficient lookup, the aircraft builder and the
% performance scripts as a plain struct.
%
%   ctx = SilverScytheContext("low")     textbook / empirical model
%   ctx = SilverScytheContext("high")    report's higher-fidelity model
%
% Fidelity summary
% ---------------------------------------------------------------------
%                low                          high
%   section      thin-airfoil 2*pi,           NACA 64-206 XFOIL polar
%   lift         report alpha_L0 = -1.05 deg  (cl_alpha.m fit)
%   profile      Cfe*S_wet/S, Raymer 12.3     FRICTION (Mason), per
%   drag                                      component, Raymer form and
%                                             interference factors, 5%
%                                             stealth leakage
%   induced      Raymer Oswald efficiency,    LIDRAG (Ives/Mason) on the
%   drag         supersonic k correlation     50-station spanload
%   wave drag    Raymer 4.5*pi/S correlation  Raymer drag-divergence curve
%                                             fits + Sears-Haack, report
%                                             Eqs 15-16
%   engine       Mattingly algebraic deck     NPSS deck, Mattingly-scaled
%                (report 13.5)                outside coverage (13.6)
%
% The longitudinal static-stability layer is selected separately, by
% options.longitudinal_source, so that the fidelity switch isolates the
% aerodynamic and propulsion methods.  See the note on that option below.
%
% Options
%   rebuild_deck   force the aero deck to be regenerated (default false)
%   verbose        print progress (default true)
%   datcom_csv     path to a Name,Value CSV overriding any derivative
%                  (default data/datcom_derivatives.csv if present)
%   engine_mode    "auto" (default), "mattingly" or "npss".  "auto" is
%                  Mattingly for low fidelity and NPSS for high.  Selecting
%                  it independently of the aerodynamic fidelity is what
%                  lets SilverScytheValidation rebuild the exact
%                  combination the report used, which was the FRICTION /
%                  LIDRAG drag deck with the Mattingly twin-engine deck.
%   longitudinal_source  "auto" (default), "textbook" or "datcom".  Selects
%                  where the static margin and pitch control power come
%                  from.  "auto" resolves to "datcom" for both fidelities
%                  so that the fidelity switch isolates the aerodynamic and
%                  propulsion methods, which is the comparison of interest.
%
%                  "textbook" reproduces report Table 8 (neutral point
%                  22.7 ft, static margin -0.23 MAC) from the senior-design
%                  AC / tail-volume chain.  That configuration does not
%                  trim in cruise: the required ruddervator angle exceeds
%                  the +/-25 deg limit of report Table 7, which is a real
%                  finding about the report's empirical S&C chain rather
%                  than a modelling artefact, and is reported by
%                  SilverScytheValidation.  Report Table 9 (DATCOM), the
%                  higher-fidelity tool the team used to validate that
%                  chain, gives a near-neutral airframe that trims with a
%                  few degrees of ruddervator.

if nargin < 1 || isempty(fidelity), fidelity = "high"; end
if nargin < 2 || isempty(options), options = struct(); end
fidelity = lower(string(fidelity));
if ~any(fidelity == ["low","high"])
    error('SilverScytheContext:UnknownFidelity', ...
        'fidelity must be "low" or "high".');
end

defaults = struct('rebuild_deck',false,'verbose',true,'datcom_csv',"", ...
    'longitudinal_source',"auto",'engine_mode',"auto",'eta_tail',0.95);
options = fill_defaults(options,defaults);
options.longitudinal_source = lower(string(options.longitudinal_source));
options.engine_mode = lower(string(options.engine_mode));

ss = SilverScytheData();

if options.longitudinal_source == "auto"
    % Both fidelities use the DATCOM longitudinal layer so that the
    % fidelity switch isolates the aerodynamic and propulsion methods.
    longitudinal = "datcom";
else
    longitudinal = options.longitudinal_source;
end
if ~any(longitudinal == ["textbook","datcom"])
    error('SilverScytheContext:UnknownLongitudinalSource', ...
        'longitudinal_source must be "auto", "textbook" or "datcom".');
end

% Aerodynamic moment reference station, matched to the selected source.
if longitudinal == "textbook"
    aero_reference_station_ft = ss.static.textbook_reference_station_ft;
else
    aero_reference_station_ft = ss.static.datcom_reference_station_ft;
end
ctx_longitudinal = longitudinal;

ctx = struct();
ctx.ss = ss;
ctx.fidelity = fidelity;
ctx.options = options;

%% ------------------------------------------------------------------
%  Section data and the 3-D lift curve
%% ------------------------------------------------------------------
ctx.airfoil = SilverScytheAirfoil(ss);

switch fidelity
    case "low"
        a0_per_rad = ss.low.airfoil_a0_per_deg*180/pi;      % = 2*pi
        alpha_L0_deg = ss.low.alpha_L0_deg;                 % report value
        lift_source = "thin-airfoil 2*pi per rad, report alpha_L0";
    case "high"
        a0_per_rad = ctx.airfoil.cl_alpha_per_rad;
        alpha_L0_deg = ctx.airfoil.alpha_L0_deg;
        lift_source = "NACA 64-206 XFOIL polar, cl_alpha.m full-range fit";
end

Lambda_mid = deg2rad(ss.midchord_sweep_deg);
CL_alpha_incompressible = helmbold(a0_per_rad,Lambda_mid,ss.AR);

ctx.lift = struct();
ctx.lift.a0_per_rad   = a0_per_rad;
ctx.lift.a0_per_deg   = a0_per_rad*pi/180;
ctx.lift.alpha_L0_deg = alpha_L0_deg;
ctx.lift.alpha_L0_rad = deg2rad(alpha_L0_deg);
ctx.lift.CL_alpha_incompressible = CL_alpha_incompressible;
ctx.lift.source = lift_source;

% CLmax cases, following the senior-design stall_solver.m usage
alpha_max = ss.airfoil.alpha_max_deg;
alpha_land = ss.airfoil.alpha_max_land_deg;
CL_clean_geometric = CL_alpha_incompressible*deg2rad(alpha_max - alpha_L0_deg);
ctx.lift.CL_max_geometric = CL_clean_geometric;
ctx.lift.CL_max_clean = CL_clean_geometric*ss.high_lift.LEX_vortex_multiplier;
ctx.lift.CL_max_takeoff = CL_clean_geometric*ss.high_lift.takeoff_multiplier;
ctx.lift.CL_landing_ref = CL_alpha_incompressible* ...
    deg2rad(alpha_land - alpha_L0_deg)*ss.high_lift.landing_multiplier;

%% ------------------------------------------------------------------
%  Drag
%% ------------------------------------------------------------------
switch fidelity
    case "low"
        ctx.drag_low = textbook_drag_constants(ss,ctx.lift);
        ctx.deck = [];
    case "high"
        ctx.deck = SilverScytheAeroDeck(ss, ...
            struct('rebuild',options.rebuild_deck,'verbose',options.verbose));
        ctx.drag_low = [];
end

%% ------------------------------------------------------------------
%  Stability and control derivatives
%% ------------------------------------------------------------------
ctx.longitudinal_source = ctx_longitudinal;
ctx.aero_reference_station_ft = aero_reference_station_ft;

ctx.sd = SilverScytheStabDerivatives(ss,ctx_longitudinal,CL_alpha_incompressible, ...
    struct('cg_station_ft',aero_reference_station_ft, ...
           'eta_tail',options.eta_tail));

datcom_path = char(options.datcom_csv);
if isempty(datcom_path)
    datcom_path = fullfile(ss.paths.data_dir,'datcom_derivatives.csv');
end
if isfile(datcom_path)
    ctx.sd = override_from_csv(ctx.sd,datcom_path);
    ctx.sd.override_file = datcom_path;
end

%% ------------------------------------------------------------------
%  Propulsion
%% ------------------------------------------------------------------
if options.engine_mode == "auto"
    if fidelity == "low"
        ctx.thrust_mode = "mattingly";
    else
        ctx.thrust_mode = "npss";
    end
else
    ctx.thrust_mode = options.engine_mode;
end
if ~any(ctx.thrust_mode == ["mattingly","npss"])
    error('SilverScytheContext:UnknownEngineMode', ...
        'engine_mode must be "auto", "mattingly" or "npss".');
end
if ctx.thrust_mode == "npss"
    ctx.engine = SilverScytheEngineDeck(ss);
else
    ctx.engine = [];
end
ctx.thrust_handle = @(thr,M,alt,V) SilverScytheThrustModel( ...
    thr,M,alt,V,ctx.thrust_mode,ss,ctx.engine);

%% ------------------------------------------------------------------
%  Takeoff / landing parameter structs handed to the analysis classes
%% ------------------------------------------------------------------
drag_ref = SilverScytheDrag(ctx,0.9,0.20,0);

ctx.takeoff_params = ss.takeoff_params;
ctx.takeoff_params.CLmax_takeoff  = ctx.lift.CL_max_takeoff;
ctx.takeoff_params.CLalpha_takeoff = CL_alpha_incompressible;
ctx.takeoff_params.CD0_takeoff    = drag_ref.CD_profile + 0.020; % gear+flap
ctx.takeoff_params.K_takeoff      = max(drag_ref.k_induced,1e-6);

ctx.landing_params = ss.landing_params;
ctx.landing_params.CLmax_landing = ctx.lift.CL_landing_ref;
ctx.landing_params.CD0_landing   = drag_ref.CD_profile + 0.045; % gear+flap
ctx.landing_params.K_landing     = max(drag_ref.k_induced,1e-6);

if options.verbose
    print_summary(ctx);
end

end

% =======================================================================
function a = helmbold(a0,Lambda,AR)
num = a0*cos(Lambda);
r = num/(pi*AR);
a = num/(sqrt(1 + r^2) + r);
end

% =======================================================================
function p = textbook_drag_constants(ss,lift)
% Constants of the senior-design textbook polar drag_calc_sol.m.

AR = ss.AR;
LE = ss.LE_sweep_deg;
TE = ss.TE_sweep_deg;

% Raymer swept-wing Oswald efficiency, drag_calc_sol.m form
e_osw = 4.6*(1 - 0.033*AR^0.53)*(cosd(LE)^0.1) - 3.3;
e_osw = min(max(e_osw,0.30),1.0);

% Average sweep at maximum thickness, drag_calc_sol.m definition
Lambda_max_t = (LE - TE)/2;
e_alpha = max(2/(2 - AR + sqrt(4 + AR^2*(1 + tand(Lambda_max_t)^2))), ...
    ss.low.e_alpha_floor);

a0_deg = lift.a0_per_deg;
CL_alpha_deg = a0_deg/(1 + (57.3*a0_deg/(pi*e_alpha*AR)));
CL_minD = CL_alpha_deg*(-ss.low.alpha_L0_deg)/2;

p = struct();
p.AR            = AR;
p.LE_sweep_deg  = LE;
p.e_osw         = e_osw;
p.e_alpha       = e_alpha;
p.CL_alpha_per_deg = CL_alpha_deg;
p.CL_minD       = CL_minD;
p.k1_sub        = 1/(pi*e_osw*AR);
p.k2_sub        = -2*(1/(pi*e_osw*AR))*CL_minD;
p.CDmin_sub     = ss.low.Cfe*(ss.S_wet_ft2/ss.S_ref_ft2);
p.M_CD0         = real(1/(cosd(LE)^0.2));
p.M_cutoff      = 1/sind(90 - LE);
p.wave_multiplier = 4.5*pi/ss.S_ref_ft2* ...
    (ss.max_cross_section_ft2/ss.wave_reference_length_ft)^2* ...
    ss.low.EWD*(0.74 + 0.37*cosd(LE));
p.blend         = ss.low.transonic_blend;
p.source        = "drag_calc_sol.m textbook polar";

end

% =======================================================================
function sd = override_from_csv(sd,path)
T = readtable(path);
names = string(T.Properties.VariableNames);
iN = find(strcmpi(names,'Name'),1);
iV = find(strcmpi(names,'Value'),1);
if isempty(iN) || isempty(iV)
    error('SilverScytheContext:BadDatcomCsv', ...
        '%s must have Name and Value columns.',path);
end
applied = strings(0,1);
for i = 1:height(T)
    key = char(string(T{i,iN}));
    value = T{i,iV};
    if iscell(value), value = value{1}; end
    if isfield(sd,key)
        sd.(key) = double(value);
        applied(end+1,1) = string(key); %#ok<AGROW>
    end
end
sd.overridden = applied;
if isfield(sd,'CL_alpha') && isfield(sd,'Cm_alpha')
    sd.Cm_CL = sd.Cm_alpha/max(sd.CL_alpha,1e-9);
    sd.static_margin = -sd.Cm_CL;
    sd.static_margin_source = "external DATCOM CSV";
end
end

% =======================================================================
function options = fill_defaults(options,defaults)
f = fieldnames(defaults);
for i = 1:numel(f)
    if ~isfield(options,f{i}) || isempty(options.(f{i}))
        options.(f{i}) = defaults.(f{i});
    end
end
end

% =======================================================================
function print_summary(ctx)
ss = ctx.ss;
sd = ctx.sd;
fprintf('\n=== SILVER SCYTHE MODEL CONTEXT: %s FIDELITY ===\n', ...
    upper(char(ctx.fidelity)));
fprintf('Lift source        : %s\n',ctx.lift.source);
fprintf('  a0               : %.5f /rad (%.6f /deg)\n', ...
    ctx.lift.a0_per_rad,ctx.lift.a0_per_deg);
fprintf('  alpha_L0         : %.4f deg\n',ctx.lift.alpha_L0_deg);
fprintf('  CL_alpha (M=0)   : %.4f /rad\n',ctx.lift.CL_alpha_incompressible);
fprintf('  CLmax clean      : %.4f (geometric %.4f x LEX %.2f)\n', ...
    ctx.lift.CL_max_clean,ctx.lift.CL_max_geometric, ...
    ss.high_lift.LEX_vortex_multiplier);
fprintf('  CL at approach   : %.4f\n',ctx.lift.CL_landing_ref);
if ctx.fidelity == "high"
    fprintf('Drag source        : %s\n',char(ctx.deck.method));
    fprintf('  e_span (LIDRAG)  : %.5f   k = %.6f\n', ...
        ctx.deck.e_span,ctx.deck.k_induced);
else
    fprintf('Drag source        : %s\n',char(ctx.drag_low.source));
    fprintf('  e_osw            : %.5f   k1 = %.6f   CL_minD = %.4f\n', ...
        ctx.drag_low.e_osw,ctx.drag_low.k1_sub,ctx.drag_low.CL_minD);
    fprintf('  CD_min           : %.6f\n',ctx.drag_low.CDmin_sub);
end
fprintf('Engine source      : %s\n',char(ctx.thrust_mode));
fprintf('Static margin      : %+.4f MAC  (%s)\n', ...
    sd.static_margin,char(sd.static_margin_source));
fprintf('  neutral point    : %.3f ft aft of nose\n',sd.neutral_point_ft);
fprintf('  moment reference : %.3f ft aft of nose\n',sd.cg_station_ft);
fprintf('  Cm_alpha         : %+.5f /rad\n',sd.Cm_alpha);
fprintf('  Cm_deltae        : %+.5f /rad (geometric %+.5f)\n', ...
    sd.Cm_deltae,sd.Cm_deltae_geometric);
fprintf('  Cm_q             : %+.5f /rad\n',sd.Cm_q);
fprintf('  Cn_beta / Cl_beta: %+.5f / %+.5f /rad\n',sd.Cn_beta,sd.Cl_beta);
fprintf('  Cl_p / Cn_r      : %+.5f / %+.5f /rad\n',sd.Cl_p,sd.Cn_r);
fprintf('  Cl_deltaa        : %+.5f /rad\n',sd.Cl_deltaa);
fprintf('  Cn_deltar        : %+.5f /rad\n',sd.Cn_deltar);
fprintf('\n');
end
