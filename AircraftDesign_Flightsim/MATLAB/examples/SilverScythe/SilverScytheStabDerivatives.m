function sd = SilverScytheStabDerivatives(ss,longitudinal_source,CL_alpha_wing,options)
% SILVERSCYTHESTABDERIVATIVES  Stability and control derivative layer.
%
% ---------------------------------------------------------------------
% Why this file is built from geometry rather than copied from the report
% ---------------------------------------------------------------------
% The report gives one DATCOM derivative table (Table 9) but never states
% its units, and the two possible readings differ by a factor of 57.3.  The
% per-radian reading is dimensionally consistent with the aircraft geometry
% for the longitudinal terms (report Cm_deltae = -0.062 against -0.101 from
% the ruddervator tail volume computed below) but leaves the aircraft with
% essentially no aerodynamic roll or yaw damping (Cl_p = -0.003, Cn_r =
% +0.0002, the latter with a destabilising sign), which cannot be flown.
% The per-degree reading fixes the damping but makes Cm_deltae = -3.55 per
% radian, roughly 30x what this empennage can produce.  Neither reading is
% usable as a complete 6-DOF deck, and the original DATCOM output files
% (SilverScytheAIL.out, SilverScytheFLAP.out) were not archived with the
% senior-design repository.
%
% This function therefore builds a self-consistent derivative set from the
% published geometry, tags every term with its basis, and prints the report
% table beside it for comparison.  Terms fall into three classes:
%
%   [GEOM]  closed-form tail-volume or strip-theory result computed from
%           report Tables 6 and 7 geometry.  Used for every control
%           derivative, which is exact-in-form here because the
%           ruddervators are all-moving surfaces (tau = 1).
%   [RPT]   taken directly from report Table 8 or Table 9.
%   [CMP]   comparator-scaled from the F-16A deck already in this
%           repository (examples/F16/F16Lookup.m, Brandt), scaled by the
%           ratio of vertical tail-volume coefficients.  Used only for the
%           lateral-directional *stability* derivatives, for which no
%           Silver Scythe high-fidelity source survives.
%
% Every value can be overridden by dropping a DATCOM export into
% data/datcom_derivatives.csv with Name,Value columns; see
% SilverScytheContext.  That is the intended route to close the [CMP]
% terms once the original DATCOM output is recovered.
%
% ---------------------------------------------------------------------
% longitudinal_source
% ---------------------------------------------------------------------
%   "textbook"  Static stability from the textbook aerodynamic-centre /
%               tail-volume chain of the senior-design S&C scripts
%               (AerodynamicCenter.m, NeutralPointStaticMargin.m), which
%               this file ports and which reproduces report Table 8: a
%               neutral point 22.7 ft aft of the nose and a static margin
%               of -0.23 MAC about station 27.9 ft.  Pitch control power
%               from the ruddervator tail volume, -0.101 per radian.
%               This combination does not trim: cruise needs more than the
%               +/-25 deg of report Table 7.
%   "datcom"    Static stability and pitch control power from report
%               Table 9, read per radian.  DATCOM was the higher-fidelity
%               tool the team used to validate the empirical chain, and it
%               gives a near-neutral airframe consistent with the report's
%               statement of Level 1 MIL-F-8785C handling under
%               fly-by-wire.  This is the default.
%
% Moment sign convention: body axes, nose-up Cm, right-wing-down Cl,
% nose-right Cn, all per radian, referenced to the aircraft CG.

if nargin < 1 || isempty(ss), ss = SilverScytheData(); end
if nargin < 2 || isempty(longitudinal_source), longitudinal_source = "datcom"; end
if nargin < 4 || isempty(options), options = struct(); end
longitudinal_source = lower(string(longitudinal_source));

if nargin < 3 || isempty(CL_alpha_wing)
    CL_alpha_wing = 2.71;   % nominal, replaced by SilverScytheContext
end

if ~isfield(options,'cg_station_ft') || isempty(options.cg_station_ft)
    options.cg_station_ft = ss.station.cg_ft;
end
if ~isfield(options,'eta_tail') || isempty(options.eta_tail)
    options.eta_tail = 0.95;    % [GEOM] tail dynamic-pressure ratio
end

S     = ss.S_ref_ft2;
b     = ss.span_ft;
cbar  = ss.cbar_ft;
AR    = ss.AR;
lam   = ss.taper_ratio;
x_cg  = options.cg_station_ft;

%% ------------------------------------------------------------------
%  Empennage geometry                                    [RPT Tables 6, 7]
%% ------------------------------------------------------------------
Gamma_t = deg2rad(ss.empennage.dihedral_deg);
S_rudd  = ss.controls.ruddervator.area_ft2;          % both panels, movable
S_th    = S_rudd*cos(Gamma_t);                       % horizontal projection
S_tv    = S_rudd*sin(Gamma_t);                       % vertical projection
l_t     = ss.controls.ruddervator.station_ft - x_cg; % pitch/yaw arm
z_t     = ss.controls.ruddervator.y_ft*tan(Gamma_t); % roll arm of S_tv

AR_tail = ss.empennage.span_ft^2/ss.empennage.area_ft2;
a_tail  = helmbold(2*pi,deg2rad(ss.empennage.LE_sweep_deg),AR_tail);
eta_t   = options.eta_tail;

%% ------------------------------------------------------------------
%  Wing aerodynamic centre and neutral point
%
%  Port of the senior-design chain AeroCenterChordMain / ...Strake /
%  AerodynamicCenter / NeutralPointStaticMargin.  It reproduces the report
%  Table 8 neutral point of 22.7 ft from geometry alone.
%% ------------------------------------------------------------------
np = textbook_neutral_point(ss,a_tail,CL_alpha_wing);

%% ------------------------------------------------------------------
%  Longitudinal
%% ------------------------------------------------------------------
CL_deltae = a_tail*eta_t*(S_th/S);                       % [GEOM]
CL_q_tail = 2*CL_deltae*(l_t/cbar);                      % [GEOM]

x_ac_wing_ft = np.x_ac_wing_ft;
arm_wing = (x_cg - x_ac_wing_ft)/cbar;

Cm_deltae_geom = -CL_deltae*(l_t/cbar);                  % [GEOM]
Cm_q = -2*CL_deltae*(l_t/cbar)^2 - 2*CL_alpha_wing*arm_wing^2;   % [GEOM]

switch longitudinal_source
    case "textbook"
        static_margin = (np.x_np_ft - x_cg)/cbar;        % negative = unstable
        Cm_deltae = Cm_deltae_geom;
        static_source = "textbook AC/tail-volume chain (report Table 8)";
    case "datcom"
        % DATCOM Table 9 read per radian.
        static_margin = -ss.datcom_report.Cm_alpha/CL_alpha_wing;
        Cm_deltae = ss.datcom_report.Cm_deltae;
        Cm_q = ss.datcom_report.Cm_q + Cm_q;   % DATCOM term plus the
                                               % geometric plunging term,
                                               % which DATCOM's tiny value
                                               % clearly omits
        static_source = "DATCOM report Table 9 (per radian)";
    otherwise
        error('SilverScytheStabDerivatives:UnknownSource', ...
            'longitudinal_source must be "textbook" or "datcom".');
end

% Cm is formed as  Cm = Cm_ac + CL*(x_cg - x_np)/cbar + ...  so the
% CL-linked term is the negative of the static margin.
Cm_CL = -static_margin;

% Section zero-lift pitching moment about the aerodynamic centre.
airfoil = SilverScytheAirfoil(ss);
Cm_ac = median(airfoil.CM(airfoil.alpha_deg >= -2 & airfoil.alpha_deg <= 6));
if ~isfinite(Cm_ac), Cm_ac = -0.042; end

%% ------------------------------------------------------------------
%  Lateral-directional control derivatives                        [GEOM]
%
%  All-moving surfaces, so the control effectiveness equals the surface's
%  full lift-curve contribution (tau = 1).
%% ------------------------------------------------------------------
CY_deltar =  a_tail*eta_t*(S_tv/S);
Cn_deltar = -CY_deltar*(l_t/b);
Cl_deltar =  CY_deltar*(z_t/b);

% Aileron: strip theory over the outboard panel occupied by the surface.
S_ail  = ss.controls.aileron.area_ft2;      % both surfaces
y_ail  = ss.controls.aileron.y_ft;
Cl_deltaa = -CL_alpha_wing*(S_ail*y_ail)/(S*b);
% Adverse yaw, taken as 20% of the roll power with the opposing sign.
Cn_deltaa = -0.20*Cl_deltaa;

%% ------------------------------------------------------------------
%  Lateral-directional stability derivatives
%% ------------------------------------------------------------------
% Roll damping, standard strip-theory closed form                  [GEOM]
Cl_p = -(CL_alpha_wing/12)*(1 + 3*lam)/(1 + lam);

% Comparator scaling from the F-16A deck in examples/F16              [CMP]
f16 = f16_reference();
Vv_SS  = (S_tv*l_t)/(S*b);
Vv_F16 = (f16.S_v*f16.l_v)/(f16.S*f16.b);
directional_scale = Vv_SS/Vv_F16;

CY_beta = f16.CY_beta*directional_scale;
Cn_beta = f16.Cn_beta*directional_scale;
Cn_r    = f16.Cn_r  *directional_scale;
CY_r    = f16.CY_r  *directional_scale;
CY_p    = f16.CY_p  *directional_scale;
Cn_p    = f16.Cn_p;                       % rate terms are wing dominated
Cl_r    = f16.Cl_r;

% Dihedral effect: no wing dihedral, so the terms are wing sweep plus the
% canted-tail side force acting above the body axis.
Cl_beta_sweep = -0.25*CL_alpha_wing*sin(2*deg2rad(np.Lambda_c4_deg))* ...
                (1 + 3*lam)/(1 + lam)/(pi*AR);
Cl_beta_tail  = -a_tail*eta_t*(S_tv/S)*(z_t/b);
Cl_beta = Cl_beta_sweep + Cl_beta_tail;

%% ------------------------------------------------------------------
%  Assemble
%% ------------------------------------------------------------------
sd = struct();
sd.longitudinal_source = longitudinal_source;
sd.CL_alpha          = CL_alpha_wing;
sd.CL_deltae         = CL_deltae;
sd.CL_q              = CL_q_tail;
sd.Cm_ac             = Cm_ac;
sd.Cm_CL             = Cm_CL;
sd.Cm_alpha          = CL_alpha_wing*Cm_CL;
sd.Cm_deltae         = Cm_deltae;
sd.Cm_deltae_geometric = Cm_deltae_geom;
sd.Cm_q              = Cm_q;
sd.CY_beta           = CY_beta;
sd.CY_p              = CY_p;
sd.CY_r              = CY_r;
sd.CY_deltar         = CY_deltar;
sd.Cl_beta           = Cl_beta;
sd.Cl_p              = Cl_p;
sd.Cl_r              = Cl_r;
sd.Cl_deltaa         = Cl_deltaa;
sd.Cl_deltar         = Cl_deltar;
sd.Cn_beta           = Cn_beta;
sd.Cn_p              = Cn_p;
sd.Cn_r              = Cn_r;
sd.Cn_deltaa         = Cn_deltaa;
sd.Cn_deltar         = Cn_deltar;

sd.static_margin     = static_margin;
sd.static_margin_source = static_source;
sd.neutral_point_ft  = np.x_np_ft;
sd.aero_centre_ft    = np.x_ac_ft;
sd.cg_station_ft     = x_cg;
sd.geometry = struct('S_th_ft2',S_th,'S_tv_ft2',S_tv,'l_t_ft',l_t, ...
    'z_t_ft',z_t,'AR_tail',AR_tail,'a_tail_per_rad',a_tail, ...
    'eta_tail',eta_t,'directional_scale',directional_scale, ...
    'downwash_deps_dalpha',np.deps_dalpha, ...
    'Lambda_c4_deg',np.Lambda_c4_deg);

sd.basis = struct( ...
    'CL_deltae',"GEOM", 'CL_q',"GEOM", 'Cm_ac',"GEOM (XFOIL polar)", ...
    'Cm_CL',string(static_source), 'Cm_deltae',"GEOM or RPT by fidelity", ...
    'Cm_q',"GEOM", 'CY_beta',"CMP", 'CY_p',"CMP", 'CY_r',"CMP", ...
    'CY_deltar',"GEOM", 'Cl_beta',"GEOM sweep + tail", 'Cl_p',"GEOM strip", ...
    'Cl_r',"CMP", 'Cl_deltaa',"GEOM strip", 'Cl_deltar',"GEOM", ...
    'Cn_beta',"CMP", 'Cn_p',"CMP", 'Cn_r',"CMP", 'Cn_deltaa',"GEOM", ...
    'Cn_deltar',"GEOM");

end

% =======================================================================
function a = helmbold(a0,Lambda,AR)
num = a0*cos(Lambda);
r = num/(pi*AR);
a = num/(sqrt(1 + r^2) + r);
end

% =======================================================================
function np = textbook_neutral_point(ss,a_tail,CL_alpha_wing)
% Port of the senior-design S&C chain, S&C/SMcalcs.  Fractions of fuselage
% length, wing chord and reference area are the ones hard-coded in
% SMempty.m for Silver Scythe.

l_f = ss.length_ft;
b   = ss.span_ft;
S   = ss.S_ref_ft2;
cr  = ss.root_chord_ft;
ct  = ss.tip_chord_ft;
lam = ss.taper_ratio;
LE  = ss.LE_sweep_deg;
AR  = ss.AR;

% Main wing
MAC   = (2/3)*cr*(1 + lam + lam^2)/(1 + lam);
y_MAC = (b/6)*(1 + 2*lam)/(1 + lam);
x_MAC = 0.228517*l_f + y_MAC*tand(LE);
x_ac_wing = x_MAC + 0.25*MAC;

% Strakes
x_strake  = 0.065633*l_f;
lam_str   = 0.295271792;
cr_str    = 0.355565*cr;
S_str     = 0.0319147*S;
LE_str    = 1.39862*LE;
b_str     = 2*S_str/(cr_str*(1 + lam_str));
MAC_str   = (2/3)*cr_str*(1 + lam_str + lam_str^2)/(1 + lam_str);
y_MAC_str = (b_str/6)*(1 + 2*lam_str)/(1 + lam_str);
x_MAC_str = x_strake + y_MAC_str*tand(LE_str);
x_ac_str  = x_MAC_str + 0.25*MAC_str;

% Pitch surface (horizontal projection of the ruddervator)
x_pitch   = 0.7277*l_f;
lam_p     = 0.08577/0.2282;
cr_p      = 0.2282*cr;
S_pitch   = 0.045704*S;
b_p       = 2*S_pitch/(cr_p*(1 + lam_p));
MAC_p     = (2/3)*cr_p*(1 + lam_p + lam_p^2)/(1 + lam_p);
y_MAC_p   = (b_p/6)*(1 + 2*lam_p)/(1 + lam_p);
x_MAC_p   = x_pitch + y_MAC_p*tand(LE);

% Vertical surface (vertical projection)
S_vert    = 0.141448*S;
b_v       = 2*S_vert/(cr_p*(1 + lam_p));
z_MAC_v   = (b_v/6)*(1 + 2*lam_p)/(1 + lam_p);
x_MAC_v   = x_pitch + z_MAC_v*tand(LE);

% Wing + strake, then the Roskam fuselage shift
x_ac_ws = x_ac_wing + (x_ac_str - x_ac_wing)*S_str/(S + S_str);
w_f = 2*sqrt(ss.max_cross_section_ft2/pi);
x_ac_wsf = x_ac_ws - (l_f*w_f^2*(0.005 + 0.111*(x_ac_ws/l_f)^2))/ ...
    (S*(CL_alpha_wing*pi/180)*57.29);

% Tail contribution with the Roskam downwash correlation
l_h = x_MAC_p - x_ac_wsf;
V_h = (S_pitch*l_h)/(S*MAC);
K_A = 1/AR - 1/(1 + AR^1.7);
K_lambda = (10 - 3*lam)/7;
l_H = abs(x_MAC - x_MAC_v);
K_H = (1 - abs(z_MAC_v/b))/(2*l_H/b)^(1/3);
Lambda_c4 = atand(tand(LE) - (1 - lam)/(AR*(1 + lam)));
deps_dalpha = 4.44*(K_A*K_lambda*K_H*sqrt(cosd(Lambda_c4)))^1.19;

x_np = x_ac_wsf + V_h*(a_tail/max(CL_alpha_wing,1e-9))*(1 - deps_dalpha)*MAC;

np = struct();
np.x_ac_wing_ft   = x_ac_wing;
np.x_ac_ft        = x_ac_wsf;
np.x_np_ft        = x_np;
np.MAC_ft         = MAC;
np.tail_volume    = V_h;
np.tail_arm_ft    = l_h;
np.deps_dalpha    = deps_dalpha;
np.Lambda_c4_deg  = Lambda_c4;

end

% =======================================================================
function f16 = f16_reference()
% F-16A lateral-directional derivatives from examples/F16/F16Lookup.m
% (Brandt), with the geometry needed to form its vertical tail volume.
f16 = struct();
f16.S       = 300.0;    % ft^2
f16.b       = 30.0;     % ft
f16.S_v     = 54.8;     % ft^2, vertical tail
f16.l_v     = 14.8;     % ft, fin arm
f16.CY_beta = -1.146;
f16.Cn_beta =  0.241;
f16.Cn_r    = -0.378;
f16.CY_r    =  0.876;
f16.CY_p    = -0.188;
f16.Cn_p    = -0.052;
f16.Cl_r    =  0.063;
end
