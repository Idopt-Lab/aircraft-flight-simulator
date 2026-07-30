function C = F16ABrandtLookup(state_vec,control_vec,geometry)
% F16ABRANDTLOOKUP
% Formula-based Brandt F-16A aerodynamic model.
%
% Input convention:
%   state_vec   = [X Y Z u v w phi theta psi p q r]'
%   control_vec = [aileron stabilator rudder ...]'
%
% Model:
%   CL and Cm use the Brandt local stability/control derivatives.
%   CD uses the Brandt Mach-dependent quadratic drag polar.
%   CY, Cl and Cn use Brandt geometry-derived lateral derivatives.
%
% This is a local, low-to-moderate-alpha conceptual model. It is not the
% high-alpha nonlinear F-16 aerodynamic database.

if ~isnumeric(state_vec) || numel(state_vec) < 12
    error("F16ABrandtLookup:InvalidState", "state_vec must be numeric and contain at least 12 elements.");
end

if ~isnumeric(control_vec)
    error("F16ABrandtLookup:InvalidControl", "control_vec must be numeric.");
end

persistent data
if isempty(data)
    data = BrandtF16AData();
end

vel = state_vec(4:6);
omega = state_vec(10:12);

altitude_m = max(-state_vec(3),0);
V_mps = max(norm(vel),1e-6);
alpha = atan2(vel(3),max(abs(vel(1)),1e-9));
beta = atan2(vel(2),max(hypot(vel(1),vel(3)),1e-9));

altitude_ft = altitude_m*3.280839895;
atm = data.functions.atmosphere(altitude_ft,data.atmos);
V_fps = V_mps*3.280839895;
Mach = V_fps/max(atm.speed_of_sound_fps,1e-9);

da = 0;
de = 0;
dr = 0;

if numel(control_vec) >= 1, da = control_vec(1); end
if numel(control_vec) >= 2, de = control_vec(2); end
if numel(control_vec) >= 3, dr = control_vec(3); end

S = geometryValue(geometry,"wing_area",data.reference.S_m2);
b = geometryValue(geometry,"wing_span",data.reference.b_m);
c = geometryValue(geometry,"mean_aerodynamic_chord",data.reference.c_m);

if ~(isfinite(S) && S > 0 && isfinite(b) && b > 0 && isfinite(c) && c > 0)
    error("F16ABrandtLookup:InvalidGeometry", "Positive wing area, span and mean chord are required.");
end

validity = data.validity;
controls = data.controls;

envelope_names = ["alpha_min";"alpha_max";"beta_min";"beta_max"; "aileron_min";"aileron_max"; "stabilator_min";"stabilator_max"; "rudder_min";"rudder_max"; "Mach_min";"Mach_max"];

envelope_c = [ ...
    validity.alpha_min_rad-alpha; ...
    alpha-validity.alpha_max_rad; ...
    -validity.beta_limit_rad-beta; ...
    beta-validity.beta_limit_rad; ...
    -controls.aileron_limit_rad-da; ...
    da-controls.aileron_limit_rad; ...
    -controls.stabilator_limit_rad-de; ...
    de-controls.stabilator_limit_rad; ...
    -controls.rudder_limit_rad-dr; ...
    dr-controls.rudder_limit_rad; ...
    validity.mach_min-Mach; ...
    Mach-validity.mach_max];

envelope_violation = max([0;envelope_c]);

% Clamping provides numerically safe coefficients only. The validity flag
% must be used to reject out-of-envelope trim/performance points.
alpha_q = clampValue(alpha,validity.alpha_min_rad,validity.alpha_max_rad);
beta_q = clampValue(beta,-validity.beta_limit_rad,validity.beta_limit_rad);
da_q = clampValue(da,-controls.aileron_limit_rad,controls.aileron_limit_rad);
de_q = clampValue(de,-controls.stabilator_limit_rad,controls.stabilator_limit_rad);
dr_q = clampValue(dr,-controls.rudder_limit_rad,controls.rudder_limit_rad);
Mach_q = clampValue(Mach,validity.mach_min,validity.mach_max);

p_hat = omega(1)*b/(2*V_mps);
q_hat = omega(2)*c/(2*V_mps);
r_hat = omega(3)*b/(2*V_mps);

dalpha = alpha_q-data.trim.alpha_ref_rad;
dde = de_q-data.trim.pitch_control_ref_rad;

CL = data.trim.CL_ref +data.stability.CLa*dalpha +data.stability.CLde*dde +data.stability.CLq*q_hat;

Cm = data.trim.Cm_ref +data.stability.Cma_takeoff*dalpha +data.stability.Cmde*dde +data.stability.Cmq*q_hat;

CY = data.stability.CYb*beta_q +data.stability.CYp*p_hat +data.stability.CYr*r_hat +data.stability.CYda*da_q +data.stability.CYdr*dr_q;

Cl = data.stability.Clb*beta_q +data.stability.Clp*p_hat +data.stability.Clr*r_hat +data.stability.Clda*da_q +data.stability.Cldr*dr_q;

Cn = data.stability.Cnb*beta_q +data.stability.Cnp*p_hat +data.stability.Cnr*r_hat +data.stability.Cnda*da_q +data.stability.Cndr*dr_q;

[CD0,k1,k2,CD_base] = data.functions.aero(Mach_q,CL,data.aero.CDx,data.aero);

CD_control = data.controls.drag.stabilator_quadratic*dde^2 +data.controls.drag.aileron_quadratic*da_q^2 +data.controls.drag.rudder_quadratic*dr_q^2;

CD = max(CD_base+CD_control,1e-4);

CL_c = [validity.CL_min-CL; CL-validity.CL_max_clean];

constraint_names = [envelope_names;"CL_min";"CL_max_clean"];
constraint_c = [envelope_c;CL_c];
constraint_violation = max([0;constraint_c]);

valid = constraint_violation <= 1e-10 && all(isfinite([CL,CD,CY,Cl,Cm,Cn]));

lookup_clamped = any(abs([alpha_q-alpha,beta_q-beta,da_q-da,de_q-de,dr_q-dr, Mach_q-Mach]) > 1e-12);

takeoff_params = struct( ...
    "mu_ground",0.02, ...
    "mu_rolling",0.02, ...
    "mu_braking",0.35, ...
    "safety_factor",1.15, ...
    "CLmax_takeoff",data.aircraft.CLmax_takeoff, ...
    "V2_to_Vs_ratio",1.20, ...
    "VR_to_Vs_ratio",1.10, ...
    "V1_to_VR_ratio",0.92, ...
    "CD0_takeoff",data.aero.CD0_sub+0.02, ...
    "K_takeoff",data.aero.k1_sub, ...
    "CLalpha_takeoff",data.stability.CLa, ...
    "screen_height_takeoff_ft",35, ...
    "rotation_alpha_deg",10, ...
    "initial_climb_angle_deg",10, ...
    "reaction_time_s",1.0);

landing_params = struct( ...
    "mu_braking",0.70, ...
    "approach_angle_deg",3, ...
    "safety_factor",1.67, ...
    "CLmax_landing",data.aircraft.CLmax_landing, ...
    "CD0_landing",data.aero.CD0_sub+0.04, ...
    "CL_landing_touchdown",0.05, ...
    "screen_height_landing_ft",50, ...
    "flare_height_m",5, ...
    "Vapp_to_Vs_ratio",1.30, ...
    "Vtd_to_Vapp_ratio",0.88, ...
    "idle_throttle",0.05, ...
    "use_idle_thrust",true);

C = struct( ...
    "CL",CL, ...
    "CD",CD, ...
    "CY",CY, ...
    "Cl",Cl, ...
    "Cm",Cm, ...
    "Cn",Cn, ...
    "CD0",CD0, ...
    "k1",k1, ...
    "k2",k2, ...
    "CD_base",CD_base, ...
    "CD_control",CD_control, ...
    "alpha_rad",alpha, ...
    "beta_rad",beta, ...
    "mach",Mach, ...
    "airspeed_mps",V_mps, ...
    "altitude_m",altitude_m, ...
    "alpha_used_rad",alpha_q, ...
    "beta_used_rad",beta_q, ...
    "aileron_used_rad",da_q, ...
    "stabilator_used_rad",de_q, ...
    "rudder_used_rad",dr_q, ...
    "valid",valid, ...
    "lookup_clamped",lookup_clamped, ...
    "constraint_names",constraint_names, ...
    "constraint_c",constraint_c, ...
    "constraint_violation",constraint_violation, ...
    "reference_state",struct( ...
        "CL_ref",data.trim.CL_ref, ...
        "Cm_ref",data.trim.Cm_ref, ...
        "alpha_ref_rad",data.trim.alpha_ref_rad, ...
        "alpha_zero_lift_rad",data.trim.alpha_zero_lift_rad, ...
        "pitch_control_ref_rad",data.trim.pitch_control_ref_rad), ...
    "derivatives",data.stability, ...
    "model_type",data.validity.model_type, ...
    "approximations",data.approximations, ...
    "takeoff_params",takeoff_params, ...
    "landing_params",landing_params);
end

function value = geometryValue(geometry,name,defaultValue)
value = defaultValue;

if isstruct(geometry)
    if isfield(geometry,name)
        candidate = geometry.(name);
        if isnumeric(candidate) && isscalar(candidate)
            value = candidate;
        end
    end
elseif isobject(geometry) && isprop(geometry,name)
    candidate = geometry.(name);
    if isnumeric(candidate) && isscalar(candidate)
        value = candidate;
    end
end
end

function y = clampValue(x,lowerBound,upperBound)
y = min(max(x,lowerBound),upperBound);
end

function data = BrandtF16AData(weight_fraction)
% BRANDTF16ADATA
% Formula-based Brandt F-16A source and derived data.
%
% This file separates:
%   1) direct workbook/source inputs,
%   2) equations reproduced from the Brandt calculator,
%   3) explicit approximations needed for a controlled 6-DOF model.
%
% The runtime lookup uses SI units. Workbook validation fields retain the
% original imperial quantities where needed for direct regression checks.

if nargin < 1 || isempty(weight_fraction)
    weight_fraction = 0.899666962714079;
end

validateattributes(weight_fraction,{'numeric'}, {'real','finite','scalar','positive','<=',1});

%% Direct source inputs

aircraft = struct();
aircraft.name = "Brandt F-16A";
aircraft.Wto_lbf = 31377.0;
aircraft.Wlanding_lbf = 20680.7005781593;
aircraft.Sref_ft2 = 300.0;
aircraft.xcg_takeoff_ft = 26.1925091690947;
aircraft.xcg_landing_ft = 26.1369136582888;
aircraft.CLmax_clean = 0.9840156989811455;
aircraft.CLmax_takeoff = 1.2756651310132752;
aircraft.CLmax_landing = 1.4259087778177055;
aircraft.CLmax_air_to_air = 1.1283380014983801;
aircraft.maximum_mach = 2.0;
aircraft.analysis_weight_fraction = weight_fraction;
aircraft.analysis_weight_lbf = weight_fraction*aircraft.Wto_lbf;

geometry = struct();

geometry.wing = struct('S_ft2',300.0, 'b_ft',30.0, 'cr_ft',16.2932790224033, 'ct_ft',3.70672097759674, 'AR',3.0, 'sweep_LE_deg',40.0, 'dihedral_deg',0.0, 'x_LE_root_ft',17.7858333333333, 'z_LE_root_ft',0.0);

geometry.pitch = struct('S_ft2',108.0, 'b_ft',18.0, 'cr_ft',9.77596741344195, 'ct_ft',2.22403258655804, 'AR',3.0, 'sweep_LE_deg',40.0, 'dihedral_deg',-10.0, 'x_LE_root_ft',36.0, 'z_LE_root_ft',0.0);

geometry.strake = struct('S_ft2',20.0, 'b_ft',5.47722557505166, 'cr_ft',7.30296743340222, 'ct_ft',0.0, 'AR',1.5, 'sweep_LE_deg',74.0, 'dihedral_deg',0.0, 'x_LE_root_ft',12.0, 'z_LE_root_ft',0.0);

geometry.vertical = struct('S_ft2',60.0, 'b_ft',9.79795897113271, 'cr_ft',8.16496580927726, 'ct_ft',4.08248290463863, 'AR',3.2, 'sweep_LE_deg',30.0, 'tilt_deg',0.0, 'x_LE_root_ft',36.0);

geometry.fuselage = struct( ...
    'length_ft',46.5, ...
    'input_max_width_ft',7.0, ...
    'input_max_height_ft',5.0, ...
    'station_max_width_ft',7.5, ...
    'station_max_height_ft',7.0, ...
    'volume_integral_ft3',1106.30614811101, ...
    'side_area_at_wing_ft2',18.199022816385, ...
    'sideforce_factor',1.0);

geometry.sweep_tmax_deg = -0.000243426089565037;

aero = struct();
aero.Mcrit = 0.8727371242684938;
aero.Mwave = 1.0547492038954984;
aero.M_k2_zero = 1.5;
aero.CD0_sub = 0.016995734443133757;
aero.CD0_wave = 0.0456687347713526;
aero.k1_sub = 0.1160314457303649;
aero.k1_wave = 0.2415276545972613;
aero.k2_sub = -0.00630192501971521;
aero.CDx = 0.0;

engine = struct();
engine.number = 1;
engine.Tsl_dry_lbf = 15000.0;
engine.Tsl_ab_lbf = 23770.0;
engine.TSFCsl_dry = 0.7;
engine.TSFCsl_ab = 2.2;
engine.TR = 1.0;
engine.dry = struct('F1',0.300,'E',1.000,'F2',1.700);
engine.ab = struct('F1',0.100,'E',0.500,'F2',2.200);
engine.tsfc = struct('mach_factor',0.350, 'dry_Mref',0.000, 'ab_Mref',0.400, 'mach_exponent',1.000, 'theta_exponent',0.500);
engine.installed_weight_lb = 0.199*engine.Tsl_ab_lbf;
engine.diameter_ft = sqrt(engine.Tsl_ab_lbf/1900);
engine.length_ft = 4.5*engine.diameter_ft;

atmos = struct();
atmos.h_tropopause_ft = 36152.0;
atmos.Tsl_R = 518.69;
atmos.Tstrat_R = 389.99;
atmos.lapse_R_per_ft = -0.00356;
atmos.rho_sl_slug_ft3 = 0.002377;
atmos.rho_trop_slug_ft3 = 0.000706;
atmos.gamma = 1.4;
atmos.R_ft_lbf_slug_R = 1716.0;
atmos.g_ft_s2 = 32.2;
atmos.Psl_psf = 2116.2;

%% Derived geometry

wing = surfaceGeometry(geometry.wing,false);
pitch = surfaceGeometry(geometry.pitch,false);
strake = surfaceGeometry(geometry.strake,false);
vertical = surfaceGeometry(geometry.vertical,true);

wing.e_stability = brandtSpanEfficiency(geometry.wing.AR,geometry.wing.sweep_LE_deg);
pitch.e_stability = brandtSpanEfficiency(geometry.pitch.AR,geometry.pitch.sweep_LE_deg);
vertical.e_stability = brandtSpanEfficiency(geometry.vertical.AR,geometry.vertical.sweep_LE_deg);

wing.CLa_rad = finiteWingCLa(geometry.wing.AR,wing.e_stability);
pitch.CLa_rad = finiteWingCLa(geometry.pitch.AR,pitch.e_stability);
vertical.CLa_rad = finiteWingCLa(geometry.vertical.AR,vertical.e_stability);

wing.e_performance = performanceSpanEfficiency(geometry.wing.AR,geometry.wing.sweep_LE_deg, geometry.sweep_tmax_deg);
pitch.e_performance = performanceSpanEfficiency(geometry.pitch.AR,geometry.pitch.sweep_LE_deg, geometry.sweep_tmax_deg);

wing.CLa_per_deg = 0.1/(1+5.73/(3.1416*wing.e_performance*geometry.wing.AR));
pitch.CLa_per_deg = 0.1/(1+5.73/(3.1416*pitch.e_performance*geometry.pitch.AR));

wing_strake = struct();
wing_strake.CLa_per_deg = wing.CLa_per_deg *(geometry.wing.S_ft2+geometry.strake.S_ft2) /geometry.wing.S_ft2;

%% Aerodynamic center, downwash, neutral point and static margin

stability = struct();

stability.xac_wing_strake_ft = wing.x_ac_ft +0.5*(strake.x_ac_ft-wing.x_ac_ft) *geometry.strake.S_ft2 /(geometry.wing.S_ft2+geometry.strake.S_ft2);

stability.xac_wing_strake_fuselage_ft = ...
    stability.xac_wing_strake_ft ...
    -(geometry.fuselage.length_ft ...
      *geometry.fuselage.input_max_width_ft^2 ...
      *(0.005+0.111 ...
      *(stability.xac_wing_strake_ft ...
      /geometry.fuselage.length_ft)^2)) ...
    /(geometry.wing.S_ft2 ...
      *wing_strake.CLa_per_deg*57.29578);

stability.xbar_ac = (stability.xac_wing_strake_fuselage_ft-wing.x_MAC_ft) /wing.MAC_ft;

downwash_distance_ft = geometry.pitch.x_LE_root_ft -0.25*geometry.wing.ct_ft -geometry.wing.x_LE_root_ft -0.25*geometry.wing.cr_ft;

if abs(downwash_distance_ft) <= eps
    error("BrandtF16AData:InvalidDownwashGeometry", "The Brandt downwash distance is zero.");
end

stability.downwash_gradient = min( ...
    sign(downwash_distance_ft) ...
    *21*wing.CLa_per_deg/sqrt(geometry.wing.AR) ...
    *(((geometry.wing.cr_ft+geometry.wing.ct_ft)/2) ...
      /abs(downwash_distance_ft))^0.25 ...
    *(10-3*geometry.wing.ct_ft/geometry.wing.cr_ft)/7 ...
    *(1-(geometry.pitch.z_LE_root_ft ...
      -geometry.wing.z_LE_root_ft)/geometry.wing.b_ft), ...
    1);

stability.tail_arm_ft = pitch.x_ac_ft-stability.xac_wing_strake_fuselage_ft;

stability.horizontal_tail_volume = geometry.pitch.S_ft2/geometry.wing.S_ft2 *stability.tail_arm_ft/wing.MAC_ft;

stability.xbar_neutral_point = stability.xbar_ac +stability.horizontal_tail_volume *pitch.CLa_per_deg/wing_strake.CLa_per_deg *(1-stability.downwash_gradient);

stability.x_neutral_point_ft = wing.x_MAC_ft+stability.xbar_neutral_point*wing.MAC_ft;

stability.static_margin_takeoff = (stability.x_neutral_point_ft-aircraft.xcg_takeoff_ft) /wing.MAC_ft;

stability.static_margin_landing = (stability.x_neutral_point_ft-aircraft.xcg_landing_ft) /wing.MAC_ft;

%% Longitudinal reference and derivatives

stability.CLa = wing.CLa_rad +pitch.CLa_rad*geometry.pitch.S_ft2/geometry.wing.S_ft2;

stability.Cma_takeoff = -stability.CLa*stability.static_margin_takeoff;
stability.Cma_landing = -stability.CLa*stability.static_margin_landing;

trim = struct();
trim.reference_density_slug_ft3 = atmos.rho_sl_slug_ft3;
trim.reference_speed_fps = 1.2*sqrt(2*aircraft.Wto_lbf /(trim.reference_density_slug_ft3 *geometry.wing.S_ft2*aircraft.CLmax_takeoff));
trim.reference_q_psf = 0.5*trim.reference_density_slug_ft3 *trim.reference_speed_fps^2;

xbar_wing_ac_absolute = wing.x_ac_ft/wing.MAC_ft;
xbar_pitch_ac_absolute = pitch.x_ac_ft/wing.MAC_ft;
xbar_takeoff_cg_absolute = aircraft.xcg_takeoff_ft/wing.MAC_ft;

trim.reference_wing_lift_lbf = aircraft.Wlanding_lbf *(xbar_takeoff_cg_absolute-xbar_pitch_ac_absolute) /((xbar_takeoff_cg_absolute-xbar_pitch_ac_absolute) -(xbar_takeoff_cg_absolute-xbar_wing_ac_absolute));

trim.reference_CL_wing = trim.reference_wing_lift_lbf /(trim.reference_q_psf*geometry.wing.S_ft2);

trim.reference_CL_wing_normal = trim.reference_CL_wing/cosd(geometry.wing.dihedral_deg);

trim.reference_pitch_lift_lbf = aircraft.Wlanding_lbf-trim.reference_wing_lift_lbf;

% These two expressions reproduce the Brandt workbook exactly.
trim.reference_CL_pitch = trim.reference_pitch_lift_lbf /(trim.reference_q_psf*(geometry.pitch.S_ft2/144));

trim.reference_CL_pitch_normal = (trim.reference_pitch_lift_lbf /(trim.reference_q_psf*(geometry.pitch.S_ft2/10))) /cosd(geometry.pitch.dihedral_deg);

trim.CL_ref = trim.reference_CL_wing_normal +trim.reference_CL_pitch_normal *geometry.pitch.S_ft2/geometry.wing.S_ft2;

trim.alpha_ref_rad = trim.reference_CL_wing_normal/wing.CLa_rad;

trim.alpha_zero_lift_rad = trim.alpha_ref_rad-trim.CL_ref/stability.CLa;

trim.pitch_surface_alpha_ref_rad = trim.reference_CL_pitch_normal/pitch.CLa_rad;

trim.pitch_control_ref_rad = trim.pitch_surface_alpha_ref_rad-trim.alpha_ref_rad;

trim.Cm_ref = 0.0;

%% Lateral-directional and short-period derivatives

stability.CYb_vertical = -vertical.CLa_rad*geometry.vertical.S_ft2/geometry.wing.S_ft2;
stability.CYb_wing_dihedral = -0.00573*abs(geometry.wing.dihedral_deg);
stability.CYb_fuselage = -2*geometry.fuselage.sideforce_factor *geometry.fuselage.side_area_at_wing_ft2/geometry.wing.S_ft2;
stability.CYb = stability.CYb_vertical +stability.CYb_wing_dihedral +stability.CYb_fuselage;

stability.Cnb_vertical = vertical.CLa_rad *(geometry.vertical.S_ft2/geometry.wing.S_ft2) *(vertical.x_ac_ft-aircraft.xcg_takeoff_ft) /geometry.wing.b_ft;
stability.Cnb_wing_dihedral = wing.CLa_rad*sind(geometry.wing.dihedral_deg)^2 *(wing.x_ac_ft-aircraft.xcg_takeoff_ft) /geometry.wing.b_ft;
stability.Cnb_pitch_dihedral = pitch.CLa_rad*sind(geometry.pitch.dihedral_deg)^2 *(geometry.pitch.S_ft2/geometry.wing.S_ft2) *(pitch.x_ac_ft-aircraft.xcg_takeoff_ft) /geometry.wing.b_ft;
stability.Cnb_fuselage = -1.3*geometry.fuselage.volume_integral_ft3 *(geometry.fuselage.station_max_width_ft /geometry.fuselage.station_max_height_ft) /(geometry.wing.S_ft2*geometry.wing.b_ft);
stability.Cnb = stability.Cnb_vertical +stability.Cnb_wing_dihedral +stability.Cnb_pitch_dihedral +stability.Cnb_fuselage;

fuselage_average_diameter_ft = (geometry.fuselage.input_max_width_ft +geometry.fuselage.input_max_height_ft)/0.7854;

stability.Clb_wing_dihedral = -0.5*wing.CLa_rad/3 *geometry.wing.dihedral_deg*pi/180 *(1+2*geometry.wing.ct_ft/geometry.wing.cr_ft) /(1+geometry.wing.ct_ft/geometry.wing.cr_ft);

stability.Clb_vertical = -vertical.CLa_rad*geometry.vertical.S_ft2 /geometry.wing.S_ft2/geometry.wing.b_ft *vertical.z_ac_ft *(1+tand(geometry.vertical.tilt_deg));

stability.Clb_wing_sweep = -2*trim.reference_CL_wing *wing.y_MAC_ft/geometry.wing.b_ft *sind(2*geometry.wing.sweep_LE_deg);

stability.Clb_wing_position = -0.042*sqrt(geometry.wing.AR) *geometry.wing.z_LE_root_ft/geometry.wing.b_ft *fuselage_average_diameter_ft/geometry.wing.b_ft *57.29578;

stability.Clb_pitch_dihedral = ...
    (-0.5*pitch.CLa_rad/3 ...
    *deg2rad(geometry.pitch.dihedral_deg) ...
    *(1+2*geometry.pitch.ct_ft/geometry.pitch.cr_ft) ...
    /(1+geometry.pitch.ct_ft/geometry.pitch.cr_ft)) ...
    *(geometry.pitch.S_ft2*geometry.pitch.b_ft ...
    /(geometry.wing.S_ft2*geometry.wing.b_ft));

stability.Clb_pitch_position = ...
    -0.042*sqrt(geometry.pitch.AR) ...
    *geometry.pitch.z_LE_root_ft/geometry.pitch.b_ft ...
    *fuselage_average_diameter_ft/geometry.pitch.b_ft ...
    *57.29578 ...
    *geometry.pitch.S_ft2*geometry.pitch.b_ft ...
    /(geometry.wing.S_ft2*geometry.wing.b_ft);

stability.Clb = stability.Clb_wing_dihedral +stability.Clb_vertical +stability.Clb_wing_sweep +stability.Clb_wing_position +stability.Clb_pitch_dihedral +stability.Clb_pitch_position;

stability.CYr = 2*vertical.CLa_rad *(geometry.vertical.S_ft2/geometry.wing.S_ft2) *(vertical.x_ac_ft-aircraft.xcg_takeoff_ft) /geometry.wing.b_ft;

stability.Cnr = -2*vertical.CLa_rad *(geometry.vertical.S_ft2/geometry.wing.S_ft2) *((vertical.x_ac_ft-aircraft.xcg_takeoff_ft) /geometry.wing.b_ft)^2;

stability.Clr = 2*vertical.CLa_rad *(geometry.vertical.S_ft2/geometry.wing.S_ft2) *(vertical.x_ac_ft-aircraft.xcg_takeoff_ft) *vertical.z_ac_ft/geometry.wing.b_ft^2;

stability.CYp = -2*vertical.CLa_rad *(geometry.vertical.S_ft2/geometry.wing.S_ft2) *vertical.z_ac_ft/geometry.wing.b_ft;

stability.Cnp = 2*vertical.CLa_rad *(geometry.vertical.S_ft2/geometry.wing.S_ft2) *vertical.z_ac_ft *(vertical.x_ac_ft-aircraft.xcg_takeoff_ft) /geometry.wing.b_ft^2;

stability.Clp = -2*vertical.CLa_rad *(geometry.vertical.S_ft2/geometry.wing.S_ft2) *(vertical.z_ac_ft/geometry.wing.b_ft)^2;

stability.pitch_control_effectiveness = 1.0;
stability.Cmq = -2*pitch.CLa_rad *stability.horizontal_tail_volume *stability.tail_arm_ft/wing.MAC_ft;
stability.Cmadot = stability.Cmq*stability.downwash_gradient;
stability.CLde = pitch.CLa_rad *(geometry.pitch.S_ft2/geometry.wing.S_ft2) *stability.pitch_control_effectiveness;
stability.Cmde = -pitch.CLa_rad *stability.horizontal_tail_volume *stability.pitch_control_effectiveness;

%% Explicit controlled-6DOF approximations

controls = struct();
controls.aileron_limit_rad = deg2rad(21.5);
controls.stabilator_limit_rad = deg2rad(25);
controls.rudder_limit_rad = deg2rad(30);

controls.aileron = struct('inboard_semispan_fraction',0.55, 'outboard_semispan_fraction',0.90, 'effectiveness',0.75);

y1 = controls.aileron.inboard_semispan_fraction *geometry.wing.b_ft/2;
y2 = controls.aileron.outboard_semispan_fraction *geometry.wing.b_ft/2;
chord_gradient = (geometry.wing.cr_ft-geometry.wing.ct_ft) *2/geometry.wing.b_ft;
integral_yc = 0.5*geometry.wing.cr_ft*(y2^2-y1^2) -(chord_gradient/3)*(y2^3-y1^3);

stability.Clda = 2*wing.CLa_rad*controls.aileron.effectiveness *integral_yc/(geometry.wing.S_ft2*geometry.wing.b_ft);
stability.CYda = 0.0;
stability.Cnda = 0.0;

controls.rudder = struct('chord_effectiveness',0.75, 'dynamic_pressure_ratio',0.90);

stability.CYdr = controls.rudder.dynamic_pressure_ratio *controls.rudder.chord_effectiveness *vertical.CLa_rad *geometry.vertical.S_ft2/geometry.wing.S_ft2;

vertical_arm_ft = vertical.x_ac_ft-aircraft.xcg_takeoff_ft;

stability.Cndr = -stability.CYdr*vertical_arm_ft/geometry.wing.b_ft;
stability.Cldr = stability.CYdr*vertical.z_ac_ft/geometry.wing.b_ft;

% Approximate lift-rate derivative from the horizontal-tail contribution.
stability.CLq = 2*pitch.CLa_rad*stability.horizontal_tail_volume;

controls.drag = struct('stabilator_quadratic',0.020, 'aileron_quadratic',0.010, 'rudder_quadratic',0.015);

%% Validity envelope

validity = struct();
validity.alpha_min_rad = deg2rad(-5);
validity.alpha_max_rad = deg2rad(20);
validity.beta_limit_rad = deg2rad(15);
validity.mach_min = 0;
validity.mach_max = aircraft.maximum_mach;
validity.CL_min = -0.80;
validity.CL_max_clean = aircraft.CLmax_clean;
validity.model_type = "local_linear_stability_plus_Brandt_Mach_drag_polar";
validity.high_alpha_model = false;

%% SI conversion and mass properties

ft_to_m = 0.3048;
ft2_to_m2 = ft_to_m^2;
lbf_to_N = 4.4482216152605;

mass = struct();
mass.analysis_mass_kg = aircraft.analysis_weight_lbf*lbf_to_N/9.80665;
mass.reference_inertia_weight_lbf = 20500.0;
mass.inertia_scale = aircraft.analysis_weight_lbf/mass.reference_inertia_weight_lbf;

slug_ft2_to_kg_m2 = 1.3558179483314;
Ixx_ref = 9496*slug_ft2_to_kg_m2;
Iyy_ref = 55814*slug_ft2_to_kg_m2;
Izz_ref = 63100*slug_ft2_to_kg_m2;
Ixz_ref = -982*slug_ft2_to_kg_m2;

mass.I_cg_kgm2 = mass.inertia_scale*[Ixx_ref,0,-Ixz_ref; 0,Iyy_ref,0; -Ixz_ref,0,Izz_ref];

mass.cg_body_m = [0;0;0];
mass.inertia_source = "benchmark_F16_tensor_scaled_by_Brandt_analysis_weight";

reference = struct();
reference.S_m2 = geometry.wing.S_ft2*ft2_to_m2;
reference.b_m = geometry.wing.b_ft*ft_to_m;
reference.c_m = wing.MAC_ft*ft_to_m;

%% Cached regression values from the calculator/workbook

expected_stability = struct( ...
    "Wing_MAC_ft",11.3201786951274, ...
    "Wing_xac_ft",25.5889532142962, ...
    "Pitch_MAC_ft",6.79210721707642, ...
    "Pitch_xac_ft",40.6818719285777, ...
    "Vertical_MAC_ft",6.3505289627712, ...
    "Vertical_xac_ft",40.1017896849116, ...
    "Wing_e",0.722735839202051, ...
    "Pitch_e",0.722735839202051, ...
    "Vertical_e",0.749125335425179, ...
    "Wing_CLa",3.26837125071049, ...
    "Pitch_CLa",3.26837125071049, ...
    "Vertical_CLa",3.42537407573747, ...
    "Xac_WingStrake_ft",25.3018209420578, ...
    "Xac_WingStrakeFuselage_ft",25.2151835119495, ...
    "Downwash",0.817518542514647, ...
    "TailArm_ft",15.4666884166283, ...
    "HorizontalTailVolume",0.491865718726053, ...
    "NeutralPoint_ft",26.1677380595486, ...
    "SM_Takeoff",-0.00218822601773327, ...
    "SM_Landing",0.00272296066077779, ...
    "CLa_Total",4.44498490096627, ...
    "Cma_Takeoff",0.00972663160872596, ...
    "CYb",-0.806401633923394, ...
    "Clb",-0.293210293100817, ...
    "Cnb",0.163551599852826, ...
    "CYr",0.635259851880531, ...
    "Cnr",-0.294533582674755, ...
    "Clr",0.0922111105886412, ...
    "CYp",-0.19888399795471, ...
    "Cnp",0.0922111105886412, ...
    "Clp",-0.0288689963255569, ...
    "Cmq",-4.39290676892833, ...
    "Cmadot",-3.59128273913702, ...
    "CLde",1.17661365025578, ...
    "Cmde",-1.60759977429429, ...
    "CLref",0.797630894781094, ...
    "AlphaRef_rad",0.171503383561248, ...
    "AlphaZeroLift_rad",-0.00794174674902601, ...
    "PitchControlRef_rad",0.0300021543852092);

calculated_stability = struct( ...
    "Wing_MAC_ft",wing.MAC_ft, ...
    "Wing_xac_ft",wing.x_ac_ft, ...
    "Pitch_MAC_ft",pitch.MAC_ft, ...
    "Pitch_xac_ft",pitch.x_ac_ft, ...
    "Vertical_MAC_ft",vertical.MAC_ft, ...
    "Vertical_xac_ft",vertical.x_ac_ft, ...
    "Wing_e",wing.e_stability, ...
    "Pitch_e",pitch.e_stability, ...
    "Vertical_e",vertical.e_stability, ...
    "Wing_CLa",wing.CLa_rad, ...
    "Pitch_CLa",pitch.CLa_rad, ...
    "Vertical_CLa",vertical.CLa_rad, ...
    "Xac_WingStrake_ft",stability.xac_wing_strake_ft, ...
    "Xac_WingStrakeFuselage_ft", ...
        stability.xac_wing_strake_fuselage_ft, ...
    "Downwash",stability.downwash_gradient, ...
    "TailArm_ft",stability.tail_arm_ft, ...
    "HorizontalTailVolume",stability.horizontal_tail_volume, ...
    "NeutralPoint_ft",stability.x_neutral_point_ft, ...
    "SM_Takeoff",stability.static_margin_takeoff, ...
    "SM_Landing",stability.static_margin_landing, ...
    "CLa_Total",stability.CLa, ...
    "Cma_Takeoff",stability.Cma_takeoff, ...
    "CYb",stability.CYb, ...
    "Clb",stability.Clb, ...
    "Cnb",stability.Cnb, ...
    "CYr",stability.CYr, ...
    "Cnr",stability.Cnr, ...
    "Clr",stability.Clr, ...
    "CYp",stability.CYp, ...
    "Cnp",stability.Cnp, ...
    "Clp",stability.Clp, ...
    "Cmq",stability.Cmq, ...
    "Cmadot",stability.Cmadot, ...
    "CLde",stability.CLde, ...
    "Cmde",stability.Cmde, ...
    "CLref",trim.CL_ref, ...
    "AlphaRef_rad",trim.alpha_ref_rad, ...
    "AlphaZeroLift_rad",trim.alpha_zero_lift_rad, ...
    "PitchControlRef_rad",trim.pitch_control_ref_rad);

validation = struct();
validation.stability_expected = expected_stability;
validation.stability_calculated = calculated_stability;

validation.performance_expected = struct( ...
    "Mach",0.524147536185738, ...
    "Temperature_R",394.09, ...
    "Density_slug_ft3",0.000735278848237216, ...
    "CD0",0.0169957344431338, ...
    "k1",0.116031445730365, ...
    "k2",-0.00630192501971521, ...
    "Sonic_fps",973.016863163224, ...
    "Speed_fps",510.004391494179, ...
    "q_psf",95.6246609964589, ...
    "Theta",0.759779444369469, ...
    "Theta0",0.801526382924393, ...
    "Delta",0.234967643357857, ...
    "Delta0",0.283344018945221, ...
    "DryLapse",0.23878979811829, ...
    "ABLapse",0.262830468959845, ...
    "CL",0.984015698981146, ...
    "CD",0.123146269696171, ...
    "Drag_lbf",3532.74608780246, ...
    "DryThrust_lbf",3581.84697177435, ...
    "ABThrust_lbf",6247.48024717551);

%% Final package

data = struct();
data.source = struct('name',"Brandt Jet Designer F-16A reconstruction", 'runtime_model',"formula_based", 'spreadsheet_runtime_dependency',false);

data.aircraft = aircraft;
data.geometry = geometry;
data.aero = aero;
data.engine = engine;
data.atmos = atmos;
data.wing = wing;
data.pitch = pitch;
data.strake = strake;
data.vertical = vertical;
data.wing_strake = wing_strake;
data.stability = stability;
data.trim = trim;
data.controls = controls;
data.validity = validity;
data.mass = mass;
data.reference = reference;
data.validation = validation;

data.approximations = struct( ...
    'CLq',"tail-volume approximation", ...
    'aileron_derivatives',"geometry-based strip integration", ...
    'rudder_derivatives',"vertical-tail effectiveness approximation", ...
    'control_drag',"quadratic assumed increments", ...
    'inertia',mass.inertia_source);

data.functions = struct('atmosphere',@brandtAtmosphere, 'aero',@brandtAero, 'engine',@brandtEngineEquations, 'longitudinal_trim',@brandtLongitudinalTrim, 'validation_table',@validationTable);
end

%% Local reusable equations

function surface = surfaceGeometry(inputSurface,isVertical)
surface = inputSurface;
surface.taper = inputSurface.ct_ft/inputSurface.cr_ft;
surface.MAC_ft = 2/3*inputSurface.cr_ft *(1+surface.taper+surface.taper^2)/(1+surface.taper);

if isVertical
    spanFactor = inputSurface.b_ft/3;
else
    spanFactor = inputSurface.b_ft/6;
end

surface.y_MAC_ft = spanFactor*(1+2*surface.taper)/(1+surface.taper);
surface.x_MAC_ft = inputSurface.x_LE_root_ft +surface.y_MAC_ft*tand(inputSurface.sweep_LE_deg);
surface.x_ac_ft = surface.x_MAC_ft+0.25*surface.MAC_ft;

if isVertical
    surface.z_ac_ft = surface.y_MAC_ft;
else
    surface.z_ac_ft = surface.y_MAC_ft*cosd(inputSurface.dihedral_deg) *tand(inputSurface.dihedral_deg);
end
end

function e = brandtSpanEfficiency(AR,sweep_LE_deg)
e = 2/(2-AR+sqrt(4+AR^2*(1+tand(sweep_LE_deg/2)^2)));
end

function e = performanceSpanEfficiency(AR,sweep_LE_deg,sweep_tmax_deg)
average_sweep_deg = (sweep_LE_deg+sweep_tmax_deg)/2;
e = max(2/(2-AR+sqrt(4+AR^2*(1+tand(average_sweep_deg)^2))), 0.6);
end

function CLa = finiteWingCLa(AR,e)
CLa = 2*pi/(1+2*pi/(pi*e*AR));
end

function atm = brandtAtmosphere(altitude_ft,atmos)
validateattributes(altitude_ft,{'numeric'}, {'real','finite','nonnegative'});

temperature_R = zeros(size(altitude_ft));
density_slug_ft3 = zeros(size(altitude_ft));

troposphere = altitude_ft < atmos.h_tropopause_ft;
stratosphere = ~troposphere;

temperature_R(troposphere) = atmos.Tsl_R+atmos.lapse_R_per_ft.*altitude_ft(troposphere);

density_slug_ft3(troposphere) = atmos.rho_sl_slug_ft3 .*(temperature_R(troposphere)/atmos.Tsl_R) .^-(1+atmos.g_ft_s2 /(atmos.lapse_R_per_ft*atmos.R_ft_lbf_slug_R));

temperature_R(stratosphere) = atmos.Tstrat_R;

density_slug_ft3(stratosphere) = atmos.rho_trop_slug_ft3 .*exp(-atmos.g_ft_s2 /(atmos.R_ft_lbf_slug_R*atmos.Tstrat_R) .*(altitude_ft(stratosphere)-atmos.h_tropopause_ft));

pressure_psf = density_slug_ft3.*atmos.R_ft_lbf_slug_R.*temperature_R;
speed_of_sound_fps = sqrt(atmos.gamma*atmos.R_ft_lbf_slug_R.*temperature_R);

atm = struct( ...
    'altitude_ft',altitude_ft, ...
    'temperature_R',temperature_R, ...
    'density_slug_ft3',density_slug_ft3, ...
    'pressure_psf',pressure_psf, ...
    'speed_of_sound_fps',speed_of_sound_fps, ...
    'theta',temperature_R/atmos.Tsl_R, ...
    'delta',pressure_psf/atmos.Psl_psf);
end

function [CD0,k1,k2,CD] = brandtAero(Mach,CL,CDx,aero)
validateattributes(Mach,{'numeric'}, {'real','finite','nonnegative'});
validateattributes(CL,{'numeric'},{'real','finite'});
validateattributes(CDx,{'numeric'},{'real','finite'});

CD0 = zeros(size(Mach));
k1 = zeros(size(Mach));
k2 = zeros(size(Mach));

subsonic = Mach <= aero.Mcrit;
transonic = Mach > aero.Mcrit & Mach < aero.Mwave;
post_wave = Mach >= aero.Mwave;

CD0(subsonic) = aero.CD0_sub;
k1(subsonic) = aero.k1_sub;

CD0(transonic) = aero.CD0_sub +(aero.CD0_wave-aero.CD0_sub) .*(Mach(transonic)-aero.Mcrit)/(aero.Mwave-aero.Mcrit);

k1(transonic) = aero.k1_sub +(aero.k1_wave-aero.k1_sub) .*(Mach(transonic)-aero.Mcrit)/(aero.Mwave-aero.Mcrit);

CD0(post_wave) = aero.CD0_sub +(aero.CD0_wave-aero.CD0_sub) .*(1-0.3*sqrt(Mach(post_wave)-aero.Mwave));

k1(post_wave) = aero.k1_sub +(aero.k1_wave-aero.k1_sub) .*(1-0.3*sqrt(Mach(post_wave)-aero.Mwave));

k2_subsonic = Mach <= aero.Mcrit;
k2_transition = Mach > aero.Mcrit & Mach < aero.M_k2_zero;
k2_zero = Mach >= aero.M_k2_zero;

k2(k2_subsonic) = aero.k2_sub;
k2(k2_transition) = aero.k2_sub +(0-aero.k2_sub) .*(Mach(k2_transition)-aero.Mcrit) /(aero.M_k2_zero-aero.Mcrit);
k2(k2_zero) = 0;

CD0 = CD0+CDx;
CD = CD0+k1.*CL.^2+k2.*CL;
end

function out = brandtEngineEquations(Mach,altitude_ft,atmos,engine)
atm = brandtAtmosphere(altitude_ft,atmos);

theta = atm.theta+zeros(size(Mach));
delta = atm.delta+zeros(size(Mach));
theta0 = theta.*(1+0.2*Mach.^2);
delta0 = delta.*(1+0.2*Mach.^2).^3.5;

dry_lapse = zeros(size(Mach));
ab_lapse = zeros(size(Mach));

dry_cool = theta0 <= engine.TR;
dry_hot = ~dry_cool;
dry_lapse(dry_cool) = delta0(dry_cool) .*(1-engine.dry.F1.*Mach(dry_cool).^engine.dry.E);
dry_lapse(dry_hot) = delta0(dry_hot) .*(1-engine.dry.F1.*Mach(dry_hot).^engine.dry.E -engine.dry.F2.*(theta0(dry_hot)-engine.TR) ./theta0(dry_hot));

ab_cool = theta0 <= engine.TR;
ab_hot = ~ab_cool;
ab_lapse(ab_cool) = delta0(ab_cool) .*(1-engine.ab.F1.*Mach(ab_cool).^engine.ab.E);
ab_lapse(ab_hot) = delta0(ab_hot) .*(1-engine.ab.F1.*Mach(ab_hot).^engine.ab.E -engine.ab.F2.*(theta0(ab_hot)-engine.TR) ./theta0(ab_hot));

dry_thrust_lbf = engine.number*engine.Tsl_dry_lbf.*dry_lapse;
ab_thrust_lbf = engine.number*engine.Tsl_ab_lbf.*ab_lapse;

TSFC_dry = engine.TSFCsl_dry .*(1+engine.tsfc.mach_factor .*abs(Mach-engine.tsfc.dry_Mref) .^engine.tsfc.mach_exponent) .*theta.^engine.tsfc.theta_exponent;

TSFC_ab = engine.TSFCsl_ab .*(1+engine.tsfc.mach_factor .*abs(Mach-engine.tsfc.ab_Mref) .^engine.tsfc.mach_exponent) .*theta.^engine.tsfc.theta_exponent;

out = struct( ...
    'theta',theta, ...
    'theta0',theta0, ...
    'delta',delta, ...
    'delta0',delta0, ...
    'dry_lapse',dry_lapse, ...
    'ab_lapse',ab_lapse, ...
    'dry_thrust_lbf',dry_thrust_lbf, ...
    'ab_thrust_lbf',ab_thrust_lbf, ...
    'TSFC_dry',TSFC_dry, ...
    'TSFC_ab',TSFC_ab, ...
    'fuel_flow_dry_lbm_hr',TSFC_dry.*dry_thrust_lbf, ...
    'fuel_flow_ab_lbm_hr',TSFC_ab.*ab_thrust_lbf);
end

function out = brandtLongitudinalTrim(CL_required,trim,stability)
A = [stability.CLa,stability.CLde; stability.Cma_takeoff,stability.Cmde];

if abs(det(A)) <= 1e-12
    error("BrandtF16AData:SingularTrimMatrix", "The local longitudinal trim matrix is singular.");
end

rhs = [CL_required(:).'-trim.CL_ref; ...
       -trim.Cm_ref*ones(1,numel(CL_required))];
increment = A\rhs;

alpha = trim.alpha_ref_rad+increment(1,:);
de = trim.pitch_control_ref_rad+increment(2,:);

CL_check = trim.CL_ref +stability.CLa.*(alpha-trim.alpha_ref_rad) +stability.CLde.*(de-trim.pitch_control_ref_rad);

Cm_residual = trim.Cm_ref +stability.Cma_takeoff.*(alpha-trim.alpha_ref_rad) +stability.Cmde.*(de-trim.pitch_control_ref_rad);

out = struct('alpha_rad',reshape(alpha,size(CL_required)), 'pitch_control_rad',reshape(de,size(CL_required)), 'CL_check',reshape(CL_check,size(CL_required)), 'Cm_residual',reshape(Cm_residual,size(CL_required)));
end

function tbl = validationTable(expected,calculated)
names = fieldnames(expected);
n = numel(names);

workbook = nan(n,1);
model = nan(n,1);
difference = nan(n,1);
relative_error_percent = nan(n,1);

for i = 1:n
    name = names{i};
    workbook(i) = expected.(name);
    model(i) = calculated.(name);
    difference(i) = model(i)-workbook(i);

    if abs(workbook(i)) > eps
        relative_error_percent(i) = 100*difference(i)/workbook(i);
    end
end

tbl = table(string(names),workbook,model,difference, relative_error_percent, 'VariableNames',{'Quantity','Workbook','Model','Difference', 'RelativeError_percent'});
end
