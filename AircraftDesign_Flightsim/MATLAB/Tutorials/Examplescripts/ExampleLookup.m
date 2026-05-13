function C = ExampleLookup(state_vec, control_vec, geometry)
persistent tab FCL FCD FCm init

if isempty(init)
    tab = load('ExampleAeroTable.mat');
    FCL = griddedInterpolant({tab.alpha_deg, tab.mach}, tab.CL, 'linear', 'nearest');
    FCD = griddedInterpolant({tab.alpha_deg, tab.mach}, tab.CD, 'linear', 'nearest');
    FCm = griddedInterpolant({tab.alpha_deg, tab.mach}, tab.Cm, 'linear', 'nearest');
    init = true;
end

vel = state_vec(4:6);
V = max(norm(vel), 1e-6);
alpha = atan2(vel(3), max(abs(vel(1)), 1e-9));
beta  = atan2(vel(2), max(sqrt(vel(1)^2 + vel(3)^2), 1e-9));
altitude = max(-state_vec(3), 0);

[~, a, ~, ~] = atmosisa(altitude);
M = max(V / max(a,1e-3), 0);

alpha_deg = rad2deg(alpha);
M = max(0, min(0.8, M));
alpha_deg = max(min(alpha_deg, 15), -10);

CL = FCL(alpha_deg, M);
CD = FCD(alpha_deg, M);
Cm = FCm(alpha_deg, M);

aileron = 0; elevator = 0; rudder = 0;
if numel(control_vec) >= 1, aileron = control_vec(1); end
if numel(control_vec) >= 2, elevator = control_vec(2); end
if numel(control_vec) >= 3, rudder  = control_vec(3); end

CL = CL + tab.derivs.dCL_de * elevator;
Cm = Cm + tab.derivs.dCm_de * elevator;

Cl = tab.derivs.dCl_da * aileron + (-0.08) * beta;
Cn = tab.derivs.dCn_dr * rudder  + ( 0.15) * beta;
CY = (-1.0) * beta;

takeoff_params = struct('mu_ground',0.02,'mu_rolling',0.02,'mu_braking',0.40, ...
    'safety_factor',1.10,'CLmax_takeoff',1.6,'V2_to_Vs_ratio',1.20,'VR_to_Vs_ratio',1.10, ...
    'V1_to_VR_ratio',0.92,'CD0_takeoff',0.030,'K_takeoff',0.055,'CLalpha_takeoff',5.0, ...
    'screen_height_takeoff_ft',35,'rotation_alpha_deg',8,'initial_climb_angle_deg',8, ...
    'reaction_time_s',1.0,'continue_time_after_VR_s',2.0,'min_accel_mps2',0.8,'min_reduced_accel_mps2',0.4, ...
    'min_brake_decel_mps2',1.2);

landing_params = struct('mu_braking',0.50,'mu_spoiler',0.00,'mu_reverser',0.00, ...
    'approach_angle_deg',3,'safety_factor',1.50,'CLmax_landing',1.8,'CD0_landing',0.040, ...
    'CD_spoiler',0.00,'CD_reverser',0.00,'CL_landing_touchdown',0.30,'screen_height_landing_ft',50, ...
    'flare_height_m',5,'Vapp_to_Vs_ratio',1.30,'Vtd_to_Vapp_ratio',0.88,'idle_throttle',0.05, ...
    'use_idle_thrust',true,'min_brake_decel_mps2',2.0);

C = struct('CL',CL,'CD',CD,'CY',CY,'Cl',Cl,'Cm',Cm,'Cn',Cn, ...
    'takeoff_params',takeoff_params,'landing_params',landing_params);
end
