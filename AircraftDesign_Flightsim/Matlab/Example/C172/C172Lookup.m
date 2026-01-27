function C = C172Lookup(state_vec, control_vec, geometry)

vel = state_vec(4:6);
omega = state_vec(10:12);

altitude = max(-state_vec(3),0);
V = max(norm(vel),1e-6);

alpha = atan2(vel(3), max(vel(1),1e-9));
beta  = atan2(vel(2), max(sqrt(vel(1)^2 + vel(3)^2),1e-9));

da = 0; de = 0; dr = 0;
if numel(control_vec) >= 1, da = control_vec(1); end
if numel(control_vec) >= 2, de = control_vec(2); end
if numel(control_vec) >= 3, dr = control_vec(3); end

S = geometry.wing_area;
b = geometry.wing_span;
if isprop(geometry,'mean_aerodynamic_chord')
    c = geometry.mean_aerodynamic_chord;
else
    c = geometry.ref_chord;
end

[~,a,~] = atmosisa(altitude);
mach = V/max(a,1e-6);

clamp = @(z,lo,hi) min(max(z,lo),hi);

CL0      = 0.28;
CL_alpha = 5.10;
CL_de    = 0.35;

CD0      = 0.030;
K        = 0.055;

CY_beta  = -0.65;
CY_dr    = 0.18;

Cl_beta  = -0.12;
Cl_da    = 0.08;
Cl_dr    = 0.02;

Cm0      = 0.02;
Cm_alpha = -1.00;
Cm_de    = -1.05;

Cn_beta  = 0.10;
Cn_da    = 0.02;
Cn_dr    = -0.08;

alpha_eff = clamp(alpha, -deg2rad(10), deg2rad(15));

CL = CL0 + CL_alpha*alpha_eff + CL_de*de;
CD = CD0 + K*CL^2;
CY = CY_beta*beta + CY_dr*dr;

Cl = Cl_beta*beta + Cl_da*da + Cl_dr*dr;
Cm = Cm0 + Cm_alpha*alpha_eff + Cm_de*de;
Cn = Cn_beta*beta + Cn_da*da + Cn_dr*dr;

mach_wave = [0.0 0.6 0.8 1.2];
CD_wave   = [0.0 0.0 0.006 0.020];
CD = CD + interp1(mach_wave, CD_wave, mach, 'linear', 'extrap');

if V > 1
    p_hat = omega(1)*b/(2*V);
    q_hat = omega(2)*c/(2*V);
    r_hat = omega(3)*b/(2*V);

    CLq = 3.5;
    Cmq = -8.0;
    Clp = -0.50;
    Cnr = -0.20;

    CL = CL + CLq*q_hat;
    Cm = Cm + Cmq*q_hat;
    Cl = Cl + Clp*p_hat;
    Cn = Cn + Cnr*r_hat;
end

CD = max(CD, 0);

takeoff_params = struct('mu_ground',0.03,'mu_rolling',0.03,'mu_braking',0.35, ...
    'safety_factor',1.10,'CLmax_takeoff',1.6,'V2_to_Vs_ratio',1.20, ...
    'VR_to_Vs_ratio',1.10,'V1_to_VR_ratio',0.92,'CD0_takeoff',0.040, ...
    'K_takeoff',0.055,'CLalpha_takeoff',5.1,'screen_height_takeoff_ft',35, ...
    'rotation_alpha_deg',8,'initial_climb_angle_deg',7,'reaction_time_s',1.0, ...
    'continue_time_after_VR_s',2.0,'min_accel_mps2',0.4,'min_reduced_accel_mps2',0.25, ...
    'min_brake_decel_mps2',1.0);

landing_params = struct('mu_braking',0.35,'mu_spoiler',0.00,'mu_reverser',0.00, ...
    'approach_angle_deg',3,'safety_factor',1.40,'CLmax_landing',2.0, ...
    'CD0_landing',0.055,'CD_spoiler',0.00,'CD_reverser',0.00, ...
    'CL_landing_touchdown',0.10,'screen_height_landing_ft',50,'flare_height_m',4, ...
    'Vapp_to_Vs_ratio',1.30,'Vtd_to_Vapp_ratio',0.90,'idle_throttle',0.05, ...
    'use_idle_thrust',true,'min_brake_decel_mps2',1.5);

C = struct('CL',CL,'CD',CD,'CY',CY,'Cl',Cl,'Cm',Cm,'Cn',Cn, ...
    'takeoff_params',takeoff_params,'landing_params',landing_params);

end
