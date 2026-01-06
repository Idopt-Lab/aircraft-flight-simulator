function C = Boeing737Lookup(state_vec, control_vec, geometry)

vel = state_vec(4:6);
V = max(norm(vel), 1e-6);

alpha = atan2(vel(3), max(abs(vel(1)), 1e-9));
beta  = atan2(vel(2), max(sqrt(vel(1)^2 + vel(3)^2), 1e-9));

altitude = max(-state_vec(3), 0);
[~, a, ~, ~] = atmosisa(altitude);
a = max(a, 1e-3);
M = V / a;

p = 0; q = 0; r = 0;
if numel(state_vec) >= 12
    p = state_vec(10);
    q = state_vec(11);
    r = state_vec(12);
end

b = geometry.wing_span;
S = geometry.wing_area;
if isprop(geometry,'mean_aerodynamic_chord')
    cbar = geometry.mean_aerodynamic_chord;
else
    cbar = S / max(b, 1e-6);
end

p_hat = p * b / (2*V);
q_hat = q * cbar / (2*V);
r_hat = r * b / (2*V);

aileron  = 0;
elevator = 0;
rudder   = 0;

if numel(control_vec) >= 1, aileron  = control_vec(1); end
if numel(control_vec) >= 2, elevator = control_vec(2); end
if numel(control_vec) >= 3, rudder   = control_vec(3); end

alpha_tab = [-0.20 0.00 0.23 0.46];
CL_tab    = [-0.68 0.20 1.20 0.20];

alpha_tab_CD0 = [-1.57 -0.26 0.00 0.26 1.57];
CD0_tab       = [ 1.50  0.042 0.021 0.042 1.50];

mach_tab_CDm = [0.00 0.79 1.10 1.80];
CDm_tab      = [0.00 0.00 0.023 0.015];

mach_tab_Cmde = [0.0 2.0];
Cmde_tab      = [-1.20 -0.30];

mach_tab_Clda = [0.0 2.0];
Clda_tab      = [0.100 0.033];

k_induced = 0.043;

CLalpha = interp1(alpha_tab, CL_tab, alpha, 'linear', 'extrap');
CLde = 0.2 * elevator;
CL = CLalpha + CLde;

CD0 = interp1(alpha_tab_CD0, CD0_tab, alpha, 'linear', 'extrap');
CDi = k_induced * (CL^2);
CDmach = interp1(mach_tab_CDm, CDm_tab, M, 'linear', 'extrap');
CDde = 0.059 * abs(elevator);
CD = CD0 + CDi + CDmach + CDde;

CYb = -1.0;
CY = CYb * beta;

Clb = -0.09;
Clp = -0.4;
Clr = 0.09;
Cldr = 0.01;
Clda = interp1(mach_tab_Clda, Clda_tab, M, 'linear', 'extrap');

Cl = (Clb * beta) + (Clp * p_hat) + (Clr * r_hat) + (Clda * aileron) + (Cldr * rudder);

Cmalpha = -0.6;
Cmq = -27.0;
Cmde = interp1(mach_tab_Cmde, Cmde_tab, M, 'linear', 'extrap') * elevator;

Cm = (Cmalpha * alpha) + Cmde + (Cmq * q_hat);

Cnb = 0.26;
Cnr = -0.35;
Cndr = -0.20;

Cn = (Cnb * beta) + (Cnr * r_hat) + (Cndr * rudder);

takeoff_params = struct('mu_ground',0.02,'mu_rolling',0.02,'mu_braking',0.40, ...
'safety_factor',1.15,'CLmax_takeoff',2.2,'V2_to_Vs_ratio',1.20, ...
'VR_to_Vs_ratio',1.10,'V1_to_VR_ratio',0.92,'CD0_takeoff',0.025, ...
'K_takeoff',0.043,'CLalpha_takeoff',5.5,'screen_height_takeoff_ft',35, ...
'rotation_alpha_deg',8,'initial_climb_angle_deg',8,'reaction_time_s',1.0, ...
'continue_time_after_VR_s',2.0,'min_accel_mps2',0.8,'min_reduced_accel_mps2',0.4, ...
'min_brake_decel_mps2',1.2);

landing_params = struct('mu_braking',0.50,'mu_spoiler',0.05,'mu_reverser',0.05, ...
'approach_angle_deg',3,'safety_factor',1.67,'CLmax_landing',2.5, ...
'CD0_landing',0.035,'CD_spoiler',0.06,'CD_reverser',0.04, ...
'CL_landing_touchdown',0.30,'screen_height_landing_ft',50,'flare_height_m',5, ...
'Vapp_to_Vs_ratio',1.30,'Vtd_to_Vapp_ratio',0.88,'idle_throttle',0.08, ...
'use_idle_thrust',true,'min_brake_decel_mps2',2.0);

C = struct('CL',CL,'CD',CD,'CY',CY,'Cl',Cl,'Cm',Cm,'Cn',Cn, ...
    'takeoff_params',takeoff_params,'landing_params',landing_params);

end
