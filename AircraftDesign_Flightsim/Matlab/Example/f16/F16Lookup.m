function C = F16Lookup(state_vec, control_vec, geometry)

alpha_range = [-0.1750, -0.0870, 0.0000, 0.0870, 0.1750, 0.2620, 0.3490, 0.4360, 0.5240, 0.6110, 0.6980, 0.7850];
elevator_range = [-0.4360, -0.2180, 0.0000, 0.2180, 0.4360];

CLDh_table = [-0.6590, -0.7090, -0.7540, -0.7920, -0.8250; 
              -0.1510, -0.1960, -0.2380, -0.2780, -0.3160; 
               0.1830,  0.1410,  0.1000,  0.0590,  0.0170; 
               0.4910,  0.4540,  0.4140,  0.3710,  0.3260; 
               0.7970,  0.7630,  0.7250,  0.6800,  0.6300; 
               1.1080,  1.0790,  1.0410,  0.9930,  0.9400; 
               1.3950,  1.3660,  1.3270,  1.2740,  1.2140; 
               1.6150,  1.5870,  1.5470,  1.4900,  1.4270; 
               1.8040,  1.7770,  1.7370,  1.6740,  1.6100; 
               1.9000,  1.8720,  1.8290,  1.7660,  1.6990; 
               1.8980,  1.8690,  1.8220,  1.7570,  1.6890; 
               1.7530,  1.7240,  1.6740,  1.6120,  1.5460];

CDDh_table = [0.2170, 0.1740, 0.1560, 0.1810, 0.2300; 
              0.0940, 0.0550, 0.0410, 0.0620, 0.1010; 
              0.0810, 0.0400, 0.0210, 0.0390, 0.0760; 
              0.1060, 0.0610, 0.0400, 0.0570, 0.1010; 
              0.1660, 0.1190, 0.0960, 0.1140, 0.1580; 
              0.2520, 0.2030, 0.1820, 0.2020, 0.2400; 
              0.4040, 0.3620, 0.3470, 0.3710, 0.4160; 
              0.6280, 0.5880, 0.5770, 0.6010, 0.6370; 
              0.8750, 0.8400, 0.8260, 0.8520, 0.8800; 
              1.1270, 1.0950, 1.0840, 1.1020, 1.1250; 
              1.3650, 1.3340, 1.3260, 1.3380, 1.3560; 
              1.5170, 1.4870, 1.4780, 1.4820, 1.4890];

CmDh_table = [ 0.2050,  0.0810, -0.0460, -0.1740, -0.2590; 
               0.1680,  0.0770, -0.0200, -0.1450, -0.2020; 
               0.1860,  0.1070, -0.0090, -0.1210, -0.1840; 
               0.1960,  0.1100, -0.0050, -0.1270, -0.1930; 
               0.2130,  0.1110, -0.0060, -0.1290, -0.1990; 
               0.2510,  0.1410,  0.0100, -0.1020, -0.1500; 
               0.2450,  0.1270,  0.0060, -0.0970, -0.1600; 
               0.2380,  0.1190, -0.0010, -0.1130, -0.1670; 
               0.2520,  0.1330,  0.0140, -0.0870, -0.1040; 
               0.2310,  0.1080,  0.0000, -0.0840, -0.0760; 
               0.1980,  0.0810, -0.0130, -0.0690, -0.0410; 
               0.1920,  0.0930,  0.0320, -0.0060, -0.0050];

CYb = -1.1460;
Clb = -0.186;
Cnb = 0.241;

Cmalpha_vs_alpha = [-0.0280, -0.0350, -0.0430, -0.0500, -0.0590, -0.0650, -0.0720, -0.0800, -0.0850, -0.0900, -0.0950, -0.1000];

CLq = 35.0;
Cmq = -10;
Clp = -0.80;
Cnr = -0.75;

Clr_vs_alpha = [-0.126, -0.026, 0.063, 0.113, 0.208, 0.230, 0.319, 0.437, 0.680, 0.100, 0.447, -0.330];
Cnp_vs_alpha = [-0.061, -0.052, -0.052, 0.012, 0.013, 0.024, -0.050, -0.150, -0.130, -0.158, -0.240, -0.150];
CYp_vs_alpha = [-0.108, -0.108, -0.188, 0.110, 0.258, 0.226, 0.344, 0.362, 0.611, 0.529, 0.298, -0.227];
CYr_vs_alpha = [0.882, 0.852, 0.876, 0.958, 0.962, 0.974, 0.819, 0.483, 0.590, 1.210, -0.493, -1.040];

mach_wave = [0.0, 0.81, 1.1, 1.8];
CD_wave = [0.0, 0.0, 0.023, 0.015];

vel = state_vec(4:6);
omega = state_vec(10:12);
altitude = max(-state_vec(3),0);
V = max(norm(vel),1e-6);
alpha = atan2(vel(3), max(abs(vel(1)), 1e-9));
beta = atan2(vel(2), max(sqrt(vel(1)^2 + vel(3)^2), 1e-9));
[~,a,~] = atmosisa(altitude);
mach = V/max(a,1e-6);

de = 0;
if numel(control_vec) >= 2
    de = control_vec(2);
end

S = geometry.wing_area;
b = geometry.wing_span;
if isprop(geometry,'mean_aerodynamic_chord')
    c = geometry.mean_aerodynamic_chord;
else
    c = geometry.wing_chord;
end

CL = interp2(elevator_range, alpha_range, CLDh_table, de, alpha, 'linear', 0);
CD = interp2(elevator_range, alpha_range, CDDh_table, de, alpha, 'linear', 0);
Cm = interp2(elevator_range, alpha_range, CmDh_table, de, alpha, 'linear', 0);

CD_mach_increment = interp1(mach_wave, CD_wave, mach, 'linear', 0.015);
CD = CD + CD_mach_increment;

Cmalpha = interp1(alpha_range, Cmalpha_vs_alpha, alpha, 'linear', 'extrap');
Cm = Cm + Cmalpha*alpha;

CY = CYb*beta;
Cl = Clb*beta;
Cn = Cnb*beta;

if V > 1
    p_hat = omega(1)*b/(2*V);
    q_hat = omega(2)*c/(2*V);
    r_hat = omega(3)*b/(2*V);

    CL = CL + CLq*q_hat;
    Cm = Cm + Cmq*q_hat;

    Clr = interp1(alpha_range, Clr_vs_alpha, alpha, 'linear', 'extrap');
    Cl = Cl + Clp*p_hat + Clr*r_hat;

    Cnp = interp1(alpha_range, Cnp_vs_alpha, alpha, 'linear', 'extrap');
    Cn = Cn + Cnp*p_hat + Cnr*r_hat;

    CYp = interp1(alpha_range, CYp_vs_alpha, alpha, 'linear', 'extrap');
    CYr = interp1(alpha_range, CYr_vs_alpha, alpha, 'linear', 'extrap');
    CY = CY + CYp*p_hat + CYr*r_hat;
end

takeoff_params = struct('mu_ground',0.02,'mu_rolling',0.02,'mu_braking',0.35, ...
    'safety_factor',1.15,'CLmax_takeoff',1.8,'V2_to_Vs_ratio',1.20, ...
    'VR_to_Vs_ratio',1.10,'V1_to_VR_ratio',0.92,'CD0_takeoff',0.02, ...
    'K_takeoff',0.05,'CLalpha_takeoff',5.0,'screen_height_takeoff_ft',35, ...
    'rotation_alpha_deg',10,'initial_climb_angle_deg',10,'reaction_time_s',1.0, ...
    'continue_time_after_VR_s',2.0,'min_accel_mps2',0.5,'min_reduced_accel_mps2',0.3, ...
    'min_brake_decel_mps2',1.0);

landing_params = struct('mu_braking',0.70,'mu_spoiler',0.00,'mu_reverser',0.00, ...
    'approach_angle_deg',3,'safety_factor',1.67,'CLmax_landing',2.0, ...
    'CD0_landing',0.035,'CD_spoiler',0.04,'CD_reverser',0.00, ...
    'CL_landing_touchdown',0.05,'screen_height_landing_ft',50,'flare_height_m',5, ...
    'Vapp_to_Vs_ratio',1.30,'Vtd_to_Vapp_ratio',0.88,'idle_throttle',0.05, ...
    'use_idle_thrust',true,'min_brake_decel_mps2',1.5);

C = struct('CL',CL,'CD',CD,'CY',CY,'Cl',Cl,'Cm',Cm,'Cn',Cn, ...
    'takeoff_params',takeoff_params,'landing_params',landing_params);
end