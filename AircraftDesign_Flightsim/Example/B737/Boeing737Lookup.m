function C = Boeing737Lookup(state_vec, control_vec, geometry)

alpha_range = deg2rad([-4 -2 0 2 4 6 8 10 12 14 16 18]);
elevator_range = deg2rad([-17 -10 0 10 17]);

CL_table = [
    -0.20 -0.08  0.00  0.08  0.12;
     0.00  0.12  0.20  0.28  0.32;
     0.20  0.32  0.40  0.48  0.52;
     0.40  0.52  0.60  0.68  0.72;
     0.60  0.72  0.80  0.88  0.92;
     0.80  0.92  1.00  1.08  1.12;
     1.00  1.12  1.20  1.28  1.32;
     1.20  1.32  1.40  1.48  1.52;
     1.38  1.50  1.58  1.66  1.70;
     1.54  1.66  1.74  1.82  1.86;
     1.68  1.80  1.88  1.96  2.00;
     1.78  1.90  1.98  2.06  2.10
];

CD_table = [
    0.028 0.024 0.022 0.024 0.028;
    0.025 0.021 0.020 0.022 0.026;
    0.024 0.020 0.019 0.021 0.025;
    0.026 0.022 0.021 0.023 0.027;
    0.030 0.026 0.025 0.027 0.031;
    0.038 0.034 0.033 0.035 0.039;
    0.050 0.046 0.045 0.047 0.051;
    0.068 0.064 0.063 0.065 0.069;
    0.093 0.089 0.088 0.090 0.094;
    0.125 0.121 0.120 0.122 0.126;
    0.165 0.161 0.160 0.162 0.166;
    0.213 0.209 0.208 0.210 0.214
];

Cm_table = [
     0.050  0.030  0.000 -0.030 -0.050;
     0.040  0.020 -0.010 -0.040 -0.060;
     0.030  0.010 -0.020 -0.050 -0.070;
     0.020  0.000 -0.030 -0.060 -0.080;
     0.010 -0.010 -0.040 -0.070 -0.090;
     0.000 -0.020 -0.050 -0.080 -0.100;
    -0.010 -0.030 -0.060 -0.090 -0.110;
    -0.020 -0.040 -0.070 -0.100 -0.120;
    -0.030 -0.050 -0.080 -0.110 -0.130;
    -0.040 -0.060 -0.090 -0.120 -0.140;
    -0.050 -0.070 -0.100 -0.130 -0.150;
    -0.060 -0.080 -0.110 -0.140 -0.160
];

CYb = -0.85;
Clb = -0.10;
Cnb = 0.15;

CLq = 5.5;
Cmq = -18.0;
Clp = -0.45;
Cnr = -0.30;
Clr = 0.15;
Cnp = -0.05;
CYp = -0.08;
CYr = 0.35;

mach_wave = [0.0 0.75 0.82 0.88 0.95];
CD_wave = [0.0 0.0 0.008 0.025 0.045];

vel = state_vec(4:6);
omega = state_vec(10:12);
altitude = max(-state_vec(3), 0);
V = max(norm(vel), 1e-6);
alpha = atan2(vel(3), max(abs(vel(1)), 1e-9));
beta = atan2(vel(2), max(sqrt(vel(1)^2 + vel(3)^2), 1e-9));
[~, a, ~] = atmosisa(altitude);
mach = V / max(a, 1e-6);

de = 0;
if numel(control_vec) >= 2
    de = control_vec(2);
end

S = geometry.wing_area;
b = geometry.wing_span;
if isprop(geometry, 'mean_aerodynamic_chord')
    c = geometry.mean_aerodynamic_chord;
else
    c = geometry.wing_chord;
end

CL = interp2(elevator_range, alpha_range, CL_table, de, alpha, 'linear', 0);
CD = interp2(elevator_range, alpha_range, CD_table, de, alpha, 'linear', 0);
Cm = interp2(elevator_range, alpha_range, Cm_table, de, alpha, 'linear', 0);

CD_mach_increment = interp1(mach_wave, CD_wave, mach, 'linear', 0.045);
CD = CD + CD_mach_increment;

CY = CYb * beta;
Cl = Clb * beta;
Cn = Cnb * beta;

if V > 1
    p_hat = omega(1) * b / (2 * V);
    q_hat = omega(2) * c / (2 * V);
    r_hat = omega(3) * b / (2 * V);

    CL = CL + CLq * q_hat;
    Cm = Cm + Cmq * q_hat;

    Cl = Cl + Clp * p_hat + Clr * r_hat;
    Cn = Cn + Cnp * p_hat + Cnr * r_hat;
    CY = CY + CYp * p_hat + CYr * r_hat;
end

takeoff_params = struct('mu_ground', 0.02, 'mu_rolling', 0.02, 'mu_braking', 0.40, ...
    'safety_factor', 1.15, 'CLmax_takeoff', 2.2, 'V2_to_Vs_ratio', 1.20, ...
    'VR_to_Vs_ratio', 1.10, 'V1_to_VR_ratio', 0.92, 'CD0_takeoff', 0.025, ...
    'K_takeoff', 0.045, 'CLalpha_takeoff', 5.0, 'screen_height_takeoff_ft', 35, ...
    'rotation_alpha_deg', 12, 'initial_climb_angle_deg', 8, 'reaction_time_s', 1.0, ...
    'continue_time_after_VR_s', 2.0, 'min_accel_mps2', 0.5, 'min_reduced_accel_mps2', 0.3, ...
    'min_brake_decel_mps2', 1.5);

landing_params = struct('mu_braking', 0.60, 'mu_spoiler', 0.00, 'mu_reverser', 0.00, ...
    'approach_angle_deg', 3, 'safety_factor', 1.67, 'CLmax_landing', 2.4, ...
    'CD0_landing', 0.035, 'CD_spoiler', 0.08, 'CD_reverser', 0.00, ...
    'CL_landing_touchdown', 0.1, 'screen_height_landing_ft', 50, 'flare_height_m', 5, ...
    'Vapp_to_Vs_ratio', 1.30, 'Vtd_to_Vapp_ratio', 0.88, 'idle_throttle', 0.05, ...
    'use_idle_thrust', true, 'min_brake_decel_mps2', 2.0);

C = struct('CL', CL, 'CD', CD, 'CY', CY, 'Cl', Cl, 'Cm', Cm, 'Cn', Cn, ...
    'takeoff_params', takeoff_params, 'landing_params', landing_params);
end