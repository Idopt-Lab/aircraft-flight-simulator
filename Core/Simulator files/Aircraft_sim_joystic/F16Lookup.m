function lookup = F16Lookup()
alpha_range    = [-0.1750, -0.0870, 0.0000, 0.0870, 0.1750, 0.2620, 0.3490, 0.4360, 0.5240, 0.6110, 0.6980, 0.7850];
elevator_range = [-0.4360, -0.2180, 0.0000, 0.2180, 0.4360];

CLDh_table = [
    -0.6590, -0.7090, -0.7540, -0.7920, -0.8250;
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

CDDh_table = [
     0.2170,  0.1740,  0.1560,  0.1810,  0.2300;
     0.0940,  0.0550,  0.0410,  0.0620,  0.1010;
     0.0810,  0.0400,  0.0210,  0.0390,  0.0760;
     0.1060,  0.0610,  0.0400,  0.0570,  0.1010;
     0.1660,  0.1190,  0.0960,  0.1140,  0.1580;
     0.2520,  0.2030,  0.1820,  0.2020,  0.2400;
     0.4040,  0.3620,  0.3470,  0.3710,  0.4160;
     0.6280,  0.5880,  0.5770,  0.6010,  0.6370;
     0.8750,  0.8400,  0.8260,  0.8520,  0.8800;
     1.1270,  1.0950,  1.0840,  1.1020,  1.1250;
     1.3650,  1.3340,  1.3260,  1.3380,  1.3560;
     1.5170,  1.4870,  1.4780,  1.4820,  1.4890];

CmDh_table = [
     0.2050,  0.0810, -0.0460, -0.1740, -0.2590;
     0.1680,  0.0770, -0.0200, -0.1450, -0.2020;
     0.1860,  0.1070, -0.0090, -0.1210, -0.1840;
     0.1960,  0.1100, -0.0050, -0.1270, -0.1930;
     0.2130,  0.1100, -0.0060, -0.1290, -0.1990;
     0.2510,  0.1410,  0.0100, -0.1020, -0.1500;
     0.2450,  0.1270,  0.0060, -0.0970, -0.1600;
     0.2380,  0.1190, -0.0010, -0.1130, -0.1670;
     0.2520,  0.1330,  0.0140, -0.0870, -0.1040;
     0.2310,  0.1080,  0.0000, -0.0840, -0.0760;
     0.1980,  0.0810, -0.0130, -0.0690, -0.0410;
     0.1920,  0.0930,  0.0320, -0.0060, -0.0050];

CYb  = -1.1460;  CYdr = 0.0860;
Clb  = -0.186;   Cnb  = 0.241;
Clda_vs_alpha = [0.0420, 0.0530, 0.0520, 0.0520, 0.0480, 0.0480, 0.0420, 0.0370, 0.0310, 0.0260, 0.0170, 0.0120];
Cndr_vs_alpha = [-0.0480, -0.0450, -0.0450, -0.0450, -0.0440, -0.0450, -0.0470, -0.0480, -0.0490, -0.0450, -0.0330, -0.0160];
CLq = 28.9; Cmq = -5.23; Clp = -0.443; Cnr = -0.378;
mach_range        = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 2.0];
altitude_range_ft = [0, 10000, 20000, 30000, 40000, 50000, 60000];
altitude_range    = altitude_range_ft * 0.3048;

mil_thrust_lbf = [
    14670, 10856,  7834,  5457, 3535, 2186, 0;
    12680,  9150,  6313,  4040, 2470, 1400, 0;
    12610,  9312,  6610,  4290, 2600, 1560, 0;
    12640,  9839,  7090,  4660, 2840, 1660, 0;
    12390, 10176,  7750,  5320, 3250, 1930, 0;
    11680,  9848,  8050,  6100, 3800, 2310, 0;
    16353, 12962,  9891,  7245, 4840, 2991, 0;
    17509, 14190, 11041,  8205, 5539, 3435, 0;
    16200, 13500, 10500,  7800, 5200, 3200, 0;
    14000, 12000,  9500,  7000, 4500, 2800, 0];

max_thrust_lbf = [
    23830, 19502, 15792, 12582,  8950,  5545, 0;
    21420, 15700, 11225,  7323,  4435,  2600, 0;
    22700, 16860, 12250,  8154,  5000,  2835, 0;
    24240, 18910, 13760,  9285,  5700,  3215, 0;
    26070, 21075, 15975, 11115,  6860,  3950, 0;
    28886, 23319, 18300, 13484,  8642,  5057, 0;
    28951, 23670, 19114, 15155, 10744,  6658, 0;
    37413, 31184, 24959, 19578, 13788,  8543, 0;
    35000, 29000, 23000, 18000, 12500,  7800, 0;
    32000, 26500, 21000, 16000, 11000,  7000, 0];

idle_thrust_lbf = [
     716,   775,  1018,  1319,  1735,  2152, 0;
     635,   425,   690,  1010,  1330,  1700, 0;
      60,    25,   345,   755,  1130,  1525, 0;
   -1020,  -710,  -300,   350,   910,  1360, 0;
   -2700, -1900, -1300,  -247,   600,  1100, 0;
   -3600, -1400,  -595,  -342,  -200,   700, 0;
   -4000, -1800,  -800,  -500,  -300,   500, 0;
   -4200, -2000,  -900,  -600,  -400,   400, 0;
   -4500, -2200, -1000,  -700,  -500,   300, 0;
   -5000, -2500, -1200,  -800,  -600,   200, 0];

mil_thrust_table = mil_thrust_lbf * 4.448;
max_thrust_table = max_thrust_lbf * 4.448;
idle_thrust_table= idle_thrust_lbf * 4.448;

lookup = @(state_vec, control_vec, geometry) interpolate_all_data( ...
    state_vec, control_vec, geometry, ...
    alpha_range, elevator_range, CLDh_table, CDDh_table, CmDh_table, ...
    CYb, CYdr, Clb, Cnb, Clda_vs_alpha, Cndr_vs_alpha, ...
    CLq, Cmq, Clp, Cnr, ...
    mach_range, altitude_range, idle_thrust_table, mil_thrust_table, max_thrust_table);
end


function C = interpolate_all_data(state_vec, control_vec, geometry, ...
    alpha_range, elevator_range, CLDh_table, CDDh_table, CmDh_table, ...
    CYb, CYdr, Clb, Cnb, Clda_vs_alpha, Cndr_vs_alpha, ...
    CLq, Cmq, Clp, Cnr, ...
    mach_range, altitude_range, idle_thrust_table, mil_thrust_table, max_thrust_table)

vel   = state_vec(4:6);
omega = state_vec(10:12);
altitude = max(-state_vec(3), 0);

V     = max(norm(vel), 1e-6);
alpha = atan2(vel(3), vel(1));
beta  = asin(max(-0.999, min(0.999, vel(2)/V)));

[~, a, ~] = atmosisa(altitude);   
mach = V / max(a,1e-6);

mach_limited     = max(0, min(2.0, mach));
altitude_limited = max(0, min(18288, altitude)); 

ua       = getControl(control_vec, 1, 0);
ue       = getControl(control_vec, 2, 0);
ur       = getControl(control_vec, 3, 0);
throttle = getControl(control_vec, 4, 0);

S = geometry.wing_area;
b = geometry.wing_span;
c = geometry.wing_chord;

CL = interp2(elevator_range, alpha_range, CLDh_table, ue, alpha, 'linear', 0);
CD = interp2(elevator_range, alpha_range, CDDh_table, ue, alpha, 'linear', 0);
Cm = interp2(elevator_range, alpha_range, CmDh_table,  ue, alpha, 'linear', 0);

if mach > 1.0
    CD = CD + 0.1 * (mach - 1.0)^2;
end

CY   = CYb  * beta + CYdr * ur;

Clda = interp1(alpha_range, Clda_vs_alpha, alpha, 'linear', 'extrap');
Cl   = Clda * ua + Clb * beta;

Cndr = interp1(alpha_range, Cndr_vs_alpha, alpha, 'linear', 'extrap');
Cn   = Cndr * ur + Cnb * beta;

if V > 1
    CL = CL + CLq * (omega(2) * c / (2*V));
    Cm = Cm + Cmq * (omega(2) * c / (2*V));
    Cl = Cl + Clp * (omega(1) * b / (2*V));
    Cn = Cn + Cnr * (omega(3) * b / (2*V));
end

if throttle <= 0.8
    idle_thrust = interp2(altitude_range, mach_range, idle_thrust_table, altitude_limited, mach_limited, 'linear', 0);
    mil_thrust  = interp2(altitude_range, mach_range,  mil_thrust_table, altitude_limited, mach_limited, 'linear', 0);
    thrust_available = idle_thrust + (throttle/0.8) * (mil_thrust - idle_thrust);
else
    mil_thrust  = interp2(altitude_range, mach_range,  mil_thrust_table, altitude_limited, mach_limited, 'linear', 0);
    max_thrust  = interp2(altitude_range, mach_range,  max_thrust_table, altitude_limited, mach_limited, 'linear', 0);
    ab_ratio = (throttle - 0.8) / 0.2;
    thrust_available = mil_thrust + ab_ratio * (max_thrust - mil_thrust);
end
thrust_available = max(0, thrust_available);

C = struct('CL', CL, 'CD', CD, 'CY', CY, 'Cl', Cl, 'Cm', Cm, 'Cn', Cn, ...
           'thrust_available', thrust_available);
end


function val = getControl(u, i, defaultVal)
if i <= numel(u)
    val = u(i);
else
    val = defaultVal;
end
end