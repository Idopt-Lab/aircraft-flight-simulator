function [T, fuel_flow] = F100ThrustModel(throttle, M, alt_m, V, rho, T_max_ref)

persistent mach_range altitude_range idle_T mil_T max_T T_ref
if isempty(mach_range)
    mach_range = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4];
    altitude_range = [0, 10000, 20000, 30000, 40000, 50000, 60000] * 0.3048;
    
    idle_thrust_lbf = [716 775 1018 1319 1735 2152 0; 
                       635 425 690 1010 1330 1700 0; 
                       60 25 345 755 1130 1525 0; 
                       -1020 -710 -300 350 910 1360 0; 
                       -2700 -1900 -1300 -247 600 1100 0; 
                       -3600 -1400 -595 -342 -200 700 0; 
                       -4000 -1800 -800 -500 -300 500 0; 
                       -4200 -2000 -900 -600 -400 400 0];
    
    mil_thrust_lbf = [14670 10856 7834 5457 3535 2186 0; 
                      12680 9150 6313 4040 2470 1400 0; 
                      12610 9312 6610 4290 2600 1560 0; 
                      12640 9839 7090 4660 2840 1660 0; 
                      12390 10176 7750 5320 3250 1930 0; 
                      11680 9848 8050 6100 3800 2310 0; 
                      16353 12962 9891 7245 4840 2991 0; 
                      17509 14190 11041 8205 5539 3435 0];
    
    max_thrust_lbf = [23830 19502 15792 12582 8950 5545 0; 
                      21420 15700 11225 7323 4435 2600 0; 
                      22700 16860 12250 8154 5000 2835 0; 
                      24240 18910 13760 9285 5700 3215 0; 
                      26070 21075 15975 11115 6860 3950 0; 
                      28886 23319 18300 13484 8642 5057 0; 
                      28951 23670 19114 15155 10744 6658 0; 
                      37413 31184 24959 19578 13788 8543 0];
    
    idle_T = idle_thrust_lbf * 4.448;
    mil_T = mil_thrust_lbf * 4.448;
    max_T = max_thrust_lbf * 4.448;
    T_ref = 128992;
end

throttle = max(0, min(1, throttle));
mach_limited = max(0, min(1.4, M));
alt_limited = max(0, min(altitude_range(end), alt_m));

idle_thrust = interp2(altitude_range, mach_range, idle_T, alt_limited, mach_limited, 'linear', 0);
mil_thrust = interp2(altitude_range, mach_range, mil_T, alt_limited, mach_limited, 'linear', 0);
max_thrust = interp2(altitude_range, mach_range, max_T, alt_limited, mach_limited, 'linear', 0);

if throttle <= 0.8
    T_abs = idle_thrust + (throttle/0.8)*(mil_thrust - idle_thrust);
else
    ab_ratio = (throttle - 0.8)/0.2;
    T_abs = mil_thrust + ab_ratio*(max_thrust - mil_thrust);
end

if nargin >= 6 && ~isempty(T_max_ref) && T_max_ref > 0
    scale = T_max_ref / T_ref;
    T_abs = T_abs * scale;
end

T = max(0, T_abs);

TSFC_idle = 0.8;
TSFC_mil = 0.7;
TSFC_max = 1.5;

if throttle <= 0.8
    TSFC = TSFC_idle + (throttle/0.8)*(TSFC_mil - TSFC_idle);
else
    ab_ratio = (throttle - 0.8)/0.2;
    TSFC = TSFC_mil + ab_ratio*(TSFC_max - TSFC_mil);
end

fuel_flow = (T / 9.81) * (TSFC / 3600);

end