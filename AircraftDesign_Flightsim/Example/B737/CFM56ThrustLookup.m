function [T, fuel_flow] = CFM56ThrustLookup(throttle, Mach, alt_m, V, rho, max_thrust_ref)
persistent init mach_idle mach_mil alt_vec IdleInterp MilInterp milThrust_N
if isempty(init)
    mach_idle = [0.0 0.2 0.4 0.6 0.8 1.0];
    alt_vec   = [-10000 0 10000 20000 30000 40000 50000 60000];

    IdleMat = [0.0420 0.0436 0.0528 0.0694 0.0899 0.1183 0.1467 0.0;
               0.0500 0.0501 0.0335 0.0544 0.0797 0.1049 0.1342 0.0;
               0.0040 0.0047 0.0020 0.0272 0.0595 0.0891 0.1203 0.0;
               0.0000 0.0000 0.0000 0.0000 0.0276 0.0718 0.1073 0.0;
               0.0000 0.0000 0.0000 0.0000 0.0174 0.0468 0.0900 0.0;
               0.0000 0.0000 0.0000 0.0000 0.0000 0.0422 0.0700 0.0];

    mach_mil = [0.0 0.2 0.4 0.6 0.8 1.0 1.2];

    MilMat = [1.2600 1.0000 0.7400 0.5340 0.3720 0.2410 0.1490 0.0;
              1.1710 0.9340 0.6970 0.5060 0.3550 0.2310 0.1430 0.0;
              1.1500 0.9210 0.6920 0.5060 0.3570 0.2330 0.1450 0.0;
              1.1810 0.9510 0.7210 0.5320 0.3780 0.2480 0.1540 0.0;
              1.2580 1.0200 0.7820 0.5820 0.4170 0.2750 0.1700 0.0;
              1.3690 1.1200 0.8710 0.6510 0.4750 0.3150 0.1950 0.0;
              0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0];

    IdleInterp = griddedInterpolant({mach_idle, alt_vec}, IdleMat, 'linear', 'nearest');
    MilInterp  = griddedInterpolant({mach_mil,  alt_vec}, MilMat,  'linear', 'nearest');

    milThrust_N = 89000.0 * 4.44822;
    init = true;
end

if nargin < 6 || isempty(max_thrust_ref)
    max_thrust_ref = milThrust_N;
end

throttle = max(0, min(1, throttle));
Mach = max(0, min(1.2, Mach));

alt_ft = alt_m / 0.3048;
alt_ft = min(max(alt_ft, alt_vec(1)), alt_vec(end));

idle_factor = IdleInterp(Mach, alt_ft);
mil_factor  = MilInterp(Mach, alt_ft);

idle_factor = idle_factor(1);
mil_factor  = mil_factor(1);

T_idle = max_thrust_ref * idle_factor;
T_mil  = max_thrust_ref * mil_factor;

T = T_idle + throttle * (T_mil - T_idle);
T = T(1);
if ~isfinite(T), T = 0; end
T = max(T, 0);

TSFC_base = 0.657;
fuel_flow = (T / 9.81) * (TSFC_base / 3600);
fuel_flow = fuel_flow(1);
if ~isfinite(fuel_flow), fuel_flow = 0; end

end
