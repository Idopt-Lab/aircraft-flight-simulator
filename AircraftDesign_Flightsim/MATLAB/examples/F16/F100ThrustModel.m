function out = F100ThrustModel(throttle, M, alt_m, V, rho, max_thrust_N)
% F100THRUSTMODEL  F-100 style thrust lookup for F-16 propulsion.
%
%   Returns a struct compatible with TurbofanPropulsion:
%     out.thrust      - scalar thrust [N]
%     out.fuel_flow   - fuel flow [kg/s]
%     out.direction   - optional body-axis thrust direction []
%     out.moment_body - optional local/reaction moment [N-m]
%
%   Lookup tables are in lbf and converted to N.
%   Thrust is interpolated over Mach number and altitude.
%   throttle:
%     0.0 to 0.8 -> idle to military power
%     0.8 to 1.0 -> military to max/afterburner

persistent mach_idle mach_mil mach_max altitude_range idle_T mil_T max_T

if isempty(mach_idle)
    altitude_range = [-10000, 0, 10000, 20000, 30000, 40000, 50000, 60000] * 0.3048;

    mach_idle = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0];
    idle_thrust_lbf = [
        766   870   941  1236  1602  2107  2612     0;
        890   893   597   969  1420  1868  2391     0;
         71    84    36   484  1060  1588  2143     0;
      -1431 -1431  -998  -422   492  1279  1911     0;
      -3790 -3790 -2668 -1825   844  1546  1603     0;
      -5057 -5057 -1966  -835  -481   983  1425     0];

    mach_mil = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4];
    mil_thrust_lbf = [
       22428 17800 13172  9505  6622  4290  2652     0;
       20844 16633 12412  9007  6321  4112  2545     0;
       20470 16394 12318  9007  6355  4149  2581     0;
       21022 16928 12834  9470  6728  4414  2742     0;
       22393 18156 13922 10360  7423  4895  3026     0;
       24368 19936 15504 11590  8455  5607  3471     0;
       26433 21894 17355 13243  9701  6479  4005     0;
       28375 23852 19331 15041 11178  7548  4682     0];

    mach_max = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6];
    max_thrust_lbf = [
       34267 29000 23732 19219 15312 10893  6748     0;
       32794 27857 22878 18578 14836 10572  6548     0;
       32335 27475 22619 18386 14703 10485  6496     0;
       32738 27809 22890 18618 14889 10617  6578     0;
       33951 28833 23709 19277 15397 10974  6800     0;
       35991 30234 24837 20106 16038 11418  7075     0;
       38373 31947 26447 21602 17220 12241  7583     0;
       41647 34879 29213 23686 18924 13453  8340     0;
       45564 38454 31344 25230 19935 14094  8732     0;
       50176 42278 34385 27585 21736 15339  9503     0;
       53131 45530 37950 30374 23827 16779 10398     0;
       57106 49002 40922 35952 26390 18441 11426     0;
       60028 52200 44370 38860 29000 20880 13340     0;
       63800 55680 47552 41760 31900 23200 15080     0];

    idle_T = idle_thrust_lbf * 4.44822;
    mil_T  = mil_thrust_lbf  * 4.44822;
    max_T  = max_thrust_lbf  * 4.44822;
end

if nargin < 6 || isempty(max_thrust_N)
    max_thrust_N = inf;
end

throttle = max(0, min(1, throttle));
alt_limited = max(altitude_range(1), min(altitude_range(end), alt_m));
M = max(0, M);

mach_idle_q = max(mach_idle(1), min(mach_idle(end), M));
mach_mil_q  = max(mach_mil(1),  min(mach_mil(end),  M));
mach_max_q  = max(mach_max(1),  min(mach_max(end),  M));

idle_thrust = interp2(altitude_range, mach_idle, idle_T, alt_limited, mach_idle_q, 'linear', 0);
mil_thrust  = interp2(altitude_range, mach_mil,  mil_T,  alt_limited, mach_mil_q,  'linear', 0);
max_thrust  = interp2(altitude_range, mach_max,  max_T,  alt_limited, mach_max_q,  'linear', 0);

if throttle <= 0.8
    tau = throttle / 0.8;
    T = idle_thrust + tau * (mil_thrust - idle_thrust);
else
    tau = (throttle - 0.8) / 0.2;
    T = mil_thrust + tau * (max_thrust - mil_thrust);
end

T = max(0, T);
T = min(T, max_thrust_N);

TSFC_idle = 0.8;
TSFC_mil  = 0.74;
TSFC_max  = 2.05;

if throttle <= 0.8
    tau = throttle / 0.8;
    TSFC = TSFC_idle + tau * (TSFC_mil - TSFC_idle);
else
    tau = (throttle - 0.8) / 0.2;
    TSFC = TSFC_mil + tau * (TSFC_max - TSFC_mil);
end

fuel_flow = (T / 9.80665) * (TSFC / 3600);

out = struct();
out.thrust      = T;
out.fuel_flow   = fuel_flow;
out.direction   = [];
out.moment_body = zeros(3,1);
end