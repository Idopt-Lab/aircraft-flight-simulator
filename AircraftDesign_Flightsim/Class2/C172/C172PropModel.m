function [thrust, mdot_fuel] = C172PropModel(throttle, M, alt, V, rho, max_thrust_sl)

throttle = max(min(throttle,1),0);
alt = max(alt,0);
V = max(V,0);

rho0 = 1.225;
sigma = max(min(rho/rho0,1.0),0.05);

V_ref = 55;
V_ratio = min(V / max(V_ref,1e-6), 3.0);

T = max_thrust_sl * sigma;
T = T * (1 - 0.55*V_ratio^1.2);
T = max(T, 0);

thrust = throttle * T;

bsfc = 0.30/3600;
mdot_fuel = bsfc * max(thrust,0);

if ~isfinite(thrust), thrust = 0; end
if ~isfinite(mdot_fuel), mdot_fuel = 0; end
