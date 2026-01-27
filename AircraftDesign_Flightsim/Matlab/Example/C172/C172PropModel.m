function [thrust, mdot_fuel] = C172PropModel(throttle, M, alt_m, V, rho)
throttle = max(0.01, min(1, throttle));
V = max(V, 0.1);

P_rated_hp = 160;
rpm_rated = 2700;
alt_ft = alt_m * 3.28084;

rpm_fraction = 0.5 + 0.5*throttle;
rpm = rpm_rated * rpm_fraction;

mp_ideal = 29.92;
if throttle >= 0.95
    mp = mp_ideal;
elseif throttle >= 0.8
    mp = 26 + (mp_ideal - 26) * (throttle - 0.8) / 0.15;
elseif throttle >= 0.6
    mp = 23 + 3 * (throttle - 0.6) / 0.2;
elseif throttle >= 0.4
    mp = 20 + 3 * (throttle - 0.4) / 0.2;
else
    mp = 18 + 2 * throttle / 0.4;
end

if alt_ft <= 0
    hp_sea = interp_power_chart(mp, rpm);
    hp_available = hp_sea;
elseif alt_ft <= 5000
    hp_sea = interp_power_chart(mp, rpm);
    hp_5k = interp_power_chart(mp - 1.5, rpm);
    frac = alt_ft / 5000;
    hp_available = hp_sea * (1 - frac) + hp_5k * frac;
elseif alt_ft <= 10000
    hp_5k = interp_power_chart(mp - 1.5, rpm);
    hp_10k = interp_power_chart(mp - 3.5, rpm);
    frac = (alt_ft - 5000) / 5000;
    hp_available = hp_5k * (1 - frac) + hp_10k * frac;
elseif alt_ft <= 15000
    hp_10k = interp_power_chart(mp - 3.5, rpm);
    hp_15k = interp_power_chart(mp - 6.0, rpm);
    frac = (alt_ft - 10000) / 5000;
    hp_available = hp_10k * (1 - frac) + hp_15k * frac;
else
    hp_15k = interp_power_chart(mp - 6.0, rpm);
    hp_available = hp_15k * exp(-(alt_ft - 15000) / 8000);
end

hp_available = max(hp_available, 5);
P_available = hp_available * 745.7;

D_prop = 1.905;
n_rps = rpm / 60;
J = V / max(n_rps * D_prop, 1e-9);

if J < 0.4
    eta = 0.35 + 0.45*J;
elseif J < 0.8
    eta = 0.53 + 0.32*(J-0.4);
elseif J < 1.2
    eta = 0.66 + 0.12*(J-0.8);
elseif J < 1.6
    eta = 0.71 - 0.05*(J-1.2);
else
    eta = 0.69 - 0.25*(J-1.6);
end
eta = max(0.30, min(eta, 0.75));

if V < 5.0
    v_eff = sqrt(2 * P_available / (rho * pi * (D_prop/2)^2));
    thrust = 2 * rho * pi * (D_prop/2)^2 * v_eff * (v_eff - V * cos(atan2(1, V/max(v_eff,1))));
    thrust = max(0, min(thrust, 2 * rho * pi * (D_prop/2)^2 * v_eff^2));
else
    thrust = eta * P_available / V;
end

thrust = max(0, thrust);

power_frac = hp_available / P_rated_hp;
if power_frac < 0.5
    bsfc = 0.52;
elseif power_frac < 0.75
    bsfc = 0.52 - 0.08 * (power_frac - 0.5) / 0.25;
else
    bsfc = 0.44;
end

mdot_fuel = bsfc * hp_available * 0.453592 / 3600;
end

function hp = interp_power_chart(mp, rpm)
mp = max(18, min(30, mp));
rpm = max(2000, min(2700, rpm));

rpm_points = [2000, 2200, 2400, 2500, 2600, 2700];
mp_points = [18, 20, 22, 24, 26, 28, 29.5];

power_data = [
    50,  60,  70,  80,  90,  100;
    60,  70,  82,  90,  100, 110;
    72,  84,  96,  106, 118, 130;
    85,  99,  113, 124, 138, 152;
    99,  115, 131, 143, 159, 172;
    114, 132, 150, 163, 178, 180;
    125, 145, 165, 175, 180, 180
];

hp = interp2(rpm_points, mp_points, power_data, rpm, mp, 'linear');
hp = max(5, min(180, hp));
end