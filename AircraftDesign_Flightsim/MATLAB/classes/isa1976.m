function [T,a,P,rho] = isa1976(h_m)

if ~isscalar(h_m) || ~isreal(h_m) || ~isfinite(h_m)
    error('isa1976:InvalidAltitude', 'Altitude must be a finite real scalar.');
end

h = max(-5000,min(h_m,47000));
g0 = 9.80665;
R = 287.05287;
gamma = 1.4;
hb = [0;11000;20000;32000;47000];
Lb = [-0.0065;0;0.0010;0.0028];
Tb = zeros(numel(hb),1);
Pb = zeros(numel(hb),1);
Tb(1) = 288.15;
Pb(1) = 101325;

for k = 1:numel(Lb)
    dh = hb(k+1)-hb(k);
    if Lb(k) == 0
        Tb(k+1) = Tb(k);
        Pb(k+1) = Pb(k)*exp(-g0*dh/(R*Tb(k)));
    else
        Tb(k+1) = Tb(k)+Lb(k)*dh;
        Pb(k+1) = Pb(k)*(Tb(k+1)/Tb(k))^(-g0/(R*Lb(k)));
    end
end

idx = find(h >= hb(1:end-1) & h <= hb(2:end),1,'last');
if isempty(idx)
    idx = 1;
end

dh = h-hb(idx);
if Lb(idx) == 0
    T = Tb(idx);
    P = Pb(idx)*exp(-g0*dh/(R*T));
else
    T = Tb(idx)+Lb(idx)*dh;
    P = Pb(idx)*(T/Tb(idx))^(-g0/(R*Lb(idx)));
end

rho = P/(R*T);
a = sqrt(gamma*R*T);
end
