function [C, takeoff_params, landing_params] = C172Lookup(x, u, geom)

ub = x(4); vb = x(5); wb = x(6);

Vh = max(sqrt(ub^2 + wb^2), 1e-6);
V  = max(sqrt(ub^2 + vb^2 + wb^2), 1e-6);

alpha = atan2(wb, max(ub, 1e-6));
beta  = atan2(vb, Vh);

da = 0; de = 0; dr = 0;
if numel(u) >= 1, da = u(1); end
if numel(u) >= 2, de = u(2); end
if numel(u) >= 3, dr = u(3); end

clamp = @(z,lo,hi) min(max(z,lo),hi);

CL0      = 0.28;
CL_alpha = 5.1;
CL_de    = 0.35;

CD0      = 0.030;
K        = 0.055;

CY_beta  = -0.65;
CY_dr    = 0.18;

Cl_beta  = -0.12;
Cl_da    = 0.08;
Cl_dr    = 0.02;

Cm0      = 0.02;
Cm_alpha = -1.00;
Cm_de    = -1.05;

Cn_beta  = 0.10;
Cn_da    = 0.02;
Cn_dr    = -0.08;

alpha_stall = deg2rad(15);
alpha_eff   = clamp(alpha, -deg2rad(10), alpha_stall);

CL = CL0 + CL_alpha*alpha_eff + CL_de*de;
CD = CD0 + K*CL^2;
CY = CY_beta*beta + CY_dr*dr;

Cl = Cl_beta*beta + Cl_da*da + Cl_dr*dr;
Cm = Cm0 + Cm_alpha*alpha_eff + Cm_de*de;
Cn = Cn_beta*beta + Cn_da*da + Cn_dr*dr;

CD = max(CD, 0);

C = struct('CL',CL,'CD',CD,'CY',CY,'Cl',Cl,'Cm',Cm,'Cn',Cn);

takeoff_params = struct();
takeoff_params.CLmax_takeoff = 1.6;
takeoff_params.CD0_takeoff   = 0.040;
takeoff_params.mu_roll       = 0.03;

landing_params = struct();
landing_params.CLmax_landing = 2.0;
landing_params.CD0_landing   = 0.055;
landing_params.mu_brake      = 0.35;
end
