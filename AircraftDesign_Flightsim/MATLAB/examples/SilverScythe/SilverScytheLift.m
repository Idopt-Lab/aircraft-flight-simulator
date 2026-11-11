function L = SilverScytheLift(ctx,alpha_rad,mach)
% SILVERSCYTHELIFT  Silver Scythe lift curve for both fidelity levels.
%
% Both fidelities use the compressible swept finite-wing relationship of the
% senior-design stall_solver.m (Nikolai / Raymer form, mid-chord sweep from
% aspect and taper ratio):
%
%   a_3D = a0*cos(L_mid) / ( sqrt(1 + (a0*cos(L_mid)/(pi*AR))^2)
%                            + a0*cos(L_mid)/(pi*AR) )
%
% The fidelities differ in where the 2-D section data come from:
%
%   low   thin-airfoil theory, a0 = 2*pi per radian, and the report's
%         published zero-lift angle of -1.05 deg
%   high  the NACA 64-206 XFOIL polar used by the report, reduced by the
%         same full-range fit as cl_alpha.m (a0 = 0.10232 /deg,
%         alpha_L0 = -1.686 deg)
%
% Compressibility is applied on top as a Prandtl-Glauert correction with
% sweep below Mach 0.9, and as supersonic linear theory (4/beta, swept)
% above Mach 1.1, blended in between.  The report's steady-state analyses
% are all at fixed conditions and did not need this; a 6-DOF model flown
% from the carrier deck to Mach 1.64 does.  Set
% ss.sim.apply_compressibility = false to recover the incompressible slope.
%
% Inputs
%   ctx        SilverScytheContext struct
%   alpha_rad  angle of attack [rad], scalar or vector
%   mach       Mach number, scalar or matching alpha_rad
%
% Output struct
%   CL, CL_alpha [1/rad], alpha_L0 [rad], CL_max, stalled (logical)

ss = ctx.ss;

alpha_rad = alpha_rad(:).';
if nargin < 3 || isempty(mach), mach = 0; end
if isscalar(mach), mach = mach*ones(size(alpha_rad)); end
mach = mach(:).';

a0_per_rad = ctx.lift.a0_per_rad;
alpha_L0   = ctx.lift.alpha_L0_rad;
AR         = ss.AR;
Lambda_mid = deg2rad(ss.midchord_sweep_deg);

a_incompressible = helmbold(a0_per_rad,Lambda_mid,AR);

if ss.sim.apply_compressibility
    factor = compressibility_factor(mach,a_incompressible,AR, ...
        deg2rad(ss.LE_sweep_deg),Lambda_mid);
else
    factor = ones(size(mach));
end

CL_alpha = a_incompressible*factor;
CL_linear = CL_alpha.*(alpha_rad - alpha_L0);

% Smooth stall limiter.  tanh keeps the trim and linearisation numerically
% well behaved while respecting the 3-D CLmax of the clean wing plus the
% LEX vortex-lift credit the report applies in its manoeuvre analysis.
CL_max = ctx.lift.CL_max_clean;
CL = CL_max*tanh(CL_linear/max(CL_max,1e-9));

L = struct();
L.CL        = CL;
L.CL_linear = CL_linear;
L.CL_alpha  = CL_alpha;
L.alpha_L0  = alpha_L0;
L.CL_max    = CL_max;
L.stalled   = abs(CL_linear) > CL_max;
L.compressibility_factor = factor;

end

% -----------------------------------------------------------------------
function a = helmbold(a0,Lambda_mid,AR)
num = a0*cos(Lambda_mid);
r = num/(pi*AR);
a = num/(sqrt(1 + r^2) + r);
end

% -----------------------------------------------------------------------
function f = compressibility_factor(mach,a_inc,AR,Lambda_LE,Lambda_mid)
% Ratio of the compressible to the incompressible 3-D lift slope.
M1 = 0.90;
M2 = 1.10;
f = ones(size(mach));

for i = 1:numel(mach)
    M = max(abs(mach(i)),0);

    if M < M1
        beta = sqrt(max(1 - M^2*cos(Lambda_mid)^2,0.05));
        f(i) = 1/beta;
    else
        % Supersonic linear theory for a swept lifting surface, capped so
        % that the slope stays finite through Mach 1.
        Msup = max(M,1.10);
        beta_sup = sqrt(max(Msup^2 - 1,0.05));
        a_sup = 4/beta_sup*cos(Lambda_LE);
        a_sup = min(a_sup,2*a_inc);
        f_sup = a_sup/a_inc;

        if M >= M2
            f(i) = f_sup;
        else
            beta = sqrt(max(1 - M1^2*cos(Lambda_mid)^2,0.05));
            f_sub = 1/beta;
            w = (M - M1)/(M2 - M1);
            f(i) = (1-w)*f_sub + w*f_sup;
        end
    end
end

f = max(f,0.2);
% Slender-wing ceiling.  Prandtl-Glauert applied to the Helmbold slope keeps
% growing towards Mach 1, but slender-wing theory caps the lift slope of a
% low-aspect-ratio delta at pi*AR/2 = 3.47 per radian independent of Mach.  The
% correction is therefore not allowed to push the slope past that value.
f = min(f,max(1,(pi*AR/2)/max(a_inc,1e-9)));

end
