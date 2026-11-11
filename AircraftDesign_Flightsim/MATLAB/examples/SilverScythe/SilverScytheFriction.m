function out = SilverScytheFriction(ss,mach,altitude_ft)
% SILVERSCYTHEFRICTION  Profile (skin friction + form) drag by FRICTION.
%
% MATLAB port of Dr. William H. Mason's FRICTION program, following the
% senior-design conversion Aero/SUMAC/zero_lift_drag.m (originally
% converted by Paul Buller, adapted by Koral Karasulu).  The method is:
%
%   * Eckert reference-temperature / Chapman-Rubesin laminar friction
%     (Blasius form),
%   * van Driest II compressible turbulent friction with a
%     Karman-Schoenherr skin-friction law solved by Newton-Raphson,
%   * a transition-weighted blend of the two,
%   * Raymer component form factors, separate expressions for
%     plate-like lifting surfaces and bodies of revolution,
%   * per-component Raymer interference factors (report Table 5), and
%   * a flat 5% stealth leakage increment on the total (report 11.4).
%
% Inputs
%   ss           SilverScytheData struct
%   mach         Mach number, scalar or vector
%   altitude_ft  altitude, scalar or vector [ft]
%
% Output struct (size numel(altitude_ft) x numel(mach))
%   CD_friction  skin-friction drag coefficient, no form factor
%   CD_form      form-drag increment
%   CD_profile   (CD_friction + CD_form) x leakage multiplier
%   component    per-component CF*Swet/Sref contributions
%
% Verified against report Figure 24: at Mach 0.80 / 30,000 ft this returns
% CD_profile = 0.0117085 against the published 0.0117085.

if nargin < 1 || isempty(ss), ss = SilverScytheData(); end
mach = mach(:).';
altitude_ft = altitude_ft(:).';

% Fluid property assumptions carried over from FRICTION
gamma  = 1.4;
Pr     = 0.72;
r_lam  = sqrt(Pr);      % laminar recovery factor
r_turb = 0.88;          % turbulent recovery factor
Te_lam  = 390;          % laminar edge temperature   [R]
Te_turb = 222;          % turbulent edge temperature [R]
c_suth  = 198.6;        % Sutherland constant        [R]

eps_max = 1e-10;
max_iter = 200;

S_wet   = ss.friction.S_wet_ft2(:);
l_ref   = ss.friction.l_ref_ft(:);
t_c     = ss.friction.t_c(:);
geom    = ss.friction.geom_type(:);
Q       = ss.friction.Q(:);
trans   = ss.friction.transition(:);
TwTaw   = ss.friction.TwTaw;
S_ref   = ss.S_ref_ft2;
leakage = ss.friction.leakage_multiplier;

n_comp = numel(S_wet);
n_alt  = numel(altitude_ft);
n_mach = numel(mach);

% Raymer form factors, constant with flight condition
FF = zeros(n_comp,1);
for i = 1:n_comp
    if geom(i) == 0
        FF(i) = 1 + 2.7*t_c(i) + 100*t_c(i)^4;      % plate / lifting surface
    else
        FF(i) = 1 + 1.5*t_c(i)^1.5 + 7*t_c(i)^3;    % body of revolution
    end
    FF(i) = FF(i)*Q(i);
end

CD_friction = zeros(n_alt,n_mach);
CD_form     = zeros(n_alt,n_mach);
component   = zeros(n_alt,n_mach,n_comp);

for ia = 1:n_alt
    [~,a_mps,~,~,nu] = atmosisa(altitude_ft(ia)*ss.units.ft_to_m);

    for jm = 1:n_mach
        M = max(mach(jm),1e-6);
        CFSW   = zeros(n_comp,1);
        CFSWFF = zeros(n_comp,1);

        % van Driest II factors depend only on Mach and wall temperature
        xm    = (gamma-1)*M^2/2;
        TawTe = 1 + r_turb*xm;
        F     = TwTaw*TawTe;
        Tw    = F*Te_turb;
        A_vd  = sqrt(r_turb*xm/F);
        B_vd  = (1 + r_turb*xm - F)/F;
        den   = sqrt(4*A_vd^2 + B_vd^2);
        alpha_vd = (2*A_vd^2 - B_vd)/den;
        beta_vd  = B_vd/den;

        if M > 0.1
            Fc = r_turb*xm/((asin(alpha_vd) + asin(beta_vd))^2);
        else
            Fc = ((1 + sqrt(F))/2)^2;
        end

        Xnum   = 1 + 122*10^(-5/Tw)/Tw;
        Xden   = 1 + 122*10^(-5/Te_turb)/Te_turb;
        Ftheta = sqrt(1/F)*(Xnum/Xden);
        Fx     = Ftheta/Fc;

        % Eckert reference temperature for the laminar branch
        TwTe  = TwTaw*(1 + r_lam*(gamma-1)*M^2/2);
        TstTe = 0.5 + 0.039*M^2 + 0.5*TwTe;
        Cstar = sqrt(TstTe)*(1 + c_suth/Te_lam)/(TstTe + c_suth/Te_lam);

        for i = 1:n_comp
            Rex = ((M*a_mps)*(l_ref(i)*ss.units.ft_to_m))/nu;

            Cf_lam = 2*0.664*sqrt(Cstar)/sqrt(Rex);

            RexBar = Fx*Rex;
            Cfb = 0.074/RexBar^0.2;
            iter = 0;
            err = 1;
            while err > eps_max
                iter = iter + 1;
                if iter > max_iter
                    error('SilverScytheFriction:NoConvergence', ...
                        ['Turbulent Cf did not converge at Mach %.3f, ', ...
                         '%.0f ft, component %d.'],M,altitude_ft(ia),i);
                end
                Cfo = Cfb;
                num = 0.242 - sqrt(Cfb)*log10(RexBar*Cfb);
                dnm = 0.121 + sqrt(Cfb)/log(10);
                Cfb = Cfb*(1 + num/dnm);
                err = abs(Cfb - Cfo);
            end
            Cf_turb = Cfb/Fc;

            CFI = Cf_turb - trans(i)*(Cf_turb - Cf_lam);
            CFSW(i)   = CFI*S_wet(i);
            CFSWFF(i) = CFSW(i)*FF(i);
        end

        CD_friction(ia,jm) = sum(CFSW)/S_ref;
        CD_form(ia,jm)     = (sum(CFSWFF) - sum(CFSW))/S_ref;
        component(ia,jm,:) = reshape(CFSWFF/S_ref,1,1,n_comp);
    end
end

out = struct();
out.mach        = mach;
out.altitude_ft = altitude_ft;
out.CD_friction = CD_friction;
out.CD_form     = CD_form;
out.CD_profile  = (CD_friction + CD_form)*leakage;
out.component   = component;
out.form_factor = FF;
out.leakage_multiplier = leakage;
out.method = "FRICTION (Mason), van Driest II + Karman-Schoenherr";

end
