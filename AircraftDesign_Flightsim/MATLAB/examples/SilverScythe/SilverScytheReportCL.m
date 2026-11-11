function out = SilverScytheReportCL(ss,mach,altitude_ft,weight_lb)
% SILVERSCYTHEREPORTCL  The report's own lift-coefficient bookkeeping.
%
% Reproduces report Equations 13 and 14 as implemented in
% Aero/SUMAC/lift_calc.m, purely so that the published drag components of
% report Figure 24 can be reproduced exactly.  It is NOT used by the flight
% model.
%
%   alpha_i   = (W0/(q*S_ref))/(pi*AR)          report Eq 13, radians
%   alpha_eff = alpha_geom + alpha_i - alpha_L0 report Eq 14
%   cl        = a0_2D * alpha_eff               constant across the span
%   ccl/ca    = c(eta)*cl/c_a  ->  LIDRAG  ->   CL
%
% Because c_a = S_ref/b is the true average chord, the spanload integral
% collapses to CL = cl, so the method returns the *two-dimensional* section
% lift coefficient as the aircraft CL.  At the report's cruise point that
% is 0.345 against a level-flight CL of 0.272, i.e. 27% high, because the
% induced-angle term of Equation 13 is added to the geometric angle instead
% of being subtracted from it, and because the 2-D rather than the 3-D lift
% slope is used.
%
% The consequence for the published numbers is confined to induced drag:
% report Figure 24 gives CD_induced = 0.01857 at Mach 0.80 and 30,000 ft,
% whereas k*CL_level^2 with the same LIDRAG k is 0.01194.  Profile drag and
% wave drag are unaffected at that condition.  The flight model uses the
% true CL from the aircraft state, so its cruise drag is lower and its
% predicted maximum Mach correspondingly higher; SilverScytheValidation
% reports both.
%
% Inputs
%   ss           SilverScytheData struct
%   mach         Mach number, scalar or vector
%   altitude_ft  altitude [ft], scalar
%   weight_lb    reference weight, default ss.weights.aero_deck_lb
%
% Output struct
%   CL_report, CL_level, alpha_i_deg, cl_section, ratio

if nargin < 1 || isempty(ss), ss = SilverScytheData(); end
if nargin < 4 || isempty(weight_lb), weight_lb = ss.weights.aero_deck_lb; end

mach = mach(:).';
alpha_geom_deg = 0;    % the deck was generated at zero geometric alpha

a0 = ss.airfoil.nominal_cl_alpha_per_deg;    % [SD] aircraft_params
aL0 = ss.airfoil.nominal_alpha_L0_deg;       % [SD] aircraft_params

[~,a_mps,~,rho] = atmosisa(altitude_ft*ss.units.ft_to_m);
a_fps = a_mps/ss.units.ft_to_m;
rho_slug = rho*0.00194032;

V = mach*a_fps;
q = 0.5*rho_slug*V.^2;

CL_level = weight_lb./(q*ss.S_ref_ft2);
alpha_i_rad = CL_level/(pi*ss.AR);
cl = a0*(alpha_geom_deg + rad2deg(alpha_i_rad) - aL0);

out = struct();
out.mach        = mach;
out.altitude_ft = altitude_ft;
out.weight_lb   = weight_lb;
out.q_psf       = q;
out.CL_level    = CL_level;
out.alpha_i_deg = rad2deg(alpha_i_rad);
out.cl_section  = cl;
out.CL_report   = cl;
out.ratio       = cl./max(CL_level,1e-12);
out.method      = "report Eqs 13-14 via lift_calc.m (bookkeeping only)";

end
