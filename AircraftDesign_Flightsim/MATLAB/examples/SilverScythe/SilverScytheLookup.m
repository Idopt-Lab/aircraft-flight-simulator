function C = SilverScytheLookup(state_vec,control_vec,geometry,ctx)
% SILVERSCYTHELOOKUP  AAMF coefficient callback for the Silver Scythe.
%
% Signature matches the framework's CoefficientAerodynamics contract, so
% the model is registered exactly like the F-16 and Boeing 737 cases:
%
%   ctx  = SilverScytheContext("high");
%   aero = CoefficientAerodynamics(@(x,u,g) SilverScytheLookup(x,u,g,ctx));
%
% CL, CD and CY are wind-axis force coefficients; Cl, Cm and Cn are
% body-axis moment coefficients about the geometry reference point, which
% SilverScytheBuild places at fuselage station 30.23 ft, the report's
% takeoff CG (report Table 13).  The framework transfers the moment to the
% instantaneous mass CG, so a CG shift from fuel burn is handled by the
% force/moment pipeline rather than by this file.
%
% Control-vector order (see SilverScytheBuild):
%   1 aileron              [rad]   differential wing trailing edge
%   2 elevator             [rad]   symmetric ruddervator
%   3 rudder               [rad]   differential ruddervator
%   4 nozzle_pitch         [rad]   thrust-vector pitch command
%   5 nozzle_yaw           [rad]   thrust-vector yaw command
%   6 throttle left engine
%   7 throttle right engine
%
% Channels 4 and 5 carry no aerodynamic derivative: thrust vectoring is
% modelled physically by rotating each engine's thrust vector inside
% SilverScytheVectoringTurbofan, which is the honest representation of a
% force that scales with thrust rather than with dynamic pressure.
%
% The pitching moment is written in the CL-linked textbook form
%
%   Cm = Cm_ac - SM*CL + Cm_deltae*delta_e + Cm_q*qhat
%
% rather than as Cm_alpha*alpha, so the static margin stays the single
% governing parameter through stall and through the compressibility
% correction on the lift slope.

ss = ctx.ss;
sd = ctx.sd;

x = state_vec(:);
u = control_vec(:);

vel   = x(4:6);
omega = x(10:12);
alt_m = max(-x(3),0);
V     = max(norm(vel),1e-6);
alpha = atan2(vel(3),max(abs(vel(1)),1e-9));
beta  = atan2(vel(2),max(hypot(vel(1),vel(3)),1e-9));

[~,a_sound,~,~] = atmosisa(alt_m);
mach = V/max(a_sound,1e-9);

[S,b,cbar] = geometry.get_reference_geometry(); %#ok<ASGLU>
if ~(b > 0), b = ss.b_ref; end
if ~(cbar > 0), cbar = ss.cbar; end

p_hat = omega(1)*b   /(2*V);
q_hat = omega(2)*cbar/(2*V);
r_hat = omega(3)*b   /(2*V);

da = get_u(u,1);
de = get_u(u,2);
dr = get_u(u,3);

%% Lift
lift = SilverScytheLift(ctx,alpha,mach);
CL = lift.CL + sd.CL_deltae*de + sd.CL_q*q_hat;

%% Drag
drag = SilverScytheDrag(ctx,CL,mach,alt_m);
CD = drag.CD_total + control_drag(sd,da,de,dr);

%% Pitching moment about the geometry reference point
Cm = sd.Cm_ac - sd.static_margin*CL + sd.Cm_deltae*de + sd.Cm_q*q_hat;

%% Lateral-directional
CY = sd.CY_beta*beta + sd.CY_p*p_hat + sd.CY_r*r_hat + sd.CY_deltar*dr;
Cl = sd.Cl_beta*beta + sd.Cl_p*p_hat + sd.Cl_r*r_hat + ...
     sd.Cl_deltaa*da + sd.Cl_deltar*dr;
Cn = sd.Cn_beta*beta + sd.Cn_p*p_hat + sd.Cn_r*r_hat + ...
     sd.Cn_deltaa*da + sd.Cn_deltar*dr;

%% Envelope bookkeeping, in the style of the F-16 lookup
beta_limit = deg2rad(ss.sim.aero_beta_limit_deg);
mach_max   = ss.sim.aero_mach_max;
constraint_names = ["mach_max";"beta_limit_pos";"beta_limit_neg"];
constraint_c = [mach - mach_max; beta - beta_limit; -beta_limit - beta];
envelope_violation = max([0; constraint_c]);

C = struct( ...
    'CL',CL,'CD',CD,'CY',CY,'Cl',Cl,'Cm',Cm,'Cn',Cn, ...
    'control_effects_included',true, ...
    'takeoff_params',ctx.takeoff_params, ...
    'landing_params',ctx.landing_params, ...
    'valid',envelope_violation <= 1e-12, ...
    'constraint_names',{constraint_names}, ...
    'constraint_c',constraint_c, ...
    'envelope_violation',envelope_violation, ...
    'lookup_clamped',mach > ctx_deck_mach_max(ctx), ...
    'fidelity',char(ctx.fidelity), ...
    'debug_alpha_rad',alpha, ...
    'debug_beta_rad',beta, ...
    'debug_mach',mach, ...
    'debug_CL_wing',lift.CL, ...
    'debug_CL_alpha',lift.CL_alpha, ...
    'debug_CL_max',lift.CL_max, ...
    'debug_stalled',lift.stalled, ...
    'debug_CD_profile',drag.CD_profile, ...
    'debug_CD_induced',drag.CD_induced, ...
    'debug_CD_wave',drag.CD_wave, ...
    'debug_e_span',drag.e_span, ...
    'debug_k_induced',drag.k_induced, ...
    'debug_static_margin',sd.static_margin, ...
    'debug_drag_source',char(drag.source));

end

% =======================================================================
function v = get_u(u,idx)
if numel(u) >= idx
    v = u(idx);
else
    v = 0;
end
end

% =======================================================================
function dCD = control_drag(sd,da,de,dr)
% Quadratic control-deflection drag.  The coefficients are the standard
% preliminary-design estimate dCD = k*delta^2 with k tied to each
% surface's own effectiveness, so a surface that produces little moment
% also produces little drag.
dCD = 0.15*abs(sd.Cl_deltaa)*da^2 + ...
      0.25*abs(sd.Cm_deltae)*de^2 + ...
      0.25*abs(sd.Cn_deltar)*dr^2;
end

% =======================================================================
function m = ctx_deck_mach_max(ctx)
if ctx.fidelity == "high" && ~isempty(ctx.deck)
    m = ctx.deck.mach_grid(end);
else
    m = ctx.ss.sim.aero_mach_max;
end
end
